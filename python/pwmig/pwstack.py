#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file contains the python wrapper version of the pwstack 
program.  pwstack originally a c++ executable.  This set of functions 
allow pwstack to be implemented in a mspass workflow.  The 
algorthms are dogmatically parallel because it is known this program 
is highly parallelizable and preferable to be  run parallel. 
"""
import dask
import math

from mspasspy.ccore.seismic import SeismogramEnsemble,TopMute
from mspasspy.ccore.utility import MsPASSError
from pwmig.ccore.pwmigcore import (RectangularSlownessGrid,
                                   DepthDependentAperture,
                                   pwstack_ensemble)
from pwmig.ccore.gclgrid import GCLscalarfield,GCLdbread

def TopMuteFromPf(pf,tag):
    """
    The C++ class TopMute does not have an AntelopePf driven constructor.
    This function is a front end to work around that in python.  
    Could have been done in C++ but a trivial constuct either way.
    
    :param pf:  AntelopePf that is assumed to have a branch with label tag
    :param tag: string defining branch name to extract
    
    :return: TopMute constructed from pf data

    """
    pfbranch=pf.get_branch(tag)
    mutetype=pfbranch.get_string('mute_type')
    t0=pfbranch.get_double('t0')
    t1=pfbranch.get_double('t1')
    return TopMute(t0,t1,mutetype)
        
class pwstack_control:
    """
    This class is used to contain all the control parameters that 
    define processing with the pwstack algorithm.   The constructor 
    is pretty much like the initialization section of the old 
    pwstack main.   
    """
    def __init__(self,db,pf,slowness_grid_tag='RectangularSlownessGrid',
            data_mute_tag='data_top_mute',
                 stack_mute_tag='stack_top_mute',
                     save_history=False,instance='undefined'): 
        # These are control parameters passed to pwstack_ensemble
        self.SlowGrid=RectangularSlownessGrid(pf,slowness_grid_tag)
        self.data_mute=TopMuteFromPf(pf,data_mute_tag)
        self.stack_mute=TopMuteFromPf(pf,stack_mute_tag)
        use_fresnel=pf.get_bool('use_fresnel_aperture')
        if use_fresnel:
            fvs=pf.get_double("fresnel_vs")
            fvp=pf.get_double("fresnel_vp")
            fcm=pf.get_double("fresnel_cutoff_multiplier")
            fperiod=pf.get_double("fresnel_period")
            fdt=pf.get_double("fresnel_lag_time_sample_interval")
            fnt=pf.get_int("fresnel_number_lag_samples")
            self.aperture=DepthDependentAperture(fvs,fvp,fperiod,fdt,fnt,fcm,True)
        else:
            self.aperture=DepthDependentAperture(pf,"depth_dependent_aperture")
        self.aperture_taper_length=pf.get_double("aperture_taper_length")
        self.centroid_cutoff=pf.get_double("centroid_cutoff")
        self.save_history=save_history 
        self.algid=instance
        # These parameters define thing things outside pwstack_ensemble 
        # but are important control
        self.db=db
        gridname=pf.get_string("pseudostation_grid_name")
        # We query by name and object type - this may require additional 
        # tags to find the unique grid requested but for how keep it simple
        query={'name' : gridname,'object_type' : 'pwmig::gclgrid::GCLgrid'}
        nfound=db.GCLfielddata.count_documents(query)
        if nfound<1:
            raise MsPASSError("GCLgrid with name="+gridname+" was not found in database",
                              "Fatal")
        elif nfound>1:
            raise MsPASSError("GCLgrid with name="+gridname+" is ambiguous - multiple documents have that name tag",
                              "Fatal")
        doc=db.GCLfielddata.find_one(query)
        self.stagrid=GCLdbread(db,doc)
        self.fold=GCLscalarfield(self.stagrid)

def site_query(db,lat,lon,ix1,ix2,cutoff):
    """
    This function is called in in parallel to do a query to the site 
    collection for stations with a center at lat, lon within 
    a distance of cutoff (in degrees). 
    
    :param db:  Database class (top level MongoDB handle)
    :param lat:  latitude (in degrees) of center point
    :param lon: longitude (in degrees) of center point
    :param cutoff:  distance in degrees of circlular region centered on 
       lat,lon to select. 
       
    :return:  python list of ObjectIds of site documents linked to 
      stations within the search area (note this includes all times so 
      a station location can appear more than once. Assume when matched 
      against wf collection this will become one-to-one.

    """
    # cutoff is assumed in degrees but mongo wants that in radians
    rcutoff=math.rad(cutoff)
    query={'coordinates' :
           {'$nearSphere' : 
            { '$geometry' : 
                { 
                    'type' : 'Point',
                    'coordinates' : [ lon , lat ]
                },
                '$maxDistance' : rcutoff
             }
           }
        }
    cursor=db.site.find(query)
    idlist=list()
    for doc in cursor:
        id=doc['_id']
        idlist.append(id)
    result=dict()
    result['idlist']=idlist
    result['lat']=lat
    result['lon']=lon
    result['ix1']=ix1
    result['ix2']=ix2

    return result
def build_wfquery(sid,ridlist):
    """
    This function is more or less a reformatter.  It expands each 
    source_id by a set of queries for a pseudostation grid 
    defined by the ridlist which a nexted data structre (details below).
    The result is a list of dict data structures that are passed down 
    this workflow chain to build ensembles with lat, lon and the 
    each pseudostation and seismogram objects linked to stations 
    inside the search radius.  Those ultimately are passed to the 
    function that actually runs the pwstack algorithm.
    
    :param sid:  expected to contain a single ObjectID the resolves 
      with the source collection and defines the seismic source 
      data to be processed (pwmig is prestack in that sense)
    :param ridlist: is a nested data structure.  It is expected to 
      contain a list of dict structures.  As a minimum each dict for 
      this function is expected to contain data with the key 'idlist'.
      The value associated with idlist is assumed to be a list of 
      ObjectIDs in the site collection that define the set of seismic 
      station data that can be stacked to define each pseudostation grid 
      position.  The function passes along the attributes that define 
      the pseudostation point - lat, lon, and two integer indicec 
      (ix1, and ix2)
    :return: a list of dict data expanded from sid with one entry 
      per pseudostation defined in ridlist.  The contents of each dict 
      are copies of lat, lon, ix1, and ix2 plus  'query' attribute that 
      should be used to query the wf_Seismogram collection to fetch 
      data to be processed.  Note an empty list will be returned immediately 
      if the ridlist received is empty.   That is common with a sparse 
      and/or irregular station geometry.
    """
    if len(ridlist)==0:
        return list()
    else:
        allqdata=list()
        for rid in ridlist:
            # Safer to make a copy of this for a tiny cost
            allq=dict(rid)
            q=dict()
            q['source_id'] = { '$eq' : sid }
            q['site_id'] = { '$in' : rid['idlist']}
            allq['query']=q
            allq.popitem['idlist']
            allqdata.append(allq)
        return allqdata

def read_ensembles(db,querydata):
    """
    Constructs a query from dict created by build_wfquery, runs it 
    on the wf_Seismogram collection of db, and then calls read_ensemble_data 
    on the cursor returned by MongoDB.  It the sets the ensemble 
    metadata for lat, lon, ix1, and ix2 before returning the ensembles
    """
    query=querydata['query']
    n=db.wf_Seismogram.count_documents(query)
    if n==0:
        d=SeismogramEnsemble()
    else:
        cursor=db.wf_Seismogram.find(query)
        d=db.read_ensemble_data(cursor,collection='wf_Seismogram')
    # We always load these into ensemble's metadata even if the ensemble 
    # has no data.  Null data are otherwise ambiguous.   With live data 
    # these are required for running the C++ code pwstack_ensemble
    d.put('pseuostation_lat',querydata['lat'])
    d.put('pseuostation_lon',querydata['lon'])   
    d.put('pseuostation_ix1',querydata['ix1']) 
    d.put('pseuostation_ix2',querydata['ix2'])
    return d
                              
            
def pwstack(db,pf,source_query=None,
        slowness_grid_tag='RectangularSlownessGrid',
            data_mute_tag='data_top_mute',
                 stack_mute_tag='stack_top_mute',
                     save_history=False,instance='undefined'):
    # the control structure pretty much encapsulates the args for 
    # this driver function
    control=pwstack_control(db,pf,slowness_grid_tag,data_mute_tag,
                    stack_mute_tag,save_history,instance)
    if source_query==None:
        base_query={}
    else:
        base_query=source_query
    base_cursor=db.source.find(base_query)
    # We make an assumption here that the array of source ids created 
    # here is tiny and not a memory problem.  Unambiguously true when 
    # source side stack was previously completed. 
    source_id_list=list()
    for doc in base_cursor:
        id=doc['id']
        source_id_list.append(id)
    #TODO  Need to add this method to api
    cutoff=control.aperture.maximum_aperture()
    staids=list()
    for i in range(control.stagrid.nx1):
        for j in range(control.stagrid.nx2):

            # We use this instead of the lat and lon methods for 
            # a minor efficiency difference.  lat and lon have to call 
            # both call this method and return only one of 3 elements
            gc=control.stagrid.geo_coordinates(i,j)
            lat=math.degrees(gc.lat)
            lon=math.degrees(gc.lon)
            # note this function returns a dict.  staids is then a 
            # list of dict 
            ids=dask.delayed(site_query)(db,lat,lon,i,j,cutoff)
            staids.append(ids)
    # I think we need to do for the steps below to work - we need data in 
    # staids
    staids.compute()
    # Now we create a dask bag with query strings for all source_ids and 
    # all staids.  We do that in parallel and put the result in bag 
    # that we will use for downstream processing of this whole workflow
    # i.e. everything from here on is a parallel workflow
    mybag=dask.bag()
    for sid in source_id_list:
        for rids in staids:
            # build_wfquery returns a dict with lat, lon, 
            # i, j, and a query string
            q=dask.delayed(build_wfquery)(sid,rids)
            mybag.append(q)
    # parallel reader - result is a bag of ensembles created from 
    # queries held in query
    mybag.map(lambda q : read_ensembles(q))
    # Now run pwstack_ensemble - it has a long arg list
    mybag.map(lambda d : pwstack_ensemble(d,control.data_mute,
                control.stack_mute,
                  control.stack_count_cutoff,
                    control.tstart,
                      control.tend, 
                        control.aperture,
                          control.aperture_taper_length,
                            control.centroid_cutoff,
                              control.mdlcopy,
                                False,'') )
    mybag.map(lambda d : db.save_ensemble_data(d))