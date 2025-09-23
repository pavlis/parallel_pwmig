#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module implements a C++ algorithm from the older pwmig set to python. 
The older C++ program was a command line tool used to build a data structure 
linking sources in an Antelope database into a set of "clusters".  A 
"cluster" in that context is a set of source located inside a radial sector 
defined by a range of distances and station-to-event azimuths.   The result 
was as still is used to drive a source-side stacking algorithm called 
RFEventStacker.  

This module contains all the functional code of telecluster.   A separate 
file in the scripts directory implements the CLI tool using command line 
arguments that make that CLI tool act much like the old C++ code.
"""

from pwmigpy.ccore.seispp import (EventCatalog,Hypocenter,RadialGrid,SectorTest)
from mspasspy.ccore.utility import (Metadata,AntelopePf)
from mspasspy.db.database import Database
from mspasspy.db.client import DBClient
import numpy as np

def dbload_EventCatalog(db,mdlist=[],query={}):
    """
    Creates the working EventCatalog object from MongoDB source 
    collection.   This function is assumed to be run in an interactive 
    environment as it has no safeties. I can throw exceptions for 
    a number of reasons.
    
    :param db:  mspass Database handle containing source collection to be 
       loaded.
    :param mdlist:  python list of keys to be loaded into Metadata 
      section of EventCatalog's internal container. Default is an 
      empty list (all that is needed for telecluster where this 
      translation from C++ is expected to be used for now.  )
    :param query:  optional query to be applied to the source collection.
      default is the find all definition of Mongodb.
    """
    evcat=EventCatalog()
    curs=db.source.find(query)
    for doc in curs:
        md=Metadata()
        # Intentionallty let this throw an exception if the
        # key is missing.  Assuming this function will always 
        # be used interactively.
        for k in mdlist:
            md[k]=doc[k]
        # Always load the source_id as the objectid of this doc
        id=doc['_id']
        md['source_id']=id
        # Load the core coordinate datta to create a hypocenter 
        # object cleanly.  This isn't elegant as we hard code
        # the keys 
        hmd=Metadata()
        hmd['source_lat']=doc['lat']
        hmd['source_lon']=doc['lon']
        hmd['source_depth']=doc['depth']
        hmd['source_time']=doc['time']
        h=Hypocenter(hmd)
        evcat.add(h,md)
    return evcat
def pfload_radial_grid(pf):
    """
    This function recreates the top section of telecluster. 
    It handles what really is an inadequacy of the RadialGrid 
    pf constructor.  There is a switch to use a uniform or nonuniform
    grid that really should just be in the c code, but the cost is 
    tiny to implement here and much faster.
    
    :param pf:  AntelopePf object to parse

    """
    # These are assumed in degrees for regular grid but are 
    # radians internally
    origin_lat=pf.get_double("origin_latitude")
    origin_lon=pf.get_double("origin_longitude")
    use_regular=pf.get_bool("use_regular_grid")
    if use_regular:
        azmin=pf.get_double("grid_minimum_azimuth")
        azmax=pf.get_double("grid_maximum_azimuth")
        naz=pf.get_long("number_grid_points_for_azimuth")
        delmin=pf.get_double("grid_minimum_delta")
        delmax=pf.get_double("grid_maximum_delta")
        ndel=pf.get_long("number_grid_points_for_delta")
        grid=RadialGrid(azmin,azmax,naz,delmin,delmax,ndel,
                                origin_lat,origin_lon)
    else:
        grid=RadialGrid(pf)
    return grid

def cat2dict(evcat):
    """
    Takes an EventCatalog input and returns a dict keyed by 
    str representation of the source_id with associated Hypocenters 
    as values.  These are conveniently written as subdocuments for 
    the cluster collection and are used to compute the hypocentroid 
    of each group.

    """
    hypos=dict()
    evcat.rewind()
    for i in range(evcat.size()):
        h=evcat.current()
        md=evcat.current_aux()
        id=str(md['source_id'])
        hypos[id]=h
        evcat.advance(1)
    return hypos
def compute_centroid(hypos):
    """
    Takes dict returned by function above and computes the hypocentroid 
    of that contents as a Hypocenter object.   Note the lat and lon 
    are in original units (radians here) and time is a (normnally) meaningless 
    mean origin time.  
    """
    centroid=Hypocenter()
    for k in hypos:
        h=hypos[k]
        centroid.lat += h.lat 
        centroid.lon += h.lon 
        centroid.depth += h.depth 
        centroid.time += h.time 
    n=len(hypos)
    centroid.lat /= n
    centroid.lon /= n
    centroid.depth /= n
    centroid.time /=n
    return centroid 

def telecluster(dbname,pfname="telecluster.pf",query={},othermd=[]):
    """
    This is a python function that does approximately the same thing 
    as a C++ program of the same name in the original pwmig package.
    See man page for telecluster in original pqmig for more details.
    
    :param dname:  mongodb database name (string) to read source 
      collection and write to special cluster collection.
    :poram pfname:  AntelopePf file name (sans .pf) that contains 
      control parameters (all control parameters must be in that 
      container or you can expect an exception. 
    :param query:  optional query to apply to the source collection 
      to pass to algorithm.   Default is to load the entire collection.
      (Note this is a pure memory algorithm and all of the collection 
       will be translated and loaded to create associations.)
    :param othermd:  optional list of keys of other md to load with the 
      hypocenter coordinates - these get posted to cluster collection.
    :return: tuple with two componens.  0 component is size of catalog 
      processed and 1 is the number of added to cluster collection 
    
    """
    # First make sure the pf file can be read and contains everything 
    # we need
    pf=AntelopePf(pfname)
    # This parameter is not used in the radial grid constructor so 
    # we fetch it immediately. We let it throw an exception and abort if 
    # it is missing
    gridname=pf.get_string('gridname')
    grid=pfload_radial_grid(pf)
    # Now attempt to load the source data
    dbclient=DBClient()
    db=Database(dbclient,dbname)
    cluster_collection=db['telecluster']
    evcat=dbload_EventCatalog(db,mdlist=othermd,query=query)
    ncluster=0
    for i in range(grid.number_azimuth_bins()):
        for j in range(grid.number_distance_bins()):
            tester=SectorTest(grid,i,j)
            catsubset=evcat.sector_subset(tester)
            nthis=catsubset.size()
            if nthis>0:
                hypos=cat2dict(catsubset)
                centroid=compute_centroid(hypos)
                doc=dict()
                doc['gridname']=gridname
                doc['hypocentroid']={'lat' : np.rad2deg(centroid.lat),
                                     'lon' : np.rad2deg(centroid.lon),
                                     'depth' : centroid.depth,
                                     'time' : centroid.time
                                     }
                # the keys method returns a "dict_keys" object which 
                # may not be groked by mongodb - didn't test but why 
                # gamble for this tiny cost.
                doc['events']=list(hypos.keys())
                celldata=grid.cell(i,j)
                cellsubdoc=dict()
                for k in celldata:
                    cellsubdoc[k]=celldata[k]
                doc['gridcell']=cellsubdoc
                cluster_collection.insert_one(doc)
                ncluster += 1
    return [evcat.size(),ncluster]

