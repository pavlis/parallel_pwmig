import os
from mspasspy.ccore.utility import (AntelopePf,
                    MsPASSError,
                    ErrorSeverity)
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth,kilometers2degrees
import pwmig.db.database as pwmigdb
class GRTray_projector:
    """
    Used to hold parameters and data for accumulation inner loop of pwmig
    done with dask fold or spark accumulate.
    """
    __init__(self,pf,svm0,Vp1d,TPptr,
       rcomp_wt,stack_only):
        self.use_depth_variable_transformation=pf.get_bool("use_depth_variable_transformation")
        self.ApplyElevationStatics=pf.get_bool("apply_elevation_statics")
        self.static_velocity=pf.get_double("elevation_static_correction_velocity")
        self.use_grt_weights=pf.get_bool("use_grt_weights")

def BuildSlownessGrid(grid,source_lat, source_lon, source_depth,model='iasp91',phase='P'):
    model=TauPyModel(model=model)
    Rearth=6378.164
    field=GCLvectorfield(grid)
    for i in range(g.n1):
        for j in range(g.n2):
            stalat=g.lat(i,j)
            stalon=g.lon(i,j)
            georesult=gps2dist_azimuth(source_lat,source_lon,stalat,stalon)
            # obspy's function we just called returns distance in m in element 0 of a tuple
            # their travel time calculator it is degrees so we need this conversion
            dist=kilometers2degrees(georesult[0]/1000.0)
            arrivals=model.get_travel_times(source_depth_in_km=depth,distance_in_degree=dist,phase_list=phase)
            ray_param=arrivals[0].ray_param
            umag=ray_param/Rearth    # need slowness in s/km but ray_param is s/radian
            baz=georesult[2]   # The obspy function seems to return back azimuth
            az=baz+180.0
            # az can be > 360 here but all known trig function algs handl this automatically
            ux=umag*sin(az)
            uy=umag*cos(az)
            # This may take some hacking - not sure this will work with current
            # pybind11 wrappers
            vector_value=[ux,uy]
            field.set_value(vector_value,i,j)

def query_by_id(gridid, db, source_id, collection='wf_Seismogram'):
    """
    Small function used in map below to parallelize mongodb query.
    Querys and returns a cursor of all wf_Seismogram entries matching
    source_id and gridid (keys are those names).  Cursors are the output
    and used in a map to parallelize the query.

    :param gridid:  integer gridid for plane wave component to be selected
    :param db:  Database handle
    :param source_id:  ObjectID of the source to set in query

    """
    query={'source_id': source_id, "gridid" : gridid}
    collection=db[collection]
    return collection.find(query)

def migrate_component(cursor,parent,VPsvm,Vs1d,zmax,tmax,dt,
        number_partitions=None):
    """
    This is an intermediate level function used in the mspass version of
    pwmig.  It will migrate one plane wave component of data created by
    pwstack and return the result as a GCLvectorfield3d object.
    The caller may either immediately stack the plane wave components
    or save them and run a seperate stacker later.  In either case the
    expectation is this function appears as the function object in a
    map operator driven by a bag of cursors. The cursors are assumed
    define a correct and consistent set of data.  No checking is made
    here to verify the consistency of the inputs.

    The code first sets up the raygrid geometry based on input slowness
    data passed through the VPsvm object (here a GCLvectorfield but originally
    this was a smaller thing read from the old pwstack to pwmig file structure).

    :param number_partitions:  Explicitly set the number of partitions.
      Default sets the number of partitions as the size of the grid in the 2
      dimension (field.n2)
    """

    # old pwmig had a relic calculation of a variable ustack and deltaslow here
    # Replacing it here with slowness vectors for each grid point in the
    # function called in the fold operation below.  Note deltaslow is
    # invariant anyway and can be extracted from the metadata of each Seismogram
    # also skipping coherence grid stuff as that was found earlier to not be
    # worth the extra compute time and it would be better done in his framework
    # as an indepndent processing step
    VPVSmax=2.0   # Intentionally large as this is just used to avoid infinite loops where used
    # This sigature of this function has changed from old pwmig - depricated
    # consant u mode
    raygrid=Build_GCLraygrid(parent,PVsvm,Vs1d,zmax,VPVSmax*tmax,dt*VPVSmax)
    seisbag=read_distributed_data(cursor,collection='wf_Seismogram')
    # Only dask fold supports binop that can accumulate a different
    # type than the inputs of the bag.  spark uses the accumulate
    # method for the same purpose.  It seems fold doesn' accept any
    # arguments so he function called here needs to be a class method that
    # allows the input parameters to be defined before beginning the reduce
    migrator=GRTray_projector(svm0,Vp1d,TPptr,ApplyElevationStatics,static_velocity
    use_depth_variable_transformation,rcomp_wt,use_grt_weights,stack_only)
    seisbag.map(migrate_one_seismogram,raygrid,migrator.TPgrid,
      migrator.ApplyElevationStatics,migrator.use_grt_weights,
        migrator.stack_only)
    # The documentation for bag accumulate implies the following will work.
    # This function might be better done as a lambda, but we'll try this initially
    pwdgrid=PWMIGfielddata(raygrid)
    # Don't assume the constructor initializes the data arrays to 0.  It does
    # now but his is a small cost for stabiliy it buys
    del raygrid
    # Set the default partitioning as n2.   We use partitioning to reduce
    # memory usage
    if number_partitions == None:
        nparitions = pwdgrid.n2
    else:
        nparitions = number_partitions
    seisbag.repartition(npartitions=npartitions)
    delayed_data = seisbag.to_delayed()
    for dgroup in delayed_data:
        dlist=dgroup.compute()
        for d in dlist:
            pwdgrid.accumulate(d)
    return pwdgrid

def migrate_event(source_id=None,db,pf='pwmig.pf',collection='GCLfielddata'):
    gclcollection=db[collection]
    # Freeze use of source collection for source_id consistent with MsPASS
    # default schema.
    if source_id=None:
        raise MsPASSError("migrate_event:  usage error.  Must specify source_id",ErrorSeverity.Fatal)
    doc=db.source.find_one({'_id' : source_id})
    if doc == None:
        # This is fatal because expected use is parallel processing of multiple event
        # if any are not defined the job should abort.  Also a preprocess
        # checker is planned to check for such problems
        raise MsPASSError("migrate_event:  source_id="+str(source_id)+" not found in database",ErrorSeverity.Fatal)
    source_lat=doc['lat']
    source_lon=doc['lon']
    source_depth=doc['depth']
    source_time=doc['time']
    control=AntelopePf(pf)
    base_message="class pwmig constructor:  "
    border_pad = control.get_int("border_padding")
    zpad = control.get_double("depth_padding_multiplier")
    if zpad>1.5 or zpad<=1.0:
        message='Illegal value for depth_padding_multiplier={zpad}\nMust be between 1 and 1.5'.format(zpad=          zpad)
        raise MsPASSError(base_message+message,ErrorSeverity.Invalid)
    fielddir=control.get_string("output_field_directory");
    if os.path.exists(fielddir):
        if not os.path.isdir():
            message='fielddir parameter defined in parameter file as {}\n'.format(          fielddir)
            message+='File exists but is not a directory as required'
            raise MsPASSError(base_message+message,ErrorSeverity.Invalid)
    else:
         os.mkdir(fielddir)
    dfilebase=control.get_string("output_filename_base");
    # Following C++ pwmig the output fieldnames, which become file names
    # defined with a dir/dfile combo i a MongoDB collection, are dfilebase+'_'+source_id
    use_depth_variable_transformation=control.get_bool("use_depth_variable_transformation")
    zmax=control.get_double("maximum_depth")
    tmax=control.get_double("maximum_time_lag")
    dux=control.get_double("slowness_grid_deltau")
    dt=control.get_double("data_sample_interval")
    zdecfac=control.get_int("Incident_TTgrid_zdecfac")
    duy=dux
    dumax=control.get_double("delta_slowness_cutoff")
    dz=control.get_double("ray_trace_depth_increment")
    rcomp_wt=control.get_bool("recompute_weight_functions")
    nwtsmooth=control.get_int("weighting_function_smoother_length")
    if nwtsmooth<=0:
        smooth_wt=False
    else:
        smooth_wt=True
    taper_length=control.get_double("taper_length_turning_rays")
    # Parameters for handling elevation statics.
    ApplyElevationStatics=control.get_bool("apply_elevation_statics")
    static_velocity=control.get_double("elevation_static_correction_velocity")
    use_grt_weights=control.get_bool("use_grt_weights")
    use_3d_vmodel=control.get_bool("use_3d_velocity_model")
    # The C++ version of pwmig loaded velocity model data and constructed
    # slowness grids from that.  Here we always use a 3D slowness model
    # loaded from the database using the new gclgrid library that
    # is deriven by mongodb.  Conversion from velocity to slowness and
    # building on from a 1d model is considered a preprocessing step
    # to assure the following loads will succeed.
    # WARNING - these names are new and not in any old pf files driving C++ version
    up3dname=control.get_string('Pslowness_model_name')
    us3dname=control.get_string('Sslowness_model_name')
    Up3d=load_3d_slowness_model(up3dname)
    Us3d=sefl.load_3d_slowness_model(us3dname)
    # Similar for 1d models.   The velocity name key is the pf here is the
    # same though since we don't convert to slowness in a 1d model
    Pmodel1d_name=control.get_string("P_velocity_model1d_name");
    Smodel1d_name=control.get_string("S_velocity_model1d_name");
    Vp1d=load_1d_velocity_model(Pmodel1d_name)
    Vs1d=load_1d_velocity_model(Smodel1d_name)
    # Now bring in the grid geometry.  First the 2d surface of pseudostation points
    parent_grid_name=control.get_string("Parent_GCLgrid_Name");
    # CAUTION:  name may not be the right key
    query={'name' : parent_grid_name}
    doc=gclcollection.find_one(query)
    parent=pwmigdb.GCLdbread(db,doc)
    # For now use defaults here iasp91 and P phase
    # Note pwmig svm0 is a different class that is an abbreviated version of
    # what this function returns.  All functions using svm0 will need to be
    # modified - TODO - WATCH
    svm0=BuildSlownessGrid(parent,source_lat, source_lon, source_depth)
    TPfield=ComputeIncidentWaveRaygrid(parent,border_pad,
       Up3d,Vp1d,svm0,zmax*zpad,tmax,dt,zdecfac,use_3d_model)
    # The loop over plane wave components is driven by a list of cursors
    # created in this MongoDB incantation
    query={'source_id' : source_id}
    gridid_list=db.wf_Seismogram.find(query).distinct('gridid')
    evbag=dask.bag.from_sequence(gridid_list)
    # This loads the bag with cursors for each plane wave component
    planewaves=evbag.map(query_by_id,db,source_id)
    # This function processes one plane wave component and returns a
    # GCLvectorfield3d object containing the migrated data in ray geometry.
    migrated_data=evbag.map(migrate_component,ARGS)
    # Think this needs to be followed by a bag.accumulate call or bag.fold
