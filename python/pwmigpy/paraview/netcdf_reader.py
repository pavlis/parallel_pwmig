from scipy.io import netcdf
import numpy as np
import sys
# We need this until I figure out how to build a setup.py file for the
# package
#sys.path.append('/home/pavlis/src/parallel_pwmig/python')
from pwmigpy.ccore.gclgrid import (GCLgrid3d,
                     GCLscalarfield3d,
                     GCLvectorfield3d,
                     r0_ellipse,
                     Geographic_point,
                     DoubleVector)
def float_key(this_float, this_format='%0.4f'):
    """" creates a key for a dictionary based on a given value

    Keyword arguments:
        this_float: given float number
        this_format: format to use to construct the key

    Return values:
        the corresponding depth_key string
    """
    return this_format % float(this_float)
def min_max(these_values):
    """find min and max of a list or numpy array

    Keyword arguments:
        these_values: values to search

    Return values:
        this_min: minimum value
        this_max: maximum value
    """
    if type(these_values) is list:
        this_min = min(these_values)
        this_max = max(these_values)
    else:
        this_min = np.min(these_values)
        this_max = np.max(these_values)
    return this_min, this_max
def lon_is_360(lon):
    """checks longitude values to determine if they represent 0/360 range

    Keyword arguments:
        lon: list of longitudes to check

    Return values:
        boolean: True or False
    """
    if isinstance(lon, float) or isinstance(lon, int) or isinstance(lon, np.float64):
        if float(lon) > 180.0:
            return True
        else:
            return False
    else:
        lon_min, lon_max = min_max(lon)
        if lon_min >= 0.0 and lon_max > 180.0:
            return True
        else:
            return False

def lon_180(lon, fix_gap=False):
    """convert longitude values from 0/360 to -180/180

    Keyword arguments:
        lon: list of longitudes to check and convert

    Return values:
        lon: updated longitudes
        lon_map: mapping between old and new longitudes
    """

    lon_360 = lon_is_360(lon)
    lon_map = {}

    lon_start, lon_end = min_max(lon)
    lon_start_index = np.argmin(lon)
    lon_end_index = np.argmax(lon)

    if isinstance(lon, float) or isinstance(lon, int) or isinstance(lon, np.float64):
        if float(lon) > 180.0:
            this_lon = float(lon) - 360.0
            lon_map[float_key(lon)] = this_lon
            lon = this_lon

    else:

        grid = abs(sorted(lon)[1] - sorted(lon)[0])

        for lon_val in lon:
            lon_map[float_key(lon_val)] = float(lon_val)

        if lon_360:
            for i, lon_value in enumerate(lon):
                this_lon = float(lon_value)
                if this_lon > 180.0:
                    lon[i] = this_lon - 360.0
                    lon_map[float_key(lon_value)] = lon[i]
                    lon_map[float_key(lon[i])] = lon[i]
                else:
                    lon[i] = this_lon
                    lon_map[float_key(lon_value)] = lon[i]

        # this logic tries to address the gap at -180/180  or 0/360 longitudes. If the start and end longitudes
        # are within one grid away from each other they are moved closer to closed the gap.
        if fix_gap:
            if 360.0 - abs(lon_end - lon_start) <= grid:
                if lon_360:
                    lon_0 = -0.0009
                    lon_1 = 0.0
                else:
                   lon_0 = -179.9999
                   lon_1 = 179.9999
                lon_map[float_key(lon[lon_start_index])] = lon_0
                lon_map[float_key(lon[lon_end_index])] = lon_1
                lon_map[float_key(lon_start)] = lon_0
                lon_map[float_key(lon_end)] = lon_1
                lon_map[float_key(lon_0)] = lon_0
                lon_map[float_key(lon_1)] = lon_1
                lon[lon_start_index] = lon_0
                lon[lon_end_index] = lon_1
    if fix_gap:
        return lon, lon_map
    else:
        return lon
def show_netcdf_variable_names(model_file):
    data = netcdf.netcdf_file(model_file, 'r')
    names=data.variables.keys()
    print('List of data variable keys in file=',model_file)
    for n in names:
        print(n)
    #data.close()
def read_netcdf_model(model_file, lat_variable='latitude',
                      lon_variable='longitude', depth_variable='depth',
                      extent=False):
    """read in an EMC Earth model in the netCDF format.

    This function is a modification of code originally in the file
    IrisEMC_Paraview_lib.py (part package EMC-paraview).   Original did
    a conversion to the emc's coordinates.  Here we return lat,lon,depth
    with the idea to convert the result to a GCLgrid in separate function.

      Keyword arguments:
      model_file: model file
      ll: lower-left coordinate
      ur: upper-right coordinate
      depth_min: minimum depth
      depth_max: maximum depth
      inc: grid sampling interval
      extend:  when true (default is false) just return the extends of
        a 3d volume - defined as the integer ranges of lat, lon, and depth

      Return values:
      x: x-coordinate  normalized to the radius of the Earth
      y: y-coordinate  normalized to the radius of the Earth
      z: z-coordinate  normalized to the radius of the Earth
      meta: file metadata information
    """


    # NetCDF files, when opened read-only, return arrays that refer directly to memory-mapped data on disk:
    data = netcdf.netcdf_file(model_file, 'r')
    variables = []
    for name in list(data.variables.keys()):
        if name not in (depth_variable, lon_variable, lat_variable):
            variables.append(name)
    attribute_names=list(variables)

    # expects variables be a function of latitude, longitude and depth, find the order
    var = variables[0]
    for i, value in enumerate(data.variables[var].dimensions):
        if value == depth_variable:
            depth_index = i
        elif value == lon_variable:
            lon_index = i
        else:
            lat_index = i

    latitude = data.variables[lat_variable][:].copy()
    longitude = data.variables[lon_variable][:].copy()
    longitude, lon_map = lon_180(longitude, fix_gap=True)

    depth = data.variables[depth_variable][:].copy()

    # select the values within the ranges (this is to get a count only)
    #latitude, longitude, depth2 = get_points_in_volume(lat, lon, depth, ll, ur, inc, depth_min, depth_max)

    # model data grid definition
    V = {}
    nx = len(longitude)
    ny = len(latitude)
    nz = len(depth)
    if extent:
        return nx - 1, ny - 1, nz - 1
    dimensions=list()
    dimensions.append(nx)
    dimensions.append(ny)
    dimensions.append(nz)
    # I think this should always work - gets number of attributes
    nv=len(variables)
    dimensions.append(nv)
    index = [-1, -1, -1]

    missing_value = None
    if hasattr(data.variables[var], 'missing_value'):
        missing_value = float(data.variables[var].missing_value)

    # The IRIS EMC code had what seems to me to be an unnecessary
    # complexity. All the data files we have for the alaska paper have
    # what paraview, at least, would call a rectilinear grid.  That means
    # we only need to store a vector of coordinates.   EMC stuff stored
    # 3D arrays that were a factor in slowing the code enormously.
    # see the EMC code if you need to restore the 3d arrays.
    # They also called these X,Y,Z.  I found the same data was
    # already store din latitude, longitude, and depth vectors

    for l, var_val in enumerate(variables):

        v = np.zeros((nx, ny, nz))
        data_in = data.variables[var_val][:].copy()

        # increment longitudes, we want to keep the first and last longitude regardless of inc
        for i, lon_val in enumerate(longitude):
            for j, lat_val in enumerate(latitude):
                for k, depth_val in enumerate(depth):


                        index[depth_index] = k
                        index[lat_index] = j
                        index[lon_index] = i
                        this_value = data_in[index[0]][index[1]][index[2]]
                        #print(x,y,z,this_value)

                        if this_value is None:
                            v[i, j, k] = None
                        elif this_value is not None:
                            if this_value == missing_value:
                                v[i, j, k] = None
                                this_value = None
                            else:
                                v[i, j, k] = this_value
                        else:
                            v[i, j, k] = this_value

        V[var_val] = v

    data.close()
    return longitude, latitude, depth, V, attribute_names, dimensions


def netcdf_arrays_to_GCLfield(lon,lat,depth,V,vkey=None,
        lat0=np.deg2rad(44.967243394444445),
            lon0=np.deg2rad(-103.77671383333333),
                r0=r0_ellipse(np.deg2rad(44.967243394444445)),
                    azimuth_y=0.0,
                        i0=0,
                            j0=0,
                                gridname='UNDEFINED',
                                    save_as_vector=False):

    """
    This function is a compansion to the read_netcdf_model.  For conversion
    to a gcl field it would normally be called immediately after
    read_netcdf_model using the set of numpy arrays it returns.
    The result can be guaranteed only if the scan order for the
    arrays is [longitude low to high][latitude low to high ] [z low to high]
    Note z as depth is reversed for the GCL structure so we fill the cells
    backward from the EMC netcdf convention.  A warning is issued if
    there are hints those conditions are not satisified.

    This function is nearly guaranteed to be super slow because of all
    the array indexing and nexted loops.

    :param lon,lat,depth: are coordinate arrays returned by read_netcdf_model.
      they are assume to be (in order) longitude, latitude, and depth
      coordinate 3D arrays.  For a uniform grid in lat,lon, depth these
      could be computed from 6 numbers, but to assure some flexibility they are not
      computed but carried as baggage.
    :param V:  value dict returned by read_netcdf_model.  This is a bit of
      a weird data structure inherited from adapting the original EMC code.
      V is return as a dict with np arrays linked to one or more keys.
      This parameter is ignored if save_as_vector is True. In the False
      case it uses this name to find the component desired if it is defined.
      The default sets it None.  In the default situation the function
      retrieves the first attribute listed if there are multiple attributes.
      For true scalar data that is the only attribute so it should work fine.
    :param lat0, lon0, r0: are the coordinates of the origin of the RegionalCoordinates
      Cartesian system used to define the local coordinate system. As
      describe in the origin Fan et al paper the origin gets translated from
      the center of earth to this point on the earth. The values passed
      should be in units of degrees for lat0 and lon0 and km for
      radius.  (Note GCL internally converts lat0 and lon0 to radians so
      be careful you poke at the attributes of what is created)
      The default is the geographical center of the US at the surface
      defined at the latitude by r0_ellipse.
    :param azimuth_y: is the rotation angle (in degrees) to apply to the
      coordinate system as an azimuth angle from north.  i.e. for 0 the
      n2 axis points north, but if azimuth_y were 10 degrees the
      n2 axis would point 10 degrees east of north at the origin.
    :param i0:  offset of the origin from coordinate center in
      grid cell units - a GCLgrid thing that is baggage here (default is 0)
    :param j0:  like i0 for the n2 axis (default 0)
    :param gridname:  Name to assign the grid geometry. (Default is undefined)
    :param save_as_vector:  When true data with multiple attributes for
      each grid cell (e.g. a P and S tomography model ).   there is no
      clean way to return the names because of a mismatch in concept of
      the GCLvectorfield and what netcdfs are commonly used for.
      If names are needed they must be retained from V and passed to
      the vtk converter where they can be displayed in paraview.
      Awkward but I can see now simple solution to resolve that problem
      without some major retooling.

    :return:  a GCLscalarfield3D f save_as_vector is false and a
      GCLvectorfield3d if True.

    """
    # this function is expected to be run interactively so these print
    # warnings using this name tag are appropriate
    prog='netcdf_arrays_to_GCLfield:  '
    nx=len(lon)
    ny=len(lat)
    nz=len(depth)
    # these cell size estimates are approximate and only close if the
    # grid is uniform in lat,lon, and depth -I think all the EMC
    # models are but unclear if all netcdfs would be
    latmin=np.amin(lat)
    latmax=np.amax(lat)
    lonmin=np.amin(lon)
    lonmax=np.amax(lon)
    zmax=np.amax(depth)
    zmin=np.amin(depth)
    delta_lon=(lonmax-lonmin)/(nx-1)
    delta_lat=(latmax-latmin)/(ny-1)
    dz=(zmax-zmin)/(nz-1)
    # gclgrid wants dx and dy in units of km so we have do some rough conversion
    center_lat=(latmax-latmin)/2
    radius=r0_ellipse(np.deg2rad(center_lat))
    dy=radius*np.deg2rad(delta_lat)
    dx=radius*np.deg2rad(delta_lon)*np.cos(abs(np.deg2rad(center_lat)))
    g=GCLgrid3d(nx,ny,nz)
    g.name=gridname
    g.r0=r0
    g.lat0=np.deg2rad(lat0)
    g.lon0=np.deg2rad(lon0)
    g.azimuth_y=np.deg2rad(azimuth_y)
    g.set_transformation_matrix()
    g.dx1_nom=dx
    g.dx2_nom=dy
    g.dx3_nom=dz
    # These don't matter much for this context but we set them as follows. 
    # This would be problematic for a global model
    i0=int(nx/2)
    j0=int(ny/2)
    k0=int(nz/2)
    g.set_lookup_origin(i0,j0,k0)
    #g=GCLgrid3d(nx,ny,nz,gridname,
    #        np.deg2rad(lat0),np.deg2rad(lon0),r0,np.deg2rad(azimuth_y),
    #           dx,dy,dz,i0,j0)
    if save_as_vector:
        # V is expected to be a dict with name keys. If we make
        # this a vector field we can set the dimensions from len like this
        nv=len(V)
        keys=list(V.keys())
        f=GCLvectorfield3d(g,nv)
        for k in range(nz):
            kk=nz-k-1
            for j in range(ny):
                for i in range(nx):
                    gp=Geographic_point()
                    gp.lon=np.deg2rad(lon[i])
                    gp.lat=np.deg2rad(lat[j])
                    r0ij=r0_ellipse(gp.lat)
                    gp.r=r0ij-depth[k]
                    cp=f.gtoc(gp)
                    f.set_coordinates(cp,i,j,kk)
                    val=DoubleVector()
                    for l in range(nv):
                        # This is obscure syntax for sure
                        keyl=keys[l]
                        val.append(V[keyl][i][j][k])
                    f.set_value(val,i,j,kk)


    else:
        if vkey:
            v=V[vkey]
        else:
            keys=list(V.keys())
            v=V[keys[0]]
        f=GCLscalarfield3d(g)
        for k in range(nz):
            kk=nz-k-1
            for j in range(ny):
                for i in range(nx):
                    gp=Geographic_point()
                    gp.lon=np.deg2rad(lon[i])
                    gp.lat=np.deg2rad(lat[j])
                    r0ij=r0_ellipse(gp.lat)
                    gp.r=r0ij-depth[k]
                    cp=f.gtoc(gp)
                    f.set_coordinates(cp,i,j,kk)
                    #print(i,j,k,kk,X[i][j][k],Y[i][j][k],Z[i][j][k],cp.x1,cp.x2,cp.x3)
                    f.set_value(v[i][j][k],i,j,kk)
    # This method is needed to set the bounding box attributes
    f.compute_extents()
    return f
