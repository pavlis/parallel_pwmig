#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 14:33:18 2021

@author: pavlis
"""
import math
from mspasspy.ccore.utility import MsPASSError
from pwmigpy.ccore.gclgrid import (GCLgrid3d,
                                   GCLscalarfield3d,
                                   GCLvectorfield3d,
                                   Geographic_point,
                                   r0_ellipse,
                                   DoubleVector)
from pwmigpy.ccore.pwmigcore import remove_mean_x3
import numpy as np

def increment_index(n,i,increment,dimension_i):
    """
    This small and rather cryptic function handles the problem of generic 
    incrementing of indices in multidimensional array data stored sequentially 
    in files. It is common practice to store such array data in linear 
    text or binary files clocked into a linear structure through nested 
    loops over each array index.  An generic problem is that the order of the 
    loops is not always specified and different authors will produce output 
    data in different scan orders = order of the loops over different 
    indices.   This small function handles that general problem but requires 
    the user to input 4 parameters that allow the algorithm to 
    properly increment a particular, single index of such an array.  
    That means in a read loop this function must be called on each 
    array index after handling each line.  
    
    The input parameter n is the file (or panda array row) counter defining 
    the current block of data for which an increment is to be computed.  
    All array indices depend upon that number.   The "i" argument is the 
    current index value to be incremented and "increment" is a number of 
    lines(rows) between which i is expected to cycle.  Calculating increment 
    for i is always the product of all the array dimensions with faster 
    varying indices.   Hence, for an array a(i,j,k) written in i,j,k 
    order (i.e. i is the fastest varying index = innermost loop) with 
    dimensions n1, n2, and n3 the increments for i, j, and k 
    would be 1, n1, and n1*n2.  dimension_i is the corresponding array 
    dimension for index i (for example above n1, n2, and n3 for i, j, k respectively)
    
    The function returns the index for i that should be assigned to the NEXT 
    datum.  Be careful for off by one errors in using this function 
    in a loop.  It should be used after setting values with the indices 
    that are the inputs to this function.
    """
    if increment == 1:
        return (n+1)%dimension_i
    else:
        if n==0:
            return i
        else:
            itest=(n+1)%increment
            if itest==0:
                j = i+1
                if j == dimension_i:
                    return 0
                else:
                    return j
            else:
                return i

def dftogcl(df,name_index=['lon','lat','depth','vp','hitcount'],
            gridscan_order=[1,0,2],n_header_lines=0,invert_z=True,
            gridname='generic_grid',lat0=60.5,lon0=-142.8,r0=6162.94,azimuth_y=20.0):
    """
    Fairly general function to convert a panda dataframe created from a text 
    file of an earth model defined on a rectilinear grid (uniform grids are subset 
    of rectilinear grids).   Requires name fields for each columns.  
    Should work for any scan order defined by scan_order.  
    
    It computes the size of each coordinate using the dataframe unique method. 
    For this reason expect trouble if the grid is not rectilinear. 
    
    The number of expected attributes in the field output is computed from 
    length of the name_index array - 3 (assumes 3d coordinates are always in 3 columns somewhere)
    
    :param df:  panda dataframe to be converted 
    :param name_index:  list of dataframe names assigned to each column
      The order of the list is CRITICAL.   0 must be the x1 coordinate 
      which here ALWAYS means longitude, 1 must be the x2 which always 
      means latitude, and 2 is expected to be depth.   Fields after that 
      are assumed to all hold data.   If there are fields you want to 
      delete from the output use awk to remove them before reading the 
      data into a panda.   
    :param gridscan_order:  Array of 3 integers defining input scan order 
      in terms of the expected generalized coordinates.   These must 
      be integers 0, 1, or 2.  The order of the array value defines 
      a sequence of decreasing speed of change in the file order.  
      Default is 1,0,2 which means 1 (coordinate x2) is the fastest 
      varying, 0 (coordinate x1) varies next, and the slowest varying 
      coordinate is 2 (coordinate x3)
    :param gridname:  name assigned to the GCLgrid3d object created
      as the frame for the output field
    :param lat0:  origin latitude (in degrees) of coordinate system
    :param lon0:  origin longitude (in degrees) of coordinate system
    :param r0:  origin radius (in km) of coordinate system
    :param azimuth_y:  axis rotation parameter of coordinate gclgrid 
      local coordinate system (in degrees)
     
    :return:  GCLscalarfield3D or GCLvectorfield3d depending on number of 
      attributes (a scalar field by definition has only one attribute)
    
    """
    nv=len(name_index)-3
    if nv<1:
        raise RuntimeError("Illegal name_index parameter input - must have at least 4 components defining dataframe name tags")
    if len(gridscan_order) != 3:
        raise RuntimeError("scan_order array has incorrect size - must be exactly 3")

    # In a panda data frame the name fields define and index.   What 
    # we require here is the names used must be the x1, x2, and x3 
    # coordinate 
    x1 = df[name_index[0]].unique()
    x2 = df[name_index[1]].unique()
    x3 = df[name_index[2]].unique()
    ngridpoints=tuple([len(x1),len(x2),len(x3)])
    n1=ngridpoints[0]
    n2=ngridpoints[1]
    n3=ngridpoints[2]
    # We use these ranges to set up an initial grid geometry that only 
    # needs to be distorted into actual grid.  That is, we use these to 
    # build a "regular" gclgrid that is then distorted into the geometry 
    # required by the panda data in df.  
    minlon=df[name_index[0]].min()
    maxlon=df[name_index[0]].max()
    minlat=df[name_index[1]].min()
    maxlat=df[name_index[1]].max()
    mindep=df[name_index[2]].min()
    maxdep=df[name_index[2]].max()
    dx1nom=np.radians(maxlon-minlon)*r0*math.sin(np.radians(lat0))/n1
    dx2nom=np.radians(maxlat-minlat)*r0/n2
    dx3nom=(maxdep-mindep)/n3
    # int calls requried because python returns a float if division has a remainder
    #grid=GCLgrid3d(n1,n2,n3,
    #        gridname,np.radians(lat0),np.radians(lon0),r0,np.radians(azimuth_y),
    #        dx1nom,dx2nom,dx3nom,int(n1/2),int(n2/2))
    # We need all this to define the coordinate system properly.  We just 
    # alloc the space and then edit the coordinate data.  Very klunky
    # I should have started over with this library
    grid=GCLgrid3d(n1,n2,n3)
    grid.name=gridname
    grid.lat0=np.radians(lat0)
    grid.lon0=np.radians(lon0)
    grid.r0=r0
    grid.azimuth_y=np.radians(azimuth_y)
    grid.set_transformation_matrix()
    grid.dx1_nom=dx1nom
    grid.dx2_nom=dx2nom
    grid.dx3_nom=dx3nom
    i0=int(n1/2)
    j0=int(n2/2)
    k0=int(n3/2)
    grid.set_lookup_origin(i0,j0,k0)
    if nv==1:
        gcf = GCLscalarfield3d(grid)
    else:
        gcf = GCLvectorfield3d(grid,nv)
    
    # Now the tricky part of this algorithm.   We use the scan_order array 
    # in combination with the sizes in the tuple n to compute the repeat
    # cycle of each coordinate.  Then when we read in we use mod arthmetic 
    # to compute the i,j,k indices from the values we set here in the 
    # array called increments. 
    increments=[1,1,1] #initialize 
    # this gridscan_order definition makes this particularly algorithm 
    # remarkably compact
    inc=1
    for j in range(3):
        increments[gridscan_order[j]]=inc
        inc *= ngridpoints[gridscan_order[j]]

    ndata_lines=len(df)
    ndata_lines -= n_header_lines
    i=0
    j=0
    k=0
    # in this loop m is the data row position and n is line number 
    # in the file.  n is used for panda indexing while m is the 
    # count used for array incrementing
    #for m in range(ndata_lines):
    m=0
    for index, row in df.iterrows():
        if m<n_header_lines:
            m+=1
            continue
        londeg = row[name_index[0]]
        latdeg = row[name_index[1]]
        z = row[name_index[2]]
        #DEBUG
        #print(i,j,k,londeg,latdeg,z)
        # GCLgrid library requires geo coordinates be in radians
        gp=Geographic_point()
        gp.lon = np.radians(londeg)
        gp.lat = np.radians(latdeg)
        gp.r = r0_ellipse(gp.lat) - z
        cp = gcf.gtoc(gp)
        if invert_z:
            kk=n3-k-1
        else:
            kk=k
        gcf.set_coordinates(cp,i,j,kk)
        #print(i,j,kk,londeg,latdeg,z,gp.r)
        if nv==1:
            val = row[name_index[3]]
            gcf.set_value(val,i,j,kk)
        else:
            val=DoubleVector()   #  pybind11 construct for std::vector
            for l in range(nv):
                thisval=row[name_index[l+3]]
                val.append(thisval)
            gcf.set_value(val,i,j,kk)
        # set values using current i,j,k - next section increments i,j, and k
        # Now increment each of the loop counters
        i=increment_index(m,i,increments[0],n1)
        j=increment_index(m,j,increments[1],n2)
        k=increment_index(m,k,increments[2],n3)
        m+=1
    gcf.compute_extents()
    return gcf
            
def Velocity3DToSlowness(Vmod3D,ConvertToPerturbation=False):
    """
    Convert a model defined by a GCLscalarfield3d object with velocities 
    as the node values to slowness.   (Note because the operator is 
    just 1/value of each grid value calling this on a model defined as 
    slowness convert it to velocity.)  Optionally remove the mean 
    value of each constant x3 slice to convert the model to an approximate 
    perturbation model.  Note this removes average slowness for a layer 
    which is NOT the same as average velocity.  The old pure C version of 
    pwmig did the opposite (remove mean velocity).
    
    Warning:  this function alters the data of Vmod3D in place to 
    save memory.   
    
    :param Vmod3D the GCLscalar3d object to be converted
    :param ConvertToPerturbation:  when true the mean value is 
      removed from each constant x3 slice to approximate a perturbation 
      model.  Note this works badly on some models near the Moho, but 
      the impact on pwmig depth projects is usually minimal.  default is false
      
    """
    if not isinstance(Vmod3D,GCLscalarfield3d):
        raise MsPASSError("Velocity3DToSlowness:  arg0 must define a GCLscalarfield3d object",
                          "Fatal")

    for i in range(Vmod3D.n1):
        for j in range(Vmod3D.n2):
            for k in range(Vmod3D.n3):
                vel = Vmod3D.get_value(i,j,k)
                u = 1/vel
                Vmod3D.set_value(u,i,j,k)
    if ConvertToPerturbation:
        remove_mean_x3(Vmod3D)
        