import os
import vtk
import sys
# We need this until I figure out how to build a setup.py file for the 
# package
sys.path.append('/home/pavlis/src/parallel_pwmig/python')
from ccore.gclgrid import (GCLscalarfield3d,GCLvectorfield3d)
def GCLfield2vtksg(field,gridname="GCLgrid",fieldnames=None):
    """
    Converts GCLfields (3d scalar of vector only) to vtkStructuredGrid
    that is returned.   Use the C++ code to convert a 2D GCL object 
    for use in paraview.  This function is used to convert GCLfield 
    data to vts (XML based legacy vtk format) file to be read into 
    paraview.  
    
    This function is a translation from a C++ function in GeoParaview
    that was depricated due to problem sorting out the changes in 
    vtk's api.   The python approach was much easier and the 
    performance loss in using python is not significant, in the authors
    opinion, as this is function is not expected to used on large 
    number of GCL objects because a human an only look at a limited 
    numbe.  
    
    :param field:  GCLscalarfield3d or GCLvectorfield3d object to 
     be converted to vtkStructuredGrid container.  
    :param gridname:  overall name to be assigned to the output
    :param fieldnames:  "scalars" names vtk assigns to each data 
      component of the field.  Should be a list or tuple (indexing
      is used so must work with integer indexing) of name strings for each 
      component of the field data at each node.  For a vector field 
      the length should be the length of the number of vector components 
      in the field (nv).  For scalars it should contain only one 
      string.  If the length is less than the field size the names 
      will be set to "UNDEFINED".  If passed as a None (the default) 
      all name tags will be set "UNDEFINED".
      
    :return:  vtkStructuredGrid container of the converted field data. 
      (normally written immediately to an output file by a vtk writer)
    

    """
    if isinstance(field,GCLscalarfield3d):
        nv=1
    elif isinstance(field,GCLvectorfield3d):
        nv=field.nv
    else:
        raise RuntimeError("GCLfield2vtksg:  cannot handle input data of type="+str(type(field)))
    # These vtk incantations are direct translations from the overloaded C++ 
    # function called convert_gcl3d_to_vtksg in GeoParaview
    grid=vtk.vtkStructuredGrid()
    points=vtk.vtkPoints()
    data_tuples=vtk.vtkDoubleArray()
    n1=field.n1
    n2=field.n2
    n3=field.n3
    grid.SetDimensions(n1,n2,n3)
    # This method is only called for vector data - seems to need to be 
    # called before SetNumberOfTuples
    if isinstance(field,GCLvectorfield3d):
        data_tuples.SetNumberOfComponents(nv)
    data_tuples.SetNumberOfTuples(n1*n2*n3)
    data_tuples.SetName(gridname)
    if fieldnames == None:
        for i in range(nv):
            data_tuples.SetComponentName(i,"UNDEFINED")
    else:
        for i in range(nv):
            if i < len(fieldnames):
                data_tuples.SetComponentName(i,fieldnames[i])
            else:
                print("Warning:  fieldnames array length=",len(fieldnames),
                      " is less input vector field vector length=",nv)
                print("Setting component ",i," to UNDEFINED")
                data_tuples.SetComponentName(i,"UNDEFINED")

    # Now translate the grid geometry and data into the vtk structure
    # The weird backward counting for the k index is to make the 
    # standard gcl grid geometry translate to top down order 
    # Translated from old code, but certain that is necessary
    k=n3-1
    for kk in range(n3):
        k=n3-kk-1
        for j in range(n2):
            for i in range(n1):
                p=field.get_coordinates(i,j,k)
                points.InsertNextPoint(p.x1,p.x2,p.x3)
                # array indexing calculation is same for vectors and 
                # scalars - works because an extra arg is present for 
                # SetComponent that must handle that indexing
                offset=(kk*n1*n2)+(j*n1)+i
                if nv==1:
                    val=field.get_value(i,j,k)
                    data_tuples.SetValue(offset,val)
                else:
                    vals=field.get_value(i,j,k)
                    for l in range(nv):
                        data_tuples.SetComponent(offset,l,vals[l])
    grid.SetPoints(points)
    # The C code had a messy double pointer that avoided the temporary
    # variable handle used here.   
    handle=grid.GetPointData()
    handle.SetScalars(data_tuples)
    return grid

def vtkFieldWriter(vsgrid,filename,format="vts",use_binary=True):
    """
    Wrapper for vtk writers for StructuredGrid data.   This function 
    is a companion to GCLfield2vtksg that is in this same module.  
    GCLfield2vtksg converts one of the 3D GCL field times (scalar or vector)
    to a vtk "structured grid".  This function takes the structured grid 
    data created that way and writes the data as an xml vts (default) or 
    legacy vtk file. For the later the default is to save the data in 
    binary, but that can be turned off to save an ascii file.  That is 
    most useful if something is corrupted in the binary file and you 
    need a readable version.
    
    :param vsgrid: is assumed to be a vtkStructedGrid object with the data to 
     be saved.
    :param filename: is the filename to which the data should be saved.  
      This function does no checking for write permission for the path 
      that filename defines.  It can be a simple file name, an absolute 
      path, or a relative path and it should still work.  One feature is 
      that paraview cares about the "extension" of the file name so the 
      function force a ".vts" for xml format files (format='vts') and 
      '.vtk' otherwise.  If the file had a previous extension 
      (something with a trail .*) it will ALWAY be stripped.  Hence you 
      MUST NOT use filenames like "AAK.BHZ.00" if you care that ".00" would
      be stripped.
    :param format:  string defining the write format.  Only two values 
      are currently accepted:  'vts' (default) writes in the xml vtk format, 
      while 'vtk' writes the data in the "legacy vtk" format.  
    :param use_binary:  boolean used only when format is 'vtk'.  When true, 
      which is the default, the data are written in binary for speed. 
      Experience has shown that the ascii version (what you get with false)
      is very slow unless the grid is small.  

    """
    s=os.path.splitext(filename)
    fnbase=s[0]   # Contains filename with extension stripped if it had one
    if(format=='vts'):
        writer=vtk.vtkXMLStructuredGridWriter()
        writer.SetFileName(fnbase+'.vts')
        #writer.SetInputGrid(vsgrid)
        writer.SetInputDataObject(vsgrid)
        writer.Write()
    elif(format=='vtk'):
        writer=vtk.vtkStructuredGridWriter()
        writer.SetFileName(fnbase+'.vtk')
        if use_binary:
            writer.SetFileTypeToBinary()
        writer.SetInputDataObject(vsgrid)
        writer.Write()
    else:
        raise RuntimeError("Cannot handle format="+format+"\nCurrently support vts or vtk")
 

