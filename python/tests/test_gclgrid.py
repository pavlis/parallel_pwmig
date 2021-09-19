#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 10:06:01 2021

@author: pavlis
"""
import pickle
from pwmigpy.db.database import (GCLdbsave, GCLdbread)
from pwmigpy.ccore.gclgrid import (GCLgrid, GCLgrid3d, GCLscalarfield,GCLscalarfield3d,
                     GCLvectorfield,GCLvectorfield3d,r0_ellipse,DoubleVector)
import numpy as np
from mspasspy.db.client import DBClient
from mspasspy.db.database import Database 
from mspasspy.util.converter import dict2Metadata
from mspasspy.ccore.utility import MsPASSError

dbclient=DBClient()
db=Database(dbclient,'gcltests')
db.drop_collection('GCLfielddata')

def printmd(md):
    for k in md:
        print(k,md[k])
def comparemd(md1,md2,ignore=['field_data_foff','_id','tag']):
    for k in md1:
        if k in ignore:
            continue
        if md1[k] != md2[k]:
            return False
    return True
def find_and_compare(db,query,md):
    """
    First object matching query in GCLfielddata collection, reads it, 
    and tests that it's attributes match pattern in md.  Used in this
    test to validate database save and readback. 
    """
    doc=db.GCLfielddata.find_one(query)
    md2=dict2Metadata(doc)
    return comparemd(md2,md)
def find_and_load(db,query):
    """
    Convenience function to apply query and return a GCL object
    found in the GCLfieldata collection.  uses find_one so if there
    are multiple copies it could give the wrong answer.

    :param db:  database handle 
    :param query:  dict mongodb query to use on GCLfieldata to 
     get desired data.
     
    Raise a RuntimeError if the query returns a None.

    """
    doc=db.GCLfielddata.find_one(query)
    if doc == None:
        raise RuntimeError("query="+str(query)+" failed and doc return was set to None")
    gclobject = GCLdbread(db,doc)
    return gclobject
    
def find_and_compare_sampledata_3d(db,query,gcldata):
    """
    Read data defined by query from Mongod and compare 
    sample data retrieved with copy in gcldata.  We have to have 
    different testers for different object types - ugly but 
    necessary.  
    
    This is the 3d version - below is the 2d version.  
    Done to keep the code more compact.
    
    Does a ton of assert calls for comparison so it will abort 
    on the first assertion failure.  Returns true if 
    there match is perfect so the call can itself be in an assert

    """
    data_from_db = find_and_load(db,query)
    assert type(data_from_db) == type(gcldata)
    return compare_sampledata_3d(data_from_db,gcldata)
    
def compare_sampledata_3d(gclcopy,gcldata):
    """
    Compares sample data of two gcl objects.  This function 
    assume the scalar attributes of the object have already 
    been tested for equality.  If not seg faults are likely 
    if the array sizes are not compatible because gcl routines 
    do not range check.  
    
    This is the 3d object version.  There is a different 2d
    version.
    
    gclcopy and gcldata are assumed to be the same type.  
    If not peculiar errors may cascade.
    """
    # All types have common coordinate arrays 
    # This assumes the method comparing metadata has already 
    # been called so we can be sure the attributes of the two 
    # objects being compare are equal
    n1=gcldata.n1
    n2=gcldata.n2 
    n3=gcldata.n3
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                # Currently we don't have subscripting for these 
                # arrays so we have to use these getters
                cp1=gcldata.get_coordinates(i,j,k)
                cp2=gclcopy.get_coordinates(i,j,k)
                assert np.isclose(cp1.x1,cp2.x1)
                assert np.isclose(cp1.x2,cp2.x2)
                assert np.isclose(cp1.x3,cp2.x3)
    # Now scan the val data for field objects
    if isinstance(gcldata,GCLscalarfield3d):
        for i in range(n1):
            for j in range(n2):
                for k in range (n3):
                    val1=gcldata.get_value(i,j,k)
                    val2=gclcopy.get_value(i,j,k)
                    assert np.isclose(val1,val2)
    elif isinstance(gcldata,GCLvectorfield3d):
        for i in range(n1):
            for j in range(n2):
                for k in range (n3):
                    val1=gcldata.get_value(i,j,k)
                    val2=gclcopy.get_value(i,j,k)
                    # now val? are vectors so we add this 
                    assert gcldata.nv == len(val1)
                    assert gclcopy.nv == len(val2)
                    assert len(val1) == len(val2)
                    for l in range(len(val1)):
                        assert np.isclose(val1[l],val2[l])
    elif isinstance(gcldata,GCLgrid3d):
        return True
    else:
        raise MsPASSError("find_and_compare_sampledata_3d:  type mismatch - this statement is a safety valve that should never be executed",
                          "Fatal")
    return True
        
def find_and_compare_sampledata_2d(db,query,gcldata):
    """
    Read data defined by query from Mongodb and compare 
    sample data retrieved with copy in gcldata.  We have to have 
    different testers for different object types - ugly but 
    necessary.  
    
    This is the version for 2 objects - the 3d version is above
    and is very similar. It is somewhat of a design flaw that
    they cannot be done easily through a single api.
    
    Does a ton of assert calls for comparison so it will abort 
    on the first assertion failure.  Returns true if 
    there match is perfect so the call can itself be in an assert

    """
    data_from_db = find_and_load(db,query)
    assert type(data_from_db) == type(gcldata)
    return compare_sampledata_2d(data_from_db,gcldata)
def compare_sampledata_2d(gclcopy,gcldata):
    """
    Compares sample data of two gcl objects.  This function 
    assume the scalar attributes of the object have already 
    been tested for equality.  If not seg faults are likely 
    if the array sizes are not compatible because gcl routines 
    do not range check.  
    
    This is the 3d object version.  There is a different 2d
    version.
    
    gclcopy and gcldata are assumed to be the same type.  
    If not peculiar errors may cascade.
    """
    # All types have common coordinate arrays 
    # This assumes the method comparing metadata has already 
    # been called so we can be sure the attributes of the two 
    # objects being compare are equal
    n1=gcldata.n1
    n2=gcldata.n2 
    for i in range(n1):
        for j in range(n2):
            # Currently we don't have subscripting for these 
            # arrays so we have to use these getters
            cp1=gcldata.get_coordinates(i,j)
            cp2=gclcopy.get_coordinates(i,j)
            assert np.isclose(cp1.x1,cp2.x1)
            assert np.isclose(cp1.x2,cp2.x2)
            assert np.isclose(cp1.x3,cp2.x3)
    # Now scan the val data for field objects
    if isinstance(gcldata,GCLscalarfield):
        for i in range(n1):
            for j in range(n2):
                val1=gcldata.get_value(i,j)
                val2=gclcopy.get_value(i,j)
                assert np.isclose(val1,val2)
    elif isinstance(gcldata,GCLvectorfield):
        for i in range(n1):
            for j in range(n2):
                val1=gcldata.get_value(i,j)
                val2=gclcopy.get_value(i,j)
                # now val? are vectors so we add this 
                assert gcldata.nv == len(val1)
                assert gclcopy.nv == len(val2)
                assert len(val1) == len(val2)
                for l in range(len(val1)):
                    assert np.isclose(val1[l],val2[l])
    elif isinstance(gcldata,GCLgrid):
        return True
    else:
        raise MsPASSError("find_and_compare_sampledata_2d:  type mismatch - this statement is a safety valve that should never be executed",
                          "Fatal")
    return True
        
    
# First test default constructors and space allocating constructors
g=GCLgrid()
g=GCLgrid(10,5)
g=GCLgrid3d()
g=GCLgrid3d(5,10,15)
lat0=37.0
lon0=-87.0
r0=r0_ellipse(np.deg2rad(lat0))
g=GCLgrid(10,5,'test2d',np.deg2rad(lat0),np.deg2rad(lon0),r0,
          0.0,10.0,10.0,2,2)
g3=GCLgrid3d(10,5,10,'test3d',np.deg2rad(lat0),np.deg2rad(lon0),r0,
          0.1,10.0,10.0,5.0,2,2)
f1=GCLscalarfield(3,3)
f3=GCLscalarfield3d(3,3,5)
f1=GCLvectorfield(3,3,2)
f3=GCLvectorfield3d(5,5,10,2)
f1=GCLscalarfield(g)
f2=GCLscalarfield3d(g3)
f3=GCLvectorfield(g,2)
f4=GCLvectorfield3d(g3,4)
# Here we fill the field data with data that are just counters
# tests the putters and will allow other followups for operators
count = 1.0
vec=DoubleVector()
# There may be a better way to do this.  DoubleVector is std::vector 
# bound by stock pybind11 code to that name.  
for i in range(f3.nv):
    vec.append(0.0)
for i in range(f1.n1):
    for j in range(f1.n2):
        f1.set_value(count,i,j)
        valtest=f1.get_value(i,j)
        assert np.isclose(count,valtest)
        for l in range(f3.nv):
            vec[l] = count*float(l+1)
        f3.set_value(vec,i,j)
        vec_retrieved=f3.get_value(i,j)
        for l in range(len(vec)):    
            assert np.isclose(vec_retrieved[l],vec[l])
        count += 1.0
vec.clear()
# similar for 3d objects - 
# There may be a better way to do this.  DoubleVector is std::vector 
# bound by stock pybind11 code to that name.  
for i in range(f4.nv):
    vec.append(0.0)
count = 1.0
for i in range(f2.n1):
    for j in range(f2.n2):
        for k in range(f2.n3):
            f2.set_value(count,i,j,k)
            valtest=f2.get_value(i,j,k)
        assert np.isclose(count,valtest)
        for l in range(f4.nv):
            vec[l] = count*float(l+1)
        f4.set_value(vec,i,j,k)
        vec_retrieved=f4.get_value(i,j,k)
        for l in range(len(vec)):    
            assert np.isclose(vec_retrieved[l],vec[l])
        count += 1.0

print('testing save and read of GCLgrid ')
tag={'tag' : 'test1'}
md=GCLdbsave(db,g,dir='/home/pavlis/tmp',auxdata=tag)
print('2d grid Doc contents saved')
printmd(md)
assert find_and_compare(db,tag,md)
assert find_and_compare_sampledata_2d(db,tag,g)
# this tests pickle.  It is a little more expensive 
# to use recover the saved version for comparison instead 
# of the one already in memory, but I was lazy and didn't want 
# to revise the compare funcions
gcopy=pickle.loads(pickle.dumps(g))
assert find_and_compare_sampledata_2d(db,tag,gcopy)


print('\n Testing save and read of GCLgrid3D ')
tag={'tag' : 'test2'}
md=GCLdbsave(db,g3,dir='/home/pavlis/tmp',auxdata=tag)
print('3d grid Doc contents saved')
printmd(md)
assert find_and_compare(db,tag,md)
assert find_and_compare_sampledata_3d(db,tag,g3)

print('\n Testing save and read of 2d scalar field ')
tag={'tag' : 'test3'}
md=GCLdbsave(db,f1,dir='/home/pavlis/tmp',auxdata=tag)
print('3d grid Doc contents saved')
printmd(md)
assert find_and_compare(db,tag,md)
assert find_and_compare_sampledata_2d(db,tag,f1)

print('\n Testing save and read of 3d scalar field ')
tag={'tag' : 'test4'}
md=GCLdbsave(db,f2,dir='/home/pavlis/tmp',auxdata=tag)
print('3d grid Doc contents saved')
printmd(md)
assert find_and_compare(db,tag,md)
assert find_and_compare_sampledata_3d(db,tag,f2)
gcopy=pickle.loads(pickle.dumps(f2))
assert find_and_compare_sampledata_3d(db,tag,gcopy)

print('\n Testing save and read of 2d vector field')
tag={'tag' : 'test5'}
md=GCLdbsave(db,f3,dir='/home/pavlis/tmp',auxdata=tag)
print('3d grid Doc contents saved')
printmd(md)
assert find_and_compare(db,tag,md)
assert find_and_compare_sampledata_2d(db,tag,f3)


print('\n Testing save and read of 3d vectorfield ')
tag={'tag' : 'test6'}
md=GCLdbsave(db,f4,dir='/home/pavlis/tmp',auxdata=tag)
print('3d grid Doc contents saved')
printmd(md)
assert find_and_compare(db,tag,md)
assert find_and_compare_sampledata_3d(db,tag,f4)
gcopy=pickle.loads(pickle.dumps(f4))
assert find_and_compare_sampledata_3d(db,tag,gcopy)

print('Testing copy constructors')
# We make copies of each of the field objects

f1copy=GCLscalarfield(f1)
f2copy=GCLscalarfield3d(f2)
f3copy=GCLvectorfield(f3)
f4copy=GCLvectorfield3d(f4)
# best use this opportunity to validate the copy constructors
assert compare_sampledata_2d(f1copy,f1)
assert compare_sampledata_2d(f3copy,f3)
assert compare_sampledata_3d(f2copy,f2)
assert compare_sampledata_3d(f4copy,f4)
# Now add two copies of each field and compare to the same data *2.0
print('testing operator+=')
f1cpy2=GCLscalarfield(f1)
f2cpy2=GCLscalarfield3d(f2)
f3cpy2=GCLvectorfield(f3)
f4cpy2=GCLvectorfield3d(f4)
f1 += f1copy
f2 += f2copy
f3 += f3copy
f4 += f4copy
print('sums succeeded - now test operator*= scalar and see if this really worked')
f1cpy2 *= 2.0
f2cpy2 *= 2.0
f3cpy2 *= 2.0
f4cpy2 *= 2.0
assert compare_sampledata_2d(f1cpy2,f1)
assert compare_sampledata_2d(f3cpy2,f3)
assert compare_sampledata_3d(f2cpy2,f2)
assert compare_sampledata_3d(f4cpy2,f4)




