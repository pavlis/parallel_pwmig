#!/usr/bin/env python3
"""
This is a python replacement for a C++ program that was part of the original 
pwmig package.   The purpose of the program is to create a 3D velocity 
model that is a projection of an input 1d velocity model.   
The optional ray parameter (-p option) can be used to scale velocities to 
convert an equivalent velocity that would be used to compute times 
using delay times for that slowness.  That is a helpful approximation 
for using a stock moveout for simple ccp stacking.   I retained that 
from the original pwmig although I question the wisdom of using that 
approximation for most types of analysis. 

I did not generalize this as much as might have been desired had I had 
any thoughts this was worth doing.  That is, for now the collection 
name for the grid data is the (default) GCLfieldata and the 
velocity model collection name is frozen as the (default) VelocityModel_1d
collection names.  
"""
import sys
import argparse
from bson import json_util
from mspasspy.db.database import Database
from mspasspy.db.client import DBClient
from pwmigpy.ccore.gclgrid import GCLscalarfield3d
from pwmigpy.db.database import (GCLdbsave,
                                 GCLdbread,
                                 vmod1d_dbread)



def main(args=None):
    # As a script that would be run from the shell we let
    # any functions below that throw exception do so and assume they
    # will write a message that can help debug what went wrong
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(prog="project1dmod",
                    usage="%(prog)s dbname gridname vmodname [-field fieldname] [-mt modtype][-p rayparameter]",
                    description="Loads a specified grid geometry and creates a field of velocities interpolated from a specified velocity model")
    parser.add_argument('dbname',
                       metavar='dbname',
                       type=str,
                       help='MongoDB database name to read model and grid geometry and to save result')
    parser.add_argument('gridname',
                        metavar='gridname',
                        type=str,
                        help='name attribute for grid to be loaded as template for 3d model field')
    parser.add_argument('modelname',
                       metavar='modelname',
                       type=str,
                       help='name attribute for the velocity model to be constructed from the database')
    parser.add_argument('-field',
                        metavar='fieldname',
                        action='store',
                        type=str,
                        default='DERIVED',
                        help='name attribute to be saved with result (default derives name as gridname+modtype with an _ seperator'
                        )
    parser.add_argument('-mt',
                        '--modtype',
                        metavar='modtype',
                        action='store',
                        type=str,
                        default='Pvelocity',
                        help='Define the model property expected as tag for the model being retrieved')
    parser.add_argument('-p',
                        '--rayparameter',
                        metavar='rayparameter',
                        type=float,
                        default=0.0,
                        help='Ray parameter to used for (optional) scaling of velocities - default 0 returns true velocity.  Use a nonzero value for delay time calculations in ccp processing'
                        )
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='Be more verbose - retrieve and print the document saved in MongoDB '
                        )
    parser.add_argument('-cs', '--connection_string',
                        metavar='connection_string',
                        type=str,
                        default='mongodb://localhost:27017',
                        help='MongoDB connection string (default mongodb://localhost:27017)'
                        )
    args = parser.parse_args(args)
    dbname=args.dbname
    dbclient=DBClient(args.connection_string)
    db=Database(dbclient,dbname)
    gridname = args.gridname
    modname = args.modelname
    modtype = args.modtype
    if args.field == 'DERIVED':
        fieldname = gridname + '_' + modtype
    else:
        fieldname = args.field
    rayp = args.rayparameter
    verbose = args.verbose
    query = {'name' : gridname,'object_type': 'pwmig::gclgrid::GCLgrid3d'}
    n = db.GCLfielddata.count_documents(query)
    if n != 1:
        if n>1:
            print('project1dmod(WARNING):  duplicate documents for for query=',query,
                  ' of GCLfielddata collection')
            print('Will use the first one found')
        else:
            print('project1dmod(Fatal):  no document found in GCLfielddata collection matching query=',query)
            exit(-1)
    doc=db.GCLfielddata.find_one(query)
    grid = GCLdbread(db,doc)
    outdir = doc['dir']
    vmquery={'name' : modname, 'property' : modtype}
    n = db.VelocityModel_1d.count_documents(vmquery)
    if n != 1:
        if n>1:
            print('project1dmod(WARNING):  duplicate documents for for query=',
                  vmquery,' of VelocityModel_1d collection')
            print('Will use the first one found')
        else:
            print('project1dmod(Fatal):  no document found in VelocityModel_1d collection matching query=',vmquery)
            exit(-1)
    vmoddoc = db.VelocityModel_1d.find_one(vmquery)
    vmod = vmod1d_dbread(db,vmoddoc)
    fmodel = GCLscalarfield3d(grid)
    del grid
    for i in range(fmodel.n1):
        for j in range(fmodel.n2):
            for k in range(fmodel.n3):
                # The old C++ code called a C++ procedure here called initialize_1Dscalar
                # It has some unnecessary complications.  Here we use use the 
                # interpolation implied in vmod.getv and set the 3d grid 
                # values with the interpolated value. 
                z = fmodel.depth(i,j,k)
                vel =vmod.getv(z)
                fmodel.set_value(vel,i,j,k)
    # this assumes the default collection in GCLdbsave is GCLfielddata
    # outdir was extracted above and will cause this result to be save 
    # in the same directory as tbhe parent grid.  The name assigned was
    # either input as an argument or derived above from the grid's name
    fmodel.name = fieldname
    doc = GCLdbsave(db,fmodel,dir=outdir,dfile=fieldname)
    print('ObjectID of doc contents saved for 2d grid')
    print(doc['_id'])
    if verbose:
        print('Document contents saved for GCLscalarfield3d object holding this model')
        saveddoc = db.GCLfielddata.find_one({'_id' : doc['_id']})
        print(json_util.dumps(saveddoc,indent=2))
                
if __name__ == "__main__":
    main()
