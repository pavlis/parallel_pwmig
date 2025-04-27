#!/usr/bin/env python3
"""
******
makegclgrid - python replacement for old c++ makegclgrid program
******

.. topic:: Usage

    makegclgrid dbname [-pf pffile] [-c collection] [-f dfile_base_name] [-v]
"""
import os
import sys
import argparse
import numpy as np
from bson import json_util
from mspasspy.db.database import Database
from mspasspy.db.client import DBClient
from mspasspy.ccore.utility import AntelopePf
from pwmigpy.ccore.gclgrid import (GCLgrid,
                                   GCLgrid3d,
                                   r0_ellipse)
from pwmigpy.db.database import GCLdbsave
from mspasspy.client import Client

def main(args=None):
    # As a script that would be run from the shell we let
    # any functions below that throw exception do so and assume they
    # will write a message that can help debug what went wrong
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(prog="makegclgrid",
                    usage="%(prog)s dbname [-pf pffile] [-c collection] [-f dfile_base_name]",
                    description="Build a gclgrid object and store in MongoDB database")
    parser.add_argument('dbname',
                       metavar='dbname',
                       type=str,
                       help='MongoDB database name to save results')
    parser.add_argument('-pf',
                        metavar='pffile',
                        action='store',
                        type=str,
                        default='makegclgrid.pf',
                        help='Set parameter file input file name (default makegclgrid.pf)'
                          )
    parser.add_argument('-c',
                        metavar='collection',
                        action='store',
                        type=str,
                        default='GCLfielddata',
                        help='Set MongoDB collection name for saving results')
    parser.add_argument('-f',
                        metavar='dfile_base',
                        action='store',
                        type=str,
                        default='USE_GRIDNAME',
                        help='Set the base file name for storing data defining this grid (default takes name from the pf file gridname)'
                        )
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='Be more verbose - retrieve and print the document in MongoDB '
                        )
    args = parser.parse_args(args)
    dbname=args.dbname
    client = Client()
    dbclient = client.get_database_client()
    db=Database(dbclient,dbname)
    collection=args.c
    verbose = args.verbose
    
    # This is a translation of the old C++ code
    # necessary for mongodb interaction
    
    pf = AntelopePf(args.pf)
    gridname = pf.get_string('gridname')
    lat0 =pf.get_double('origin_latitude')
    lon0 = pf.get_double("origin_longitude")
    lat0 = np.radians(lat0)
    lon0 = np.radians(lon0)
    r0 = r0_ellipse(lat0)
    azimuth_x = pf.get_double("x_axis_azimuth")
    rotation_angle = np.radians(90.0 - azimuth_x)
    azimuth_y = -rotation_angle
    dx1 = pf.get_double("delta_x1")
    dx2 = pf.get_double("delta_x2")
    n1 = pf.get_long("n1")
    n2 = pf.get_long("n2")
    i0 = pf.get_long("origin_offset_x1")
    j0 = pf.get_long("origin_offset_x2")
    create_3d = pf.get_bool("build_3d_grid")
    dir = pf.get_string("gridfile_directory")
    if not os.path.exists(dir):
        os.makedirs(dir)
    if args.f == 'USE_GRIDNAME':
        dfile_base=gridname 
    else:
        dfile_base=args.f
    if create_3d:
        dx3 = pf.get_double("delta_x3")
        n3 = pf.get_long("n3")
    # this is the 2d grid
    g = GCLgrid(n1, n2, gridname,lat0, lon0, r0, azimuth_y,
                        dx1, dx2, i0, j0)
    doc = GCLdbsave(db,g,collection=collection,dir=dir,dfile=dfile_base+'_2d')
    #GCLdbsave returns the dict assembled to post to the database not MongoDB
    # doc.  
    print('ObjectID of MongoDB doc contents saved for 2d grid')
    print(doc['_id'])
    if verbose:
        print('Document contents saved for 2d grid')
        saveddoc = db[collection].find_one({'_id' : doc['_id']})
        print(json_util.dumps(saveddoc,indent=2))
    if create_3d:
        g3d = GCLgrid3d(n1, n2, n3, gridname, lat0, lon0,r0, azimuth_y,
                        dx1, dx2, dx3, i0, j0)
        doc = GCLdbsave(db,g3d,collection=collection,dir=dir,dfile=dfile_base+'_3d')
        print('ObjectID of doc contents saved for 2d grid')
        print(doc['_id'])
        if verbose:
            print('Document contents saved for 3d grid')
            saveddoc = db[collection].find_one({'_id' : doc['_id']})
            print(json_util.dumps(saveddoc,indent=2))
        

if __name__ == "__main__":
    main()


