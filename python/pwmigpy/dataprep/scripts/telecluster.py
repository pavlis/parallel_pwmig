#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Command line tool to build cluster collection used for stacking converted 
wave data from common common source area.  

This file handles command line argument options to assembled the 
arguments used to run the telecluster python function that does the 
actual work.  i.e. this particular file just the CLI interface.  
To see how telecluster works see python module telecluster.

Usage for CLI tool is:
    
    telecluster dbname [-pf pffile -q source_query -o keylst]

Created on Tue Sep 23 05:29:20 2025

@author: pavlis
"""

import sys
import argparse
import json
from pwmigpy.dataprep.telecluster import telecluster

def main(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = argparse.ArgumentParser(
        prog = "telecluster",
        usage="%(prog)s dbname [-pf pffile -q source_query -o keylist]",
        description="Populate cluster collection for source data",
    )
    parser.add_argument(
        "dbname",
        metavar="dbname",
        type=str,
        help="MongoDB database name containg source collection to be processed",
    )
    parser.add_argument(
        "-pf",
        "pffile",
        action="store",
        type=str,
        default="Telecluster.pf",
        help="Set parameter file that sets up clustering geometry ",
    )
    parser.add_argument(
        "-q",
        "--query",
        action="store",
        type=str,
        default=None,
        help="Specify optional query (json format) to apply to source collection",
    )
    parser.add_argument(
        "-o",
        "--otherdata",
        action="store",
        type=str,
        default=None,
        help="Use with comma separated list of keys for other data to be loaded from source and stored in cluster collection"
    )
    args = parser.parse_args(args)
    # maybe should alter the function to take the mongodb handle not just 
    # the name.   Only issue is exception handling but for now leave  
    # as is - if it ain't broken don't fix it
    dbname = args.dbname
    # The telecluster functionis implemented to use the pfname not the 
    # object constructed from it
    pffile = args.pffile
    # depend on None as default
    if args.query is None:
        query = dict()
    else:
        try:
            query = json.loads(args.query)
        except json.decoder.JSONDecodeError as e:
            print("telecluster:  Error in json specification of query operator with -q option")
            print("string received=",args.query)
            print("Exception message from json module:")
            print(e)
            exit(-1)
    if args.otherdata is None:
        otherdata = []
    else:
        otherdata = args.otherdata.split(",")
    # this main function could throw exceptions but we don't try to 
    # handlle them.  If I (glp) find them too cryptic could put an error 
    # hander here
    Ntotal, Ncluster = telecluster(dbname,
                                   pfname=pffile,
                                   query=query,
                                   othermd=otherdata)
    print("telecluster processed ",Ntotal,
          " source records.  ",Ncluster," of those were added to cluster collection")
    
    
if __name__ == "__main__":
    main()

