#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 11:40:25 2021

@author: pavlis
"""
import os
from bson.objectid import ObjectId

from mspasspy.ccore.utility import (MsPASSError, Metadata)
from mspasspy.util.converter import Metadata2dict
# This is a temporary for testing
from pwmigpy.ccore.gclgrid import (GCLgrid, GCLgrid3d, GCLscalarfield, GCLscalarfield3d,
                     GCLvectorfield, GCLvectorfield3d)


def GCLdbsave_attributes(db, md, collection="GCLfielddata"):
    """
    This function can be used to save attributes of a GCL object
    stored in a Metadata container (md arg).  The save is
    what we called promiscuous in MsPASS because none of the
    attributes are given a type check.  That works here because
    this should only be used in conjunction with the C++ code
    that insulates this python function from type collisions.

    :db: MongoDB database handle
    :md: Metatdata container to be inserted
    :collection:  name of MongoDB collection into which md should
      be added.  (default GCLfieldata)

    :return: id of inserted doc
    """
    dbh = db[collection]
    # Here we always just blindly save all metadata attributes
    # assuming this is used only for GCL objects so there is not
    # real need for schema enforcement.
    if isinstance(md, Metadata):
        doc_to_save = Metadata2dict(md)
        insertion_output = dbh.insert_one(doc_to_save)
        return insertion_output.inserted_id
    else:
        raise MsPASSError(
            "GCLdbsave_attributes:  arg 1 must be a Metadata container", 'Invalid')


def GCLdbsave(db, obj, collection="GCLfielddata",
              dir=None, dfile=None, auxdata=None):
    """
    This procedure saves a GCLgrid object as a file with a parallel save
    of the data attributes to MongoDB collection.

    The family of data objects in the GCLgrid library all have a save
    method that saves the data to a pair of files.  The (usually large)
    arrays used to defined the grid and field data are saved to one
    file (name+".dat") with binary fwrite in C++ and the attributes that
    define the grid and field properties are stored in an text file with the
    ending ".pf".  This function adds a few additional parameters to those
    stored in pf (notably the file names and foff data) and writes those
    to a MongoDB document.   Optionally auxdata can be used to add
    additional attributes to the mongodb document. That is often needed,
    for example, to tag a particular data set or pieces of a data set
    that need to be sorted out later.

    :param db:  required MongoDB database handle. MUST be the
      mspasspy Database child of MongoDB's database handle as we call
      a mspass method in this function.
    :param obj:  GCLgrid related object to be saved (required)
    :param collection:  MongoDB collection name to save the parametric
      data to (default is GCLfieldata).
    :param dir:  file system directory name to write files containing
      the data image to save.  Good practice is to use a fully qualified
      path as the string is copied verbatim to mongodb as the dir
      attribute.  The default is None.  When None the current working
      directory will be used and dir will be set to the full qualified
      path for the current directory.
    :param dfile:  name to use for the file name component of the output
      files.  Note two files will be created with name dfile+".dat" and
      difle+".pf" containing binary data and parametric data respectively.
      The default is None, and is highly recommended.  When set as None
      the file name is created as the grid "name" attribute combined with
      a unique objectid string using Mongodb's ObjectId generator.
      WARNING, however, is that the objectid created will NOT
      match the "_id" of the document created by this call.
      We are just using the ObjectId as a convenience fway to create a uuid.
      Specifically if, for example, the "name" tag in obj was "usarray"
      the dfile name created with default would be something like this:
          usarray_60781ece54ce3ee01bab9449.dat
          and
          usarray_60781ece54ce3ee01bab9449.pf
      Be warned that if you use custom dfile name inputs you must be
      careful to keep the names unique.  The writer will abort
      with RuntimeError exception if the pf file created as dfile+".pf"
      already exists.

    :param auxdata:  optional metadata (or dict - any dict like container
      should work) of additional attributes to be saved in the output
      mongodb document.

    :return:  Metadata container copy of saved attributes with the id of
    the saved record added with key "_id"
    """
    # This use of an exception mechanism is a brutal way to validate that
    # obj is a GCL related object.  If name doesn't resolve it will
    # throw an exception
    try:
        objname = obj.name
    except:
        raise MsPASSError("GCLdbsave:  Invalid data - required arg 2 must by a child of BasicGCLgrid",
                          "Invalid")
    dbh = db[collection]


    # Default to current directory if dir is not defined
    if dir == None:
        cwd = os.getcwd()
        outdir = cwd
    else:
        outdir = os.path.abspath(dir)
    if dfile == None:
        id = ObjectId()
        outfile = objname+"_"+str(id)
    else:
        outfile = dfile

    # The save methods can throw an exception for a number of
    # reasons.  for now we just let them do so and let the caller
    # handle the error if it happens.
    # If successful md will be a metadata container whose contents
    # we will write to mongodb.
    md = obj.save(outfile, outdir, "pfhdr")
    # Now we have to add these attributes that are not part of the
    # object's internal attributes returned in md OR the
    # foff parameters that can only be known during the save
    # that is grid_data_foff and field_data_foff
    md['dir']=outdir
    md['grid_data_file']=outfile
    md['field_data_file']=outfile
    # these are a bit of a kludge.  They are frozen in the C++ library
    # but we need them here to mesh with the C++  api without changing it
    md['field_data_file_extension']='dat'
    md['grid_data_file_extension']='dat'
    # WE add this as a useful more generic name  - it is ignored by reader
    md['dfile']=outfile

    if auxdata != None:
        for k in auxdata:
            # We cannot allow auxdata to set any of the parameters set
            # previously or mysterious results could fullow, most of them bad
            if k in md:
                message = 'GCLdbsave:  key=',k,' is defined in auxdata dict\nNot allowed because that is a keyword in the required parameters'
                raise MsPASSError(message,'Fatal')
            md[k]=auxdata[k]
    id = GCLdbsave_attributes(db, md)
    md['_id']=id
    return md


def GCLdbread(db, id_or_doc, collection="GCLfielddata"):
    """
    These need to be defined in the database:
    dir
    grid_data_file_extension
    grid_data_file
    grid_data_foff
    field_data_file
    field_data_file_extension
    field_data_foff


    save method set only the following (other than object attributes):
    grid_data_foff
    field_data_foff

    """
    dbcol = db[collection]
    if isinstance(id_or_doc, ObjectId):
        doc = dbcol.find_one({"_id": id_or_doc})
    elif isinstance(id_or_doc, dict):
        doc = id_or_doc
    else:
        raise MsPASSError("GCLdbread:  Arg 2 has invalid type - must be objectid or doc",
                          "Invalid")
    md = Metadata(doc)

    # This attribute is is created by calling typeid in the C++ code.
    # The names are demangled to make them readable.  Demangling
    # restores the full namespace typing of the symbols so we
    # get the following somewhat weird incantation.  This would
    # be really ugly and nonportable if we didn't demangle the names
    objtype = md['object_type']
    if objtype == 'pwmig::gclgrid::GCLgrid':
        return GCLgrid(md)
    elif objtype == 'pwmig::gclgrid::GCLgrid3d':
        return GCLgrid3d(md)
    elif objtype == 'pwmig::gclgrid::GCLscalarfield':
        return GCLscalarfield(md)
    elif objtype == 'pwmig::gclgrid::GCLscalarfield3d':
        return GCLscalarfield3d(md)
    elif objtype == 'pwmig::gclgrid::GCLvectorfield':
        return GCLvectorfield(md)
    elif objtype == 'pwmig::gclgrid::GCLvectorfield3d':
        return GCLvectorfield3d(md)
