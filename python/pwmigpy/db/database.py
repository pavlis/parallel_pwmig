
import os
from bson.objectid import ObjectId

from mspasspy.ccore.utility import (MsPASSError, Metadata,AntelopePf)
from mspasspy.util.converter import Metadata2dict
from pwmigpy.ccore.gclgrid import (GCLgrid, GCLgrid3d, GCLscalarfield, GCLscalarfield3d,
                     GCLvectorfield, GCLvectorfield3d)
from pwmigpy.ccore.seispp import VelocityModel_1d
import pandas as pd


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
    # These attributes are not needed when we only store the grid data
    # Hence we load them only for field data.  Note a not of testing 
    # for grid doesn't work because of inheritance 
    if isinstance(obj,GCLscalarfield) or isinstance(obj,GCLscalarfield3d) or isinstance(obj,GCLvectorfield) or isinstance(obj,GCLvectorfield3d):
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

# This set of procedures are used for managing layered earth models
# encapsulated in the seispp object called VelocityModel_1d.
def read_1d_model_file(fname,format='plain',property='P',
     depth_key='depth',gradient_key=None,model='ak135'):
    """
    Use this function to read a text file defining a layered earth
    model of any kind.   The result will be returned as a VelocityModel_1d
    object but any scalar property could be return in the v vector of the object.
    The format argument defines the format of the text file to be read.
    The idea is this function can be readily hacked to support new formats.
    The currently supported list of formats is defined below:

    :param fname:  string defining the file name containing the data to be read.
    :param format:  string defining the format of the file.  The following is a
      list of supported formats:
        plain - text file of just z, velocity, and grad values separated by
          white space.  This format is a legacy format form seispp and the
          original pwmig
        rbh - Herrmann's velocity model format.  Similar to plain but with a
          lengthy header
        csv - a standard csv file BUT this file must contain header lines
          that are consistent with the values passed through the arguments
          property, depth_key, and gradient_key.
    :param property: physical property to be loaded.  For plain or rbh format
      must be either "P" or "S".  For csvit must match a column heading name
      for the column of data to be loaded.  for mod1d must be the name used 
      in the property field for the antelope mod1d table.  
      ('Pvelocity' or 'Svelocity' for all examples I know)
    :param depth_key:  header name for depth field in a csv file.  Ignore for
      rbh and plain formats.  (default "P")
    :param gradient_key:  some models use gradients between grid points in the
      model definition.   This parameter is referenced for csv file and when
      defined the column with which it is associated is loaded as the gradient
      for the model BELOW the comparable depth point.   Default is None which
      means grad will be set zero (constant velocity layer model)
    :param model:  model name - used only for mod1d format to select 
      subset of mod1d table for a particular named model 

    """
    if format == 'plain' or format == 'rbh':
        # These two formats are handled by the VelocityModel_1d constructor
        # That contructor can throw a MsPASSError exception but we assume
        # a usage where normally this function would not be called deep
        # within  large parallel job
        vmod = VelocityModel_1d(fname,format,property)
    elif format == 'csv':
        # This example illustrates a generic way to read a text file and
        # create a VelocityModel_1d object.  This case uses pandas and
        # requires a header with a label matching property
        df = pd.read_csv(fname)
        z = df[depth_key]
        v = df[property]
        if gradient_key != None:
            grad = df[gradient_key]
        else:
            grad = []
            for k in range(len(v)):
                grad.append(0.0)
        # In a dataframe v an z can always assumed to be the same length
        vmod = VelocityModel_1d(len(z))
        vmod.nlayers = len(z)
        for i in range(len(z)):
            vmod.z.append(z[i])
            vmod.v.append(v[i])
            vmod.grad.append(grad[i])
    elif format == 'mod1d':
        # These names and columns are locked to Antelope schema css3.0
        df=pd.read_csv(fname,sep='\s+',
                names=['model_name','property','depth','value','gradient','units','auth','lddate'])
        # panda subset with match to model_name and property columns
        df = df.loc[lambda dfl: dfl['model_name'] == model]
        df = df.loc[lambda dfl: dfl['property'] == property]
        vmod = VelocityModel_1d()
        for index, row in df.iterrows():
            vmod.z.append(row['depth'])
            vmod.v.append(row['value'])
            vmod.grad.append(row['gradient'])
        vmod.nlayers = len(vmod.z)
    else:
        raise MsPASSError('read_1d_model_file:  unsupported format name with tag='+format,'fatal')

    return vmod

def vmod1d_dbsave(db,vmod1d,model_name,property='Pvelocity',collection="VelocityModel_1d"):
    dbcol = db[collection]
    doc={'property' : property, 'name' : model_name}
    doc['nlayers'] = vmod1d.nlayers
    # the C++ object here uses an std::vector container.  
    # mongodb doesn't support that native so we have convert the 
    # vectors to python arrays
    z=[]
    v=[]
    grad=[]
    for i in range(vmod1d.nlayers):
        z.append(vmod1d.z[i])
        v.append(vmod1d.v[i])
        grad.append(vmod1d.grad[i])
    doc['depth'] = z
    doc['velocity'] = v
    doc['gradient'] = grad
    ret = dbcol.insert_one(doc)
    return ret.inserted_id

def vmod1d_dbread(db,id_or_doc,collection='VelocityModel_1d'):
    dbcol = db[collection]
    if isinstance(id_or_doc, ObjectId):
        doc = dbcol.find_one({"_id": id_or_doc})
    elif isinstance(id_or_doc, dict):
        doc = id_or_doc
    else:
        raise MsPASSError("vmod1d_dbread:  Arg 2 has invalid type - must be objectid or doc",
                          "Invalid") 
    vmod = VelocityModel_1d()
    z=doc['depth']
    v=doc['velocity']
    grad=doc['gradient']
    # Assume z,v,and grad are same length.  Let python throw 
    # an exception if they aren't - should not happen anyway
    vmod.nlayers = doc['nlayers']
    for i in range(vmod.nlayers):
        vmod.z.append(z[i])
        vmod.v.append(v[i])
        vmod.grad.append(grad[i])
    return vmod
    

# This set of functions are used to implement a verify procedure on
# parametric data that can be cast into the tree structure of an AntelopePf.
# With these procedures pf data can be stored and retrieved from
# MongoDB and verified against a set of test restrictions.
def pf_dbsave(db,pf,name_tag, collection='AntelopePf'):
    """
    Save an AntelopePf object to MongoDB as document. If

    The data structure of an AntelopePf is actually not dependent upon the
    format invented by Dan Quinlan of BRTT they call a "pf file".
    A Pf basically defines a tree structured set of data indexed with key-value
    pairs.  That maps exactly into a MongoDB document when subdocuments of
    any level are allowed.  The one item in a Pf that does not mesh perfectly
    is a Tbl which in the API is converted to an std::list of strings.  The
    consumer of the data in a Tbl is responsible for parsing the data in that
    list of strings.

    This function saves an AntelopePf as a document in a collection defined by
    the collection argument.   No schema is enforced on what is saved
    the AntelopePf is strongly typed.

    :param db:  MsPASS database handle (a parent MongoDB database handle should work also)
    :param pf:  The AntelopePf object that is to be saved.
    :param name_tag:   This is required to be a string that is saved with
      the frozen key name of 'pfname'.   (e.g. if name_tag is 'pwmig' if
      if successful this function will create an entry {'pfname' : 'pwmig'}
      in the document stored by that call to this function).  The AntelopePf
      must not contain the keyword pfname.
    :param collection:  MongoDB collection to which the data should be stored.
      (default is 'AntelopePf')

    :return: ObjectID of the document saved.
    """
    # First we do a few sanity checks.
    if not isinstance(name_tag,str):
        raise MsPASSError('pf_dbsave:  name_tag argument must be a string','Fatal')
    # All the work is done in this function.  It returns a dict but we call it
    # doc only because we will treat it that way when we save it to Mongodb
    doc = pf_to_dict(pf)
    # We can't allow doc to contain the keyword pfname
    if 'pfname' in doc.keys():
        raise MsPASSError('pfsave:  pf object contains a value for the keyword=pfname.  Not allowed','Fatal')
    doc['pfname'] = name_tag
    col = db[collection]
    result = col.insert_one(doc)
    return result.inserted_id

def pf_dbread(db,id_or_doc,collection='AntelopePf'):
    """
    Creates and returns an AntelopePf from a MongoDB database. 
    
    This function is the inverse of pf_dbsave.  It reads the data stored in 
    MongoDb as a set of documents and subdocuments and reconstructs an 
    AntelopePf object from that data.   The types of the returned data 
    are determined by how MongoDB defined them, which should be the 
    same as the original pf used to create the db record.   

    """
    dbcol = db[collection]
    if isinstance(id_or_doc, ObjectId):
        doc = dbcol.find_one({"_id": id_or_doc})
    elif isinstance(id_or_doc, dict):
        doc = id_or_doc
    else:
        raise MsPASSError("pf_dbread:  Arg 2 has invalid type - must be objectid or doc",
                          "Invalid")
    result = dict_to_pf(doc)
    return result
    
#def pf_verify(db,pfname):
def pf_to_dict(pf):
    """
    Returns a dict representation of an AntelopePf object.  Simple parameters
    map directly to key-value pairs.  Branches (Arr& in the pf format) are
    converted to dict values keyed to the Arr key name.   (Note the conversion
    is recursive so branches of artibrarily deep level can be converted.)
    This function can be called standalone for a conversion but is used
    by the pf_dbsave function as a first step because a dict can be saved
    directly in MongoDB with pymongo.

    :param pf:  AntelopePf object to be converted.
    :return: dict that is an alternate representation of pf
    """
    # This is a trick to extract only the simple name-value pairs.
    # It can work because Arr and Tbl data are stored in separate maps
    # in the pf object
    # TODO:  that is theoretical - test me
    md = Metadata(pf)
    result = Metadata2dict(md)
    # Tbls are stored as an array of strings 
    for tag in pf.tbl_keys():
        tbl = pf.get_tbl(tag)
        result[tag] = tbl
    for tag in pf.arr_keys():
        branchpf = pf.get_branch(tag)
        # WARNING:  this is a recursion
        branch_dict = pf_to_dict(branchpf)
        result[tag] = branch_dict
    return result

def dict_to_pf(doc):
    """
    Inverse of pf_to_dict.
    """
    pf = AntelopePf()
    # For simple types we can use the Metadata put method 
    # and it will resolve the type correctly (VERIFY)
    for key in doc.keys():
        val = doc[key]
        if isinstance(val,str) or isinstance(val,float) or isinstance(val,bool) or isinstance(val,int):
            pf.put(key,val)
        if isinstance(val,list):
            pf.put(key,val)
            
    # needs some testing externally - next handle Tbl data at this level 
    #  Think we can do that with a test for python array type 
    # Then call this function on each entry found as a python dict - also need to verify before continuing
            
