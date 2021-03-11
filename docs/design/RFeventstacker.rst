RFeventstacker
===================
Overview
~~~~~~~~~~~~
This program can definitely exploit mspass effectively.  The processing
flow is fundamentally a loop over gathers defined by a combination of
site_id and a list of source_id's (output of telecluster).   MongoDB
sort of a wf collection find can make that process quite simple.  The
existing C code does that with antelope, but the revision should be much
easier.  Parallelization should be by ensembles with common site_id and
with suite of source_ids defined by clustering.  The first task for the
design of this revision is how to do that.

First, worth reviewing the current algorithm.  This is it is pseudocode:

#. Read a bunch of control parameters from a pf (current program has an
   output file read by pwmig - needs to change to ensembles returned by a
   function.)
#. Long string of antelope db operations culminating in a dbgroup to
   assemble a view grouped to define individual ensembles (need to
   translate this to mongodb)
#. loop over ensemble groups
   #. Load ensemble
   #. resample if necessary (we can drop this for a mspass version and
      throw an exception if a mismatch is found)
   #. Load previously computed hypocentroid of group
   #. Compute stack (uses Stack class - may want to convert that for mspass)
   #. posts a bunch of stuff to metadata of each computed stack
   #. writes result using dbprocess approach in antelope.   (This can be
      simplified a lot for mspass)

The main problems to solve are:  (1) grouping equivalent for a parallel operation
with mspass and (2) database schema needs this algorithm will likely required.
A couple brief subsections follow on each topic.

Grouping to form ensemble groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
I think the key to this is that the telecluster program will need to build a
collection where each document maps to a single ensemble that needs to be
stacked.   That means, I think, that document written by telecluster should
be something like the following::

     sta AAK
     grid_id XXXXXXXX    # object id of a point in an "image" collection
     sources (XXXXX,YYYYYY,....)   # list of source_ids  (may not be needed)
     centroid_lat   22.455  # this and following are hypocentroid of the group
     centroid_lon 122.44
     centroid_depth 44.5
     wfids (list of linking wf_Seismogram ObjectIds)

I am pretty certain if documents are put into a special collection, which
I will tentatively call :code:`cluster`, we can parallelize the processing on
a cursor created from a find query on the cluster collection.

The stacker would then need to a python function that would take one of these
documents and a mspass Database as args.  The algorithm would do the steps
in the pseudocode above following "loop over ensemble groups" (step 3)

Schema issues
~~~~~~~~~~~~~~~~~~
Think we need the following new collections:

#.  :code:`image`  collection.  This is a bit like source, site, and channel
    in that the core data are spatial coordinates.   Provides a geneal way to
    define a spatial set of points.  Requires some thought to make it generic
    but relatable to more efficient storage structures like the gclgrid or
    (more so) regular spatial grids in lat-lon or cartesian space.   This
    collection definitely will need some kind of name tag for queries.
#.  :code:`cluster` source cluster collection used as described in the previous
    section but a generic method to link sources to any image point.  Should
    work, for example, just as well for Shawn's PP precursor bounce points or
    something like PcP or ScP studies.

Actual attributes for each of these collections are better sorted out when
I start implementing this.
