Overview
===============
This document an overview of key initial ideas how to use MsPASS to
parallelize the pwmig package.  The focus will be on the pwmig_core
programs as they are the most time consuming algorithms.

Initial Steps
~~~~~~~~~~~~~~
Before attacking the pwmig_core programs there are at least two major
initial steps for libraries and data prep programs.

1.  :code:`libgclgrid` is a core library for pwmig.  That code is pretty
    standard C++ so building a shared library of the library should not
    be hard.  The more difficult design decision is how to interact with
    python.  A closely related issue is if to utilize mongodb and develop
    constructors that would save and read gclgrid data from mongodb.
2.  :code:`telecluster` is currently an antelope-centric program to assemble
    sources in a radial grid.  It provides key data used as input to
    do source side stacks which proved essential to make improve the performance
    and overall quality of pwmig images.   The new version will almost certainly
    need to be a hybrid C++ and python program using pybing11 and pymongo.
    Basically the new code will need to interact with the mspass source
    collection to produce an output similar to the current program but
    translated to pymongo.

Unambiguous MsPASS Conversions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These programs will definitely profit from being adapted to MsPASS:

1.  :code:`RFeventstacker` existing program that uses telecluster's output to
    do vertical stacks of deconvolved signals.   Likely will become a hybrid
    C++ and python application.
2.  :code:`pwstack`  can be extensively reworked with MSPASS and could
    provide a major speedup.   How to define the parallelism is to be
    determined.  Input from mongodb is an obvious choice.   The existing
    function to compute slant stack over a specified slowness grid
    should be mapable into a map function that has an ensemble in and a
    and ensemble out.  A complication is if the granularity should be like
    the current main program and have each call process one ensemble
    for all image points and all plane waves or make a major change.
    The alternative would be to parallelize at the level of pwmig_process.
    That would require returning the data currenly saved to a scratch
    file passed to pwmig and saving those data to mongodb.  That will, however,
    require creating a python interaction with pwmig.

Less clear
~~~~~~~~~~~~
:code:`pwmig` and :code:`gridstacker` are potenially more problematic.   Both
currently utlize files for input and output.   The files are huge but
are known to be read very quickly.  Not so sure mongodb gridfs would be
advised in this case.  I suspect strongly will want to modify the gclgrid
library to interact cleanly with python.  I think it can be done with a variant
of the older antelope db readers.    
