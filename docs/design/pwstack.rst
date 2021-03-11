pwstack
============
This program probably takes the most thought to get right.   I think I need to
actually get the gclgrid library and RFeventstacker converted before take
this one on seriously.  Nonetheless, I think a few things are immediately
clear:

#. The existing core function called :code:`pwstack_process` can likely be
   adapted to the new structure the core algorithm run in parallel.
   i.e. the new basic structure is pwstack would be replaced by a
   parallel structure driving pwstack_process as the core function.
   How that links to ensemble grouping is a key issue that will need some
   work.   Will plan to extend that below.
#. The GCLgrid in 2d used to drive pwstack can probably be well handled with the
   image collection we are going to need for the RFeventcluster algorithm
   anyway.
#. The slowness grid that drives pwstack can be defined by input parameters
   but it will need a common method to distribute it to all the workers.
   I think that just means the pybind11 code will need pickle serialization of
   the slowness grid.   A bit depends on what we aim to parallelize.  I am
   pretty sure the python wrapper should define the loop over slowness calling
   pwstack_process for each input ensemble.
#. The ensemble inputs should use a spatial query using the outer distance
   aperture cutoff.
#. Pretty sure on second pass that parallelization should be at the pseudostation
   and event ensemble level.  i.e. the basic grouping for handling is an
   ensemble formed from a spatial query around the image point combined with
   grouping by cluster source point.  
