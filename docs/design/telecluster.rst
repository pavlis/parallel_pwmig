telecluster and EventCatalog
==============================
It turns out that telecluster depends heavily on the memory intensive
EventCatalog class.   The reason is it uses a predicate function object
to do the selection efficiently.   It also uses a pure memory algorithm
but for the expected use that should not be an issue.

Here are design points:

*  The easiest path to solution is to use a file or a fifo to create a
   simplified image of the source collection.   I think a python main
   could be created that would (a) automatically create a tmp file,
   (2) write coordinate data to the scratch file, and (3) run a revision
   of the EventCatalog with a constructor with a file name that defines
   the data file to read.   Actually, when I look at this I think the
   right solution is to build a RadialGridClustering class that either inherits or
   "has a" EventCatalog.  It would be driven by a pf file that could contain
   all the stuff needed including a dir and dfile used to load the source
   data.  That greatly simplifies the C++ bindings because it is reduced to
   a single class.
*  The main of telecluster can and should be converted to python.
   Needs to implement the following:

   1.  A python main or perhaps just a driver function - actually I think
       a function is better as a query option will be essential.
   2.  Main would need to open a database and create a tmp file of source
       coordinates to be the tmp file.
   3.  It next call the C++ constructor to build the new class (tentatively
       called RadialGridClustering).
   4.  Should delete the tmp file when the class is constructed.
   5.  After that the algorithm of the C++ main would be converted to python.
       A loop over the radial grid saving cluster data to mongodb as completed.
       Needed there would be a new set of output function(s) to write the
       results.
