pwmig
============
As the end of the chain this is the most uncertain design.   A key initial
design is if we retain the current structure and build an input file to
drive a single event imaging step or redesign pwmig as a reduce function.
A key point for the later that I'm not sure about is at what level the
imaging step is commutative - required by reduce.

A tentative thought about this program is that the initial step would be
to use the existing file structure to build the input to pwmig.   There
are two very good reasons for this:

#.  pwmig is the most complicated program of the group and the fewer
    changes the better.  That would also complicate the interfacing to
    python.
#.  The approach is known to be fast.  io time has always been a small fraction
    of the wall time for this program.

Perhaps a better reason is I think I have a solution for using this model
for pwmig.   Two things would be required besides getting pwmig to compile
within the mspass framework:

#. Write a python function that would assemble a scratch file using the
   existing file format.   The function would assemble the file from pieces
   written by pwstack into MongoDB.   This function would effectively be
   using the database to reorder the data and assemble it into the
   special file format.  This function would be parallel to write and rdd
   of scratch file names (maybe file handles).
#. A small driver function that would be driven by the rdd of file names
   (or file handles) and run a minor modification of the current pwmig
   program.  
