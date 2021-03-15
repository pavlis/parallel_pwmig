makegclgrid
==================
This is an essential command line tool for creating image grid geometries
for pwmig.  The most elaborate solution would be to convert it to a python
function that would have the same functionality.  However, I think that
what I'll do is keep the existing program and have it write a file as it
does now (outside antelope).  I will then write a simple little python
function that will do a form of loading the data into mongodb.  Since the
plan initially is to just use file storage for grids that function will
be a useful feature anyway.   i.e. the sequence to build the image grids
in pwmig will be:

1.  Run makegclgrid and save to files.
2.  Run the (as yet unnamed) function to read the files and put an entry in
    a collection (likely called the gclgrid collection).
