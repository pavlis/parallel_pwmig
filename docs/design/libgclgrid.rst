Changes for libglcgrid
==========================
1.  Will need to wrap this with a multiple level namespace as in mspass.
    Probably something like mspass::gclgrid.
2.  Fix the relic formatting of the include file to work better with doxygen.
3.  Need a constructor that can be called from python to load all types of
    grids and fields from mongodb.   I think a simple way to do that is
    with a Metadata constructor.   The metadata container would contain a
    dir and dfile entries that would be used by the constructor
4.  dir and dfiles need also support foff.   With hat I think a pybind11
    lambda can be used to have a python constructor that includes only an
    ObjectID.  If that doesn't work otehr code could just use the dir dfile
    constructors.
5.  Revision should only throw MsPASSError objects.
6.  May want to add an ErrorLogger to each of the data objects.
