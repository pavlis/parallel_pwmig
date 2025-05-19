#ifndef _GCLGRID2_H_
#define _GCLGRID2_H_
#include <stdlib.h>
#include <string>
#include <vector>
//#include "pwmig/dsap/stock.h"
//#include "pwmig/dsap/coords.h"
#include "mspass/utility/dmatrix.h"
#include "mspass/utility/Metadata.h"
#include "pwmig/gclgrid/GCLgridError.h"
#include "pwmig/gclgrid/swapbytes_pwmig.h"
namespace pwmig::gclgrid
{
//using namespace pwmig::gclgrid;

//==================================================================
/*!
 Geographical Curvilinear Grid library.  This is an object oriented
 library for distorted grid objects.  These are common in 3d visualization
 but the distinguishing factor here is that a GCLgrid is
 geographically referenced.  Points in a grid have both a Cartesian and
 Geographical equivalent.  The data are actually stored in a Cartesian
 format to allow common distorted grid algorithms to work.  Geographical
 elements are added through member functions that convert to a
 geographical reference frame.  The concepts of the library are described in
 a paper by Fan and Pavlis (2006)

 The version here is several generations removedf rom the original Fan and Pavlis
 paper version.  It probably should have been burned and rewritten as it is
 now heavily patched.  Two large changes have happened since the original:
 1.  Addition of file save and file-based constructors.  Created to
     support first full open source version in the mid 2010s.
 2.  Revisions to assimilate the library into mspass including python bindings.
     That was done in 2021 and led to an even larger set of changes.

 For there are some common oddities needed to mesh with mspass and the
 database save method implemented in python.
 1.  The file saves are used as the only method to save the object for now.
     I used the pfhdr format functions.  Something important added was
		 an error trap that forbids the ".pf" part of the file created from
		 being reused.  The original file format supported multiple data in
		 the same file, but anachronisms make that not really functional in the
		 mspass framework.  Beware that each save must have a unique name.
		 This is guaranteed in the python wrapper code by using a mongodb objectid
		 in the file name, but if you use this version in C be care of that
		 oddity.
	2. There are some frozen metadata attributes the mongodb saves have to
	   coordinate with this C++ code.  For the documentation here are
		 the keys that need are set in the pf that are not part of the object
		 attributes:

		 dir
     grid_data_file_extension
     grid_data_file
     grid_data_foff
     field_data_file
     field_data_file_extension
     field_data_foff

	3.  Curently the pfhdr format has two frozen "file extensions".  The
	    attributes are stored as a pf file with the extension ".pf".  The
			internal array data are stored in contigous blocks read with fread
			for the large double array into a buffer.  This helps make the
			readers very fast.  These file have the frozen extension ".dat"

 \author Gary L. Pavlis
*/
//==================================================================

/*! Standardize name for output file format keyword.*/
const std::string default_output_format("pfhdr");
/*! Standardize name for output file data component extention. */
const std::string dfileext("dat");
/*! Standardize name for the parameter file extention added to a dfile base*/
const std:: string pfext("pf");

/*!
 This data structure is used to encapsulate data to describe a point on
 or inside the earth using a geographical (spherical geometry) reference
 system.  Note that geographical angles are always assumed to be radians
 in any internal use of this package.

 These were previously C struct.  for mspass python conversion
 changed to class and added default constructor with initalization.
*/

class  Geographic_point
{
public:
/*!
  Latitude of the point (radians).
*/
	double lat;
/*!
 Longitude of point (radians).
*/
	double lon;
/*!
 Radius of point from Earth's center (kilometers).
*/
	double r;
	Geographic_point()
	{
		lat=0.0; lon=0.0; r=0.0;
	};
	/*! Fully parameterized constructor.

	This constructor was added 2022 to simply pickle interface
	in pybind11 */
	Geographic_point(const double lat_in, const double lon_in, const double r_in)
	{
		lat = lat_in;
		lon = lon_in;
		r = r_in;
	};
	Geographic_point(const Geographic_point& parent)
	{
		lat=parent.lat;
		lon=parent.lon;
		r=parent.r;
	};
};
/*!
 GCLgrid objects hold points internally in a Cartesian reference fram.
 This data structure encapsulates such a coordinate.  It perhaps should
 be a class with a member to return a 3 vector alternative to the
 verbose naming.
*/
class Cartesian_point
{
public:
/*!
 The coordinate in the internal x1 direction.
*/
	double x1;
/*!
 The coordinate in the internal x2 direction.
*/
	double x2;
/*!
 The coordinate in the internal x3 direction.
*/
	double x3;
	Cartesian_point()
	{
		x1=0.0; x2=-0.0; x3=0.0;
	};
	/*! Fully parameterized constructor.

	This constructor was added 2022 to simply pickle interface
	in pybind11 */
	Cartesian_point(const double x1_in, const double x2_in, const double x3_in)
	{
		x1 = x1_in;
		x2 = x2_in;
		x3 = x3_in;
	};
	Cartesian_point(const Cartesian_point& parent)
	{
		x1=parent.x1;
		x2=parent.x2;
		x3=parent.x3;
	};
	std::vector<double> coordinates()
	{
		std::vector<double> v;
		v.push_back(x1);
		v.push_back(x2);
		v.push_back(x3);
		return v;
	}
};

/*!
 This is a base class that contains common attributes and virtual
 declarations for higher level objects that are derived from this base.
 Note that both 2d and 3d grid and field objects are derived from here.
 \author Gary L. Pavlis
*/
class BasicGCLgrid
{
public:
/*!
 Name assigned to this object
*/
	std::string name;
/*!
 Latitude (radians) of origin of the grid's Cartesian coordinate system.
*/
	double lat0;
/*!
 Longitude (radians) of origin of the grid's Cartesian coordinate system.
*/
	double lon0;
/*!
 Radial distance of origin of the grid's Cartesian coordinate system from Earth's center.
*/
	double r0;
/*!
 Nominal azimuth (radians) of positive y axis of grid line at origin location.
*/
	double azimuth_y;
/*!
 Nominal spacing (km) of grid lines along the 1 position gridlines.
*/
	double dx1_nom;
/*!
 Nominal spacing (km) of grid lines along the 2 position gridlines.

*/
	double dx2_nom;
/*!
 Number of grid points in generalized coordinate index 1 direction (first array index).
*/
	int n1;
/*!
 Number of grid points in generalized coordinate index 2 direction (second array index).
*/
	int n2;
/*!
 Offset in first index to origin grid point.  i.e. origin is at first coordinate index i0.
*/
	int i0;
/*!
 Offset in second index to origin grid point.  i.e. origin is at second coordinate index j0. ([i0][j0])
*/
	int j0;
/*!
 Cartesian x1 coordinate of lower limit of bounding box for grid.
 This edge of the  bounding box is defined as the smallest x1 value minus dx1_nom.
*/
	double x1low;
/*!
 Cartesian x1 coordinate of upper limit of bounding box for grid.
 This edge of the  bounding box is defined as the largest x1 value plus dx1_nom.
*/
	double x1high;
/*!
 Cartesian x2 coordinate of lower limit of bounding box for grid.
 This edge of the  bounding box is defined as the smallest x2 value minus dx2_nom.
*/
	double x2low;
/*!
 Cartesian x2 coordinate of upper limit of bounding box for grid.
 This edge of the  bounding box is defined as the largest x2 value plus dx2_nom.
*/
	double x2high;
/*!
 Cartesian x3 coordinate of lower limit of bounding box for grid.
 This edge of the  bounding box is defined as the smallest x3 value minus dx3_nom.
*/
	double x3low;
/*!
 Cartesian x3 coordinate of upper limit of bounding box for grid.
 This edge of the  bounding box is defined as the largest x3 value plus dx3_nom.
*/
	double x3high;
/*!
 Default Constructor.

 Implemented to initialize all base attributes explicitly.
*/
	BasicGCLgrid();
	/*! Need explicit virtual destructor for this base class.
	An oddity of inheritance discussed in many C++ textbooks.
	*/
  virtual ~BasicGCLgrid(){}

	;
/*! Copy constructor.*/
	BasicGCLgrid(const BasicGCLgrid& old);
/*!
 This is a rotation matrix that defines the transformation from standard spherical
 coordinates (the geographical reference frame) to the local Cartesian coordinate
 system used in a GCLgrid.  Was private in an earlier incarnation, but this is
 is messy when we depend on inheritance so it is public.  Users should not manipulate
 this directly, however, but if it is desired they should use the fetch_transformation_matrix
 member function.
*/
	double gtoc_rmatrix[3][3];
/*!
 This is a close companion to the gtoc_rmatrix. The full transformation used to
 define the Cartesian system in a GCLgrid is a translation from the earth's center
 to the grid coordinate system origin.  This is followed by a rotation by gtoc_rmatrix.
 The translation_vector attribute defines what the name implies.
*/
	double translation_vector[3];
/*!
 Set or reset the transformation properties for the grid.  The transformation properties
 are uniquely defined by the coordinate system origin and azimuth_y so if these are
 all that are known this low level member function can be called to set the transformation
 properties.  It is of minimal use to most users and should be used cautiously and only
 if you thoroughly understand the way this all works.

*/
	void set_transformation_matrix();
/*!
 Returns the transformation matrix for this grid as a 3x3 dmatrix object.
*/
	mspass::utility::dmatrix fetch_transformation_matrix() const;
/*!
 Returns a newly allocated 3 vector of double containing a copy of the translation
 vector defining the GCLgrid transformation property.  The user must be sure to
 call delete [] after using this vector to avoid a memory leak.
*/
	double *fetch_translation_vector() const;
/*!
 Convert Cartesian coordinates to geographical coordinates.
 \return Geographic_point data structure

 \param x1p - Cartesian x1 coordinate of point to convert
 \param x2p - Cartesian x2 coordinate of point to convert
 \param x3p - Cartesian x3 coordinate of point to convert
*/
	Geographic_point ctog(const double x1p,const double x2p, const double x3p) const;
/*!
 Convert from Cartesian coordinates to geographical coordinates.
 \return Geographic_point data structure
 \param p point to convert stored in a Cartesian_point data structure
*/
	Geographic_point ctog(const Cartesian_point p) const;
/*!
 Convert from geographical to Cartesian coordinates in the GCLgrid
 coordinate system.

 \return Cartesian_point data structure
 \param plat Latitude (radians) of point to convert.
 \param plon Longitude (radians) of point to convert.
 \param pr Earth radius (km) of point to convert.
*/
	Cartesian_point gtoc(const double plat, const double plon, const double pr) const;
/*!
 Convert from geographical to Cartesian coordinates in the GCLgrid
 coordinate system.

 \return Cartesian_point data structure
 \param p point to convert stored in a Geographic_point data structure.
*/
	Cartesian_point gtoc(const Geographic_point p) const;
	double depth(const Cartesian_point p) const;
	double depth(const Geographic_point p) const;
/*!
 This member function sets the x1min, x1max, x2min, x2max, x3min, and x3max bounding
 box attribute.  These define the "extents" of the name.  Should be needed only if one
 builds a custom grid from user defined coordinates.  Said a different way if you build a
 GCLgrid from lower level routines be sure to call this function to set the bounding box
 correctly.
*/
	virtual void compute_extents()=0;
/*!
Reset lookup index to origin.
The lookup functions used in the GCLgrid library keeps an index of
the previous lookup results under an assumption the next point
requested will be nearby.  This can cause convergence problems in
some grids, however, if that position is a poor place to start
(e.g. sudden jump to a completely different location).  This function
should be called if a lookup fails to try to recover.  Internal methods
like the += operator does this automatically.

This method may be depricated but is retained for now.
*/
	virtual void reset_index()=0;
/*!
Query for the current lookup index position.  The lookup functions used
in the GCLgrid library keep and index of the previous lookup result under
an assumption that the next point requested will be nearby.  This method
is used to ask what the current index position.

This method may be depricated. the new parallel_lookup does not use the
internal index preserving the previous state.

\param ind vector of ints of length sufficient to hold the index
(2 for 2d and 3 for 3d grids).  This is very dangerous and there is no
bound checking.  Usual construct would be to declare a local int array
of the right length and pass the array name to this method.
*/
	virtual void get_index(int *ind)=0;

/*! \brief Return all attributes of a GCL object in a Metadata container.

This virtual method was added for the mspass conversion.  It is a core
need to save the attributes of any of the GCl objects with MongoDB.  In the
implementation we always save the grid and field data to files, but the
attributes are saved in MongoDB.  Because Metadata maps cleanly to MongoDB's
name-value pair paradigm that works as it does in MsPASS for seismic data.
*/
  virtual mspass::utility::Metadata get_attributes() const = 0;
/*!
 Comparison of two grids for equality.  Equality in this context is NOT the obvious.
 Equality means the transformation properties of the two grids being compared are the same.
 That is, the operator tests only the transformation matrix and translation vector for
 equality.  This is a necessary condition to allow to grids to be mapped onto each other
 for higher level operations like +=.
*/
	bool operator==(const BasicGCLgrid&);
/*!
 Comparison of two grids for inequality.  Equality in this context is NOT the obvious.
 Equality means the transformation properties of the two grids being compared are the same.
 That is, the operator tests only the transformation matrix and translation vector for
 equality.  This is a necessary condition to allow to grids to be mapped onto each other
 for higher level operations like +=.  Returns true if the transformation properties of
 two matrices do not match.
*/
	bool operator!=(const BasicGCLgrid&);
};

/*!
 This is the working two-dimensional version of a GCLgrid.  A GCLgrid defines a two-dimensional
 surface in Earth coordinates.  This could be any surface, but the lookup method used in
 this version will not work if the surface is multivalued or even strongly warped.
 See the algorithm description in the paper by Fan and Pavlis (in review) for details.

 Positions of grid points are stored in three, two-dimensional arrays. Each of these
 arrays define cartesian coordinates for a grid point at a particular index position.
 The topology is that the location in the Cartesian coordinate system of the grid
 point defined by index positions i and j are (x1[i][j], x2[i][j], x3[i][j]).
 The coordinates are thus matrices indexing each position in grid.

 Note the object inherits most data attributes from the BasicGCLgrid object.
*/


class GCLgrid : public BasicGCLgrid
{
public:
/*!
 n1 x n2 matrix of the x1 coordinate values for
 grid elements.
*/
	double **x1;
/*!
 n1 x n2 matrix of the x2 coordinate values for
 grid elements.
*/
	double **x2;
/*!
 n1 x n2 matrix of the x3 coordinate values for
 grid elements.
*/
	double **x3;
/*!
 Default constructor.  Note it sets the x1,x2, and x3 pointers to NULL
 which is used as a test to avoid duplicate free calls on these potentially
 large arrays.
*/
	GCLgrid(){
		n1=0;n2=0;x1=NULL;x2=NULL;x3=NULL;
	};
/*!
 Simple constructor.  Allocates space for x1, x2, and x3 arrays and initializes
 other data attributes to zero.

  \param n1 number of grid points on generalized coordinate axis 1.
  \param n2 number of grid points on generalized coordinate axis 2.
*/
	GCLgrid(const int n1size, const int n2size);
/*!
  Constructor for what we call a "regular" GCLgrid in the Fan and Pavlis (in review) paper.
  The object this constructs is a constant geoid elevation surface (Follows the reference
  ellipsoid at elevation set by the origin radius, r0.) with approximately regular spacing
  between grid points (within the limits of what is possible on a spherical surface).
  Note that the makegclgrid program is little more than a wrapper around this and the 3d
  version of this constructor.

  \param n1 number of grid points on generalized coordinate axis 1.
  \param n2 number of grid points on generalized coordinate axis 2.
  \param n name to assign this grid.
  \param la0 latitude to use for origin.
  \param lo0 longitude to use for origin.
  \param radius0 Earth radius to use for origin.  If 0 or negative will use
		r0_ellipse at la0.
  \param az azimuth of great circle path through the origin that defines the generalized coordinate 2 axis direction.
  \param dx1n nominal spacing of grid points on the generalized coordinate 1 axis.
  \param dx2n nominal spacing of grid points on the generalized coordinate 2 axis.
  \param iorigin 1 axis grid index of the origin in generalized coordinate grid frame.
  \param jorigin 2 axis grid index of the origin in generalized coordinate grid frame.

*/
	GCLgrid(const int n1size, const int n2size, const std::string n,
		const double la0, const double lo0,
		const double radius0, const double az, double dx1n, const double dx2n,
		const int iorigin, const int jorigin);
        /*! \brief Constuct from a file.

          A standard simple way to build any object is from a data file.
          This provides that interface.  The interface is relatively
          generic, however, as the format string provides a mechanism
          to specify variations in possible file format.  In the
          current implementation the file name should not contain
          a .ext string (e.g. mygrid.dat) as this may be used to
          split a header section from a data section.

          \param fname is the file name to read (may be just a root
          name depending on format option)
          \param format is a string describing allowed data format.
          Current default is split structure with one file used to
          hold an antelope pf to store scalar attributes and a second
          file used to store the actual data. Most applications will
          likely want to force default and use only one argument to
          this constructor.

          */
  GCLgrid(const std::string fname, const std::string format=default_output_format);
	/*! Metadata driven constructor.

	This constructor was added to allow this package to work with MongoDB.
	It acts very similar to the file-based constructors but instead of
	reading a pf representation of the attributes it assumes the same and a
	bit more are found in a Metadata container.   The idea of this constructor
	is that it will be used with pymongo.  The pymongo function reads the
	Metadata attributes from a MongoDB document as in MsPASS.   The grid
	data, however, is always assumed to be present in an external file
	accessible from the compute node calling this constructor.  I say that
	because that model may fail in cloud systems.

	One difference with this constructor compared to the file based constructor
	is that the following attributes are required that are NOT written with
	a file-based save to the pf file:

	dir - directory name string
	grid_data_file - file base name (with extension)
	grid_data_file_extension - file name "extension for data portion of file extension
	  (normally should use keyword used internally as ".dat" but constructor
	  allows others)
  grid_data_foff - integer file offset to locate start of block of data for
	 grid coordinates of this object.

	\param md is the Metadata container normally assumed produced by cracking a
	  MongoDB document.

	\throw GCLgridError for a long list of possible failuers.   Any
	  exceptions mean the object was not constucted.
	 */
	GCLgrid(const mspass::utility::Metadata& md);
/*! Standard copy constructor.*/
	GCLgrid(const GCLgrid&);  //standard copy constructor
/*!
 Destructor.  Nontrivial destructor has to destroy the coordinate arrays correctly
 and handle case when they are never defined.  Handles this by checking for
 NULL pointers for these arrays.  If the pointers are NULL the free routines
 are not called.  This is important to know if you try to create a GCLgrid
 object by the default constructor.
*/
	~GCLgrid();
/*!
 Standard assignment operator.
*/
	GCLgrid& operator=(const GCLgrid& );
    /*! \brief Save to a file.

      Many applications will find it more convenient to write
      this object to standard data files.  This abstracts the
      write process putting result in one or more formats.
      The concept here is the most applications would use
      the defaults for the the format specifications.
      Alternate formats are defined by using a different keyword
      to this method and implementing that alternative format
      in the code for this method.  Similarly alternate formats
      should implement a parallel code change for the file-based
      constructor which has a similar format argument.

      \param fname is the name of the output file to save data
      \param dir is a directory name to store data.  If zero
        length assume current directory.
      \param format is a keyword used to describe the format.
        (current default is a two file output separating the
        header and data.  Hence users should avoid .ext
        file names like myfile.dat.)


			\return copy of Metadata container with attributes defining the
			 object sans grid data.

      \exception GCLgridError is throw if save fails.
    */
    mspass::utility::Metadata save(const std::string fname, const std::string dir,
            const std::string format=default_output_format);
	/*!
	 Find the index position of a point in a GCLgrid.
	 This is a low level function to find the location of a point
	 specified as the Cartesian, ordered pair (x1p,x2p) in a grid.
	 It does not return the actual index positions, but only sets the
	 internal index.  The routine is very procedural returning an
	 integer code (see below) indicating success or failure of the lookup.
	 This was intentional as this routine is called millions of times in
	 some contexts so efficiency is critical.  The alternative would be
	 to throw an exception when a lookup failed, but since this is viewed
	 as a common problem that could happen millions of times this was
	 a potential efficiency problems (the books all say throwing exceptions
	 is an expensive operation).
	 \return 2 when point is in gray area within on nominal grid spacing of the edge
	 \return 1 when the point is outside the bounding box.
	 \return 0 on success.\
	 \return -1 if the lookup function did not converge.

	 \param x1p - Cartesian x1 coordinate of point to find within the grid
	 \param x2p - Cartesian x2 coordinate of point to find within the grid
	*/
	int lookup(const double x1p, const double x2p);
	void reset_index() {ix1=i0; ix2=j0;};
	void get_index(int *ind) {ind[0]=ix1; ind[1]=ix2;};
	/*!
	 Returns the geographical coordinates of a point in the grid specified by grid index positions.

	 If you need the actual coordinates of the points that define the grid
	 converted to geographic coordinates use this function.  Use the
	 ctog function to convert an arbitrary ordered triplet.

	 \return grid point requested in an Geographic_point data structure.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	*/
	Geographic_point geo_coordinates(const int i1,const int i2) const;
	/*!
	 Get the latitude (in radians) of a grid point.

	 \return latitude in radians of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	*/
	double lat(const int i1, const int i2) const;
	/*!
	 Get the longitude (in radians) of a grid point.

	 \return longitude in radians of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	*/
	double lon(const int i1, const int i2) const;
	/*!
	 Get the radius from the center of the Earth of a grid point.

	 \return Earth radius in kilometers of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	*/
	double r(const int i1, const int i2) const;
	/*!
	 Get the depth below the standard reference ellipsoid of a grid point.

	 \return depth in kilometers of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	*/
	double depth(const int i1, const int i2) const;

	/*! Returns all scalar attributes of the object in Metadata container. */
	virtual mspass::utility::Metadata get_attributes() const;

	/*! Sets extents attributes based on min and max values */
	void compute_extents();
	/*! Set coordinates at a grid point.

	Sets the coordinates of a point at index i,j.  The only checks are
	that i and j received are in bounds.  It will throw a GCLgridError if the
	points given are invalid.  This is largely a convenience function for python
	as in C++ the arrays are public.

	\param p Geogrpahic_point definition of the point to be set.
	\param i first grid index value to set (i.e i of  [i][j])
	\param j second grid index value to set
	*/
	void set_coordinates(const Cartesian_point& p, const int i, const int j);
	/*! Set coordinates at a grid point.

	Sets the coordinates of a point at index i,j.  The only checks are
	that i and j received are in bounds.  It will throw a GCLgridError if the
	points given are invalid.  This is largely a convenience function for python
	as in C++ the arrays are public.  This is an overloaded version for
	input a Geographic coordinates versus Cartesian

	\param p Geogrpahic_point definition of the point to be set.
	\param i first grid index value to set (i.e i of  [i][j])
	\param j second grid index value to set
	*/
	void set_coordinates(const Geographic_point& p, const int i, const int j);
	/*! Getter for grid point */
	Cartesian_point get_coordinates(const int i, const int j) const;
	friend class GCLscalarfield;
	friend class GCLvectorfield;
private:
	int ix1, ix2;
};
/*!
 Three-dimensional version of a GCLgrid object.

 Adds additional attributes required to deal with an added dimension and
 uses a higher dimensional array to contain points.
*/
class GCLgrid3d : public BasicGCLgrid
{
public:
/*!
 Nominal spacing (km) of grid lines along the 3 position gridlines.

*/
	double dx3_nom;
/*!
 Number of grid points in generalized coordinate index 3 direction (second array index).
*/
	int n3;
/*!
 Offset in third index to origin grid point.  i.e. origin is at third coordinate index k0. ([i0][j0][k0])
*/
	int k0;
/*!
 n1 x n2 x n3 Three-dimensional array that holds the Cartesian
 x1 component of the positions of each grid point.
*/
	double ***x1;
/*!
 n1 x n2 x n3 Three-dimensional array that holds the Cartesian
 x2 component of the positions of each grid point.
*/
	double ***x2;
/*!
 n1 x n2 x n3 Three-dimensional array that holds the Cartesian
 x3 component of the positions of each grid point.
*/
	double ***x3;

	/*!
	 Default constructor.
	 Note sets pointers to NULL to make destructor work correctly
	 when a grid is created by this constructor.
	*/
	GCLgrid3d(){
		n1=0;n2=0;n3=0;
		x1=NULL;x2=NULL;x3=NULL;
                fast_lookup=true;
	};
/*!
 Simple constructor.  Allocates space for x1, x2, and x3 arrays and initializes
 other data attributes to zero.

  \param n1size number of grid points on generalized coordinate axis 1.
  \param n2size number of grid points on generalized coordinate axis 2.
  \param n3size number of grid points on generalized coordinate axis 3.
  \parm fl when true use fast lookup method.  False is slower and more sure.
*/
	GCLgrid3d(const int n1size, const int n2size, const int n3size, const bool fl=true);
/*!
  Constructor for what we call a "regular" GCLgrid in 3D in the Fan and Pavlis (in review) paper.
  The object this constructs is spherical shell, boxlike object built up of elemental spherical
  shell cube-like grid components.  The top surface of the box is defined by the reference
  ellipsoid.  The bottom is a constant depth below this.  Note that the grid this creates
  is oriented with the third generalized coordinate index running from the bottom (deeper
  inside the Earth) upward as the grid index increases.  Note no index is allowed for
  index 3 for the origin.  This constuctor always puts the origin at the bottom of the grid.
  This isn't necessary, but something that is just frozen into this implementation.
  For the same reason r0 is ignored and just set internally by this constructor.
  Note the fast lookup method is always on for this constructor because the
  grid is always regular making this extra overhead of the slow method
  totally unnecessary.

  Note that the makegclgrid program is little more than a wrapper around this and the 2d
  version of this constructor.

  \param n1 number of grid points on generalized coordinate axis 1.
  \param n2 number of grid points on generalized coordinate axis 2.
  \param n3 number of grid points on generalized coordinate axis 3.
  \param n name to assign this grid.
  \param la0 latitude to use for origin.
  \param lo0 longitude to use for origin.
  \param radius0 Earth radius to use for origin.  If 0 or negative will use
		r0_ellipse at la0.
  \param az azimuth of great circle path through the origin that defines the generalized coordinate 2 axis direction.
  \param dx1n nominal spacing of grid points on the generalized coordinate 1 axis.
  \param dx2n nominal spacing of grid points on the generalized coordinate 2 axis.
  \param dx3n nominal spacing of grid points on the generalized coordinate 3 axis.
  \param iorigin 1 axis grid index of the origin in generalized coordinate grid frame.
  \param jorigin 2 axis grid index of the origin in generalized coordinate grid frame.


*/
	GCLgrid3d(const int n1size, const int n2size, const int n3size, const std::string n,
		const double la0, const double lo0, const double radius0, const double az,
		const double dx1n, const double dx2n, const double dx3n,
		const int iorigin, const int jorigin);

        /*! \brief Constuct from a file.

          A standard simple way to build any object is from a data file.
          This provides that interface.  The interface is relatively
          generic, however, as the format string provides a mechanism
          to specify variations in possible file format.  In the
          current implementation the file name should not contain
          a .ext string (e.g. mygrid.dat) as this may be used to
          split a header section from a data section.

          \param fname is the file name to read (may be just a root
          name depending on format option)
          \param format is a string describing allowed data format.
          Current default is split structure with one file used to
          hold an antelope pf to store scalar attributes and a second
          file used to store the actual data. Most applications will
          likely want to force default and use only one argument to
          this constructor.
          */
        GCLgrid3d(const std::string fname, const std::string format=default_output_format,
                const bool fl=true);

	/*! Metadata driven constructor.

	This constructor was added to allow this package to work with MongoDB.
	It acts very similar to the file-based constructors but instead of
	reading a pf representation of the attributes it assumes the same and a
	bit more are found in a Metadata container.   The idea of this constructor
	is that it will be used with pymongo.  The pymongo function reads the
	Metadata attributes from a MongoDB document as in MsPASS.   The grid
	data, however, is always assumed to be present in an external file
	accessible from the compute node calling this constructor.  I say that
	because that model may fail in cloud systems.

	One difference with this constructor compared to the file based constructor
	is that the following attributes are required that are NOT written with
	a file-based save to the pf file:

	dir - directory name string
	grid_data_file - file base name (with extension)
	grid_data_file_extension - file name "extension for data portion of file extension
		(normally should use keyword used internally as ".dat" but constructor
		allows others)
	grid_data_foff - integer file offset to locate start of block of data for
	 grid coordinates of this object.

	\param md is the Metadata container normally assumed produced by cracking a
		MongoDB document.

	\throw GCLgridError for a long list of possible failuers.   Any
		exceptions mean the object was not constucted.
	 */
	GCLgrid3d(const mspass::utility::Metadata& md);
	/** Standard copy constructor. */
	GCLgrid3d(const GCLgrid3d&);
	/** Standard assignment operator. */
	GCLgrid3d& operator=(const GCLgrid3d& );
        /*! \brief Save to a file.

          Many applications will find it more convenient to write
          this object to standard data files.  This abstracts the
          write process putting result in one or more formats.
          The concept here is the most applications would use
          the defaults for the the format specifications.
          Alternate formats are defined by using a different keyword
          to this method and implementing that alternative format
          in the code for this method.  Similarly alternate formats
          should implement a parallel code change for the file-based
          constructor which has a similar format argument.

          \param fname is the name of the output file to save data
          \param dir is a directory name to store data.  If zero
            length assume current directory.
          \param format is a keyword used to describe the format.
            (current default is a two file output separating the
            header and data.  Hence users should avoid .ext
            file names like myfile.dat.)
					\return copy of Metadata container with attributes defining the
						 object sans grid data.
          \exception GCLgridError is throw if save fails.
        */
        mspass::utility::Metadata save(const std::string fname, const std::string dir,
                const std::string format=default_output_format);
	/*!
	 Find the index position of a point in a GCLgrid3d object.
	 This is a low level function to find the location of a point
	 specified as the Cartesian, ordered triplet (x1p,x2p,x3p) in a grid.
	 It does not return the actual index positions, but only sets the
	 internal index.  The routine is very procedural returning an
	 integer code (see below) indicating success or failure of the lookup.
	 This was intentional as this routine is called millions of times in
	 some contexts so efficiency is critical.  The alternative would be
	 to throw an exception when a lookup failed, but since this is viewed
	 as a common problem that could happen millions of times this was
	 a potential efficiency problems (the books all say throwing exceptions
	 is an expensive operation).

	 This function was the original api for this library.  To adapt the code
	 to mspass it was quickly realized the algorithm was not thread safe because
	 of the internally maintained index it depended up.  This method
	 was retained for backward compatibility but users are warned it may be
	 depricated in the future and should not be used.

	 \return 2 when point is in gray area within on nominal grid spacing of the edge
	 \return 1 when the point is outside the bounding box.
	 \return 0 on success.\
	 \return -1 if the lookup function did not converge.

	 \param x1p - Cartesian x1 coordinate of point to find within the grid
	 \param x2p - Cartesian x2 coordinate of point to find within the grid
	 \param x3p - Cartesian x3 coordinate of point to find within the grid
	*/
	int lookup(const double x1p, const double x2p, const double x3p);
	/*! \brief Thread safe lookup method.

	 Finds the index position of a point in a GCLgrid3d object.
	 This is a low level function to find the location of a point
	 specified as the Cartesian, ordered triplet (x1p,x2p,x3p) in a grid.
	 The routine is very procedural returning an
	 integer code (see below) indicating success or failure of the lookup.
	 This was intentional as this routine is called millions of times in
	 some contexts so efficiency is critical.  The alternative would be
	 to throw an exception when a lookup failed, but since this is viewed
	 as a common problem that could happen millions of times this was
	 a potential efficiency problems (the books all say throwing exceptions
	 is an expensive operation).

	 This method uses the same algorithm as the origina lookup method
	 (see above) but the user is responsible for maintaining the lookup
	 index through the arguments ix1_0, ix2_0, and ix3_0.   This makes
	 the method stateless which the old code was note.   The decision about
	 what value to use for the starting index is dependent on data and
	 the use.  Random lookups should normally alway just use the lookup origin
	 and ignore what is returned.  Lookups working along a coordinate axis,
	 which is the norm for most grid processing algorithms, should maintain
	 the index within the loop and use the previous value as a starting point
	 for the next point.  A key point of this method is that type of lookup
	 in distorted grids can be an order of magnitude faster than a search
	 from the origin.

	 \return 1 when the point is outside the bounding box.
	 \return 0 on success.\
	 \return -1 if the lookup function did not converge.

	 \param x1p - Cartesian x1 coordinate of point to find within the grid
	 \param x2p - Cartesian x2 coordinate of point to find within the grid
	 \param x3p - Cartesian x3 coordinate of point to find within the grid
	 \param ix1_0 - first grid index value from which the iteration should start.
	 \param ix2_0 - second grid index value from which the iteration should start.
	 \param ix3_0 - thrid grid index value from which the iteration should start.

	*/
	int parallel_lookup(const double x1p, const double x2p, const double x3p,
	    int& ix1_0, int& ix2_0, int& ix3_0) const;
	/*! Overloaded convenience method to use an std::vector to hold the index. */
	int parallel_lookup(const double x1p, const double x2p, const double x3p,
	    std::vector<int>& index) const;
	/*! \brief Test if a point is inside a grid.

	It is not generally trivial to determine if a point is inside a
	distorted box (the concept a GCLgrid3d implements).  This method can
	be used a convenience to test for whether or not a point is inside the
	grid.  It is a useful alternative to calling the low level parallel_lookup
	method.   This function is only appropriate for random tests or where
	readability is more important than speed because it always uses a search
	starting from the lookup origin. */
	bool point_is_inside_grid(const Geographic_point gp);
	/*! Companion to the old lookup method.  It should not be used as it
	will may be depricated in the fugure. */
	void reset_index() {ix1=i0; ix2=j0; ix3=k0;};
	/*! Companion to the old lookup method.  It should not be used as it
	will may be depricated in the fugure. */
	void get_index(int *ind) {ind[0]=ix1; ind[1]=ix2; ind[2]=ix3;};
	/*!
	 Returns the geographical coordinates of a point in the grid specified by grid index positions.

	 If you need the actual coordinates of the points that define the grid
	 converted to geographic coordinates use this function.  Use the
	 ctog function to convert an arbitrary ordered triplet.

	 \return grid point requested in an Geographic_point data structure.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	 \param i3 Value of grid index 3 for point desired.
	*/
	Geographic_point geo_coordinates(const int i1,const int i2,const int i3) const;
	/*!
	 Get the latitude (in radians) of a grid point.

	 \return latitude in radians of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	 \param i3 Value of grid index 3 for point desired.
	*/
	double lat(const int i1,const int i2,const int i3) const;
	/*!
	 Get the longitude (in radians) of a grid point.

	 \return latitude in radians of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	 \param i3 Value of grid index 3 for point desired.
	*/
	double lon(const int i1,const int i2,const int i3) const;
	/*!
	 Get the radius from the center of the Earth of a grid point.

	 \return Earth radius in kilometers of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	 \param i3 Value of grid index 3 for point desired.
	*/
	double r(const int i1,const int i2,const int i3) const;
	/*!
	 Get the depth below the standard reference ellipsoid of a grid point.

	 \return depth in kilometers of the requested grid point.
	 \param i1 Value of grid index 1 for point desired.
	 \param i2 Value of grid index 2 for point desired.
	 \param i3 Value of grid index 3 for point desired.
	*/
	double depth(const int i1,const int i2,const int i3) const;
	/*! Returns all scalar attributes of the object in Metadata container. */
	virtual mspass::utility::Metadata get_attributes() const;

	/*! Sets extents attributes based on min and max values */
	void compute_extents();
        /*! \brief Enable high accuracy lookup.

          Default lookup method can sometimes find the wrong cell when
          a grid is extremely distorted.  The direction set method uses
          a Jacobian computed with vectors computed in the direction of
          increasing grid index values.   If the point of a distorted
          cell box not touched by the set of edge vectors defined by this
          operation differs strongly from a projection from the opposite
          edges, the lookup will find the wrong cell.  The high accuracy
          method checks for this condition and attempts to recover using
          a search in adjacent cells.  This is, however, a very expensive
          operation and should normally be avoided unless your application
          is prone to large, local distortions.
          */
        void enable_high_accuracy()
        {
            fast_lookup=false;
        };
        /*! Inverse of enable_high_accuracy method - turn on fast method.*/
        void enable_fast_lookup()
        {
            fast_lookup=true;
        };

/*!  Destructor.
 Nontrivial destructor has to destroy the coordinate arrays correctly
 and handle case when they are never defined.  Handles this by checking for
 NULL pointers for these arrays.  If the pointers are NULL the free routines
 are not called.  This is important to know if you try to create a GCLgrid
 object by the default constructor.
*/
	~GCLgrid3d();
	/*! Set coordinates at a grid point.

	Sets the coordinates of a point at index i,j.  The only checks are
	that i and j received are in bounds.  It will throw a GCLgridError if the
	points given are invalid.  This is largely a convenience function for python
	as in C++ the arrays are public.

	\param p Geogrpahic_point definition of the point to be set.
	\param i first grid index value to set (i.e i of  [i][j][k])
	\param j second grid index value to set
  \param k third grid index value to set
	*/
	void set_coordinates(const Cartesian_point& p,
		const int i, const int j,const int k);
	/*! Set coordinates at a grid point.

	Sets the coordinates of a point at index i,j.  The only checks are
	that i and j received are in bounds.  It will throw a GCLgridError if the
	points given are invalid.  This is largely a convenience function for python
	as in C++ the arrays are public.  This is an overloaded version for
	input a Geographic coordinates versus Cartesian

	\param p Geogrpahic_point definition of the point to be set.
	\param i first grid index value to set (i.e i of  [i][j][k])
	\param j second grid index value to set
  \param k third grid index value to set
	*/
	void set_coordinates(const Geographic_point& p,
		const int i, const int j,const int k);
	/*! Getter for grid point */
	Cartesian_point get_coordinates(const int i, const int j, const int k) const;
	/*! Get the lookup origin.

	The parallel and single threaded lookup functions both utilize an iterative
	algorithm to find the cell a specified point is within.  The iteration
	can sometimes fail and the algorithm the always resets to the "lookup origin"
	and tries again before giving up.  The parallel_lookup function uses this
	method to do that reset.   The original version that was single threaded
	had the index hard coded.
	*/
	std::vector<int> get_lookup_origin() const;
	/*! \brief Set the lookup origin to the default value.

	The default lookup origin is the grid center computed from n1/2, n2/2, and
	n3/2.
	*/
	void set_lookup_origin();
	/*! \brief Set the lookup origin to a user specified point.

	The grid origin is used for recovery if the initial iteration fails.
	Use this method to set the origin explicitly.  The lookup origin is
	set to the grid index [i][j][k].  An exception will be thrown if any of
	the indices given are outside the array dimensions.
	*/
	void set_lookup_origin(const int i, const int j, const int k);

private:
	int ix1, ix2, ix3;
        bool fast_lookup;
};
/** Two-dimensional scalar field defined on a GCLgrid framework. */
class GCLscalarfield :  public GCLgrid
{
public:
	/*!
	 Scalar field variable stored in a two-dimensional C array.
	 Index of the field is parallel with coordinate arrays x1, x2, and x3
	 that define grid point positions in space.
	*/
	double **val;

	/** Default constructor */
	GCLscalarfield();
	/*!
	 Simple constructor.  Allocates space for all arrays but loads nothing.
	 Do not assume anything is initialized.

	  \param n1size number of grid points on generalized coordinate axis 1.
	  \param n2size number of grid points on generalized coordinate axis 2.
	*/
	GCLscalarfield(const int n1size, const int n2size);
	/** Standard copy constructor */
	GCLscalarfield(const GCLscalarfield& parent);
	/*!
	 Partial copy constructor cloning grid but not setting field variable.

	 A common situation is to have a grid that is already defined that
	 needs to be cloned and have field variables set through some other
	 mechanism.  For example, one might create a standard grid and then
	 plan to load it with values from a different grid or compute values
	 of a field variable at the grid points.

	 \param g Grid to be cloned.
	*/
	GCLscalarfield(const GCLgrid& g);
        /*! \brief Constuct from a file.

          A standard simple way to build any object is from a data file.
          This provides that interface.  The interface is relatively
          generic, however, as the format string provides a mechanism
          to specify variations in possible file format.  In the
          current implementation the file name should not contain
          a .ext string (e.g. mygrid.dat) as this may be used to
          split a header section from a data section.

          \param fname is the file name to read (may be just a root
          name depending on format option)
          \param format is a string describing allowed data format.
          Current default is split structure with one file used to
          hold an antelope pf to store scalar attributes and a second
          file used to store the actual data. Most applications will
          likely want to force default and use only one argument to
          this constructor.
          \parm enforce_object_type is used to turn type checking
          on an off.   When true (default) file object type must
          match the type of the constructor.   Useful for subclasses
          in some situations and to read a grid linked to a field
          ignoring the field data.
          */
  GCLscalarfield(const std::string fname, const std::string format=default_output_format,
                const bool enforce_object_type=true);

	/*! Metadata driven constructor.

	This constructor was added to allow this package to work with MongoDB.
	It acts very similar to the file-based constructors but instead of
	reading a pf representation of the attributes it assumes the same and a
	bit more are found in a Metadata container.   The idea of this constructor
	is that it will be used with pymongo.  The pymongo function reads the
	Metadata attributes from a MongoDB document as in MsPASS.   The grid
	data, however, is always assumed to be present in an external file
	accessible from the compute node calling this constructor.  I say that
	because that model may fail in cloud systems.

	One difference with this constructor compared to the file based constructor
	is that the following attributes are required that are NOT written with
	a file-based save to the pf file:

	dir - directory name string
	grid_data_file - file base name (with extension)
	grid_data_file_extension - file name "extension for data portion of file extension
		(normally should use keyword used internally as ".dat" but constructor
		allows others)
	grid_data_foff - integer file offset to locate start of block of data for
	 grid coordinates of this object.
	field_data_file - file name where field data is stored (can and often is
   same as grid_data_file with different foff)
	field_data_extension - as in grid_data_extension and usually the same
	field_data_foff - integer (long) file offset in bytes of block of data
	  that defines the field data.

	\param md is the Metadata container normally assumed produced by cracking a
		MongoDB document.

	\throw GCLgridError for a long list of possible failuers.   Any
		exceptions mean the object was not constucted.
	 */
  GCLscalarfield(const mspass::utility::Metadata& md);
	/*!
	 Destructor.
	 Nontrivial destructor as it has to destroy components
	 stored in plain C arrays.
	*/
	~GCLscalarfield();
	/** Zeros the field variable */
	void zero();
	/** Standard assignment operator */
	GCLscalarfield& operator=(const GCLscalarfield&);
        /*! \brief Save to a file.

          Many applications will find it more convenient to write
          this object to standard data files.  This abstracts the
          write process putting result in one or more formats.
          The concept here is the most applications would use
          the defaults for the the format specifications.
          Alternate formats are defined by using a different keyword
          to this method and implementing that alternative format
          in the code for this method.  Similarly alternate formats
          should implement a parallel code change for the file-based
          constructor which has a similar format argument.

          \param fname is the name of the output file to save data
          \param dir is a directory name to store data.  If zero
            length assume current directory.
          \param format is a keyword used to describe the format.
            (current default is a two file output separating the
            header and data.  Hence users should avoid .ext
            file names like myfile.dat.)
					\return copy of Metadata container with attributes defining the
						 object sans grid data.

          \exception GCLgridError is throw if save fails.
        */
  mspass::utility::Metadata save(const std::string fname, const std::string dir,
                const std::string format=default_output_format);
	/*! Returns all scalar attributes of the object in Metadata container. */
	mspass::utility::Metadata get_attributes() const;
	/*! Get the field value at specified grid cell */
	double get_value(const int i, const int j) const;


	/*! Sets the value of the field variable at grid cell i,j */
	void set_value(const double newvalue, const int i, const int j);
	/*! \brief Add one field to another.

	This acts like the += operator for simple types, but does so
	for an entire field defined on a grid.  As for simple types the
	field on the right hand side is accumulated to the field on
	the left had side of the operator.  This works even if the
	the two fields have different transformation properties and
	do not match in positions.  The most important limitation to
	understand, however, is basic sampling.  The field to contain
	the accumulation needs to be as denser or more densely sampled
	as the field on the right hand side to avoid aliasing effects.
	This is exactly like time series sampling -- upsampling can just be
	done by interpolation but downsampling requires an antialiasing
	(smoother) filter.
	*/
	GCLscalarfield& operator+=(GCLscalarfield& other);
	/*!
	Multiply all field values by a constant scale factor.

	\param c constant by which the field is to be scaled.
	*/
	GCLscalarfield& operator*=(const double c);
	/*!
	Linear interpolation function for a field.

	This is one of the core algorithms described in detail in the
	Fan and Pavlis (in review) paper.  It is intended as a low level
	routine accessible to the user if desired.  Most applications, however,
	will likely prefer to use a higher level application of this method
	through something like the += operator.  It is CRITICAL to recognize
	that this function must be called AFTER a previous call to lookup on
	the same point.  This routine blindly uses the index stored in the
	object and will return total garbage if the lookup was not called or
	returned an error condition that was not handled.  i.e. for efficieny
	this function simply assumes the interpolation requested will work.
	For this reason it will never throw an exception.
	\return interpolated value of the field at the requested point.
	\param x1p - Cartesian x1 coordinate of point to where the field is to be interpolated
	\param x2p - Cartesian x2 coordinate of point to where the field is to be interpolated
	\param x3p - Cartesian x3 coordinate of point to where the field is to be interpolated
	*/
	double interpolate(const double x1p, const double x2p, const double x3p);
	/*!
	 stream output operator for a field.
	 Format is:
	 <pre>
	@@line 1:  n1 n2
	@@line 2+:  x1, x2, x3, lat(deg), lon(deg), r(km), val in C output order
	@@		(index 2 most rapidly varying).
	</pre>
	*/
	friend ostream& operator << (ostream&,GCLscalarfield&);
};

/** Two-dimensional vector field defined on a GCLgrid framework. */
class GCLvectorfield : public GCLgrid
{
public:
	/** Number of components to vectors stored in this field */
	int nv;
	/*!
	Vector field variable stored in a three-dimensional C array.
	Index of the field is parallel with coordinate arrays x1, x2, and x3
	that define grid point positions in space.  The last index is the
	vector component.
	*/
	double ***val;

	/** Default constructor.  */
	GCLvectorfield();
	/*!
	Simple constructor.  Allocates space for x1, x2, and x3 arrays and initializes
	other data attributes to zero.

	 \param n1size number of grid points on generalized coordinate axis 1.
	 \param n2size number of grid points on generalized coordinate axis 2.
	 \param n3size number of grid points on generalized coordinate axis 3.
	*/
	GCLvectorfield(const int,const int,const int);
	/** Standard copy constructor. */
	GCLvectorfield(const GCLvectorfield& parent);
	/*!  \brief Construct a vector field using a parent grid geometry.

  Sometimes we have a grid geometry that we want to attach field data to.
	This constructor does that initializing the vector components to all 0s.

	\param g is the 2D grid to use as the basis for the field to be constructed
	\param nv is the number of components to create for each grid point (vector size)
	*/
	GCLvectorfield(const GCLgrid &g,const int nv);
        /*! \brief Construct from a file.

          A standard simple way to build any object is from a data file.
          This provides that interface.  The interface is relatively
          generic, however, as the format string provides a mechanism
          to specify variations in possible file format.  In the
          current implementation the file name should not contain
          a .ext string (e.g. mygrid.dat) as this may be used to
          split a header section from a data section.

          \param fname is the file name to read (may be just a root
          name depending on format option)
          \param format is a string describing allowed data format.
          Current default is split structure with one file used to
          hold an antelope pf to store scalar attributes and a second
          file used to store the actual data. Most applications will
          likely want to force default and use only one argument to
          this constructor.
          \parm enforce_object_type is used to turn type checking
          on an off.   When true (default) file object type must
          match the type of the constructor.   Useful for subclasses
          in some situations and to read a grid linked to a field
          ignoring the field data.
          */
  GCLvectorfield(const std::string fname, const std::string format=default_output_format,
                const bool enforce_object_type=true);
	/*! Metadata driven constructor.

	This constructor was added to allow this package to work with MongoDB.
	It acts very similar to the file-based constructors but instead of
	reading a pf representation of the attributes it assumes the same and a
	bit more are found in a Metadata container.   The idea of this constructor
	is that it will be used with pymongo.  The pymongo function reads the
	Metadata attributes from a MongoDB document as in MsPASS.   The grid
	data, however, is always assumed to be present in an external file
	accessible from the compute node calling this constructor.  I say that
	because that model may fail in cloud systems.

	One difference with this constructor compared to the file based constructor
	is that the following attributes are required that are NOT written with
	a file-based save to the pf file:

	dir - directory name string
	grid_data_file - file base name (with extension)
	grid_data_file_extension - file name "extension for data portion of file extension
		(normally should use keyword used internally as ".dat" but constructor
		allows others)
	grid_data_foff - integer file offset to locate start of block of data for
	 grid coordinates of this object.
	field_data_file - file name where field data is stored (can and often is
   same as grid_data_file with different foff)
	field_data_extension - as in grid_data_extension and usually the same
	field_data_foff - integer (long) file offset in bytes of block of data
	  that defines the field data.

	\param md is the Metadata container normally assumed produced by cracking a
		MongoDB document.

	\throw GCLgridError for a long list of possible failuers.   Any
		exceptions mean the object was not constucted.
	 */
	GCLvectorfield(const mspass::utility::Metadata& md);
	/** Standard assignment operator. */
	GCLvectorfield& operator=(const GCLvectorfield&);
	/*!
	Destructor.
	Note the same precautions about application of the default constructor as noted
	in GCLgrid::~GCLgrid
	*/
	~GCLvectorfield();
	/** Zeros the field variable */
	void zero();
        /*! \brief Save to a file.

          Many applications will find it more convenient to write
          this object to standard data files.  This abstracts the
          write process putting result in one or more formats.
          The concept here is the most applications would use
          the defaults for the the format specifications.
          Alternate formats are defined by using a different keyword
          to this method and implementing that alternative format
          in the code for this method.  Similarly alternate formats
          should implement a parallel code change for the file-based
          constructor which has a similar format argument.

          \param fname is the name of the output file to save data
          \param dir is a directory name to store data.  If zero
            length assume current directory.
          \param format is a keyword used to describe the format.
            (current default is a two file output separating the
            header and data.  Hence users should avoid .ext
            file names like myfile.dat.)
					\return copy of Metadata container with attributes defining the
						 object sans grid data.
          \exception GCLgridError is throw if save fails.
        */

  mspass::utility::Metadata save(const std::string fname, const std::string dir,
                const std::string format=default_output_format);
	/*! Get the field value at specified grid cell */
	std::vector<double> get_value(const int i, const int j) const;

	/*! Set the field value at grid point i,j.

        The input vector, newvals, size is tested against this->nv.
        This simple method will throw an exception if the sizes do
        not match.  Note for this function to work in python the
        pybind11 code needs to define how an std::vector<double>
        is bound to python.  Following mspass I map this to a
        symbol there called DoubleVector.
        */
	void set_value(const std::vector<double>newvals, const int i, const int j);
	/*! Add one field to another.

	This acts like the += operator for simple types, but does so
	for an entire field defined on a grid.  As for simple types the
	field on the right hand side is accumulated to the field on
	the left had side of the operator.  This works even if the
	the two fields have different transformation properties and
	do not match in positions.  The most important limitation to
	understand, however, is basic sampling.  The field to contain
	the accumulation needs to be as denser or more densely sampled
	as the field on the right hand side to avoid aliasing effects.
	This is exactly like time series sampling -- upsampling can just be
	done by interpolation but downsampling requires an antialiasing
	(smoother) filter.
	*/
	void operator+=(GCLvectorfield& other);
	/*!
	Multiply all field values by a constant scale factor.

	\param c constant by which the field is to be scaled.
	*/
	GCLvectorfield& operator*=(const double c);
	/*!
	Linear interpolation function for a field.

	This is one of the core algorithms described in detail in the
	Fan and Pavlis (in review) paper.  It is intended as a low level
	routine accessible to the user if desired.  Most applications, however,
	will likely prefer to use a higher level application of this method
	through something like the += operator.  It is CRITICAL to recognize
	that this function must be called AFTER a previous call to lookup on
	the same point.  This routine blindly uses the index stored in the
	object and will return total garbage if the lookup was not called or
	returned an error condition that was not handled.  i.e. for efficieny
	this function simply assumes the interpolation requested will work.
	For this reason it will never throw an exception.
	\return pointer to (freshly allocated with new) interpolated value of the
	vector field at the requested point.  The user should call delete [] to
	release the memory allocated for this vector to avoid a memory leak.

	\param x1p - Cartesian x1 coordinate of point to where the field is to be interpolated
	\param x2p - Cartesian x2 coordinate of point to where the field is to be interpolated
	\param x3p - Cartesian x3 coordinate of point to where the field is to be interpolated
	*/
	double *interpolate(const double x1p, const double x2p, const double x3p);
	/*! Returns all scalar attributes of the object in Metadata container. */
	mspass::utility::Metadata get_attributes() const;
	/*!
	 stream output operator for a vector field.
	 Format is:
	 <pre>
	@@line 1:  n1 n2 nv
	@@line 2+:  x1, x2, x3, lat(deg), lon(deg), r(km), val[0..nv-1] in C output order
	@@		(index 2 most rapidly varying).
	</pre>
	*/
	friend ostream& operator << (ostream&,GCLvectorfield&);
};
/** Three-dimensional scalar field defined on a GCLgrid framework. */
class GCLscalarfield3d : public GCLgrid3d
{
public:
	/*!
	Scalar field variable stored in a three-dimensional C array.
	Index of the field is parallel with coordinate arrays x1, x2, and x3
	that define grid point positions in space.
	*/
	double ***val;

	/** Default constructor. */
	GCLscalarfield3d();
	/*!
	Simple constructor.  Allocates space for x1, x2, x3, and val arrays and initializes
	object data attributes to zero.

	 \param n1size number of grid points on generalized coordinate axis 1.
	 \param n2size number of grid points on generalized coordinate axis 2.
	 \param n3size number of grid points on generalized coordinate axis 3.
	*/
	GCLscalarfield3d(const int n1size, const int n2size, const int n3size);
	/** Standard copy constructor. */
	GCLscalarfield3d(const GCLscalarfield3d&);
	/*!
	Partial copy constructor cloning grid but not setting field variable.

	A common situation is to have a grid that is already defined that
	needs to be cloned and have field variables set through some other
	mechanism.  For example, one might create a standard grid and then
	plan to load it with values from a different grid or compute values
	of a field variable at the grid points.

	\param g Grid to be cloned.
	*/
	GCLscalarfield3d(const GCLgrid3d &g);
        /*! \brief Constuct from a file.

          A standard simple way to build any object is from a data file.
          This provides that interface.  The interface is relatively
          generic, however, as the format string provides a mechanism
          to specify variations in possible file format.  In the
          current implementation the file name should not contain
          a .ext string (e.g. mygrid.dat) as this may be used to
          split a header section from a data section.

          \param fname is the file name to read (may be just a root
          name depending on format option)
          \param format is a string describing allowed data format.
          Current default is split structure with one file used to
          hold an antelope pf to store scalar attributes and a second
          file used to store the actual data. Most applications will
          likely want to force default and use only one argument to
          this constructor.
          */
        GCLscalarfield3d(const std::string fname,
					const std::string format=default_output_format);

	/*! Metadata driven constructor.

	This constructor was added to allow this package to work with MongoDB.
	It acts very similar to the file-based constructors but instead of
	reading a pf representation of the attributes it assumes the same and a
	bit more are found in a Metadata container.   The idea of this constructor
	is that it will be used with pymongo.  The pymongo function reads the
	Metadata attributes from a MongoDB document as in MsPASS.   The grid
	data, however, is always assumed to be present in an external file
	accessible from the compute node calling this constructor.  I say that
	because that model may fail in cloud systems.

	One difference with this constructor compared to the file based constructor
	is that the following attributes are required that are NOT written with
	a file-based save to the pf file:

	dir - directory name string
	grid_data_file - file base name (with extension)
	grid_data_file_extension - file name "extension for data portion of file extension
		(normally should use keyword used internally as ".dat" but constructor
		allows others)
	grid_data_foff - integer file offset to locate start of block of data for
	 grid coordinates of this object.
	field_data_file - file name where field data is stored (can and often is
   same as grid_data_file with different foff)
	field_data_extension - as in grid_data_extension and usually the same
	field_data_foff - integer (long) file offset in bytes of block of data
	  that defines the field data.


	\param md is the Metadata container normally assumed produced by cracking a
		MongoDB document.

	\throw GCLgridError for a long list of possible failuers.   Any
		exceptions mean the object was not constucted.
	 */
	GCLscalarfield3d(const mspass::utility::Metadata& md);
	/*!
	Destructor.
	Note the same precautions about application of the default constructor as noted
	in GCLgrid3d::~GCLgrid3d
	*/
	~GCLscalarfield3d();
	/** Zeros field variable array */
	void zero();
	/*! Get the field value at specified grid cell */
	double get_value(const int i, const int j, const int k) const;
	/*! \brief Get the field value at a specified point.

	This overloaded method does a lookup and interpolation before returning
	a value.  It should normally always be proceeded by a call to
	GCLgrid3d::point_is_inside_grid to make sure the lookup wil look
	algorithm confirms the point can be found and is inside the grid.
	If lookup fails this function will throw an exception.  That means
	any processing function that uses this method in a batch process must
	either have an error handler or precede the call by one to
	point_is_inside_grid.

	\param gp is the point to retieve a value.
	*/
	double get_value(const Geographic_point gp) const;
	/*! Sets the value of the field variable at grid cell i,j */
	void set_value(const double newvalue, const int i, const int j, const int k);
	/** Standard assignment operator. */
	GCLscalarfield3d& operator=(const GCLscalarfield3d& parent);

        /*! \brief Save to a file.

          Many applications will find it more convenient to write
          this object to standard data files.  This abstracts the
          write process putting result in one or more formats.
          The concept here is the most applications would use
          the defaults for the the format specifications.
          Alternate formats are defined by using a different keyword
          to this method and implementing that alternative format
          in the code for this method.  Similarly alternate formats
          should implement a parallel code change for the file-based
          constructor which has a similar format argument.

          \param fname is the name of the output file to save data
          \param dir is a directory name to store data.  If zero
            length assume current directory.
          \param format is a keyword used to describe the format.
            (current default is a two file output separating the
            header and data.  Hence users should avoid .ext
            file names like myfile.dat.)
					\return copy of Metadata container with attributes defining the
						 object sans grid data.
          \exception GCLgridError is throw if save fails.
        */
        mspass::utility::Metadata save(const std::string fname, const std::string dir,
                const std::string format=default_output_format);
	/*!
	Add one field to another.

	This acts like the += operator for simple types, but does so
	for an entire field defined on a grid.  As for simple types the
	field on the right hand side is accumulated to the field on
	the left had side of the operator.  This works even if the
	the two fields have different transformation properties and
	do not match in positions.  The most important limitation to
	understand, however, is basic sampling.  The field to contain
	the accumulation needs to be as denser or more densely sampled
	as the field on the right hand side to avoid aliasing effects.
	This is exactly like time series sampling -- upsampling can just be
	done by interpolation but downsampling requires an antialiasing
	(smoother) filter.
	*/
	GCLscalarfield3d& operator+=(const GCLscalarfield3d& other);
	/*!
	Multiply all field values by a constant scale factor.

	\param c constant by which the field is to be scaled.
	*/
	GCLscalarfield3d& operator*=(const double c);
	/*! Interpolate a 3d scalar field.

	Usage and caveats are the same as described in GCLscalarfield.
	*/
	double interpolate(const double,const double,const double);
	double parallel_interpolate(const double x1p, const double x2p, const double x3p,
	                const int ix1_0, const int ix2_0, const int ix3_0) const;
	/*! Returns all scalar attributes of the object in Metadata container. */
	mspass::utility::Metadata get_attributes() const;
	/*!
	 stream output operator for a 3d scalar field.
	 Format is:
	 <pre>
	@@line 1:  n1 n2 n3
	@@line 2+:  x1, x2, x3, lat(deg), lon(deg), r(km), val in C output order
	@@		(index 3 most rapidly varying).
	</pre>
	*/
	friend ostream& operator << (ostream&,GCLscalarfield3d&);
};
/** Three-dimensional vector field defined on a GCLgrid framework. */
class GCLvectorfield3d : public GCLgrid3d
{
public:
	/** Number of components to each vector of the field. */
	int nv;
	/*!
	Scalar field variable stored in a four-dimensional C array.
	Index of the field is parallel with coordinate arrays x1, x2, and x3
	that define grid point positions in space.  The fourth index defines
	the component of each vector of the field.
	*/
	double ****val;

	/** Default constructor. */
	GCLvectorfield3d();
	/*!
	Simple constructor.  Allocates space for x1, x2, x3, and val arrays and initializes
	object data attributes to zero.

	 \param n1size number of grid points on generalized coordinate axis 1.
	 \param n2size number of grid points on generalized coordinate axis 2.
	 \param n3size number of grid points on generalized coordinate axis 3.
	 \param nvsize number of components of vectors in this field
	*/
	GCLvectorfield3d(const int n1size, const int n2size,
		const int n3size, const int nvsize);
	/** Standard copy constructor. */
	GCLvectorfield3d(const GCLvectorfield3d& parent);
	/*!
	Partial copy constructor cloning grid but not setting field variable.

	A common situation is to have a grid that is already defined that
	needs to be cloned and have field variables set through some other
	mechanism.  For example, one might create a standard grid and then
	plan to load it with values from a different grid or compute values
	of a field variable at the grid points.

	\param g Grid to be cloned.
	\param nvsize number of vector components for vector field val array.
	*/
	GCLvectorfield3d(const GCLgrid3d &g,const int nvsize);
        /*! \brief Constuct from a file.

          A standard simple way to build any object is from a data file.
          This provides that interface.  The interface is relatively
          generic, however, as the format string provides a mechanism
          to specify variations in possible file format.  In the
          current implementation the file name should not contain
          a .ext string (e.g. mygrid.dat) as this may be used to
          split a header section from a data section.

          \param fname is the file name to read (may be just a root
          name depending on format option)
          \param format is a string describing allowed data format.
          Current default is split structure with one file used to
          hold an antelope pf to store scalar attributes and a second
          file used to store the actual data. Most applications will
          likely want to force default and use only one argument to
          this constructor.
          */
  GCLvectorfield3d(const std::string fname, const std::string format=default_output_format);
	/*! Metadata driven constructor.

	This constructor was added to allow this package to work with MongoDB.
	It acts very similar to the file-based constructors but instead of
	reading a pf representation of the attributes it assumes the same and a
	bit more are found in a Metadata container.   The idea of this constructor
	is that it will be used with pymongo.  The pymongo function reads the
	Metadata attributes from a MongoDB document as in MsPASS.   The grid
	data, however, is always assumed to be present in an external file
	accessible from the compute node calling this constructor.  I say that
	because that model may fail in cloud systems.

	One difference with this constructor compared to the file based constructor
	is that the following attributes are required that are NOT written with
	a file-based save to the pf file:

	dir - directory name string
	grid_data_file - file base name (with extension)
	grid_data_file_extension - file name "extension for data portion of file extension
		(normally should use keyword used internally as ".dat" but constructor
		allows others)
	grid_data_foff - integer file offset to locate start of block of data for
	 grid coordinates of this object.
	field_data_file - file name where field data is stored (can and often is
   same as grid_data_file with different foff)
	field_data_extension - as in grid_data_extension and usually the same
	field_data_foff - integer (long) file offset in bytes of block of data
	  that defines the field data.


	\param md is the Metadata container normally assumed produced by cracking a
		MongoDB document.

	\throw GCLgridError for a long list of possible failuers.   Any
		exceptions mean the object was not constucted.
	 */
	GCLvectorfield3d(const mspass::utility::Metadata& md);
	/*!
	Destructor.
	Note the same precautions about application of the default constructor as noted
	in GCLgrid3d::~GCLgrid3d
	*/
	~GCLvectorfield3d();
	/** Zeros the field variable */
	void zero();
	/*! Get the field value at specified grid cell */
	std::vector<double> get_value(const int i, const int j, const int k) const;
	/*! \brief Get the field value at a specified point.

	This overloaded method does a lookup and interpolation before returning
	a value.  It should normally always be proceeded by a call to
	GCLgrid3d::point_is_inside_grid to make sure the lookup wil look
	algorithm confirms the point can be found and is inside the grid.
	If lookup fails this function will throw an exception.  That means
	any processing function that uses this method in a batch process must
	either have an error handler or precede the call by one to
	point_is_inside_grid.

	\param gp is the point to retieve a value.
	*/
	std::vector<double> get_value(const Geographic_point gp) const;
	/*! Set the field value at grid point i,j,k.

        The input vector, newvals, size is tested against this->nv.
        This simple method will throw an exception if the sizes do
        not match.  Note for this function to work in python the
        pybind11 code needs to define how an std::vector<double>
        is bound to python.  Following mspass I map this to a
        symbol there called DoubleVector.
        */
	void set_value(const std::vector<double>newvals,
           const int i, const int j, const int k);
	/** Standard assignment operator. */
	GCLvectorfield3d& operator=(const GCLvectorfield3d&);
        /*! \brief Save to a file.

          Many applications will find it more convenient to write
          this object to standard data files.  This abstracts the
          write process putting result in one or more formats.
          The concept here is the most applications would use
          the defaults for the the format specifications.
          Alternate formats are defined by using a different keyword
          to this method and implementing that alternative format
          in the code for this method.  Similarly alternate formats
          should implement a parallel code change for the file-based
          constructor which has a similar format argument.

          \param fname is the name of the output file to save data
          \param dir is a directory name to store data.  If zero
            length assume current directory.
          \param format is a keyword used to describe the format.
            (current default is a two file output separating the
            header and data.  Hence users should avoid .ext
            file names like myfile.dat.)
					\return copy of Metadata container with attributes defining the
						 object sans grid data.
          \exception GCLgridError is throw if save fails.
        */
  mspass::utility::Metadata save(const std::string fname, const std::string dir,
                const std::string format=default_output_format);
	/*!
	Add one field to another.

	This acts like the += operator for simple types, but does so
	for an entire field defined on a grid.  As for simple types the
	field on the right hand side is accumulated to the field on
	the left had side of the operator.  This works even if the
	the two fields have different transformation properties and
	do not match in positions.  The most important limitation to
	understand, however, is basic sampling.  The field to contain
	the accumulation needs to be as denser or more densely sampled
	as the field on the right hand side to avoid aliasing effects.
	This is exactly like time series sampling -- upsampling can just be
	done by interpolation but downsampling requires an antialiasing
	(smoother) filter.
	*/
	GCLvectorfield3d& operator+=(const GCLvectorfield3d& other);
	GCLvectorfield3d operator+(const GCLvectorfield3d& other) const;
	/*!
	Multiply all field values by a constant scale factor.

	\param c constant by which the field is to be scaled.
	*/
	GCLvectorfield3d& operator*=(const double c);
	/*! Returns all scalar attributes of the object in Metadata container. */
	mspass::utility::Metadata get_attributes() const;
	/*! Interpolate a 3d scalar field.

	Usage and caveats are the same as described in
	    GCLvectorfield::interpolate
	 \return pointer to vector of doubles allocated with new [].
	 User must take caution to free this array to avoid memory leaks.
	*/
	double *interpolate(const double,const double,const double);
	double *parallel_interpolate(const double x1p, const double x2p, const double x3p,
	                const int i0, const int j0, const int k0) const;
	/*!
	 stream output operator for a 3d scalar field.
	 Format is:
	 <pre>
	@@line 1:  n1 n2 n3 nv
	@@line 2+:  x1, x2, x3, lat(deg), lon(deg), r(km), val[0..nv-1] in C output order
	@@		(index 3 most rapidly varying).
	</pre>
	*/
	friend ostream& operator << (ostream&,GCLvectorfield3d&);
};

// Function prototypes and templates from here to end of file.

/*! \brief Returns distance from the center of the Earth (km) of the standard ellipsoid at a specified latitude.

The reference ellipsoid depends only on latitude.  A GCLgrid uses the reference ellipsoid
as the reference datum to define a depth.

\return distance (in kilometers) from the center of the Earth to the sea level geoid surface at the requested latitude.
\param lat latitude (in radians) for which the standard ellipse radius is requested.
*/
double r0_ellipse(const double);
/*!
 retrieves a path along a gridline of a GCLgrid3d object.

 see man(3) extract_gridline.
*/
mspass::utility::dmatrix *extract_gridline(const GCLgrid3d& grid, const int ix1,
	const int ix2, const int ix3, const int comp, const bool reverse);
/*!
 Integrate a 3D field variable along a predefined path.

 see man(3) pathintegral.
*/
vector<double> pathintegral(const GCLscalarfield3d& field,const mspass::utility::dmatrix& path);
/*!
Transformation from standard spherical to local coordinates.

 see man(3) ustrans.
*/
mspass::utility::dmatrix ustrans(const GCLgrid& g, const double lat, const double lon);
/*!
Initialize a field with a layered structure.

It is often useful to initialize a 3d field to a depth dependent
function.  This function uses 1d interpolation to take a 1d function
specified in a depth format and project it into all parts of a field.
A typical example would be initializing a 3d velocity model in tomography
to a 1d starting model.

\param field field to to initialized.
\param val1d vector containing field variables at 1d grid points.
\param z1d parallel vector to val1d containing depths (NOT Radius) of the
   grid points with values stored in the val1d vector.
\param grad vector of gradients (forward looking) of field between 1d grid points.
   (Use 0.0 for all values for constant val layer models).
*/
void initialize_1Dscalar(const GCLscalarfield3d& field,
  const vector<double> val1d,
	  const vector<double> z1d,
		  const vector<double>grad);
/*!
Overloaded version of function by same name.

See long form for full description.  This form is used for block models
where the gradient is forced to zero between 1d grid points.

@see initialize_1Dscalar[0]
*/
void initialize_1Dscalar(const GCLscalarfield3d& field,
  const vector<double> val1d,const vector<double> z1d);
/*!
Map a path from one grid coordinate system to another.

Sometimes one has a line (path) in one grid that one wants to map into
another.  If the two grids have identical transformation properties
(can be determined with == or != operators) this is not necessary, but
if they are not congruent the points need to be converted between the
two coordinate systems.  This function simplifies that process.

\param parentgrid grid the path to be converted was originally defined in.
\param path path defined as an 3 x n array of points in the coordinate
	system of parentgrid.
\param newpathgrid grid onto which path is to be mapped.

\return a 3 X n matrix of points in the newpathgrid coordinate system.
*/
mspass::utility::dmatrix remap_path(const GCLgrid3d& parentgrid,
	const mspass::utility::dmatrix& path, const GCLgrid3d& newpathgrid);
/*!
Saves a 3d scalarfield to a stream in Data Explorer's native forma (dx).

\param g grid to be written to output stream.
\param out output stream.
*/
void dx_output(GCLscalarfield3d& g, ostream& out);

/*!
Allocate memory for a four dimensional array.

The GCLgrid library uses a contiguous block of memory to hold
arrays and uses C pointers to create the indexing needed for
standard C style subscripting.
A contiguous block is used in preference to a potentially
millions of malloc calls if the array were built in pieces.
This function creates a 4d array in this manner.
\param n1 size of index 1 of array.
\param n2 size of index 2 of array.
\param n3 size of index 3 of array.
\param n4 size of index 4 of array.
*/
double ****create_4dgrid_contiguous(const int n1, const int n2, const int n3, const int n4);
/*!
Allocate memory for a three dimensional array.

The GCLgrid library uses a contiguous block of memory to hold
arrays and uses C pointers to create the indexing needed for
standard C style subscripting.
A contiguous block is used in preference to a potentially
millions of malloc calls if the array were built in pieces.
This function creates a 3d array in this manner.
\param n1 size of index 1 of array.
\param n2 size of index 2 of array.
\param n3 size of index 3 of array.
*/
double ***create_3dgrid_contiguous(const int n1, const int n2, const int n3);
/*!
Allocate memory for a two dimensional array.

The GCLgrid library uses a contiguous block of memory to hold
arrays and uses C pointers to create the indexing needed for
standard C style subscripting.
This follows the pattern for 3 and 4d arrays for consistency,
although the reasons for needing a contiguous block in this
case are not as important.
This function creates a 2d array in this manner.
\param n1 size of index 1 of array.
\param n2 size of index 2 of array.
*/
double **create_2dgrid_contiguous(const int n1, const int n2);
/*!
Plain C destructor for a four-dimensional array.

This a companion free function to destroy an array
created by create_4dgrid_contiguous  .
The the GCLgrid library it is hidden in the destructor
for different objects.

\param n1 size of index 1 of array.
\param n2 size of index 2 of array.
\param n3 size of index 3 of array.
*/
void free_4dgrid_contiguous(double ****array,const int n1, const int n2, const int n3);
/*!
Plain C destructor for a three-dimensional array.

This a companion free function to destroy an array
created by create_3dgrid_contiguous  .
The the GCLgrid library it is hidden in the destructor
for different objects.

\param n1 size of index 1 of array.
\param n2 size of index 2 of array.
*/
void free_3dgrid_contiguous(double ***array,const int n1, const int n2);
/*!
Plain C destructor for a two-dimensional array.

This a companion free function to destroy an array
created by create_2dgrid_contiguous  .
The the GCLgrid library it is hidden in the destructor
for different objects.

\param n1 size of index 1 of array.
*/
void free_2dgrid_contiguous(double **array,const int n1);
/*!
FORTRAN function that does interpolation for a distorted cube.

This is a core function for interpolation of a field using shape
functions for a distorted cube element.  Computes a Jacobian to map
actual geometry into a standard space and computes weights using the
standard 8 point (cube) serendipity shape functions.

\author Kagan Tuncay
\param xx 3 vector defining point in space interpolation is requested
\param coord  array holding actual coordinates of 8 corner points
		of distorted cube.  xx assumed to be inside this element.
\param fun  return vector of length 8.  On return holds weights of
		to used in sum of values of function on 8 corners.
*/
void fme_weights_ (double *xx, double *coord, double *fun);
/*!
Converts a true depth to depth in a flattened coordinate system.

Some packages use the flattening transformation as an approximation
map spherical shells into a Cartesian reference frame.  This function
converts depths to flattened depths.

\author Steve Roecker translation to C by Pavlis
\param z - actual depth
\return depth with the flattening transformation applied to z
*/
double flatz(const double z);
/*!
Inverse flattening transformation of depth.

Some packages use the flattening transformation as an approximation
map spherical shells into a Cartesian reference frame.  This function
depths in the flattened earth coordinate system back to true depth.

\author Steve Roecker translation to C by Pavlis
\param z - flattened coordinate depth
\return actual depth in the earth
*/
double uflatz(const double z);
/*!
Applies flattening transformation to a velocity.

Some packages use the flattening transformation as an approximation
map spherical shells into a Cartesian reference frame.  This function
converts velocity at a given depth to the value with the flattening
transformation applied.

\author Steve Roecker translation to C by Pavlis
\param v - actual velocity at depth z
\param z - actual depth
\return velocity with the flattening transformation applied to v at z
*/
double flatvel(const double v,const double z);
/*!
Inverse flattening transformation of a velocity.

Some packages use the flattening transformation as an approximation
map spherical shells into a Cartesian reference frame.  This function
converts velocity in a flattened coordinate system back to the true
value at a depth that has to be determined from applying uflatz to z.

\author Steve Roecker translation to C by Pavlis
\param v - flattened coordinate velocity
\param z - flattened coordinate depth
\return true velocity with flattening correction removed
*/
double uflatvel(const double v, const double z);

/*!
Extract one component from a 2d vector field into a parallel scalar.

It is often useful to extract one component from a vector field
and treat that quantity as a scalar field.  This function
encapsulates that idea.  Note that it returns a pointer to a newly
allocated scalarfield object that the caller must deal with.
Use an shared_ptr to contain this pointer of just deal with the
usual rules of handling dynamically allocated objects to avoid
memory leaks.

\author Gary L. Pavlis
\exception int exception and posts and error to elog if the requested
     component is outside the range of the field.

\param f input vector field to be converted.
\param component component number to extract.  Assumes C convention
     of 0 being the first component.
*/
GCLscalarfield *extract_component(const GCLvectorfield& f,const int component);
/*!
Extract one component from a 2d vector field into a parallel scalar.

It is often useful to extract one component from a vector field
and treat that quantity as a scalar field.  This function
encapsulates that idea.  Note that it returns a pointer to a newly
allocated scalarfield object that the caller must deal with.
Use an shared_ptr to contain this pointer of just deal with the
usual rules of handling dynamically allocated objects to avoid
memory leaks.

\author Gary L. Pavlis
\exception int exception and posts and error to elog if the requested
     component is outside the range of the field.

\param f input vector field to be converted.
\param component component number to extract.  Assumes C convention
     of 0 being the first component.
*/
GCLscalarfield3d *extract_component(const GCLvectorfield3d& f,const int component);
/*!
Remap one grid to coordinate system of another.

Sometimes it is useful when dealing with multiple grids to
allow an algorithm to make an assumption that all grids in the
set have a common coordinate system.  This can avoid the
overhead of conversion of points to and from geographical
coordinates.  Experience has shown this is a nontrivial
calculation and needs to be minimized for algorithms that
might require large numbers of such conversions.
This function maps one grid to the coordinate system of the
one passed as "pattern".  This transformation is defined
in the BasicGCLgrid, lowest member of the heirarchy so it
higher levels should be cast to a BasicGCLgrid to be used
as a pattern.

Note that remap_grid can be called on field objects derived
from this one with no effect as the grid geometry is not
altered.  Only the coordinate system changes.
Note the grid is altered in place so the grid object passed will
be modified after this function completes.  The exception is if
the grids are already congruent in which case it silently does nothing.

\author Gary L. Pavlis

\param g grid to be remapped.
\param pattern grid whose coordinate system is to be used for
   new version of grid.
*/
void remap_grid(GCLgrid& g, const BasicGCLgrid& pattern);
/*!
Remap one grid to coordinate system of another.

Sometimes it is useful when dealing with multiple grids to
allow an algorithm to make an assumption that all grids in the
set have a common coordinate system.  This can avoid the
overhead of conversion of points to and from geographical
coordinates.  Experience has shown this is a nontrivial
calculation and needs to be minimized for algorithms that
might require large numbers of such conversions.
This function maps one grid to the coordinate system of the
one passed as "pattern".  This transformation is defined
in the BasicGCLgrid, lowest member of the heirarchy so it
higher levels should be cast to a BasicGCLgrid to be used
as a pattern.

Note that remap_grid can be called on field objects derived
from this one with no effect as the grid geometry is not
altered.  Only the coordinate system changes.
Note the grid is altered in place so the grid object passed will
be modified after this function completes.  The exception is if
the grids are already congruent in which case it silently does nothing.

\author Gary L. Pavlis

\param g grid to be remapped.
\param pattern grid whose coordinate system is to be used for
   new version of grid.
*/
void remap_grid(GCLgrid3d& g, const BasicGCLgrid& pattern);

/*!
 \brief Remap one grid to coordinate system of another.

   Sometimes it is useful when dealing with multiple grids to
   allow an algorithm to make an assumption that all grids in the
   set have a common coordinate system.  This can avoid the
   overhead of conversion of points to and from geographical
   coordinates.  Experience has shown this is a nontrivial
   calculation and needs to be minimized for algorithms that
   might require large numbers of such conversions.
   This function uses the set of parameters required to
   define this transformation as arguments.  It is
   a close cousin to two related procedures that use a
   grid object.  Use this one when it is inconvenient to
   have to actually load another grid as a pattern.

   Note that remap_grid can be called on field objects derived
   from this one with no effect as the grid geometry is not
   altered.  Only the coordinate system changes.
   Note the grid is altered in place so the grid object passed will
   be modified after this function completes.  The exception is if
   the grids are already congruent in which case it silently does nothing.

   \param g grid to be remapped (note it is a pointer).
   \param olat latitude (radians) of new origin
   \param olon longitude (radians) of new origin
   \param oradius radius (km) of new origin
   \param azn azimuth of the y axis in the new coordinate system.

   \author Gary L. Pavlis
*/
void remap_grid(BasicGCLgrid *g,
  const double olat, const double olon, const double oradius, const double azn);
/*! \brief Decimate a GCLgrid3d.

We sometimes want to decimate a grid. This procedure does this for
a 3D grid with variable decimation for each generalized coordinate
axis.  An added compliation is the fact that because GCLgrids hae
the x3 direction directed upward in earth coordinates and the free
surface is so special, we work the x3 coordinate backward compared
to the others.  i.e. we for the n3-1 point to be the n3-1 point in
the result, not the 0 points.

\param g parent grid that is to be decimated
\param dec1 decimation factor for x1
\param dec2 decimation factor for x2
\param dec3 decimamtion factor for x3

\return pointer to decimated grid object
*/
GCLgrid3d *decimate(const GCLgrid3d& g,const int dec1, const int dec2, const int dec3);

/*! Establish if binary data needs to be type swapped.

This function tries to fetch the name "datatype" from the received Metadata
container.   If it is defined it fetches the value.  If defined it must be
one of two values:  "u8" or "t8".   If it is neither it will throw a
GCLgridError exception.  It then determines the byte order of the
machine this code is running on.  "t8" is little endian and "u8" is
assumed big endian.  If the datatype and the actual do not match
the function returns true.   If they match it returns false.  Note the
default if datatype is not defined is little endian because Intel won the
byte order wars a while back. */
bool byte_swap_is_needed(const mspass::utility::Metadata& md);


}  // End of namespace
#endif
