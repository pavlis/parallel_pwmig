#ifndef _RECTANGULAR_SLOWNESS_GRID_H_
#define RECTANGULAR_SLOWNESS_GRID_H_
#include <string>

#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "mspass/seismic/SlownessVector.h"
#include "mspass/utility/AntelopePf.h"
using mspass::seismic::SlownessVector;
namespace pwmig{
namespace pwmigcore{
/*! \brief This object defines a uniform grid of points in slowness space.

 In array processing it is common to need a grid of feasible slowness
 vectors and data are commonly stacked to define all of this family
 of plane wave stacks.  The object used here grids up slowness space
 in uniform spacing in ux and uy components.
 The gridding is defined by the lower left corner (uxlow, uylow),
 spacing in EW and NS directions (dux, duy), and the number of
 points in the x (EW) and y (NS) directions.
 Note the grid and all methods assume units of s/km, although
 this application is largely blind to units.
\author Gary L. Pavlis
**/
class RectangularSlownessGrid
{
public:

/*!
 Default constructor.
**/
	RectangularSlownessGrid();  // generic default is defined
/*!
 Fully parameterized constructor.

\param nm - name to be assigned to grid.
\param uxl - lower left corner ux.
\param uyl - lower left corner uy.
\param dux - dux (increment in x direction)
\param duy - duy (increment in y direction)
\param nx - number of grid points in x direction.
\param ny - number of grid points in y direction.
**/
	RectangularSlownessGrid(const std::string nm, const double uxl, const double uyl,
		const double dux,const double duy,const int nx, const int ny);

/* \brief Metadata driven constructor.

This is similar to the Pf driven method except the data is passed
through a more generic object in SEISPP called Metadata.  Keywords
defined in the Metadata object are identical to those in the pf version.
As for the pf constructor the required parameters are:
  "Slowness_Grid_Name" (string), "uxlow" (real), "uylow" (real),
  "nux" (int), "nuy" (int), "dux" (real), and "duy" (real).
  They define an nux by nuy regular grid in slowness space with
  spacing dux and duy respectively with the lower left corner of the
  grid at (uxlow,uylow).  The name is just a tag.

\param mdin is the Metadata object with attributes set to build this object.
\param tag - name to search in pf to describe this grid object.`
              The parameters to describe the object are assumed encased in an
              &Arr{ } construct with this tag.  This allows multiple grids to
              be defined in a single parameter file with different tags.
\exception MetadataGetError (child of SeisppError) is thrown if the
  required attributes are not defined in the Metadata object passed.
  */
  RectangularSlownessGrid(const mspass::utility::AntelopePf& mdin,
		const std::string tag);
/*!
 Standard copy constructor.
**/
	RectangularSlownessGrid(const  RectangularSlownessGrid&);
/*!
 Returns x component of slowness grid at index position i.
**/
	double ux(int i) const {return(uxlow+i*dux);};
/*!
 Returns y component of slowness grid at index position j.
**/
	double uy(int i) const  {return(uylow+i*duy);};
/*!
 Returns a slowness grid object for grid position (i,j).
\exception SeisppError object if i and j are outside range.
**/
	SlownessVector slow(int i, int j) const;
	/*! Return number of grid points in x direction*/
	int ux_size() const {return nux;};
	/*! Return number of grid points in y direction */
	int uy_size() const {return nuy;};
	/*! Standard assignment operator. */
	RectangularSlownessGrid& operator=(const RectangularSlownessGrid& parent);
	/*! Offset slowness grid data by a specified vector in slowness space. */
	RectangularSlownessGrid& operator+=(const mspass::seismic::SlownessVector& u0);
private:
	/*!
	 Name assigned to this grid object.
	**/
		std::string name;
	/*!
	 Minimum ux (East-west) component of grid.
	 The location of the lower left corner of the grid defined by this object is
	 at coordinates (uxlow,uylow).
	**/
		double uxlow;
	/*!
	 Minimum uy (North-south) component of grid.
	 The location of the lower left corner of the grid defined by this object is
	 at coordinates (uxlow,uylow).
	**/
		double uylow;
	/*!
	 Grid size (s/km) of slowness grid in EW component.
	**/
		double dux;
	/*!
	 Grid size (s/km) of slowness grid in NS component.
	**/
		double duy;
	/*!
	 Number of grid points in x (EW) direction of slowness grid.
	**/
		int nux;
	/*!
	 Number of grid points in y (NW) direction of slowness grid.
	**/
		int nuy;
  friend class boost::serialization::access;
  template<class Archive>
       void serialize(Archive & ar, const unsigned int version)
  {
    ar & name;
    ar & uxlow;
    ar & uylow;
    ar & dux;
    ar & duy;
    ar & nux;
    ar & nuy;
	};
};
} // end pwmigcore namespace
} // end pwmig namespace

#endif
