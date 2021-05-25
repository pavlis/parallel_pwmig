#ifndef _DEPTH_DEPENDENT_APERTURE_H_
#define _DEPTH_DEPENDENT_APERTURE_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "mspass/seismic/SlownessVector.h"
#include "mspass/utility/AntelopePf.h"
namespace pwmig::pwmigcore
{
class DepthDependentAperture
{
public:

	DepthDependentAperture(){};
        /*! Construct using linear segments specified by parallel vectors*/
	DepthDependentAperture(const int npt){
		t.reserve(npt);
		aperture.reserve(npt);
		cutoff.reserve(npt);
	};
  /*! Use Fresnel zone formula to construct.
  *
  * A useful depth dependent aperture formula is to use the
  * Fresnel zone size for vertical incident S in a constant
  * velocity medium.  This constructor implements that formula.
  * The aperture is parameterized internally by two parallel
  * vectors just like the manual constructors.  The dtau and ntau
  * variable set the uniform grid interval and number of points
  * respectively.
  *
  *\param vs is S wave velocity.
  *\param vp is P wave velocity.
  *\param period is wave period to compute the Fresnel zone size.
  *     Be aware this is effectively a width parameter because at
  *     zero lag the width becomes this value and increases with
  *     lag.
  *\param dtau is lag sample interval for formula
  *\param ntau is the number of lags to compute (Note that the
  *      interpolation method continues the value of (ntau-1)*dtau
  *      to infinity.
  *\param cutoffmultiplier is a scale factor used to compute cutoff
  *     parameter.  For this constructor cutoff is computed as this
  *     number times the aperture.
  *\param echoresult if set true will print computed aperture sizes
  *     using cout.
  */
  DepthDependentAperture(const double vs, const double vp,const double period,
            const double dtau,const int ntau,const double cutoff_multiplier,
						 const bool echoresult);
  /*! Parameter file driven constructor.

  Uses a pf format file through a nonantelope interface called
  a PfStyleMetadata object.

    \param md is the object created from a pf file.
    \param tag is an &Arr tag bounding he parameters for this object.
  */
	DepthDependentAperture(const mspass::utility::AntelopePf& md,
		  const std::string tag);
  /*! Copy constructor. */
	DepthDependentAperture(const DepthDependentAperture& parent);
	/*! Return aperture value for specified time lag.*/
	double get_aperture(double tnow) const;
	/*! Return cutoff distance for lag specified. */
	double get_cutoff(double tnow)const;
	/*! Standard assignment operator. */
  DepthDependentAperture& operator=(const DepthDependentAperture& parent);
	/*! Return number of points used to define the geometry of this aperture object. */
	size_t npoints() const{return t.size();};
	/*! Return the largest defined cutoff distance. */
	double maximum_cutoff() const;
	/*! Return the largest defined aperture size. */
	double maximum_aperture() const;
private:
  std::vector<double> t,aperture,cutoff;
	friend class boost::serialization::access;
	template<class Archive>
			 void serialize(Archive & ar, const unsigned int version)
	{
		ar & t;
		ar & aperture;
		ar & cutoff;
	};
};
}  // end namespace pwmig::pwmigcore
#endif
