#ifndef _HYPOCENTER_H_
#define _HYPOCENTER_H_
#include <string>
#include "mspass/utility/Metadata.h"
namespace pwmig::utility
{
/*! \brief Defines a source position and a set of useful methods with which a source
position can be associated.

Seismic processing is commonly revolved around geometry defined by source and
receiver positions.  This object encapsulates information about a source and
provides a set of methods useful for waveform processing that depends on
source information.   Note that all angle quantities are assumed defined in
radians.
\author Gary L. Pavlis
**/
class Hypocenter
{
public:
	/*!
	   Latitude of the source in radians.
	**/
	double lat;
	/*!
	   Longitude of the source in radians.
	**/
	double lon;
	/*!
	   Depth below sea level of the source in kilometers.
	**/
	double z;
	/*!
	   Origin time of the source in seconds.  If absolute time is being used
	   this is an epoch time.
	**/
	double time;
	/*! Default constructor.  Initializes all to zero*/
	Hypocenter(){lat=0.0; lon=0.0; z=0.0; time=0.0;};  // default
	/*!
	     Metadata object driven constructor.  Looks for these keywords:
	     source_lat, source_lon, source_depth, and source_time.
	     Note it is assumed the latitude and longitude values contained in
	     the parent metadata are given already in radians.
	**/
	Hypocenter(const mspass::utility::Metadata& );
	/*!
	     Fully parameterized constructor fill in all data members of
	     Hypocenter object.

	    \param lat0 latitude of the source (in radians)
	    \param lon0 longitude of the source (in radians)
	    \param z0 depth of the source in km
	    \param t0 orign time of the source (epoch time).

	**/
	Hypocenter(const double lat0, const double lon0,
		       const double z0, const double t0);
	/*!   Standard copy constructor */
	Hypocenter(const Hypocenter&);
	/*!   Standard assignment operator.*/
	Hypocenter& operator=(const Hypocenter&);
	/*!   Compute distance from source to surface position (lat0,lon0).*/
	double distance(const double lat0, const double lon0) const;
	/*!
	     Compute event to station azimuth.
	     Event to source azimuth is defined as the azimuth of the great
	     circle path joining event and station at the position of the
	     source and looking toward the station.
	    \param lat0 latitude of the station
	    \param lon0 longitude of the station
	**/
	double esaz(const double lat0, const double lon0) const;
	/*!
	     Compute station to event azimuth.
	     Station to event azimuth is defined as the azimuth of the
	     great circle path joining event and station but evaluated
	     at the position of the station and looking back toward
	     the source.
	    \param lat0 latitude of the station
	    \param lon0 longitude of the station
	**/
	double seaz(const double lat0, const double lon0) const;
};

} // End SEISPP namespace declaration
#endif
