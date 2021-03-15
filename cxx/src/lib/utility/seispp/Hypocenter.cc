// These are methods for the Hypocenter object to compute common standard things.`
// These are mostly interface routines to Antelopes coords library.
// Note all returns are distances and angles are returned in units of radians.

#include "pwmig/utility/coords.h"
#include "pwmig/utility/stock.h"
#include "pwmig/utility/Hypocenter.h"
#include "mspass/seismic/keywords.h"
#include "mspass/utility/MsPASSError.h"
/* needed for rad function for conversion */
#include "mspass/utility/SphericalCoordinates.h"
using namespace std;
using namespace pwmig::utility;
using namespace mspass::seismic;

namespace SEISPP
{

// copy constructor and = operators commonly are implicitly needed
Hypocenter::Hypocenter(const Hypocenter& h0)
{
	lat = h0.lat;
	lon = h0.lon;
	z = h0.z;
	time = h0.time;
}
Hypocenter& Hypocenter::operator=(const Hypocenter& h0)
{
    if(this!=&h0)
    {
	lat = h0.lat;
	lon = h0.lon;
	z = h0.z;
	time = h0.time;
    }
    return(*this);
}
// This constructor creates a Hypocenter object using parameters
// read from a Metadata object.  It would be useful to have a
// mechanism to not have a frozen set of names used to fetch
// the required parameters but I'm not sure how to do that
// without making the interface clumsy.
//
// Note that things like TimeSeries objects can be used to create
// a Hypocenter through this mechanism through the use of a dynamic_cast
//
Hypocenter::Hypocenter(const Metadata& md)
{
	try {
		lat=md.get_double(SEISMICMD_slat);
		lon=md.get_double(SEISMICMD_slon);
		z=md.get_double(SEISMICMD_sdepth);
		time=md.get_double(SEISMICMD_stime);
		/* Assume lat and lon are in degrees so we have to convert them
		to radians.  Assume depth is already in km*/
		lat=mspass::utility::rad(lat);
		lon=mspass::utility::rad(lon);
	} catch (MsPASSError& mderr) {throw mderr;};
}
// the obvious fully parameterized constructor
Hypocenter::Hypocenter(const double lat0, const double lon0,
	const double z0, const double t0)
{
	lat=lat0;
	lon=lon0;
	z=z0;
	time=t0;
}



double Hypocenter::distance(const double lat0, const double lon0) const
{
	double epidist, az;

	dist(lat,lon,lat0,lon0,&epidist, &az);
	return(epidist);
}

double Hypocenter::esaz(const double lat0, const double lon0) const
{
	double epidist, az;

	dist(lat,lon,lat0,lon0,&epidist, &az);
	return(az);
}
double Hypocenter::seaz(const double lat0, const double lon0) const
{
	double epidist, az;

	dist(lat0,lon0,lat,lon,&epidist, &az);
	return(az);
}

} // Termination of namespace SEISPP definitions
