#include "pwmig/gclgrid/gclgrid.h"
using namespace pwmig::gclgrid;
namespace pwmig::gclgrid
{
/*  forward and backward flattening transformation routines
are in this file.

It is common practice in regional scale studies to use the
flattening transform as an approximation for converting from
spherical to cartesian coordinates.  flatz converts from
depth to an equivalent flattened depth.  uflatz does
the inverse.  flatvel converts velocity and uflatvel
does the reciprocal.

This is a C translation of routines by the same name in
FORTRAN from Steve Roecker.
*/
const double REARTH=6371.0;
double flatz(const double z)
{
	return(REARTH*log(REARTH/(REARTH-z)));
}
double uflatz(const double z)
{
	return(REARTH*(1.0-exp(-z/REARTH)));
}
double flatvel(const double v, const double z)
{
	return(v*exp(z/REARTH));
}
double uflatvel(const double v, const double z)
{
	return(v*exp(-z/REARTH));
}
} //end namespace
