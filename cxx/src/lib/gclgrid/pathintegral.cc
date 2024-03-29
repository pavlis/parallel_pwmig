#include <math.h>
#include <vector>
#include "pwmig/gclgrid/gclgrid.h"
#include "mspass/utility/dmatrix.h"
using namespace std;
using namespace pwmig::gclgrid;
using mspass::utility::dmatrix;
namespace pwmig::gclgrid
{
using mspass::utility::dmatrix;
/* General purpose utility to integrate a scalar field variable
along a path defined by the input matrix path.  The prototype
example of this for us seismologists is integration of slowness
along a ray path.  Only the cartesian system of the field is
used. The components of path are assumed stored in a 3xnc
matrix whose components use the same basis as the field object.
The algorithm works along the path.  If the lookup function
ever fails for any reason the path will be truncated at the
point where the lookup failed.  The caller should test for
this condition by verifying that the length of the output
vector of doubles is the same as the number of columns in
the path matrix.

The function returns a vector of doubles containing integral
along path of the field variable.

The routine throws an exception when the input matrix is not
correctly dimension.  This would normally be a programming error
and the caller may choose to not attempt to catch the error
except in debugging.

Author:  Gary L. Pavlis
Written:  June 2003
*/
vector <double> pathintegral(GCLscalarfield3d& field,const dmatrix& path)
{
	int sz;
	int npts;
	double outval,outval_last;
	double val1,val2;  // averaged to get field value for interval
	double dx;  // path increment
	vector<double> outvec;
	int i;

	sz=path.rows();
	npts = path.columns();
	if(sz!=3)
	{
	  throw GCLgridError("pathintegral:  input matrix of path coordinates has incorrect dimensions");
	}
	outvec.reserve(npts);
	// push 0 to the first point
	outvec.push_back(0.0);
	// This was found to alway be prudent
	// Sept 2012 changes:  this doesn't really do anything anymore but will
	// retain out of caution as the cost is tiny
	field.reset_index();
	/* New code for parallel_lookup.  This is a bit awkward because the std::vector
	container return doesn't match argument approach I used for holding the
	last integer position managed by parallel_lookup. */
	std::vector<int> lookup_origin;
	lookup_origin = field.get_lookup_origin();
	int ix1_0, ix2_0, ix3_0;
	ix1_0 = lookup_origin[0];
	ix2_0 = lookup_origin[1];
	ix3_0 = lookup_origin[2];
	for( i=1,outval=0.0,outval_last=0.0;i<npts;++i)
	{
		double dx1,dx2,dx3;
		/* the if is used to do a safe return if the ray traverses outside the grid */
		if( field.parallel_lookup(path(0,i-1),path(1,i-1),path(2,i-1),
	     ix1_0,ix2_0,ix3_0) ) break;
		val1=field.interpolate(path(0,i),path(1,i),path(2,i));
		if( field.parallel_lookup(path(0,i),path(1,i),path(2,i),
	     ix1_0,ix2_0,ix3_0) ) break;
		val2=field.interpolate(path(0,i),path(1,i),path(2,i));
		// This could be functionized, but I'll make it inline
		dx1 = path(0,i)-path(0,i-1);
		dx2 = path(1,i)-path(1,i-1);
		dx3 = path(2,i)-path(2,i-1);
		dx = sqrt(dx1*dx1+dx2*dx2+dx3*dx3);
		outval=outval_last+(val1+val2)*dx/2.0;
		outval_last = outval;
		outvec.push_back(outval);
	}
	return(outvec);
}
/* A path defined by a 3xn dmatrix is defined by the Cartesian reference frame in the
grid from which it is derived.  If one wants to use this path inside another grid,
which does not necessarily have the same Cartesian transformation, the path has
to be converted to the new reference frame.  This function does this.

pathgrid is the grid in which the curve defined by the dmatrix path was originally
defined.  othergrid is the new grid into which the path is to be mapped.  The
returned result is a new dmatrix of the same size as path, but defined in the
Cartesian reference frame for othergrid instead of pathgrid.
*/

dmatrix remap_path(const GCLgrid3d& pathgrid, const dmatrix& path, const GCLgrid3d& othergrid)
{
	int i;
	int m=path.columns();
	if( (path.rows()!=3) || (m<=0) )
		throw GCLgridError("remap_path:  input path matrix dimensions are invalid");
	dmatrix newpath(3,m);
	Geographic_point geo;
	Cartesian_point p;
	for(i=0;i<path.columns();++i)
	{
		geo=pathgrid.ctog(path(0,i),path(1,i),path(2,i));
		p=othergrid.gtoc(geo);
		newpath(0,i)=p.x1;
		newpath(1,i)=p.x2;
		newpath(2,i)=p.x3;
	}
	return(newpath);
}

} //end nameespace
