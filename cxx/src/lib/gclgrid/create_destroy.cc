#include <typeinfo>
#include <string.h>
#include <math.h>
#include <boost/core/demangle.hpp>
#include "mspass/utility/AntelopePf.h"
#include "mspass/utility/SphericalCoordinate.h"
#include "pwmig/dsap/stock.h"
#include "pwmig/dsap/coords.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "pwmig/gclgrid/gclgrid_subs.h"
#include "pwmig/gclgrid/swapbytes_pwmig.h"

namespace pwmig::gclgrid
{
using namespace std;
using namespace pwmig::gclgrid;
using mspass::utility::Metadata;
using mspass::utility::AntelopePf;
using mspass::utility::MsPASSError;
using mspass::utility::rad;


/* These are create and free routines for 2d, 3d, and 4d arrays in plain C.
They are used below for C++ constructors and destructors.  The MOST
IMPORTANT thing to note about the constructions used here is that I have
implemented these built on contiguous blocks of memory.  This is NOT necessarily
the way one would always build multidimension arrays in C.  This is probably a bad
idea for future maintainability of this code as someone might be tempted to
write their own destructor or add a constructor that used a different memory model.
Future users beware on this point.  I did this on purpose, however, under
and assumption it would improve performance.*/


/* Many gcl functions work with 3d C grids.  C implements multidimensioned
arrays pointer arrays, which makes this kind of ugly.  Anyway, this
function will create the F90 equivalent of  x[0:n1-1][0:n2-1][0:n3-1]
and return a ***double pointer for the 3d array.  The approach
used here allocates the full block of memory in one single call to
guarantee the result is in a contiguous block of memory.
This is useful for reading and writing this type of entity because
then a single call to fread and fwrite can be used.
Note, however, that the actual order in memory remains the C convention
not the FOTRAN form.  i.e. the rightmost index defines members contingous
in memory.
Author:  GAry Pavlis
Written:  June 2000
*/

double ***create_3dgrid_contiguous(const int n1, const int n2, const int n3)
{
	double ***ptr3d;
	double **ptr2ptr;
	double *ptr;
	int i,j;

	/* Previous version used an datascope/antelope macro allot
	that handled malloc errors automatically. Made manual here.
	*/
	//allot(double *,ptr,n1*n2*n3);
	//allot(double ***,ptr3d,n1);
	ptr=(double *)calloc(n1*n2*n3,sizeof(double));
	ptr3d=(double ***)calloc(n1,sizeof(double **));
	if( ptr==NULL || ptr3d==NULL)
	{
	  throw GCLgridError("create_3dgrid_contiguous:  calloc failure");
	}
	for(i=0;i<n1;++i)
	{
		//allot(double **,ptr2ptr,n2);
		ptr2ptr=(double **)calloc(n2,sizeof(double *));
		if(ptr2ptr==NULL)
		{
		  throw GCLgridError("create_3dgrid_contingous:  calloc failure while creating pointer to pointer indexing array");
		}
		ptr3d[i] = ptr2ptr;
	}
	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			ptr3d[i][j] = ptr + n2*n3*i + n3*j;
		}
	}
	return(ptr3d);
}
/* Free routine for a 3d array.  pointers arrays require extra attention
to avoid a memory leak.   This function assumes standard C indexing
*/
void free_3dgrid_contiguous(double ***x,const int n1, const int n2)
{
	int i;
	double *ptr;

	/* this clears the work space*/
	ptr = x[0][0];
	free((void *)ptr);

	/* The pointer arrays are still then and have to be freed
	seperately */
	for(i=0;i<n1;++i)  free((void *)x[i]);
	free((void *)x);
}
/* parallel routines to the above for 2d */
double **create_2dgrid_contiguous(const int n1, const int n2)
{
	double **ptr2ptr;
	double *ptr;
	int i;

	/* As above this used to used these antelope allot functions from
	stock.h*/
	//allot(double *,ptr,n1*n2);
	//allot(double **,ptr2ptr,n1);
	ptr=(double *)calloc(n1*n2,sizeof(double));
	ptr2ptr=(double **)calloc(n1,sizeof(double *));
	if(ptr==NULL || ptr2ptr==NULL)
	{
	  throw GCLgridError("create_2dgrid_contiguous:  malloc failure creating 2d grid workspace");
	}
	for(i=0;i<n1;++i)
		ptr2ptr[i] = ptr + n2*i;
	return(ptr2ptr);
}
void free_2dgrid_contiguous(double **x,const int n1)
{
	double *ptr;

	ptr = x[0];
	free((void *)ptr);
	free((void *)x);
}
//
// A bit out of order, but now similar construction for 4D
//
double ****create_4dgrid_contiguous(const int n1, const int n2, const int n3, const int n4)
{
	const string base_error("create_4dgrid_continuous:  ");
	double ****ptr4d;
	double ***ptr3d;
	double **ptr2ptr;
	double *ptr;
	int i,j,k;

	/* As above these previously used antelope's allot macro*/
	//allot(double *,ptr,n1*n2*n3*n4);
	//allot(double ****,ptr4d,n1);
	ptr=(double *)calloc(n1*n2*n3*n4,sizeof(double));
	if(ptr==NULL)
	{
	  stringstream ss;
	  ss << base_error
	    << "memory allocation failed for 4d array of size="
	    << n1<<"X"<<n2<<"X"<<n3<<"X"<<n4<<endl;
	  throw GCLgridError(ss.str());
	}
	ptr4d=(double ****)calloc(n1,sizeof(double ***));
	if(ptr4d==NULL)
	{
	  stringstream ss;
	  ss<< base_error
	    << "memory allocation failed for top level point array of length="
	    << n1<<endl;
	  throw GCLgridError(ss.str());
	}
	for(i=0;i<n1;++i)
	{
		//allot(double ***,ptr3d,n2);
		ptr3d=(double ***)calloc(n2,sizeof(double **));
		if(ptr3d==NULL)
		{
		  stringstream ss;
		  ss << base_error
		    << "memory allocation failed for 2nd level pointer array"
		    << " of length="<<n2<<endl;
		  throw GCLgridError(ss.str());
		}
		ptr4d[i] = ptr3d;
		for(j=0;j<n2;++j)
		{
			//allot(double **,ptr2ptr,n3);
			ptr2ptr=(double **)calloc(n3,sizeof(double *));
			if(ptr2ptr==NULL)
			{
			  stringstream ss;
			  ss<<base_error
			    << "memory allocation failed for 3rd level pointer array of length="
			    << n3<<endl;
			  throw GCLgridError(ss.str());
			}
			ptr4d[i][j]=ptr2ptr;
		}
	}
	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			for(k=0;k<n3;++k)
			{
				ptr4d[i][j][k] = ptr + n2*n3*n4*i
						+ n3*n4*j + n4*k;
			}
		}
	}
	return(ptr4d);
}
void free_4dgrid_contiguous(double ****x,const int n1, const int n2, const int n3)
{
	int i,j;
	double *ptr;

	/* this clears the work space*/
	ptr = x[0][0][0];
	free((void *)ptr);

	/* unwind and release the pointer arrays -- pretty ugly this deep */
	for(i=0;i<n1;++i)
	{
		for(j=0;j<n2;++j)
		{
			free((void *)x[i][j]);
		}
		free((void *)x[i]);
	}
	free((void *)x);
}
//
// C++ constructors
//
BasicGCLgrid::BasicGCLgrid()
{
	name=string("");
	lat0=0.0; lon0=0.0; r0=0.0;
	azimuth_y=0.0;  dx1_nom=0.0;  dx2_nom=0.0;
	n1=0;  n2=0;
	i0=0;  j0=0;
	x1low=0.0;  x1high=0.0;
	x2low=0.0;  x2high=0.0;
	x3low=0.0;  x3high=0.0;
}
// This could be defaulted as scalars are copied, I believe, but
// best to be explicit.
BasicGCLgrid::BasicGCLgrid(const BasicGCLgrid& g)
{
	name=g.name;
	lat0=g.lat0;
	lon0=g.lon0;
	r0=g.r0;
	azimuth_y=g.azimuth_y;
	dx1_nom=g.dx1_nom;
	dx2_nom=g.dx2_nom;
	n1=g.n1;
	n2=g.n2;
	i0=g.i0;
	j0=g.j0;
	x1low=g.x1low;
	x1high=g.x1high;
	x2low=g.x2low;
	x2high=g.x2high;
	x3low=g.x3low;
	x3high=g.x3high;
	for(int i=0;i<3;++i)
	{
		for(int j=0;j<3;++j) gtoc_rmatrix[i][j]=g.gtoc_rmatrix[i][j];
		translation_vector[i]=g.translation_vector[i];
	}
}

GCLgrid::GCLgrid (const int n1size, const int n2size) : BasicGCLgrid()
{
	n1=n1size;
	n2=n2size;
	x1=create_2dgrid_contiguous(n1size,n2size);
	x2=create_2dgrid_contiguous(n1size,n2size);
	x3=create_2dgrid_contiguous(n1size,n2size);
}

GCLgrid3d::GCLgrid3d (const int n1size, const int n2size,
       const int n3size, const bool fl) : BasicGCLgrid()
{
	n1=n1size;
	n2=n2size;
	n3=n3size;
	x1=create_3dgrid_contiguous(n1size,n2size,n3size);
	x2=create_3dgrid_contiguous(n1size,n2size,n3size);
	x3=create_3dgrid_contiguous(n1size,n2size,n3size);
  fast_lookup=fl;
}
/* This routine is used in all file-based constructors using a
   pf to store attributes. */
Metadata pfload_GCLmetadata(const string fname)
{
    try {
        /* fname is expected to be a root name with out and extension. 
        We thus need to extend the name with a pfext defined in gclgrid.h*/
        string fullname=fname+"."+pfext;
        AntelopePf md(fullname);
        return(dynamic_cast<Metadata&>(md));
    }catch(...){throw;};
}

GCLgrid::GCLgrid(const string fname, const string format)
{
    const string base_error("GCLgrid file-based constructor:  ");
    if(format==default_output_format)
    {
        try{
            Metadata params=pfload_GCLmetadata(fname);
            /* This would need to be a private method if
               attributes were not public.  Warning if interface
               is changed.*/
            pwmig::gclgrid::pfload_common_GCL_attributes<pwmig::gclgrid::GCLgrid>(*this,params);
            /* Intentionally do not check for object_type to allow
             field constructors to use this */
           //string otype=params.get_string("object_type");
            /* This and other similar routines need to get info about
               the byte order of the data.  Done with this attribute.*/
            string datatype=params.get_string("datatype");
            if( !((datatype=="u8") || (datatype=="t8") ) )
                throw GCLgridError(base_error
                        + "Do not know how to handle datatype="
                        + datatype
                        +"\nMust be u8 or t8");
            bool little_endian=IntelByteOrder();
            bool need_to_swap_bytes;
            if( (datatype=="t8") && little_endian)
                need_to_swap_bytes=true;
            else if( (datatype=="u8") && !little_endian)
                need_to_swap_bytes=true;
            else
                need_to_swap_bytes=false;
            string dfile=fname+"."+dfileext;
            FILE *fp = fopen(dfile.c_str(),"r");
            if(fp == NULL)
                throw GCLgridError(base_error
                        + "fopen failed on file "
                        + dfile);
            /* These are part of the object*/
            x1 = create_2dgrid_contiguous(n1,n2);
            x2 = create_2dgrid_contiguous(n1,n2);
            x3 = create_2dgrid_contiguous(n1,n2);
            /* Database stores data in geographic coordinates.
               File stores Cartesian form.  A bit inconsistent,
               but a design choice.  Storing geo coordinates
               would be a future format choice.*/
            int gridsize = n1*n2;
            if(fread(x1[0],sizeof(double),gridsize,fp) != gridsize)
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fread failed on file reading x1 coordinate  array"
                        + dfile);
            }
            if(fread(x2[0],sizeof(double),gridsize,fp) != gridsize)
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fread failed on file reading x2 coordinate  array"
                        + dfile);
            }
            if(fread(x3[0],sizeof(double),gridsize,fp) != gridsize)
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fread failed on file reading x3 coordinate array"
                        + dfile);
            }
            fclose(fp);
            if(need_to_swap_bytes)
            {
                swapdvec(x1[0],gridsize);
                swapdvec(x2[0],gridsize);
                swapdvec(x3[0],gridsize);
            }
            /* database version calls set_transformation_matrix()
               method here, but not needed because we called
               it earlier */
        } catch(MsPASSError& serr)
        {
            /* Translate message to common error for this package.*/
            throw GCLgridError(base_error
                    + "Problems loading attributes from pf\nDetail:  "
                    + serr.what());
        }
    }
    else
        throw GCLgridError(base_error
                + "unknown input format with tag="
                + format
                + "\nCurrently only accept format="
                +default_output_format);
}
GCLgrid::GCLgrid(const Metadata& md)
{
  try{
    pfload_common_GCL_attributes<GCLgrid>(*this,md);
    read_GCL2d_coord_arrays<GCLgrid>(*this,md);
  }catch(...){throw;};
}

// copy constructors here use inheritance of the BasicGCLgrid
// to reduce redundant code.
GCLgrid::GCLgrid(const GCLgrid& g)
	: BasicGCLgrid(dynamic_cast<const BasicGCLgrid&>(g))
{
	int i,j;
	ix1=g.ix1;
	ix2=g.ix2;
	x1=create_2dgrid_contiguous(n1,n2);
	x2=create_2dgrid_contiguous(n1,n2);
	x3=create_2dgrid_contiguous(n1,n2);
	//
	//I use separate loops for each array here as this is highly
	//optimized on most compilers
	//
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) x1[i][j]=g.x1[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) x2[i][j]=g.x2[i][j];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j) x3[i][j]=g.x3[i][j];

}
GCLgrid3d::GCLgrid3d(const string fname, const string format,const bool fl)
{
    /* This code is painfully similar to GCLgrid constructor
       with same arguments */
    const string base_error("GCLgrid3d file-based constructor:  ");
    fast_lookup=fl;
    if(format==default_output_format)
    {
        try{
            Metadata params=pfload_GCLmetadata(fname);
            /* This would need to be a private method if
               attributes were not public.  Warning if interface
               is changed.*/
            pfload_common_GCL_attributes<GCLgrid3d>(*this,params);
            pfload_3dgrid_attributes<GCLgrid3d>(*this,params);
            /* This and other similar routines need to get info about
               the byte order of the data.  Done with this attribute.*/
            string datatype=params.get_string("datatype");
            if( !((datatype=="u8") || (datatype=="t8") ) )
                throw GCLgridError(base_error
                        + "Do not know how to handle datatype="
                        + datatype
                        +"\nMust be u8 or t8");
            bool little_endian=IntelByteOrder();
            bool need_to_swap_bytes;
            if( (datatype=="t8") && little_endian)
                need_to_swap_bytes=true;
            else if( (datatype=="u8") && !little_endian)
                need_to_swap_bytes=true;
            else
                need_to_swap_bytes=false;
            string dfile=fname+"."+dfileext;
            FILE *fp = fopen(dfile.c_str(),"r");
            if(fp == NULL)
                throw GCLgridError(base_error
                        + "fopen failed on file "
                        + dfile);
            /* These are part of the object.  The file format
             stored the data in the Cartesian system for good
             or bad.  Note this is different from the db data
             where the storage is geographic`*/
	    x1=create_3dgrid_contiguous(n1,n2,n3);
	    x2=create_3dgrid_contiguous(n1,n2,n3);
	    x3=create_3dgrid_contiguous(n1,n2,n3);
            int gridsize = n1*n2*n3;
            if(fread(x1[0][0],sizeof(double),gridsize,fp) != gridsize)
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fread failed on file reading x1 coordinate array"
                        + dfile);
            }
            if(fread(x2[0][0],sizeof(double),gridsize,fp) != gridsize)
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fread failed on file reading x2 coordinate array"
                        + dfile);
            }
            if(fread(x3[0][0],sizeof(double),gridsize,fp) != gridsize)
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fread failed on file reading x3 coordinate array"
                        + dfile);
            }
            fclose(fp);
            if(need_to_swap_bytes)
            {
                swapdvec(x1[0][0],gridsize);
                swapdvec(x2[0][0],gridsize);
                swapdvec(x3[0][0],gridsize);
            }
        } catch(MsPASSError& serr)
        {
            /* Translate message to common error for this package.*/
            throw GCLgridError(base_error
                    + "Problems loading attributes from pf\nDetail:  "
                    + serr.what());
        }
    }
    else
        throw GCLgridError(base_error
                + "unknown input format with tag="
                + format
                + "\nCurrently only accept format="
                +default_output_format);
}
GCLgrid3d::GCLgrid3d(const Metadata& md)
{
  try{
    pfload_common_GCL_attributes<GCLgrid3d>(*this,md);
    pfload_3dgrid_attributes<GCLgrid3d>(*this,md);
    read_GCL3d_coord_arrays<GCLgrid3d>(*this,md);
  }catch(...){throw;};
}
GCLgrid3d::GCLgrid3d(const GCLgrid3d& g)
	: BasicGCLgrid(dynamic_cast<const BasicGCLgrid&>(g))
{
	int i,j,k;
	n3=g.n3;
	dx3_nom=g.dx3_nom;
	k0=g.k0;
	ix1=g.ix1; ix2=g.ix2; ix3=g.ix3;
	x1=create_3dgrid_contiguous(n1,n2,n3);
	x2=create_3dgrid_contiguous(n1,n2,n3);
	x3=create_3dgrid_contiguous(n1,n2,n3);
	// could use an memcpy call here and it might be faster
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k) x1[i][j][k]=g.x1[i][j][k];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k) x2[i][j][k]=g.x2[i][j][k];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k) x3[i][j][k]=g.x3[i][j][k];

  fast_lookup=g.fast_lookup;
}

/* This small function builds a baseline vector of latitude and
longitude values along a specified azimuth with a prescribed
origin.

Arguments:
	lat0, lon0 - origin in radians
	phi - angle to rotate the x axis (+ east) baseline.  That is
		the baseline is established as a great circle passing
		through the origin point with an angle phi measured positive
		northward from east.  (counterclock looking from above)
		Note that away from the origin this angle relative to north
		will change.
	dx - distance spacing between points on the baseline (in radians)
	n - number of points along the baseline.
	i0 - position of the origin point along the baseline vector.
		NOTE:  we use C convention so 0 is the first point in
		the vector.
	lat, lon - hold n vectors of latitude, longitude pairs of the
		baseline points.  Note that lat[i0], lon[i0] should
		equal lat0, lon0 respectively.
Author:  G Pavlis
Written:  Aug 2000
*/

void build_baseline(const double lat0, const double lon0, const double phi,
			const double dx,const int n,const int i0,
			  double *lat, double *lon)
{
	double azimuth;
	int i;
	double delta;

	/* We want azimuth between 0 and 2*pi */
	azimuth = M_PI_2 - phi;
	if(azimuth < 0.0 ) azimuth += (2.0*M_PI);

	for(i=0,delta=dx*((double)-i0);i<n;++i,delta+=dx)
	{
		if(i<i0)
			latlon(lat0,lon0,fabs(delta),azimuth-M_PI,lat+i,lon+i);
		else if(i==i0)
		{
			lat[i] = lat0;
			lon[i] = lon0;
		}
		else
		{
			latlon(lat0,lon0,delta,azimuth,lat+i,lon+i);
		}
	}
}
//
// This is the C++ constructor that uses a pf and coord routines to compute
// a GCLgrid on the fly using the description coming from pf
//
GCLgrid::GCLgrid(const int n1in, const int n2in,
	const string namein,
	const double lat0in, const double lon0in, const double r0in,
	const double azimuth_yin, const double dx1_nomin, const double dx2_nomin,
	const int i0in, const int j0in)

{

	/* pole to baseline */
	double pole_lat, pole_lon;
	int i,j;
	double delta;
	double rotation_angle,azimuth_i;
	double dx1_rad,dx2_rad;
	// temporary grids to hold lat-lon-r
	double **plat, **plon, **pr;

	name=namein;
	lat0=lat0in;
	lon0=lon0in;
	if(r0in>0.0)
		r0=r0in;
	else
		r0 = r0_ellipse(lat0);
	azimuth_y=azimuth_yin;
	rotation_angle=-azimuth_y;
	dx1_nom=dx1_nomin;
	dx2_nom=dx2_nomin;
	dx1_rad = dx1_nom/r0;
	dx2_rad = dx2_nom/r0;
	n1=n1in;
	n2=n2in;
	i0=i0in;
	j0=j0in;

	/* These are the vector of lat and lons used to define the
	baseline the effectively forms the equator of of the rotated
	coordinates passing through the origin*/
	double *baseline_lat=new double[n1in];
	double *baseline_lon=new double[n1in];

	build_baseline(lat0,lon0,
		rotation_angle,dx1_rad,n1,i0,baseline_lat,baseline_lon);
	/* We need to compute the pole to our base line to use as a
	target for grid lines that are locally perpendicular along
	the grid lines -- like longitude lines at the equator */
	latlon(lat0,lon0,M_PI_2,-rotation_angle,&pole_lat,&pole_lon);

	/* We now allocate memory for the large coordinate arrays
	themselves and set the elements of the structure */
	plat = create_2dgrid_contiguous(n1,n2);
	plon = create_2dgrid_contiguous(n1,n2);
	pr = create_2dgrid_contiguous(n1,n2);
	x1 = create_2dgrid_contiguous(n1,n2);
	x2 = create_2dgrid_contiguous(n1,n2);
	x3 = create_2dgrid_contiguous(n1,n2);

	/* We now complete the latitude/longitude grid by projecting
	lines from the baseline toward the computed pole.  Grid points
	will be equally spaced along the gcp toward the pole, but the
	lines they form will converge.  This is complicated greatly by
	the grid origin parameters.  We could make this a function,
	but it is simpler to just write it inline.*/
	for(i=0;i<n1;++i)
	{
		dist(baseline_lat[i],baseline_lon[i],pole_lat, pole_lon,
			&delta,&azimuth_i);
		if(azimuth_i < 0.0) azimuth_i += (2.0*M_PI);
		delta = dx2_rad*(-(double)j0);
		for(j=0;j<n2;++j,delta+=dx2_rad)
		{
			if(j<j0)
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i-M_PI,
					plat[i]+j,plon[i]+j);
			else if(j==j0)
			{
				plat[i][j] = baseline_lat[i];
				plon[i][j] = baseline_lon[i];
			}
			else
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i,
					plat[i]+j,plon[i]+j);
			pr[i][j] = r0_ellipse(plat[i][j]);
		}
	}
	//
	// This function standarizes setting the tranformation from geo to cartesian
	// This then sets the Cartesian section
	//
	set_transformation_matrix();
	Cartesian_point cpt;
	for(i=0;i<n1;++i)
	    for(j=0;j<n2;++j)
	    {
			cpt = gtoc(plat[i][j],plon[i][j],
					pr[i][j]);

			x1[i][j] = cpt.x1;
			x2[i][j] = cpt.x2;
			x3[i][j] = cpt.x3;
	    }
	/* We have to compute the extents parameters as the minimum
	and maximum in each cartesian direction */
	this->compute_extents();

	delete [] baseline_lat;
	delete [] baseline_lon;
	free_2dgrid_contiguous(plat,n1);
	free_2dgrid_contiguous(plon,n1);
	free_2dgrid_contiguous(pr,n1);
}
//
// This is the C++ constructor that uses a pf and coord routines to compute
// a GCLgrid on the fly using the description coming from pf
//
GCLgrid3d::GCLgrid3d(const int n1in, const int n2in, const int n3in,
	const string namein,
	  const double lat0in, const double lon0in, const double r0in,
	    const double azimuth_yin, const double dx1_nomin, const double dx2_nomin,
        const double dx3_nomin, const int i0in, const int j0in)
{
        /* Fast lookup method is always enabled for a grid created this
           way.   Reason is this grid is alway very close to regular
           and the speed loss for lookup is ill advised. */
  fast_lookup=false;
	/* pole to baseline */
	double pole_lat, pole_lon;
	int i,j,k;
	double delta;
	double z0,z;
	double dx1_rad,dx2_rad;
	double rotation_angle, azimuth_i;
	double ***plat, ***plon, ***pr;

	name=namein;
	lat0=lat0in;
	lon0=lon0in;
	if(r0in>0.0)
		r0=r0in;
	else
		r0 = r0_ellipse(lat0);
	azimuth_y=azimuth_yin;
	rotation_angle=-azimuth_y;
	dx1_nom=dx1_nomin;
	dx2_nom=dx2_nomin;
	dx3_nom=dx3_nomin;
	dx1_rad = dx1_nom/r0;
	dx2_rad = dx2_nom/r0;
	n1=n1in;
	n2=n2in;
	n3=n3in;
	i0=i0in;
	j0=j0in;
//
// We force the origin of the r coordinate to be at the bottom of the grid
//
	k0=0;
	z0 = ((double)(n3-1))*dx3_nom;
	/* These are the vector of lat and lons used to define the
	baseline the effectively forms the equator of of the rotated
	coordinates passing through the origin*/
	double *baseline_lat=new double[n1in];
	double *baseline_lon=new double[n1in];

	build_baseline(lat0,lon0,
		rotation_angle,dx1_rad,n1,i0,baseline_lat,baseline_lon);
	/* We need to compute the pole to our base line to use as a
	target for grid lines that are locally perpendicular along
	the grid lines -- like longitude lines at the equator */
	latlon(lat0,lon0,M_PI_2,-rotation_angle,&pole_lat,&pole_lon);

	/* We now allocate memory for the large coordinate arrays
	themselves and set the elements of the structure */
	plat = create_3dgrid_contiguous(n1,n2,n3);
	plon = create_3dgrid_contiguous(n1,n2,n3);
	pr = create_3dgrid_contiguous(n1,n2,n3);
	x1 = create_3dgrid_contiguous(n1,n2,n3);
	x2 = create_3dgrid_contiguous(n1,n2,n3);
	x3 = create_3dgrid_contiguous(n1,n2,n3);

	r0 = r0-z0;
	azimuth_y = -rotation_angle;

	/* We now complete the latitude/longitude grid by projecting
	lines from the baseline toward the computed pole.  Grid points
	will be equally spaced along the gcp toward the pole, but the
	lines they form will converge.  This is complicated greatly by
	the grid origin parameters.  We could make this a function,
	but it is simpler to just write it inline. Because this function
	descended from the a C program call makegclgrid, which built both
	a 2d and 3d grid in one pass, I'm stealing that algorithm.  That is
	the top of the grid (at n3-1) is built first, then the remainder is
	filled in with same lat lons but with varying depth.  So, first this
	loop creates the top sheet of the grid */
	int top_of_grid=n3-1;
	for(i=0;i<n1;++i)
	{
		dist(baseline_lat[i],baseline_lon[i],pole_lat, pole_lon,
			&delta,&azimuth_i);
		if(azimuth_i < 0.0) azimuth_i += (2.0*M_PI);
		delta = dx2_rad*(-(double)j0);
		for(j=0;j<n2;++j,delta+=dx2_rad)
		{
			if(j<j0)
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i-M_PI,
					&(plat[i][j][top_of_grid]),&(plon[i][j][top_of_grid]));
			else if(j==j0)
			{
				plat[i][j][top_of_grid] = baseline_lat[i];
				plon[i][j][top_of_grid] = baseline_lon[i];
			}
			else
				latlon(baseline_lat[i],baseline_lon[i],
					fabs(delta),azimuth_i,
					&(plat[i][j][top_of_grid]),&(plon[i][j][top_of_grid]));
			pr[i][j][top_of_grid] = r0_ellipse(plat[i][j][top_of_grid]);
		}
	}
//
//Now we fill in the rest of the grid keeping the lat lon of the top for all
//levels but varying r -- i.e. these are spherical shells.
//
	for(k=0,z=z0;k<n3-1;++k,z-=dx3_nom)
	{
		for(i=0;i<n1;++i)
			for(j=0;j<n2;++j)
			{
				plat[i][j][k]=plat[i][j][top_of_grid];
				plon[i][j][k]=plon[i][j][top_of_grid];
				pr[i][j][k] = r0_ellipse(plat[i][j][top_of_grid]) - z;
			}
	}
	//
	// This function standarizes setting the tranformation from geo to cartesian
	//
	set_transformation_matrix();
	Cartesian_point cpt;
	for(i=0;i<n1;++i)
	    for(j=0;j<n2;++j)
		for(k=0;k<n3;++k)
		{
			cpt = gtoc(plat[i][j][k],plon[i][j][k],
					pr[i][j][k]);

			x1[i][j][k] = cpt.x1;
			x2[i][j][k] = cpt.x2;
			x3[i][j][k] = cpt.x3;
		}
	/* We have to compute the extents parameters as the minimum
	and maximum in each cartesian direction */
	this->compute_extents();

	delete [] baseline_lat;
	delete [] baseline_lon;
	free_3dgrid_contiguous(plat,n1,n2);
	free_3dgrid_contiguous(plon,n1,n2);
	free_3dgrid_contiguous(pr,n1,n2);
}
//
// constructors for derived types for scalar and vector fields defined on
// a gcl grid.  Makes extensive use of inheritance.
//
GCLscalarfield::GCLscalarfield() : GCLgrid()
{
	val=NULL;
}

GCLscalarfield::GCLscalarfield(const int n1size, const int n2size)
	: GCLgrid(n1size, n2size)
{
	val=create_2dgrid_contiguous(n1size, n2size);
}
GCLvectorfield::GCLvectorfield() : GCLgrid()
{
	val=NULL;
}

GCLvectorfield::GCLvectorfield(const int n1size, const int n2size, const int n3size)
	: GCLgrid(n1size, n2size)
{
	nv=n3size;
	val=create_3dgrid_contiguous(n1size, n2size, n3size);
}
GCLscalarfield::GCLscalarfield(const GCLscalarfield& g)
	 : GCLgrid(dynamic_cast<const GCLgrid&>(g))
{
	int i,j;
	val = create_2dgrid_contiguous(g.n1, g.n2);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			val[i][j]=g.val[i][j];
}
GCLscalarfield::GCLscalarfield(const GCLgrid& g)
	: GCLgrid(g)
{
	// Used to copy scalar attributes and grid data here.
	// Use of call to grid constructor removes this need
	int i,j;
	val=create_2dgrid_contiguous(g.n1, g.n2);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)val[i][j]=0.0;
}
GCLvectorfield::GCLvectorfield(const GCLvectorfield& g)
	: GCLgrid(dynamic_cast<const GCLgrid&>(g))
{
	int i,j,k;
	nv = g.nv;
	val = create_3dgrid_contiguous(g.n1, g.n2, nv);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<nv;++k)
				val[i][j][k]=g.val[i][j][k];
}
GCLvectorfield::GCLvectorfield(const GCLgrid& g, const int n3) : GCLgrid(g)
{
	int i,j,k;
	nv=n3;
	val=create_3dgrid_contiguous(g.n1, g.n2,nv);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<nv;++k)val[i][j][k]=0.0;
}
//
//3d versions
//
GCLscalarfield3d::GCLscalarfield3d() : GCLgrid3d()
{
	val=NULL;
}

GCLscalarfield3d::GCLscalarfield3d(const int n1size, const int n2size, const int n3size)
	: GCLgrid3d(n1size, n2size, n3size)
{
	val=create_3dgrid_contiguous(n1size, n2size, n3size);
}
GCLvectorfield3d::GCLvectorfield3d() : GCLgrid3d()
{
	val=NULL;
}

GCLvectorfield3d::GCLvectorfield3d(const int n1size, const int n2size,
  const int n3size, const int n4) : GCLgrid3d(n1size, n2size, n3size)
{
	nv=n4;
	val=create_4dgrid_contiguous(n1size, n2size, n3size, n4);
}
GCLscalarfield3d::GCLscalarfield3d(const GCLscalarfield3d& g)
	: GCLgrid3d(dynamic_cast<const GCLgrid3d&>(g))
{
	int i,j,k;
	val = create_3dgrid_contiguous(n1,n2,n3);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k)
				val[i][j][k]=g.val[i][j][k];
}
GCLscalarfield3d::GCLscalarfield3d(const GCLgrid3d& g)
	: GCLgrid3d(g)
{
	int i,j,k;
	val=create_3dgrid_contiguous(g.n1, g.n2, g.n3);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k)
				val[i][j][k]=0.0;
}
GCLvectorfield3d::GCLvectorfield3d(const GCLvectorfield3d& g)
	: GCLgrid3d(dynamic_cast<const GCLgrid3d&>(g))
{
	int i,j,k,l;
	nv = g.nv;
	val = create_4dgrid_contiguous(n1,n2,n3,nv);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k)
				for(l=0;l<nv;++l)
					val[i][j][k][l]=g.val[i][j][k][l];
}
GCLvectorfield3d::GCLvectorfield3d(const GCLgrid3d& g, const int n4)
	: GCLgrid3d(g)
{
	int i,j,k,l;
	nv=n4;
	val=create_4dgrid_contiguous(g.n1, g.n2, g.n3,n4);
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
			for(k=0;k<n3;++k)
				for(l=0;l<nv;++l)
					val[i][j][k][l]=0.0;
}
/* File based constructors for field objects */
GCLscalarfield::GCLscalarfield(const string fname, const string format,
        const bool enforce_object_type)
    : GCLgrid(fname,format)
{
    const string base_error("GCLscalarfield file constructor:  ");
    try {
        if(format==default_output_format)
        {
            /* Necesary complication to recreate the Metadata object.
               Inefficent but assumption is large number of these
               objects are not going to be created by this mechanism. */
            Metadata params=pfload_GCLmetadata(fname);
            if(enforce_object_type)
            {
            string otype=params.get_string("object_type");
            string otypethis=string(typeid(*this).name());
            if(otype  != boost::core::demangle(otypethis.c_str()))
                throw GCLgridError(base_error
               + "Object type mismatch.  "
               + "Called GCLscalarfield constructor on a file with object_type="
               + otype);
            }
            string dfile=fname+"."+dfileext;
            FILE *fp=fopen(dfile.c_str(),"r");
            if(fp==NULL)
                throw GCLgridError(base_error
                        + "fopen failed for file="
                        + dfile);
            long fvalstart;
            // 3 is because there are 3 coordinates for each grid point
            fvalstart=sizeof(double)*n1*n2*3;
            if(fseek(fp,fvalstart,SEEK_SET))
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fseek to start of field data area of file failed");
            }
            size_t npts=n1*n2;
            val=create_2dgrid_contiguous(n1,n2);
            size_t nvread;
            nvread=fread((void *)(&(val[0][0])),sizeof(double),npts,fp);
            fclose(fp);
            if(nvread!=npts)
                throw GCLgridError(base_error
                        + "fread error reading field data area of file"
                        + dfile);
        }
        else
        {
            throw GCLgridError(base_error
                    + "Do not know how to handle format with name="
                    +format);
        }

    }catch(...){throw;};
}
GCLscalarfield3d::GCLscalarfield3d(const string fname, const string format)
    : GCLgrid3d(fname,format)
{
    const string base_error("GCLscalarfield3d file constructor:  ");
    try {
        if(format==default_output_format)
        {
            Metadata params=pfload_GCLmetadata(fname);
            string otype=params.get_string("object_type");
            string otypethis=string(typeid(*this).name());
            if(otype  != boost::core::demangle(otypethis.c_str()))
                throw GCLgridError(base_error
               + "Object type mismatch.  "
               + "Called GCLscalarfield3d constructor on a file with object_type="
               + otype);
            string dfile=fname+"."+dfileext;
            FILE *fp=fopen(dfile.c_str(),"r");
            if(fp==NULL)
                throw GCLgridError(base_error
                        + "fopen failed for file="
                        + dfile);
            long fvalstart;
            fvalstart=sizeof(double)*n1*n2*n3*3;
            if(fseek(fp,fvalstart,SEEK_SET))
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fseek to start of field data area of file failed");
            }
            size_t npts=n1*n2*n3;
            val=create_3dgrid_contiguous(n1,n2,n3);
            size_t nvread;
            nvread=fread((void *)(&(val[0][0][0])),sizeof(double),npts,fp);
            fclose(fp);
            if(nvread!=npts)
                throw GCLgridError(base_error
                        + "fread error reading field data area of file"
                        + dfile);
        }
        else
        {
            throw GCLgridError(base_error
                    + "Do not know how to handle format with name="
                    +format);
        }

    }catch(...){throw;};
}
GCLvectorfield::GCLvectorfield(const string fname, const string format,
        const bool enforce_object_type)
    : GCLgrid(fname,format)
{
    const string base_error("GCLvectorfield file constructor:  ");
    try {
        if(format==default_output_format)
        {
            Metadata params=pfload_GCLmetadata(fname);
            if(enforce_object_type)
            {
            string otype=params.get_string("object_type");
            string otypethis=string(typeid(*this).name());
            if(otype  != boost::core::demangle(otypethis.c_str()))
                throw GCLgridError(base_error
               + "Object type mismatch.  "
               + "Called GCLvectorfield constructor on a file with object_type="
               + otype);
            }
            /* nv has to be handled specially.  attribute name maintenance
               issue here */
            nv=params.get_int("nv");
            string dfile=fname+"."+dfileext;
            FILE *fp=fopen(dfile.c_str(),"r");
            if(fp==NULL)
                throw GCLgridError(base_error
                        + "fopen failed for file="
                        + dfile);
            long fvalstart;
            fvalstart=sizeof(double)*n1*n2*3;
            if(fseek(fp,fvalstart,SEEK_SET))
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fseek to start of field data area of file failed");
            }
            size_t npts=n1*n2*nv;
            val=create_3dgrid_contiguous(n1,n2,nv);
            size_t nvread;
            nvread=fread((void *)(&(val[0][0][0])),sizeof(double),npts,fp);
            fclose(fp);
            if(nvread!=npts)
                throw GCLgridError(base_error
                        + "fread error reading field data area of file"
                        + dfile);
        }
        else
        {
            throw GCLgridError(base_error
                    + "Do not know how to handle format with name="
                    +format);
        }
    }catch(...){throw;};
}
GCLvectorfield3d::GCLvectorfield3d(const string fname, const string format)
    : GCLgrid3d(fname,format)
{
    const string base_error("GCLvectorfield3d file constructor:  ");
    try {
        if(format==default_output_format)
        {
            Metadata params=pfload_GCLmetadata(fname);
            string otype=params.get_string("object_type");
            string otypethis=string(typeid(*this).name());
            if(otype  != boost::core::demangle(otypethis.c_str()))
                throw GCLgridError(base_error
               + "Object type mismatch.  "
               + "Called GCLvectorfield3d constructor on a file with object_type="
               + otype);
            /* nv has to be handled specially.  attribute name maintenance
               issue here */
            nv=params.get_int("nv");
            string dfile=fname+"."+dfileext;
            FILE *fp=fopen(dfile.c_str(),"r");
            if(fp==NULL)
                throw GCLgridError(base_error
                        + "fopen failed for file="
                        + dfile);
            long fvalstart;
            fvalstart=sizeof(double)*n1*n2*n3*3;
            if(fseek(fp,fvalstart,SEEK_SET))
            {
                fclose(fp);
                throw GCLgridError(base_error
                        + "fseek to start of field data area of file failed");
            }
            size_t npts=n1*n2*n3*nv;
            val=create_4dgrid_contiguous(n1,n2,n3,nv);
            size_t nvread;
            nvread=fread((void *)(&(val[0][0][0][0])),sizeof(double),npts,fp);
            fclose(fp);
            if(nvread!=npts)
                throw GCLgridError(base_error
                        + "fread error reading field data area of file"
                        + dfile);
        }
        else
        {
            throw GCLgridError(base_error
                    + "Do not know how to handle format with name="
                    +format);
        }

    }catch(...){throw;};
}
GCLscalarfield::GCLscalarfield(const Metadata& md) : GCLgrid(md)
{
  try{
    const string base_error("GCLscalarfield Metadata contructor:  ");
    if(this->n1 <= 0 || this->n2 <=0)
    {
      stringstream ss;
      ss << base_error
         << "Illegal dimensions defined in Metadata received"
         << "Parsed n1="<<this->n1<<" n2="<<this->n2<<endl
         << "All must be nonzero, positive numbers"<<endl;
      throw GCLgridError(ss.str());
    }
    this->val=(double **)create_2dgrid_contiguous(this->n1,this->n2);
    double *buffer=(*this->val);
    size_t bufsize=(this->n1)*(this->n2);
    read_fielddata<GCLscalarfield>(md,buffer,bufsize);
  }catch(...){throw;};
}
GCLscalarfield3d::GCLscalarfield3d(const Metadata& md) : GCLgrid3d(md)
{
  try{
    const string base_error("GCLscalarfield3d Metadata contructor:  ");
    if(this->n1 <= 0 || this->n2 <=0 || this->n3 <=0)
    {
      stringstream ss;
      ss << base_error
         << "Illegal dimensions defined in Metadata received"
         << "Parsed n1="<<this->n1<<" n2="<<this->n2<<" n3="<<this->n3<<endl
         << "All must be nonzero, positive numbers"<<endl;
      throw GCLgridError(ss.str());
    }
    this->val=(double ***)create_3dgrid_contiguous(this->n1,this->n2,this->n3);
    double *buffer=(**this->val);
    size_t bufsize=(this->n1)*(this->n2)*(this->n3);
    read_fielddata<GCLscalarfield3d>(md,buffer,bufsize);
  }catch(...){throw;};
}
GCLvectorfield::GCLvectorfield(const Metadata& md) : GCLgrid(md)
{
  try{
    const string base_error("GCLvectorfield Metadata contructor:  ");
    this->nv = md.get_int("nv");   // frozen name
    if(this->n1 <= 0 || this->n2 <=0 || this->nv <=0)
    {
      stringstream ss;
      ss << base_error
         << "Illegal dimensions defined in Metadata received"
         << "Parsed n1="<<this->n1<<" n2="<<this->n2<<" vector size nv="<<this->nv<<endl
         << "All must be nonzero, positive numbers"<<endl;
      throw GCLgridError(ss.str());
    }
    this->val=(double ***)create_3dgrid_contiguous(this->n1,this->n2,this->nv);
    double *buffer=(**this->val);
    size_t bufsize=(this->n1)*(this->n2)*(this->nv);
    read_fielddata<GCLvectorfield>(md,buffer,bufsize);
  }catch(...){throw;};
}
GCLvectorfield3d::GCLvectorfield3d(const Metadata& md) : GCLgrid3d(md)
{
  try{
    const string base_error("GCLvectorfield3d Metadata contructor:  ");
    this->nv = md.get_int("nv");   // frozen name
    if(this->n1 <= 0 || this->n2 <= 0 || this->n3 <= 0 || this->nv <= 0)
    {
      stringstream ss;
      ss << base_error
         << "Illegal dimensions defined in Metadata received"
         << "Parsed n1="<<this->n1<<" n2="<<this->n2<<" n3="<<this->n3
         << " vector size nv="<<this->nv<<endl
         << "All must be nonzero, positive numbers"<<endl;
      throw GCLgridError(ss.str());
    }
    this->val=(double ****)create_4dgrid_contiguous(this->n1,this->n2,this->n3,this->nv);
    double *buffer=(***this->val);
    size_t bufsize=(this->n1)*(this->n2)*(this->n3)*(this->nv);
    read_fielddata<GCLvectorfield3d>(md,buffer,bufsize);
  }catch(...){throw;};
}

//
//C++ destructors
//
GCLgrid::~GCLgrid()
{
	if(x1!=NULL) free_2dgrid_contiguous(x1,n1);
	if(x2!=NULL) free_2dgrid_contiguous(x2,n1);
	if(x3!=NULL) free_2dgrid_contiguous(x3,n1);
}
GCLgrid3d::~GCLgrid3d()
{
	if(x1!=NULL) free_3dgrid_contiguous(x1,n1,n2);
	if(x2!=NULL) free_3dgrid_contiguous(x2,n1,n2);
	if(x3!=NULL) free_3dgrid_contiguous(x3,n1,n2);
}
//
//Note standard rule that a derived class utilizes the base class destructor
//after calling these
//
GCLscalarfield::~GCLscalarfield()
{
	if(val!=NULL) free_2dgrid_contiguous(val,n1);
}
GCLvectorfield::~GCLvectorfield()
{
	if(val!=NULL) free_3dgrid_contiguous(val,n1,n2);
}
GCLscalarfield3d::~GCLscalarfield3d()
{
	if(val!=NULL) free_3dgrid_contiguous(val,n1,n2);
}
GCLvectorfield3d::~GCLvectorfield3d()
{
	if(val!=NULL) free_4dgrid_contiguous(val,n1,n2,n3);
}
void GCLgrid::compute_extents()
{
	int i,j;
	double extents_dx;  // this becomes max of grid sizes
	extents_dx=MAX(dx1_nom,dx2_nom);
	x1low=x1[0][0];
	x1high=x1[0][0];
	x2low=x2[0][0];
	x2high=x2[0][0];
	x3low=x3[0][0];
	x3high=x3[0][0];
	for(i=0;i<n1;++i)
		for(j=0;j<n2;++j)
		{
			x1low = MIN(x1[i][j],x1low);
			x1high = MAX(x1[i][j],x1high);
			x2low = MIN(x2[i][j],x2low);
			x2high = MAX(x2[i][j],x2high);
			x3low = MIN(x3[i][j],x3low);
			x3high = MAX(x3[i][j],x3high);
		}
	x1low -= extents_dx;
	x2low -= extents_dx;
	x3low -= extents_dx;
	x1high += extents_dx;
	x2high += extents_dx;
	x3high += extents_dx;
}
void GCLgrid3d::compute_extents()
{
	int i,j,k;
	double extents_dx;  // this becomes max of grid sizes
	extents_dx=MAX(dx1_nom,dx2_nom);
	extents_dx=MAX(extents_dx,dx3_nom);
	x1low=x1[0][0][0];
	x1high=x1[0][0][0];
	x2low=x2[0][0][0];
	x2high=x2[0][0][0];
	x3low=x3[0][0][0];
	x3high=x3[0][0][0];

	for(i=0;i<n1;++i)
	    for(j=0;j<n2;++j)
		for(k=0;k<n3;++k)
		{
			x1low = MIN(x1[i][j][k],x1low);
			x1high = MAX(x1[i][j][k],x1high);
			x2low = MIN(x2[i][j][k],x2low);
			x2high = MAX(x2[i][j][k],x2high);
			x3low = MIN(x3[i][j][k],x3low);
			x3high = MAX(x3[i][j][k],x3high);
		}
	x1low -= extents_dx;
	x2low -= extents_dx;
	x3low -= extents_dx;
	x1high += extents_dx;
	x2high += extents_dx;
	x3high += extents_dx;

}
} //end namespace
