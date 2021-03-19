#include <stdio.h>
#include <math.h>
#include "mspass/utility/dmatrix.h"
#include "pwmig/gclgrid/gclgrid.h"
#include "misc/blas.h"
namespace pwmig::gclgrid
{
using namespace std;
using mspass::utility::dmatrix;
void BasicGCLgrid::set_transformation_matrix()
{
	double pole_lat, pole_lon;
	double x0[3], xpole[3], xcros[3];

	latlon(lat0, lon0, M_PI_2,azimuth_y, &pole_lat, &pole_lon);

	dsphcar(lon0,lat0,x0);
	dsphcar(pole_lon, pole_lat, xpole);
	dr3cros(xpole, x0, xcros);

	//
	//these unit vectors are stored in rows of transformation matrix
	//yielding the matrix from geo to rotated cartesian of GCLgrid object
	//
	dcopy(3,x0,1,gtoc_rmatrix[2],1);
	dcopy(3,xpole,1,gtoc_rmatrix[1],1);
	dcopy(3,xcros,1,gtoc_rmatrix[0],1);
	dcopy(3,x0,1,translation_vector,1);
	dscal(3,r0,translation_vector,1);
}
dmatrix BasicGCLgrid::fetch_transformation_matrix() const
{
	dmatrix U(3,3);
	int i,j;
	for(i=0;i<3;++i)
		for(j=0;j<3;++j)
			U(i,j) = gtoc_rmatrix[i][j];
	return(U);
}
double *BasicGCLgrid::fetch_translation_vector() const
{
	int i;
	double *t = new double[3];
	for(i=0;i<3;++i) t[i]=translation_vector[i];
	return(t);
}
Cartesian_point BasicGCLgrid::gtoc(const double plat,
	const double plon, const double pr) const
{
	Cartesian_point p;
	double xp[3],dxp[3];
	int i;
	//dsphcar computes a unit vector.  We then scale it by radius
	dsphcar(plon, plat, xp);
	for(i=0;i<3;++i) xp[i]*=pr;
	//
	//Remove the origin translation vector
	//
	/* Weird const cast is necessary because const tag for this method
	collides with the api for d3rsub which is an old C function that does not
	have any const declarations */
	dr3sub(xp,const_cast<double *>(translation_vector),dxp);
	//
	//Rotate by transformation matrix and that is it
	//
	p.x1=ddot(3,gtoc_rmatrix[0],1,dxp,1);
	p.x2=ddot(3,gtoc_rmatrix[1],1,dxp,1);
	p.x3=ddot(3,gtoc_rmatrix[2],1,dxp,1);
	return(p);
}
Cartesian_point BasicGCLgrid::gtoc(const Geographic_point p) const
{
	return(this->gtoc(p.lat,p.lon,p.r));
}
Geographic_point BasicGCLgrid::ctog(const double x1p, const double x2p,
	const double x3p)const
{
	Geographic_point p;
	double dxp[3], x[3];
	int i;

	//
	//First apply rotation matrix to put this point in the
	//standard geographic reference frame.  We use a safe
	//BLAS algorithm for the multiply of the transpose
	//
	for(i=0;i<3;++i) dxp[i]=0.0;
	daxpy(3,x1p,gtoc_rmatrix[0],1,dxp,1);
	daxpy(3,x2p,gtoc_rmatrix[1],1,dxp,1);
	daxpy(3,x3p,gtoc_rmatrix[2],1,dxp,1);
	//
	//dxp now contains the relative location vector in standard
	//coordinates.  We have to apply a translation vector
	//relative to the center of the earth before we can convert
	//to correct latitude and longitude
	//
	for(i=0;i<3;++i) x[i]=dxp[i]+translation_vector[i];
	dcarsph(x,&(p.lon),&(p.lat));
	p.r = dnrm2(3,x,1);
	return(p);
}
Geographic_point BasicGCLgrid::ctog(const Cartesian_point p) const
{
	return(this->ctog(p.x1,p.x2,p.x3));
}
// The following should probably be in a separate file as
// they aren't the BasicGCLgrid type.
Geographic_point GCLgrid::geo_coordinates(const int i, const int j) const
{
	return(ctog(x1[i][j],x2[i][j],x3[i][j]));
}
Geographic_point GCLgrid3d::geo_coordinates(const int i, const int j,
	const int k) const
{
	return(ctog(x1[i][j][k],x2[i][j][k],x3[i][j][k]));
}
double GCLgrid::lat(const int i, const int j) const
{
	Geographic_point p;
	p = ctog(x1[i][j],x2[i][j],x3[i][j]);
	return(p.lat);
}
double GCLgrid::lon(const int i, const int j) const
{
	Geographic_point p;
	p = ctog(x1[i][j],x2[i][j],x3[i][j]);
	return(p.lon);
}
double GCLgrid::r(const int i, const int j) const
{
	Geographic_point p;
	p = ctog(x1[i][j],x2[i][j],x3[i][j]);
	return(p.r);
}
double GCLgrid3d::lat(const int i, const int j, const int k) const
{
	Geographic_point p;
	p = ctog(x1[i][j][k],x2[i][j][k],x3[i][j][k]);
	return(p.lat);
}
double GCLgrid3d::lon(const int i, const int j, const int k) const
{
	Geographic_point p;
	p = ctog(x1[i][j][k],x2[i][j][k],x3[i][j][k]);
	return(p.lon);
}
double GCLgrid3d::r(const int i, const int j, const int k) const
{
	Geographic_point p;
	p = ctog(x1[i][j][k],x2[i][j][k],x3[i][j][k]);
	return(p.r);
}

} // End namespace enscapsulation
