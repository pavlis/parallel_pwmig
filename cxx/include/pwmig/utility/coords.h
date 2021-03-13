
#ifndef __COORDS__
#define __COORDS__

#include <math.h>

/* Physics Vade Mecum, Herbert Anderson, pp 67 */
#define EQUATORIAL_EARTH_RADIUS 6378.164 

#define deg(r)    ((r) * 180.0/M_PI)
#define rad(d)    ((d) * M_PI/180.0)


/* 3D geometry and Spherical Geometry */

#ifdef	__cplusplus
extern "C" {
#endif

extern double km2deg ( double km );
extern double deg2km ( double degrees );

extern void carsph ( float x[3], float *xlong, float *xlat );
extern void dcarsph ( double x[3], double *xlong, double *xlat );
extern void sphcar ( double lng, double lat, float x[3] );
extern void dsphcar ( double lng, double lat, double x[3] );

extern double ddistance ( double xlat1, double xlon1, double xlat2, double xlon2 );
extern void dist ( double xlat1, double xlong1, double xlat2, double xlong2, double *del, double *az );
extern void latlon ( double lat1, double lon1, double delta, double azi, double *lat2, double *lon2 );
extern void sphdelta ( double xlat1, double xlon1, double xlat2, double xlon2, double * eastdel, double * northdel );
extern void deltasph ( double xlat1, double xlon1, double eastdel, double northdel, double * xlat2, double * xlon2 );

extern double dr3dot ( double x[3], double y[3] );
extern double dr3mag ( double x[3] );
extern void dr3add ( double v1[3], double v2[3], double v3[3] );
extern void dr3cros ( double x[3], double y[3], double z[3] );
extern void dr3mov ( double a[3], double b[3] );
extern void dr3mxm ( double a[9], double b[9], double c[9] );
extern void dr3mxv ( double a[9], double b[3], double c[3] );
extern void dr3norm ( double v[3] );
extern void dr3norm ( double x[3] );
extern void dr3ortho ( double a[3], double b[3] );
extern void dr3sub ( double v1[3], double v2[3], double v3[3] );
extern void dr3sxv ( double a, double b[3], double c[3] );
extern void dr3tran ( double a[9], double b[9] );
extern void drotmat ( double x[3], double y[9], double theta );
extern void xrotate(double x[3],double y[3], double z[3], double theta);

#ifdef	__cplusplus
}
#endif

#endif 

