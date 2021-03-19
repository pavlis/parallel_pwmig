#include <stdio.h>
#include <math.h>
#include "pwmig/dsap/coords.h"
/*
#include "pwmig/dsap/stock.h"
*/
#define deg(r)    ((r) * 180.0/M_PI)
#define rad(d)    ((d) * M_PI/180.0)

double ddistance ( xlat1, xlon1, xlat2, xlon2 )
double xlat1, xlon1, xlat2, xlon2 ;
{
    double x1[3], x2[3] ;
    dsphcar(rad(xlon1), rad(xlat1), x1 ) ;
    dsphcar(rad(xlon2), rad(xlat2), x2 ) ;
    return deg(acos(dr3dot(x1, x2))) ;
}

void
dist (xlat1, xlong1, xlat2, xlong2, del, az)

double xlat1, xlong1, xlat2, xlong2;
double *del, *az;

/*
 *
 *   Subroutine dist will compute the great circle distance in radians
 *   and the relative azimuth angle from north between two points
 *   specified as latitude-longitude coordinates.
 *
 *   Inputs  -	xlat1, xlong1	= The latitude-longitude in radians of the
 *			  	  first point. Latitude is positive north
 *			  	  from the equator and longitude is positive
 *			  	  east from Greenwich, England.
 *		xlat2, xlong2	= Same for the second point.
 *
 *   Outputs -	del		= The shortest angular great circle distance
 *			  	  between the two points.
 *		az		= The relative azimuth in radians positive
 *			  	  east from north of point 2 relative to
 *			  	  point 1 as viewed from point 1.
 *
 */

{
      double clat, slat, clat2, slat2, s;
      double x1, z1, x2, y2, z2, x, z, xpp, ypp, zpp;

      slat = sin(xlat1);
      clat = cos(xlat1);
      x1 = clat;
      z1 = slat;

      slat2 = sin(xlat2);
      clat2 = cos(xlat2);
      x2 = clat2*cos(xlong2-xlong1);
      y2 = clat2*sin(xlong2-xlong1);
      z2 = slat2;

      x = x2 - x1;
      z = z2 - z1;
      xpp = x*clat + z*slat;
      ypp = y2;
      zpp = -x*slat + z*clat;

      s = sqrt ( xpp*xpp + ypp*ypp + zpp*zpp );
      *del = 2.0 * asin ( 0.5*s );
      if ( ypp == 0.0 && zpp == 0.0 )
	  *az = 0.0 ;
      else
	  *az = atan2 ( ypp, zpp );
      if (*az < 0.0) *az = *az + 2.0*M_PI;
}


void dist_ ( f_xlat1, f_xlong1, f_xlat2, f_xlong2, del, az )
double *f_xlat1 ;
double *f_xlong1 ;
double *f_xlat2 ;
double *f_xlong2 ;
double *del ;
double *az ;
{
    dist(*f_xlat1, *f_xlong1, *f_xlat2, *f_xlong2, del, az );
}

void latlon_ ( f_xlat1, f_xlong1, f_del, f_az, xlat2, xlong2 )
double *f_xlat1 ;
double *f_xlong1 ;
double *f_del ;
double *f_az ;
double *xlat2 ;
double *xlong2 ;
{
    latlon(*f_xlat1, *f_xlong1, *f_del, *f_az, xlat2, xlong2 );
}


void latlon (lat1,lon1, delta, azi, lat2, lon2)
double lat1, lon1, delta, azi, *lat2, *lon2 ;
{
      double new_rotated[3], axis[3], new[3] ;

      /* The method here is
      * 1) use a coordinate system which puts the north pole
      *    at lat1, lon1.
      * 2) rotate this north pole an angle delta about the axis
      *    defined by the great circle longitude=lon1-azi.
      * 3) rotate back to the original coordinate system.
      */

      if ( fabs(delta) > 0.0 )  {
	  dsphcar (lon1-azi+M_PI, M_PI_2-delta, new_rotated ) ;
	  dsphcar (lon1+M_PI_2, 0.0, axis ) ;
	  xrotate (axis, new_rotated, new, M_PI_2 - lat1) ;
	  dcarsph (new, lon2, lat2) ;
      } else {
	  *lat2 = lat1 ;
	  *lon2 = lon1 ;
      }
}

void
sphdelta (xlat1, xlon1, xlat2, xlon2,
              eastdel, northdel)

double        xlat1, xlon1, xlat2, xlon2;
double *      eastdel;
double *            northdel;

{
      double xlt1, xln1, xlt2, xln2, del, az ;

      if (xlat1 == xlat2 && xlon1 == xlon2) {
        *eastdel = 0.0;
        *northdel = 0.0;
        return;
      } else {
        xlt1 = rad(xlat1) ;
        xln1 = rad(xlon1) ;
        xlt2 = rad(xlat2) ;
        xln2 = rad(xlon2) ;
        dist (xlt1, xln1, xlt2, xln2, &del, &az);
      }
      del = deg(del) ;
      *eastdel = del * sin (az);
      *northdel = del * cos (az);

      return;
}

void
deltasph (xlat1, xlon1, eastdel, northdel,
              xlat2, xlon2)

double        xlat1, xlon1, eastdel, northdel;
double *      xlat2;
double *             xlon2;

{
      double xlt1, xln1, xlt2, xln2, del, az;

      if (eastdel == 0.0 && northdel == 0.0) {
        *xlat2 = xlat1;
        *xlon2 = xlon1;
        return;
      } else {
        xlt1 = rad(xlat1) ;
        xln1 = rad(xlon1) ;
        del = sqrt(eastdel*eastdel + northdel*northdel);
        del = rad(del) ;
	if ( eastdel == 0.0 && northdel == 0.0 )
	    az = 0.0 ;
	else
	    az = atan2 ( eastdel, northdel );
        latlon (xlt1, xln1, del, az, &xlt2, &xln2);
      }
      *xlat2 = deg(xlt2) ;
      *xlon2 = deg(xln2) ;

      return;
}

/* Physics Vade Mecum, Herbert Anderson, pp 67 */
#define EQUATORIAL_EARTH_RADIUS 6378.164

double km2deg ( km )
double km ;
{
    return deg(km/EQUATORIAL_EARTH_RADIUS) ;
}

double deg2km ( degrees )
double degrees ;
{
    return (rad(degrees) * EQUATORIAL_EARTH_RADIUS) ;
}


/* $Id: dist.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
