#include <math.h>
//#include "stock.h"
#include "coords.h"
#define SIGN(d)         (((d)>0)? 1.0 : ((d)<0) ? -1.0 : 0.0)


void   dcarsph(x, xlong, xlat)
double x[3];
double *xlong;
double *xlat;
{
  if ( x[0] == 0.0 && x[1] == 0.0 ) {
    *xlong = 0.0 ; 
    *xlat = SIGN(x[2]) * M_PI_2 ; 
  } else {
      *xlong = atan2(x[1], x[0]);
      *xlat = atan2( x[2], sqrt( (x[0] * x[0] + x[1] * x[1])));
  }
}



/* $Id: dcarsph.c,v 1.2 1997/04/15 21:45:16 danny Exp $ */
