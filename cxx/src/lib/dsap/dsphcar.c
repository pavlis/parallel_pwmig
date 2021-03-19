#include <math.h>
#include "pwmig/dsap/coords.h"

void dsphcar(lng, lat, x)
double  lng, lat;
double  x[3];
{
	double coslat = cos(lat);
	x[0] = coslat * cos(lng);
	x[1] = coslat * sin(lng);
	x[2] = sin(lat);
}

/* $Id: dsphcar.c,v 1.1.1.1 1997/04/12 04:18:50 danq Exp $ */
