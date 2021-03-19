/*H+
*	Title	:  dr3mag.c
*	Author	:  Jim Lewis (adapted from r3mag.c)
*	Date	:  19 April, 1991
*	Synopsis:  Magnitude of a 3-vector
*	Keywords:  ees, matvec, r3mag, vector
*	Revisions:
*	mm/dd/yy   name		description
*U+
*	Usage	:  double dr3mag(x)
*		   double x[3];
*	Input	:  x:  3-vector
	Output	:  returns magnitude (length) of x
*U-								*/

/* ************************************************************ */


#include <math.h>

double dr3mag(x)
double x[3];

{
	double r3 = 0;
	int i;

	for (i = 0; i < 3; i++)
	   r3 += x[i] * x[i];

	return(sqrt(r3));
}

/* $Id: dr3mag.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
