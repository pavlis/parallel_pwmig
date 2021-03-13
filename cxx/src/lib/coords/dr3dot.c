/*H+
*	Title	:  dr3dot.c
*	Author	:  Jim Lewis  (adapted from r3dot.c)
*	Date	:  19 April, 1991
*	Synopsis:  Dot product of 2 3-vectors
*	Keywords:  ees, matvec, vector
*	Revisions:
*	mm/dd/yy   name		description
*H-								*/
/*U+
*	Usage	:  double dr3dot(x, y)
*		   double  x[3], y[3];
*	Input	:  x, y:  3-vectors
*	Output	:  returns dot product of x and y
*U-								*/

/* ************************************************************ */


double dr3dot(x, y)
double x[3], y[3];

{
	double r3 = 0;
	int i;

	for (i = 0; i < 3; i++)
	   r3 += x[i] * y[i];

	return(r3);
}

/* $Id: dr3dot.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
