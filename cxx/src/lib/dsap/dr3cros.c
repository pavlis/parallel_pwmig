
/*H+
*	Title	:  dr3cros.c
*	Author	:  Jim Lewis (adapted from r3cros_.f)
*	Date	:  19 April, 1991
*	Synopsis:  Cross product of two 3-vectors
*	Keywords:  ees, matvec, dr3cros, vector
*	Revisions:
*	mm/dd/yy   name		description
*	Usage	:  void dr3cros(a, b, c)
*		   double a[3], b[3], c[3];
*	Input	:  a:  3-vector
*		   b:  3-vector
*	Output	:  c:  3-vector cross product of a times b
*U-								*/
/*	Cross product of two 3-vectors
*
*	dr3cros(a,b,c)
*	a is a 3-vector
*	b is a 3-vector
*	c is the vector result (c may be the same as a or b)
*D-								*/

/* ************************************************************ */

/* SCCS id ---------------------------------------------------- */

void dr3cros(a,b,c)
double a[3],b[3],c[3];
{
   double d[3];

   d[0] = a[1]*b[2] - a[2]*b[1];
   d[1] = a[2]*b[0] - a[0]*b[2];
   d[2] = a[0]*b[1] - a[1]*b[0];

   c[0] = d[0];
   c[1] = d[1];
   c[2] = d[2];

   return;
}

/* $Id: dr3cros.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
