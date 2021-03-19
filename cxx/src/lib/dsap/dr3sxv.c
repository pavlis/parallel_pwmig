/*H+
*	Title	:  r3sxv.c
*	Author	:  Jim Lewis (adapted from r3sxv.c)
*	Date	:  19 April, 1991
*	Synopsis:  Multiply scalar times 3-vector
*	Keywords:  ees, matvec, r3sxv, scalar, vector
*	Revisions:
*	mm/dd/yy   name		description
*U+
*	Usage	:  void r3sxv(a, b, c)
*		   double a, b[3], c[3];
*	Input	:  a:  scalar
*		   b:  3-vector
*	Output	:  c:  3-vector product of a times b
*U-								*/
/*D+
*	Multiply a scalar times a 3-vector
*
*	dr3sxv(a,b,c)
*	a is the scalar
*	b is the vector
*	c is the vector result (c may be the same as b)
*D-								*/

/* ************************************************************ */


void  dr3sxv(a, b, c)
double a, b[3];
double c[3];
{
	c[0] = a * b[0];
	c[1] = a * b[1];
	c[2] = a * b[2];
}

/* $Id: dr3sxv.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
