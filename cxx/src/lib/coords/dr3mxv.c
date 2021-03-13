/*H+
*	Title	:  dr3mxv.c
*	Author	:  George Kaplan, adapted from MIT matvec library
*	Date	:  01 Mar 1990
*	Synopsis:  Multiply 3x3 matrix by a 3-vector
*	Keywords:  ees, matvec, r3mxv, matrix, vector
*	Revisions:
*	mm/dd/yy   name		description
*U+
*	Usage	:  void dr3mxv(a, b, c)
*		   double a[9], b[3], c[3];
*	Input	:  a:  3x3 matrix
*		   b:  3 vector
*	Output	:  c:  product 3-vector a times b
*U-								*/
/*D+
*	Subroutine to multiply a 3 by 3 matrix by a 3-vector
*	dr3mxv(a,b,c)
*
*	a is the input matrix
*	b is the input vector
*	c contains, on return, the vector result a x b
*		(c may be the same as b)
*D-								*/

/* ************************************************************ */


void dr3mxv(a, b, c)
double  a[9], b[3], c[3];
{
	double d[3];
	int i;

	d[0] = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	d[1] = a[3] * b[0] + a[4] * b[1] + a[5] * b[2];
	d[2] = a[6] * b[0] + a[7] * b[1] + a[8] * b[2];
	for (i = 0; i < 3; i++)  c[i] = d[i];
}

/* $Id: dr3mxv.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
