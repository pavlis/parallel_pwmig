/*H+
 *	Title	:  dr3sub.c
 *	Author	:  George Kaplan, from MIT matvec library
 *	Date	:  02 Mar 1990
 *	Synopsis:  Subtract a vector from another
 *	Keywords:  ees, matvec, r3sub, vector
 *	@(#)r3sub_.f	1.2 3/2/90  UCB SSL
 *	Revisions:
 *	mm/dd/yy   name		description
 *       7/19/91  danq  double precision version
 */
/*U+
 *	Usage	:  dr3sub(a, b, c)
 *	Input	:  a, b:  3-vectors
 *	Output	:  c:  3-vector; result of a - b
 *U-
 */


void            dr3sub(v1, v2, v3)
    double          v1[3], v2[3], v3[3];
{
    int             i;

    for (i = 0; i < 3; i++)
	v3[i] = v1[i] - v2[i];
}

/* $Id: dr3sub.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
