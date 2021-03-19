/*H+
 * 	Title	:  dr3add.c
 * 	Author	:  George Kaplan, from MIT matvec library
 * 	Date	:  01 Mar 1990
 * 	Synopsis:  Add  two 3 vectors
 * 	Keywords:  ees, matvec, r3add, vector
 * 	Revisions:
 * 	mm/dd/yy   name		description
 *       7/19/91   danq   double precision version in c
 *H-
 */
/*U+
 * 	Usage	:  r3add(v1, v2, v3)
 * 	Input	:  v1, v2:  3-vectors
 * 	Output	:  v3:  vector sump of v1 and v2
 *U-
 */


void            dr3add(v1, v2, v3)
    double          v1[3], v2[3], v3[3];
{
    int             i;

    for (i = 0; i < 3; i++)
	v3[i] = v1[i] + v2[i];
}

/* $Id: dr3add.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
