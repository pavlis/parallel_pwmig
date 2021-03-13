#include <math.h>

void xrotate(double x[3],double y[3], double z[3], double theta)
{
	double st = sin((double)theta);
	double ct = cos((double)theta);
	double ct1 = 1 - ct;

	double exp1 = x[0] * x[1] * ct1;
	double exp2 = x[2] * st;
	double exp3 = x[0] * x[2] * ct1;
	double exp4 = x[1] * st;
	double exp5 = x[1] * x[2] * ct1;
	double exp6 = x[0] * st;

	double xs0 = x[0] * x[0];
	double xs1 = x[1] * x[1];
	double xs2 = x[2] * x[2];

	z[0] = ((1 - xs0) * ct + xs0) * y[0] + 
		(exp1 - exp2) * y[1] +
		(exp3 + exp4) * y[2];
	z[1] = (exp1 + exp2) * y[0] + 
		((1 - xs1) * ct + xs1) * y[1] +
		(exp5 - exp6) * y[2];
	z[2] = (exp3 - exp4) * y[0] + 
		(exp5 + exp6) * y[1] + 
		((1 - xs2) * ct + xs2) * y[2];
}

/* $Id: rotate.c,v 1.1.1.1 1997/04/12 04:18:51 danq Exp $ */
