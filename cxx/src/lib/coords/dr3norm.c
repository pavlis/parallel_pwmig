/*H+    
*	Title	: dr3norm     *******
*	Author	: Fil Machi
*	Date	: 02/04/91
*	Synopsis: normalize vector
*	Keywords: vector, norm
*	Revisions:
*	mm/dd/yy	name	description
*H-
*U+
* 	Usage	: dr3norm(double v[3])
*	Input	: 
*	Output	: 
*U-
* ****************************************************************** */

#include <math.h>

void dr3norm(v)
double v[3];
{
  register double norm;
  
  norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (norm>0) {
    norm=1.0/norm;
    v[0] *= norm;  
    v[1] *= norm;  
    v[2] *= norm;  
    }
  else
    v[0]=v[1]=v[2]=0.0;

  return;  
}

/* $Id: dr3norm.c,v 1.1.1.1 1997/04/12 04:18:49 danq Exp $ */
