#include <time.h>
#include "stock.h"
 
double
now ()
{
 double rightnow ;
 
#ifdef __SVR4
 struct timespec tp ;
 clock_gettime(CLOCK_REALTIME, &tp);
 rightnow = tp.tv_sec*1.0 + tp.tv_nsec/1000000000. ;
#else
 rightnow = (double) time(0) ; 
#endif

 return rightnow ;
}
 

