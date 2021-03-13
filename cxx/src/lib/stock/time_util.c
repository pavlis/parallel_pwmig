/* $Name $Revision: 1.1.1.1 $ $Date: 1997/04/12 04:18:52 $  */
/************************************************************************
 *
 *  bcdtime.c
 *
 *  Decodes BCD time into secs since 00:00 Jan 1, 1970.  Will not work
 *  for any earlier times.
 *
 ***********************************************************************/
#include <time.h>
#include "coords.h"

void 
bcdtime (bcd, tstr)
unsigned char  *bcd;		       /* pointer to first byte of BCD time
				        * stamp */
struct date_time *tstr;
{
    char            tim_str[26];
    long            loc_t;
    struct tm      *tim;
    int             sec,
                    msec;

    /* Find current Year  */

    loc_t = time (0);
    tim = gmtime (&loc_t);

    tstr->year = tim->tm_year + 1900;

    sprintf (tim_str, "%02x%02x%02x%02x%02x%02x\0",
	     bcd[0], bcd[1], bcd[2], bcd[3], bcd[4], bcd[5]);
    sscanf (tim_str, "%3d%2d%2d%2d%3d", 
	    &tstr->doy, &tstr->hour, &tstr->minute, &sec, &msec);
    tstr->second = (float) sec + (float) msec / 100.0;
    tstr->date = 1000 * tstr->year + tstr->doy;

    tstr->epoch = dtoepoch (tstr->date) +
	3600.0 * tstr->hour +
	60.0 * tstr->minute +
	tstr->second;
}
 
 /* conver epoch time to BCD stamp -  YYDDDHHMMSS */
 
void epoch2bcd(bcd, epoch)
unsigned char *bcd;    /* pointer to first byte of BCD time stamp */
double epoch;
{
     double sec;
     int yr, dd, hh, min,isec, msec;
     int val[6], i;
     char tim_str[26];
     
     e2h(epoch, &yr, &dd, &hh, &min, &sec);
     isec = sec;
     msec = ( sec - ( double ) isec) * 100;
     sprintf(tim_str, "%03d%02d%02d%02d%03d\0",
                      dd, hh, min, isec, msec);
     sscanf(tim_str, "%02x%02x%02x%02x%02x%02x",
     &val[0],&val[1],&val[2],&val[3],&val[4],&val[5]);

     for( i = 0; i < 6; i++ )     
        bcd[i] = val[i];
 
}
 

/* $Id: time_util.c,v 1.1.1.1 1997/04/12 04:18:52 danq Exp $ */
