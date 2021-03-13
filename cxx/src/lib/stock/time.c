#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/param.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "stock.h"

#ifndef MAX
#define MAX(a,b)        ((a)>(b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)        ((a)<(b) ? (a) : (b))
#endif

char *dbstrf2c();
char *DBL_compose_filename();

static int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31,31};

int htoe(struct date_time *dt)
{
	dt->epoch = 
	dtoepoch(dt->date) + 
	dt->hour * 3600. + 
	dt->minute * 60. +
	dt->second;
	return(0);
}

int timeprint(struct date_time *dt)
{
	printf("%15.3f %8ld %s %2d,%4d %2d:%02d:%02.3f\n",
	dt->epoch,
	dt->date,
	dt->mname,
	dt->day,
	dt->year,
	dt->hour,
	dt->minute,
	dt->second);
	return(0);
}

int zh_today(struct date_time *dt)
{
	dt->epoch = dtoepoch( todaysdate());
	etoh(dt);
	return(0);
}

int mdtodate(struct date_time *dt)
{
	int i,dim;
	dt->doy = 0;
	for( i = 0 ; i < dt->month - 1 ; i++ ){
		dim = days_in_month[i];
		if( i == 1 && isleap(dt->year) ) dim++;
		dt->doy += dim;
	}
	dt->doy += dt->day;
	dt->date = 1000 * dt->year + dt->doy;
	return(0);
}

int time_string2epoch(char *string,double *time)
/*
 *	Accepts time as a string in one of three forms:
 *		"557656335.496"			- epoch time.
 *		"1988:feb:15:22:35:43.5"	- year:month:day:hour:min:sec
 *		"1990078:03:10:15.6"		- juliandate:hour:min:sec
 */

{
	*time = str2epoch(string) ;
	return (1);
}

/*
int
my_strncasecmp(str1, str2, len)

char *      str1;
char *            str2;
int                     len;

{
	int i1, i2, i;

	i1 = strlen(str1);
	i2 = strlen(str2);
	for (i=0; i<i1; i++) str1[i] = tolower(str1[i]);
	for (i=0; i<i1; i++) str2[i] = tolower(str2[i]);
	return (strncmp(str1, str2, len));
}
*/

/* collection of time conversion utility subroutines */
	/* convert julian date to epoch time */

#define SECONDS_PER_DAY (24*60*60)

static double year2days(double year)
{
    double leapdays ; 
    leapdays =   floor((year - 1969) /   4.0)  
	       - floor((year - 1901) / 100.0) 
	       + floor((year - 1601) / 400.0);
    return  (year - 1970) * 365 + leapdays ;
}


double dtoepoch(long date)
{
        int year, day ;
	double days, leapdays ; 

        year = date / 1000;
        day  = date % 1000;
 
	days  = day-1 + (year - 1970) * 365 ;
	leapdays =   floor((year - 1969) /   4.0)  
		   - floor((year - 1901) / 100.0) 
		   + floor((year - 1601) / 400.0);
	days += leapdays ;

        return days*SECONDS_PER_DAY ;
}

/* return true if leap year else false */
int isleap(int year)
{
	return(year % 4 == 0 && year % 100 != 0 || year % 400 == 0);
}
	/* time conversion stolen from original archive */
double timecon(char *timestr)
{
	double tnum ;
	char con[20];
	long len,i,j;

	strcpy(con,timestr);
	len = strlen(con);

	for( i=0 ; isdigit(con[i]) && i < len ; i++ );
	con[i] = '\0';		/* cut it off before first non-digit */
	tnum = atoi(con) * 3600.;
	if( i >= len ) return( tnum );

	for( j = ++i ; isdigit(con[j]) && j < len ; j++ );
	con[j] = '\0';
	tnum += atoi(&con[i]) * 60.;
	if( i >= len ) return( tnum );

	for( i= ++j ; (isdigit(con[i]) || con[i] == '.') && i < len ; i++ );
	con[i] = '\0';
	tnum += atof(&con[j]);

	return(tnum);
}
/* return todays date in a long (epoch = Jan 1,1970) */
long todaysdate()
{
	long now;
	struct tm *tsp ;

	now = time(0);		/* get epoch time (in GMT) */
	tsp = gmtime(&now);	/* disect it */
	now = 1000 * ( tsp->tm_year ) + ( tsp->tm_yday ) + 1;
	return( now + 1900000 );
}

static char *month_name[] =
{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

int month_day(struct date_time *dt)
{
	int i,dim,leap;

	leap = isleap(dt->year);
	dt->day = dt->doy;
	for( i = 0 ; i < 12 ; i ++ ){
		dim = days_in_month[i];
		if( leap && i == 1 ) dim++;
		if( dt->day <= dim ) break;
		dt->day -= dim;
	}
	dt->month = i + 1;
	strcpy(dt->mname,month_name[i]);
	return 0 ;
}

int etoh(struct date_time *dt)
{
	double doy, year, days ;
	double remainder, rdays ;

	rdays = floor(dt->epoch / 86400.);
	remainder = dt->epoch - rdays * 86400. ; 
	dt->hour = floor(remainder/3600.) ; 
	remainder -= dt->hour*3600. ; 
	dt->minute = floor(remainder/60.) ; 
	remainder -= dt->minute*60. ; 
	dt->second = remainder ;

	/* given a number of days, we want to find the year and day of year */
	year = floor(rdays/365.) + 1970 ; 
	days = year2days ( year ) ; 
	doy = rdays - days ; 
	if ( doy < -366 ) { /* epoch way out of range */
	    doy = 0 ;
	} if ( doy < 0 ) {
	    year-- ; 
	    if ( isleap ((int) year) ) 
		doy += 366 ; 
	    else
		doy += 365 ; 
	} else if ( doy >= 730 ) {
	    doy = 0 ; 
	} else if ( isleap((int) year) ) {
	    if ( doy > 365 ) {
		year++; 
		doy -= 366 ; 
	    }
	} else if ( doy > 364 ) {
	    year++ ; 
	    doy -= 365 ; 
	}

	/* All this nonsense is to keep returning somewhat reasonable
	   numbers even when the input epoch time is way out of range,
	   like -1e219.  
        */
	dt->doy = doy+1 ;
	dt->year = MAX(MIN(year,LONG_MAX/10000.),LONG_MIN/10000.) ;

	dt->date = dt->year * 1000 + dt->doy;
	month_day(dt);
	return 0;
}

char *dbstrf2c(char *string, int len)
{
        static char out[512];
        int i, j, k;

        out[0] = '\0';
        if (len < 1) return (out);
        for (i=0; i<len; i++) if (string[i] != ' ') break;
        if (i == len) return (out);
        for (j=len; j>0; j--) if (string[j-1] != ' ') break;
        if (j == 0) return (out);
        for (k=0; i<j; i++) out[k++] = string[i];
        out[k] = '\0';
        return (out);
}

/*
 ******************************************************************************
 *
 *	get_nice_times()
 *
 ******************************************************************************
 */


int get_nice_times(double tstart, double tend, int maxincs,
        double *nstart, double *ninc)
/*
 *	get_nice_times will produce "nice" timing values which are used mainly
 *	for generating timing lines on a seismic trace plot.
 *
 *	Inputs  -	tstart	= An epochal time in seconds corresponding
 *				  to the start of the time window within which
 *				  the "nice" time values will be determined.
 *	         	tend	= An epochal time in seconds corresponding
 *				  to the end of the time window within which
 *				  the "nice" time values will be determined.
 *			maxincs	= An integer which is the maximum number of
 *				  "nice" time values which are to be found
 *				  within the specified time window. The actual
 *				  number of "nice" time values is guaranteed
 *				  to be less than or equal to this number.
 *
 *	Outputs -	nstart	= A pointer to an epochal time in seconds which
 *				  corresponds to the first "nice" time value.
 *				  Note: If there are no "nice" time values
 *				  within the specified window, then *nstart
 *				  is set to zero.
 *			ninc	= A pointer to a time interval in seconds which
 *				  corresponds to the "nice" time increment.
 */

{
	double window;
	double time;
	int itime;
	double exp;
	double ref;
	struct date_time date_time;

/*
 *	Find the "nice" time increment
 */
	window = fabs(tend - tstart);
	*nstart = 0.;
	if (window == 0. || maxincs <= 0) return 0;
	date_time.epoch = tstart;
	etoh(&date_time);
	window = window / maxincs;
/*
 *	sub-second marks
 */
	if (window <= 1.0) {
		exp = 1.0;
		while (1) {
			if (window > 0.5) {
				*ninc = exp;
				time = date_time.second/exp;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*exp;
				*nstart = tstart + time - date_time.second;
				break;
			} else if (window > 0.2) {
				time = 0.5*exp;
				*ninc = time;
				time = 2.0*date_time.second/exp;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*0.5*exp;
				*nstart = tstart + time - date_time.second;
				break;
			} else if (window > 0.1) {
				time = 0.2*exp;
				*ninc = time;
				time = 5.0*date_time.second/exp;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*0.2*exp;
				*nstart = tstart + time - date_time.second;
				break;
			}
			exp *= 0.1;
			window *= 10.0;
		}
	} else {
/*
 *	second marks
 */
		if (window <= 2.0) {
			*ninc = 2.0;
			time = date_time.second/2.0;
			itime = time - 0.0001;
			itime++;
			time = itime;
			time = time*2.0;
			*nstart = tstart + time - date_time.second;
		} else if (window <= 5.0) {
			*ninc = 5.0;
			time = date_time.second/5.0;
			itime = time - 0.0001;
			itime++;
			time = itime;
			time = time*5.0;
			*nstart = tstart + time - date_time.second;
		} else if (window <= 10.0) {
			*ninc = 10.0;
			time = date_time.second/10.0;
			itime = time - 0.0001;
			itime++;
			time = itime;
			time = time*10.0;
			*nstart = tstart + time - date_time.second;
		} else if (window <= 15.0) {
			*ninc = 15.0;
			time = date_time.second/15.0;
			itime = time - 0.0001;
			itime++;
			time = itime;
			time = time*15.0;
			*nstart = tstart + time - date_time.second;
		} else if (window <= 20.0) {
			*ninc = 20.0;
			time = date_time.second/20.0;
			itime = time - 0.0001;
			itime++;
			time = itime;
			time = time*20.0;
			*nstart = tstart + time - date_time.second;
		} else if (window <= 30.0) {
			*ninc = 30.0;
			time = date_time.second/30.0;
			itime = time - 0.0001;
			itime++;
			time = itime;
			time = time*30.0;
			*nstart = tstart + time - date_time.second;
		} else if (window <= 60.0) {
			*ninc = 60.0;
			time = 60.0;
			if (date_time.second == 0.0) time = 0.0;
			*nstart = tstart + time - date_time.second;
		} else {
/*
 *	minute marks
 */
			if (window <= 120.0) {
				*ninc = 120.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = ref/120.0;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*120.0;
				*nstart = tstart + time - ref;
			} else if (window <= 300.0) {
				*ninc = 300.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = ref/300.0;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*300.0;
				*nstart = tstart + time - ref;
			} else if (window <= 600.0) {
				*ninc = 600.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = ref/600.0;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*600.0;
				*nstart = tstart + time - ref;
			} else if (window <= 900.0) {
				*ninc = 900.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = ref/900.0;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*900.0;
				*nstart = tstart + time - ref;
			} else if (window <= 1200.0) {
				*ninc = 1200.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = ref/1200.0;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*1200.0;
				*nstart = tstart + time - ref;
			} else if (window <= 1800.0) {
				*ninc = 1800.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = ref/1800.0;
				itime = time - 0.0001;
				itime++;
				time = itime;
				time = time*1800.0;
				*nstart = tstart + time - ref;
			} else if (window <= 3600.0) {
				*ninc = 3600.0;
				ref = date_time.minute;
				ref = ref*60.0 + date_time.second;
				time = 3600.0;
				if (ref == 0.0) time = 0.0;
				*nstart = tstart + time - ref;
			} else {
/*
 *	hour marks 
 */
				if (window <= 2.0*3600.0) {
					*ninc = 2.0*3600.0;
					ref = date_time.hour;
					ref = ref*60.0 + date_time.minute;
					ref = ref*60.0 + date_time.second;
					time = ref/(*ninc);
					itime = time - 0.0001;
					itime++;
					time = itime;
					time = time*(*ninc);
					*nstart = tstart + time - ref;
				} else if (window <= 4.0*3600.0) {
					*ninc = 4.0*3600.0;
					ref = date_time.hour;
					ref = ref*60.0 + date_time.minute;
					ref = ref*60.0 + date_time.second;
					time = ref/(*ninc);
					itime = time - 0.0001;
					itime++;
					time = itime;
					time = time*(*ninc);
					*nstart = tstart + time - ref;
				} else if (window <= 6.0*3600.0) {
					*ninc = 6.0*3600.0;
					ref = date_time.hour;
					ref = ref*60.0 + date_time.minute;
					ref = ref*60.0 + date_time.second;
					time = ref/(*ninc);
					itime = time - 0.0001;
					itime++;
					time = itime;
					time = time*(*ninc);
					*nstart = tstart + time - ref;
				} else if (window <= 8.0*3600.0) {
					*ninc = 8.0*3600.0;
					ref = date_time.hour;
					ref = ref*60.0 + date_time.minute;
					ref = ref*60.0 + date_time.second;
					time = ref/(*ninc);
					itime = time - 0.0001;
					itime++;
					time = itime;
					time = time*(*ninc);
					*nstart = tstart + time - ref;
				} else if (window <= 12.0*3600.0) {
					*ninc = 12.0*3600.0;
					ref = date_time.hour;
					ref = ref*60.0 + date_time.minute;
					ref = ref*60.0 + date_time.second;
					time = ref/(*ninc);
					itime = time - 0.0001;
					itime++;
					time = itime;
					time = time*(*ninc);
					*nstart = tstart + time - ref;
				} else if (window <= 24.0*3600.0) {
					*ninc = 24.0*3600.0;
					ref = date_time.hour;
					ref = ref*60.0 + date_time.minute;
					ref = ref*60.0 + date_time.second;
					time = ref/(*ninc);
					itime = time - 0.0001;
					itime++;
					time = itime;
					time = time*(*ninc);
					*nstart = tstart + time - ref;
				} else {
/*
 *	day marks
 */
					exp = 1.0;
					window = window / 86400.0;
					while (1) {
						if (window <= 2.0) {
							time = 2.0*exp;
							*ninc = time*86400.0;
							ref = date_time.doy;
							ref = ref*24.0 + date_time.hour;
							ref = ref*60.0 + date_time.minute;
							ref = ref*60.0 + date_time.second;
							time = ref/(*ninc);
							itime = time - 0.0001;
							itime++;
							time = itime;
							time = time*(*ninc);
							*nstart = tstart + time - ref;
							break;
						} else if (window <= 5.0) {
							time = 5.0*exp;
							*ninc = time*86400.0;
							ref = date_time.doy;
							ref = ref*24.0 + date_time.hour;
							ref = ref*60.0 + date_time.minute;
							ref = ref*60.0 + date_time.second;
							time = ref/(*ninc);
							itime = time - 0.0001;
							itime++;
							time = itime;
							time = time*(*ninc);
							*nstart = tstart + time - ref;
							break;
						} else if (window <= 10.0) {
							time = 10.0*exp;
							*ninc = time*86400.0;
							ref = date_time.doy;
							ref = ref*24.0 + date_time.hour;
							ref = ref*60.0 + date_time.minute;
							ref = ref*60.0 + date_time.second;
							time = ref/(*ninc);
							itime = time - 0.0001;
							itime++;
							time = itime;
							time = time*(*ninc);
							*nstart = tstart + time - ref;
							break;
						}
						exp *= 10.0;
						window *= 0.1;
					}
				}
			}
		}
	}
    return 0 ;
}

/* $Id: time.c,v 1.2 1997/04/21 14:31:21 danny Exp $ */
