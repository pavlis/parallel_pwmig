
#include <stdio.h>

double          dtoepoch();

struct date_time
{
    double          epoch;
    long            date;
    int             year;
    int             month;
    char            mname[4];
    int             day;
    int             doy;
    int             hour;
    int             minute;
    float           second;
};
void month_day(struct date_time *dt);
void mdtodate(struct date_time *dt);
void etoh(struct date_time *dt);
void 
doy2mday ( doy, year, month, day ) 
int year, doy ;
int *month, *day ; 
{
    struct date_time dt ; 

    dt.year = year ; 
    dt.doy = doy ; 
    month_day ( &dt ) ; 
    *day = dt.day ; 
    *month = dt.month ; 
}


int 
mday2doy ( year, month, day ) 
int month, day, year ; 
{
    struct date_time dt ; 

    dt.year = year ; 
    dt.month = month ;
    dt.day = day ; 

    mdtodate ( &dt ) ; 

    return dt.doy ;
}


double
h2e(iyear, iday, ihour, imin, sec)

int            iyear;
int            iday;
int            ihour;
int            imin;
double         sec;

/* htoe will compute a double precision epoch time in seconds from 00:00:00.0
 * 1 Jan 1970 GMT. */

{
    int             jdate;
    double          epoch;

    jdate = 1000 * (iyear) + (iday);
    epoch = dtoepoch(jdate);
    epoch += 3600.0 * (ihour) + 60.0 * (imin) + sec;
    return epoch ; 
}



void
e2h(epoch, iyear, iday, ihour, imin, sec)

double         epoch;
int            *iyear;
int            *iday;
int            *ihour;
int            *imin;
double         *sec;

/* e2h will compute human readable time from a double precision epoch time
 * in seconds from 00:00:00.0 1 Jan 1970 GMT. */

{
    struct date_time dt;

    dt.epoch = epoch;
    etoh(&dt);
    *iyear = dt.year;
    *iday = dt.doy;
    *ihour = dt.hour;
    *imin = dt.minute;
    *sec = dt.second;
}

#ifdef DEBUG

static char *months[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
			"Sep", "Oct", "Nov", "Dec" } ; 

main ()
{
    char aline[1024] ; 
    int year, year2, day, day2, hour, hour2, min, min2 ; 
    double second, second2 ; 
    double epoch ;
    int month, mday ; 
    int month2, mday2 ; 
    int i, y, doy ; 

    for ( y=1988 ; y < 1992 ; y++ ) 
    for ( i=0 ; i<366 ; i++ ) 
	{
	doy2mday ( i, y, &month, &day ) ; 
	doy = mday2doy ( y, month, day ) ; 
	if ( doy != i ) 
	    fprintf ( stderr, "day of year conversion fails for doy=%d, %d\n", i, y ) ; 
	}

    printf ( "\nEnter year day hour minute second: " ) ; 
    while ( gets(aline) ) 
	{
	sscanf ( aline, "%d %d %d %d %lf", &year, &day, &hour, &min, &second ) ; 
	doy2mday ( day, year, &month, &mday ) ; 
	printf ( "%2d:%02d:%05.3lf %s %d, %d\n", hour, min, second, 
			months[month-1], mday, year ) ; 

	epoch = h2e ( year, day, hour, min, second ) ; 
	printf ( "Epoch time: %10.3lf\n", epoch ) ; 

	e2h ( epoch, &year2, &day2, &hour2, &min2, &second2 ) ; 
	doy2mday ( day2, year2, &month2, &mday2 ) ; 
	printf ( "%2d:%02d:%05.3lf %s %d, %d\n", hour2, min2, second2, 
			months[month2-1], mday2, year2 ) ; 
	printf ( "\nEnter year day hour minute second: " ) ; 
	}
}

#endif

/* $Id: h2e.c,v 1.1.1.1 1997/04/12 04:18:50 danq Exp $ */
