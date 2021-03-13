#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include "coords.h"
#include "stock.h"

double epoch ( int yd )
{
    return h2e ( yd/1000, yd % 1000, 0, 0, 0.0 ) ; 
}

int yearday ( double e ) 
{
    int year, doy, hour, minute ; 
    double second ; 

    e2h ( e, &year, &doy, &hour, &minute, &second ) ; 
    return year *1000 + doy ; 
}

char *strydtime ( double e )
{
    int year, doy, hour, minute, month, day ;
    double second ;
    char s [STRSZ] ;
 
    /* round to the display precision */
    e *= 1000.0 ; 
    e += 0.5 ; 
    e = floor(e)/1000. ; 

    e2h ( e, &year, &doy, &hour, &minute, &second ) ;

    doy2mday ( doy, year, &month, &day ) ;
 
    sprintf ( s, "%2d/%02d/%4d (%03d) %2d:%02d:%06.3f",
        month, day, year, doy, hour, minute, second ) ;
    return strdup(s) ;
}

char *strtime ( double e ) 
{
    int year, doy, hour, minute, month, day ; 
    double second ; 
    char s [STRSZ] ; 
 
    /* round to the display precision */
    e *= 1000.0 ; 
    e += 0.5 ; 
    e = floor(e)/1000. ; 

    e2h ( e, &year, &doy, &hour, &minute, &second ) ; 
    doy2mday ( doy, year, &month, &day ) ; 

    sprintf ( s, "%2d/%02d/%4d  %2d:%02d:%06.3f", 
	month, day, year, hour, minute, second ) ; 
    return strdup(s) ; 
}

char *strdate ( double e )
{
    int year, doy, hour, minute, month, day ; 
    double second ; 
    char s [STRSZ] ; 

    e2h ( e, &year, &doy, &hour, &minute, &second ) ; 
    doy2mday ( doy, year, &month, &day ) ; 

    sprintf ( s, "%2d/%02d/%4d", month, day, year ) ; 
    return strdup(s) ; 
}

static char *Month[] = {
    "January", 
    "February", 
    "March", 
    "April", 
    "May", 
    "June", 
    "July",
    "August",
    "September",
    "October", 
    "November", 
    "December"
    } ;

static int
e2hz ( double epoch, char *timezone, int *year, int *doy, int *month, int *mday, int *hour, int *minute, double *second, char **zonename ) 
{
    int retcode = 0 ; 

    if ( timezone == NULL ) {
	e2h ( epoch, year, doy, hour, minute, second ) ; 
	doy2mday ( *doy, *year, month, mday ) ; 
	*zonename = "UTC" ;
    } else {
	time_t	    ut ; 
	char	    s[STRSZ] ;
	char 	   *oldTZ ; 
	struct tm       utm ; 
	double msec ;

	ut = floor ( epoch ) ; 
	msec = epoch - ut ;
	if ( *timezone != '\0' ) {
	    oldTZ = strdup(getenv("TZ")) ; 
	    sprintf(s, "TZ=%s", timezone) ; 
	    putenv ( s ) ; 
	    errno = 0 ; 
	    tzset() ; 
	    /* show_tz ( "e2hz before" ) ; */
	    if ( errno != 0 ) {
		errno = 0 ; 
		elog_log (0, "Unrecognized timezone '%s'\n", timezone ) ; 
		retcode = -1 ;
		putenv( "TZ=UTC") ;
		tzset() ; 
		/* show_tz ( "e2hz ->UTC" ) ; */
	    }
	}
	localtime_r(&ut, &utm ) ; 
	if ( *timezone != '\0' ) {
	    sprintf(s, "TZ=%s", oldTZ) ; 
	    putenv ( s ) ; 
	    free(oldTZ) ;
	    tzset() ; 
	    /* show_tz ( "e2hz restoring" ) ; */
	}
	*month = utm.tm_mon+1 ; 
	*mday = utm.tm_mday ; 
	*year = utm.tm_year+1900 ;
	*hour = utm.tm_hour ; 
	*minute = utm.tm_min ; 
	*second = (double) utm.tm_sec + msec ;
	*doy = utm.tm_yday+1 ; 
	*zonename = (utm.tm_isdst ? tzname[1] : tzname[0] ) ; 
    }
    return retcode ;
}

static char *
ztime2str ( char *f, double e, int year, int doy, int month, int mday, int hour, int minute, double second, char *zonename ) 
{
    char s [STRSZ] ;
    char *in, *out ;
    char c ;

    in = f ; 
    out = s; 
    *out = 0 ; 

    while ( (c = *in++) != 0 )  {
	if ( c != '%' )
	    *out++ = c ; 
	else
	    {
	    c = *in++ ; 
	    if ( c == 0 ) break ; 
	    switch (c) {

	    case 'b' :
		strncpy ( out, Month[month-1], 3 ) ; 
		out += 3; 
		*out = 0 ; 
		break ;

	    case 'B' :
		strcpy ( out, Month[month-1] ) ; 
		break ;

	    case 'D' : 
		sprintf ( out, "%2d/%2d/%2d", month, mday, year-1900 ) ; 
		break ;

	    case 'd' : 
		sprintf ( out, "%02d", mday ) ; 
		break ;
	    
	    case 'E' :
		sprintf ( out, "%15.3f", e ) ; 
		break ;

	    case 'e' :
		sprintf ( out, "%2d", mday ) ; 
		break ;
   
	    case 'H' : 
		sprintf ( out, "%02d", hour ) ; 
		break ;

	    case 'I' : 
		sprintf ( out, "%02d", hour % 12 ) ; 
		break ;

	    case 'g':
		sprintf ( out, "%04d-%03d", year, doy ) ; 
		break ;

	    case 'G':
		sprintf ( out, "%04d-%02d-%02d", year, month, mday ) ; 
		break ;

	    case 'j' : 
		sprintf ( out, "%03d", doy ) ; 
		break ;
	    
	    case 'k' : 
		sprintf ( out, "%2d", hour ) ; 
		break ;

	    case 'l' : 
		sprintf ( out, "%2d", hour % 12 ) ; 
		break ;

	    case 'm' :
		sprintf ( out, "%02d", month ) ; 
		break ;
	    
	    case 'M' : 
		sprintf ( out, "%02d", minute ) ; 
		break ;
	    
	    case 'n' : 
		*out++ = '\n' ; 
		*out = 0 ; 
		break ;
	    
	    case 'p' :
		sprintf ( out, "%s", hour >= 12 ? "PM" : "AM" ) ; 
		break ;
	    
	    case 'S' :
		sprintf ( out, "%02d", (int) second ) ; 
		break ;
	    
	    case 's' :
		sprintf ( out, "%03d", (int) rint(second * 1000) % 1000 ) ; 
		break ;
	    
	    case 't' :
		*out++ = '\t' ; 
		*out = 0 ; 
		break ; 

	    case 'T' :
		sprintf ( out, "%02d:%02d:%06.3f", hour, minute, second ) ; 
		break ;
	    
	    case 'y' :
		sprintf ( out, "%2d", year-1900 ) ; 
		break ; 
	    
	    case 'Y' :
		sprintf ( out, "%4d", year ) ; 
		break ; 

	    case 'z':
		sprintf ( out, "%s", zonename ) ; 
		break ;

	    default : 
		*out++ = c ; 
		*out = 0 ; 
		break ;
	    }

	    while ( *out != 0 ) out++ ; 
	}
    }
    *out = 0 ;
    return strdup(s) ; 
}

char *epoch2str ( double epoch, char *f )
{
    int year, doy, hour, minute, month, mday ; 
    double second ; 
    char *zonename ;

    e2hz ( epoch, 0, &year, &doy, &month, &mday, &hour, &minute, &second, &zonename ) ;
    return ztime2str ( f, epoch, year, doy, month, mday, hour, minute, second, zonename ) ;
}

char *zepoch2str ( double epoch, char *f, char *timezone ) 
{
    int year, doy, hour, minute, month, mday ; 
    double second ; 
    char *zonename ;

    e2hz ( epoch, timezone, &year, &doy, &month, &mday, &hour, &minute, &second, &zonename ) ;
    return ztime2str ( f, epoch, year, doy, month, mday, hour, minute, second, zonename ) ;
}

int fmttime ( double epoch, char *result, int fmt, char *timezone )
{
    int year, doy, hour, minute, month, mday ; 
    double second ; 
    int retcode = 0 ;
    char *zonename = 0 ;
 
    /* round to the display precision */
    epoch *= 1000.0 ; 
    epoch += 0.5 ; 
    epoch = floor(epoch)/1000. ; 

    e2hz ( epoch, timezone, &year, &doy, &month, &mday, &hour, &minute, &second, &zonename ) ;

    switch (fmt) {

	case STRYDTIME:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%2d/%02d/%4d (%03d) %2d:%02d:%06.3f %s",
		    month, mday, year, doy, hour, minute, second, zonename ) ;
	    else
		sprintf ( result, "%2d/%02d/%4d (%03d) %2d:%02d:%06.3f",
		    month, mday, year, doy, hour, minute, second ) ;
	    break ;

	case STRTIME:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%2d/%02d/%4d  %2d:%02d:%06.3f %s", 
		    month, mday, year, hour, minute, second, zonename ) ; 
	    else
		sprintf ( result, "%2d/%02d/%4d  %2d:%02d:%06.3f", 
		    month, mday, year, hour, minute, second ) ; 
	    break ;

	case STRDATE:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%2d/%02d/%4d %s", 
		    month, mday, year, zonename ) ;
	    else
		sprintf ( result, "%2d/%02d/%4d", 
		    month, mday, year ) ;
	    break ;

	case MDTIME:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%04d-%02d-%02d  %2d:%02d:%06.3f %s", 
		    year, month, mday, hour, minute, second, zonename ) ; 
	    else
		sprintf ( result, "%04d-%02d-%02d  %2d:%02d:%06.3f", 
		    year, month, mday, hour, minute, second ) ; 
	    break ;

	default :
	case YDTIME:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%04d-%03d %2d:%02d:%06.3f %s", 
		    year, doy, hour, minute, second, zonename ) ; 
	    else
		sprintf ( result, "%04d-%03d %2d:%02d:%06.3f", 
		    year, doy, hour, minute, second ) ; 
	    break ;

	case MD:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%04d-%02d-%02d %s", 
		    year, month, mday, zonename ) ; 
	    else
		sprintf ( result, "%04d-%02d-%02d", 
		    year, month, mday ) ; 
	    break ;

	case YD:
	    if ( zonename && strcmp(zonename, "UTC") != 0 ) 
		sprintf ( result, "%04d-%03d %s", 
		    year, doy, zonename ) ; 
	    else
		sprintf ( result, "%04d-%03d", 
		    year, doy ) ; 
	    break ;
    }

    return retcode ;
}

#define SECONDS_PER_YEAR (365.*24.*3600.)
#define SECONDS_PER_DAY       (24.*3600.)
#define SECONDS_PER_HOUR          (3600.)
#define SECONDS_PER_MINUTE          (60.)

char *strdelta(double e)
{
    char string[50], *cp ; 
    double y, d, h, m, s ; 
 
    /* round to the display precision */
    e *= 1000.0 ; 
    e += 0.5 ; 
    e = floor(e)/1000. ; 

    cp = string ; 
    if ( e < 0 ) {
	*cp++ = '-' ; 
	e = -e ;
    } else {
	*cp++ = ' ' ; 
    }

    if ( e > SECONDS_PER_YEAR ) {
	y = floor(e/SECONDS_PER_YEAR) ; 
	e -= y * SECONDS_PER_YEAR ; 
	d = floor(e/SECONDS_PER_DAY) ; 
	e -= d * SECONDS_PER_DAY ; 
	sprintf ( cp, "%.0f years %.0f days", y, d ) ; 
    } else if ( e > SECONDS_PER_DAY ) {
	d = floor(e/SECONDS_PER_DAY) ; 
	e -= d * SECONDS_PER_DAY ; 
	h = floor(e*10./SECONDS_PER_HOUR) /10. ;
	e -= h * SECONDS_PER_HOUR ;
	sprintf ( cp, "%.0f days %.1f hours", d, h ) ; 
    } else if ( e > SECONDS_PER_HOUR ) {
	h = floor(e/SECONDS_PER_HOUR) ;
	e -= h * SECONDS_PER_HOUR ;
	m = floor(e/SECONDS_PER_MINUTE) ;
	e -= m * SECONDS_PER_MINUTE ;
	sprintf ( cp, "%2.0f:%02.0f hours", h, m ) ; 
    } else if ( e > SECONDS_PER_MINUTE ) {
	m = floor(e/SECONDS_PER_MINUTE) ;
	e -= m * SECONDS_PER_MINUTE ;
	s = floor(e+.5) ;
	e -= s ;
	sprintf ( cp, "%2.0f:%02.0f minutes", m, s ) ; 
    } else {
	sprintf ( cp, "%06.3f seconds", e ) ; 
    } 

    return strdup(string) ;
}


/* $Id: epoch.c,v 1.3 1997/08/17 15:47:30 danq Exp $ */
