#ifndef _stock_
#define _stock_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define MAX_STA_SIZE	32
#define MAX_CHAN_SIZE	32

#define MAX_SIGNED_INT   (0x7fffffff)
#define MAX_UNSIGNED_INT (0xffffffffU)

#include "elog.h"
#define BITS_PER_ELEMENT (8*sizeof(char *))
#define bitnew() inittbl(0L,1L,0,0L,(int)BITS_PER_ELEMENT/8)

#define bitfree(B) freetbl(B,0)
#define bitempty(B) clrtbl(B,0) 
/*
#define HOOK_MAGIC 814615491

typedef struct Hook {
    int magic ; 
    void (*free)(void *) ;
    // These are not in old DSAP so comment out for now
    //int type ;
    void *p ;
    //pthread_mutex_t lock ;
} Hook ; 
*/

typedef struct Xlat {
        char *name ;
        int num ;
} Xlat ;


#ifdef  __cplusplus
extern "C" {
#endif

//extern Hook *new_hook ( void (*hookfree)(void *p) );
//extern void free_hook ( Hook **hookp );

//int make_pathdirs ( char *filename );
int makedir ( const char *dir );
//int mkfile ( char *name, unsigned long nbytes, char *init, int initsize );
char *datapath (const char *envname,const char *dirname,const char *filename,
        const char *suffix);
char *datafile (const char *envname,const char *filename);

int blank ( char *s );
int whitespace ( char *s );
void sncopy(char *dest,const char *source,int n);
void szcopy(char *dest,const char *source,int n);
void copystrip(char *dest,const char *source,int n);
void *memdup ( void *a, int n );
//Tbl *split (char *s,char *c );

/* These are time functions converted from dsap */
struct date_time{  /* date/time structure from old time conversion routines */
        double epoch;
        long date;
        int year;
        int month;
        char mname[4];
        int day;
        int doy;
        int hour;
        int minute;
        float second;
};

#define MDTIME                  1
#define YDTIME                  2
#define MD                      3
#define YD                      4
#define STRTIME                 5
#define STRDATE                 6
#define STRYDTIME               7
double now();
double epoch(int yd);
int yearday ( double e );
char *strydtime ( double e );
char *strtime(double e);
char *strdate(double e);
//static int e2hz ( double epoch, char *timezone, int *year, int *doy, int *month, int *mday, int *hour, int *minute, double *second, char **zonename );
//static char* ztime2str ( char *f, double e, int year, int doy, int month, int mday, int hour, int minute, double second, char *zonename );
char *epoch2str ( double epoch, char *f );
char *zepoch2str ( double epoch, char *f, char *timezone );
int fmttime ( double epoch, char *result, int fmt, char *timezone );
char *strdelta(double e);
int htoe(struct date_time *dt);
int timeprint(struct date_time *dt);
int zh_today(struct date_time *dt);
int mdtodate(struct date_time *dt);
int time_string2epoch(char *s, double *t);
//static double year2days(double year);
double dtoepoch(long date);
int isleap(int year);
double timecon(char *timestr);
long todaysdate();
int month_day(struct date_time *dt);
int etoh(struct date_time *dt);
char *dbstrf2c(char *s,int len);
int get_nice_times(double tstart, double tend, int maxincs, 
        double *nstart, double *ninc);
void doy2mday(int doy,int year, int *month, int *day);
int mday2doy(int yr,int mon, int day);
double h2e(int yr, int day, int hr, int min, double sec);
void e2h(double e, int *yr, int *day, int *hr, int *min, double *sec);
//static int check_tz();
int zstr2epoch(char *, double *e);
double str2epoch(char *);
//static int my_yyinput(char *buf,int ms);
#ifdef	__cplusplus
}
#endif

/* When a macro generates multiple output statements, there is a danger
 * that these statements might lead to odd results in a context where a
 * single statement is expected, eg. after an if -- STMT forces the macro
 * output to be a single statement, and the following stuff with ZERO is
 * to make lint work properly. */

#ifdef lint
extern int      ZERO;
#else
#define ZERO 0
#endif

#define STMT(stuff) do { stuff } while (ZERO)

/* The following macros simplify memory allocation and testing
 *
 * allot and reallot are convenient interfaces to malloc and realloc which
 * perform testing on the result without cluttering the code with
 * branches and messages which are seldom if ever used.
 *
 * ALLOTERROR specifies what happens when a malloc or realloc fails -- it
 * may be overridden with a special macro within the file */

#ifdef __STDC__
#define ALLOTDIE(ptr,size)	elog_die(1,"Out of memory in %s at line %d of '%s'\n\tCan't malloc %lu bytes for " #ptr "\n", __func__, __LINE__, __FILE__, (long)size)
#else
#define ALLOTDIE(ptr,size)	elog_die(1,"Out of memory in %s at line %d of '%s'\n\tCan't malloc %lu bytes for ptr\n", __func__, __LINE__, __FILE__, (long)size)
#endif

#ifndef BADALLOT
#define BADALLOT -1
#endif

#ifdef lint
extern int ReTuRn(int) ;
#define BADRETURN ReTuRn(BADALLOT)
#else
#define BADRETURN return BADALLOT
#endif

#ifdef __STDC__
#define ALLOTREGISTER_ERROR(ptr,size) \
    STMT(elog_log(1,"Out of memory in %s at line %d of '%s'\n\tCan't malloc %lu bytes for " #ptr "\n", __func__, __LINE__, __FILE__, (long)size); BADRETURN ;)
#else
#define ALLOTREGISTER_ERROR(ptr,size) \
    STMT(elog_log(1,"Out of memory in %s at line %d of '%s'\n\tCan't malloc %lu bytes for ptr\n", __func__, __LINE__, __FILE__, (long)size); BADRETURN ;)
#endif

#ifdef __STDC__
#define ALLOTCOMPLAIN(ptr,size) \
    STMT({elog_complain(1,"Out of memory in %s at line %d of '%s'\n\tCan't malloc %lu bytes for " #ptr "\n", __func__, __LINE__, __FILE__, (long)size); BADRETURN ;})
#else
#define ALLOTCOMPLAIN(ptr,size) \
    STMT({elog_complain(1,"Out of memory in %s at line %d of '%s'\n\tCan't malloc %lu bytes for ptr\n", __func__, __LINE__, __FILE__, (long)size); BADRETURN ;})
#endif

#ifndef ALLOTERROR
#define ALLOTERROR ALLOTDIE
#endif

#define SIZE_BUFFER(TYPE,BUFFER,BUFSZ,NEEDED) \
    if ((NEEDED)>(BUFSZ)) { \
	if ((BUFSZ) != 0) { \
	    free(BUFFER) ; \
	} \
	callot(TYPE,BUFFER,NEEDED) ; \
	BUFSZ = NEEDED ; \
    } else { ; }

#define CSIZE_BUFFER(TYPE,BUFFER,BUFSZ,NEEDED) \
    if ((NEEDED)>(BUFSZ)) { \
	if ((BUFSZ) != 0) { \
	    free(BUFFER) ; \
	} \
	callot(TYPE,BUFFER,NEEDED) ; \
	BUFSZ = NEEDED ; \
    } else { ; }

#define RESIZE_BUFFER(TYPE,BUFFER,BUFSZ,NEEDED) \
    if ((NEEDED)>(BUFSZ)) { \
	if ((BUFSZ) != 0) { \
	    reallot(TYPE,BUFFER,NEEDED) ; \
	} else { \
	    callot(TYPE,BUFFER,NEEDED) ; \
	} \
	BUFSZ = NEEDED ; \
    } else { ; }


/* allotted can be used in an if statement when the action in allot is
 * inadequate */
#define allotted(type,ptr,size) \
  ( (ptr=(type) malloc((size_t)((size)*sizeof(*(ptr))))) != NULL )
#define callotted(type,ptr,size) \
  ( (ptr=(type) calloc((size_t)(size),sizeof(*(ptr))) ) != NULL )

/* allot checks the results of malloc and generates an error message in
 * the case of failures */
#define allot(type,ptr,size)    \
  STMT( if (!callotted(type,ptr,size)) ALLOTERROR(ptr,size);)
#define callot(type,ptr,size)    \
  STMT( if (!callotted(type,ptr,size)) ALLOTERROR(ptr,size);)

/* reallot and reallotted correspond to allot and allotted and are
 * interfaces to realloc. */
#define reallotted(type,ptr,size) \
  ( (ptr=(type) realloc((char *)ptr,(size_t)((size)*sizeof(*ptr)))) != NULL )
#define reallot(type,ptr,size)    \
  STMT( if (!reallotted(type,ptr,size)) ALLOTERROR(ptr,size);)


/* insist provides a concise method for consistency checks which are
 * never expect to fail, avoiding cluttering the code with branches which
 * are seldom if ever taken.
 * 
 * INSISTFAIL specifies what happens when an insist fails -- it may be
 * overridden with a special macro within a file */

/*
#ifndef INSISTFAIL
#define INSISTFAIL die
#endif

#define insist(ex)	STMT( if (!(ex)){ INSISTFAIL(1, "*stack*Unexpected failure: %s in file %s, line %d\n", __func__, __FILE__, __LINE__ ) ;} )
#define btwarn(ex)	STMT( if (!(ex)){ elog_debug(1, "*stack*Unexpected failure: %s in file %s, line %d\n", __func__, __FILE__, __LINE__ ) ;} )

#define insistmsg(ex,msg)	STMT( if (!(ex)){ INSISTFAIL(1, "%s: %s in file %s, line %d\n", msg, __func__, __FILE__, __LINE__ ) ;} )
#define btwarnmsg(ex,msg)	STMT( if (!(ex)){ elog_debug(1, "%s: file %s, %s in line %d\n", msg, __func__, __FILE__, __LINE__ ) ;} )
*/

/* The meaning of the following macros is fairly obvious; however, they
 * are all dangerous in the sense that their arguments will be evaluated
 * multiple times, which can lead to erroneous results. for example,
 * sqr(a++) will give an undefined result which is probably not what you
 * wanted and may differ among compilers. */

#define sqr(d) 		((d)*(d))

#ifndef SQR
#define SQR(d) 		((d)*(d))
#endif

#ifndef ABS
#define ABS(d) 		((d)< 0  ? -(d) : (d))
#endif

#ifndef MAX
#define MAX(a,b)        ((a)>(b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)        ((a)<(b) ? (a) : (b))
#endif

#ifndef SIGN
#define SIGN(d) 	(((d)>0)? 1.0 : ((d)<0) ? -1.0 : 0.0)
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* The following size is judged to be large enough to avoid problems with
 * overflowing string buffers, so that the various string routines which
 * do no checking can be used without much fear. */
#define STRSZ 1024

#endif


