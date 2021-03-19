#ifndef __deviants__
#define __deviants__

#include "deviants_include.h"

#include <string.h>
#include <stdio.h>
#include <netdb.h>

typedef double STDARG_REAL ;

typedef struct Std_statvfs {
     unsigned long f_bsize;   /* preferred file system block size */
     unsigned long f_frsize;  /* fundamental filesystem block(size if supported) */
     unsigned long  f_blocks;  /* total # of blocks on file system in units of f_frsize */
     unsigned long  f_bfree;   /* total # of free blocks */
     unsigned long  f_bavail;  /* # of free blocks avail to non-super-user */
     unsigned long  f_files;   /* total # of file nodes (inodes) */
     unsigned long  f_ffree;   /* total # of free file nodes */
     unsigned long  f_favail;  /* # of inodes avail to non-super-user*/
     unsigned long f_fsid;    /* file system id (dev for now) */
} Std_statvfs ;

#define CTIME_R(t,b) 	ctime_r((t),(b))

#define btsort	psort_r		// use Apple built in threaded sort

#ifdef __APPLE_API_PRIVATE

/*
 * Structure per mounted file system.  Each mounted file system has an
 * array of operations and an instance record.  The file systems are
 * put on a doubly linked list.
*/

#endif

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */

#ifdef __BIG_ENDIAN__
#define WORDS_BIGENDIAN 1
#else /* !__BIG_ENDIAN__ */
   /* #undef WORDS_BIGENDIAN */
#endif /* __BIG_ENDIAN__ */


#include <sys/param.h>

#if 0
#define SIGPWR 30	/* actually USR1 for Darwin -- SIGPWR is undefined */
#define SIGPOLL SIGEMT  
#define SIGCLD SIGCHLD
#endif

#define sigset signal

#define fdatasync(FD)  (0)

#define MAX_IPS 32
#define now std_now

#ifdef  __cplusplus
extern "C" {
#endif

// extern void btsort (void *bot, size_t nmemb, size_t size, void *pvt, int (*compar) (void *, const void *, const void *)) ;

extern char *cuserid(char *name) ;

extern char * getTZ ( void );
extern void setTZ() ;
extern char *strsignal(int sig) ;

extern int sig2str ( int sig, char *name );
extern int str2sig ( char *str, int *sig );

extern off_t tell ( int fd );
extern char * getexecname ( void );

extern FILE *ZOPEN ( char *filename );
#define zopen ZOPEN

extern double htond ( double val );
extern double ntohd ( double val );
extern void htondp ( double * valp_from, double * valp_to );
extern void ntohdp ( double * valp_from, double * valp_to );
extern float htonf ( float val );
extern float ntohf ( float val );
extern void htonfp ( float * valp_from, float * valp_to );
extern void ntohfp ( float * valp_from, float * valp_to );

extern void swap2 ( void *from, void *to, int n );
extern void swap4 ( void *from, void *to, int n );
extern void swap8 ( void *from, void *to, int n );
extern void N32H4 ( int *to, char *from, int n );
extern void H42N3 ( char *to, int *from, int n );
extern void NF2HD ( double *to, char *from, int n );
extern void HD2NF ( char *to, double *from, int n );
extern void HI2NC ( char *to, int *from, int n );
extern void HI2NS ( char *to, int *from, int n );
extern void NC2HI ( int *to, char *from, int n );
extern void NS2HI ( int *to, char *from, int n );
extern void HI2NF ( char *to, int *from, int n );
extern void NF2HI ( int *to, char *from, int n );

#define DLUNIT  3
#define init_Dlclock() _init_Dlclock(DLUNIT)
extern void _init_Dlclock(int unit)  ;
extern void set_Dlclock(double time) ; 

extern double std_now ( void );
extern double now_ ( void );

extern int is_nfs ( char *path );

extern void * nomem ( void *a, int n );

#include <time.h>
extern int std_set_timer ( void (*alarm_handler)(void), double how_long );

extern int my_hostname ( char *name );
extern int my_os ( char *name );
extern int my_hardware ( char *name );
extern int my_ip ( char *hostname, char *ipc, int *ip );
extern int ip2name ( int addr, char *name );
extern int my_username ( char *name );

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

extern int name2ip ( char *name, struct in_addr *addr, char *ipc );

extern char *mcinfo (char **);
extern int my_ips ( char *ips, int size );
extern int my_ipv6s ( struct sockaddr_storage *sa, int size );
extern int my_ether (char *info) ;

extern void reserve_FILE(void) ;
extern void release_FILE(void) ;

extern int gxpipe ( char *cmd, char **env, int *cmdin );
extern int gxpipex ( char *cmd, int *cmdin );
extern int gxopen ( char *filename );
extern int gxline ( int fd, char *aline, int maxchars );
extern int gxtell ( int fd );
extern int gxseek ( int fd, int p );
extern int gxclose ( int fd );

extern int std_statvfs(const char *path, Std_statvfs *buf);

extern int pid_exists ( int pid );
extern int is_timezone ( char *name );

void * std_dlopen (const char *lib) ;

extern void show_stack_trace ( void );
extern void show_stack_trace_string ( int strsize, char *str );
extern void show_fault ( int sig, siginfo_t *vsi, void *vuap );

#ifdef  __cplusplus
}
#endif

#endif
