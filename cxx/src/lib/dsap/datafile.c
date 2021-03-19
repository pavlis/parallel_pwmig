/* This file is a heavily modified version of a file from the open source
   datascope library.  Early in BRTT's history the released this open
   source version.  Here is their copyright notice:  */
/*  Copyright (c) 1997 Boulder Real Time Technologies, Inc.           */
/*                                                                    */
/*  This software module is wholly owned by Boulder Real Time         */
/*  Technologies, Inc. Any use of this software module without        */
/*  express written permission from Boulder Real Time Technologies,   */
/*  Inc. is prohibited.                                               */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "pwmig/dsap/stock.h"

//static char PWMIGpath[FILENAME_MAX] ;
/* This was the old version of this.   Changed to a home directory to be
less restrictive.
static char *Default_pwmigpath = "/opt/pwmig" ;
*/

char *getpwmig()
{
    char *pwmigpath ;
    const char *defaultpath="~/pwmig";

    if ( (pwmigpath = getenv ( "PWMIG" )) != 0 ) {
	return pwmigpath ;
    }
    else
        return strdup(defaultpath);
}
    /*Original had all this stuff - removed that functionality to get rid of elog
      I am retaining it during debug in case I need to restore some of this
      functionality.  Eventually this debris should be removed. */

//    elog_query ( ELOG_ARGV, 0, (void *) &argv ) ;
//    if ( argv == 0 ) { /* elog not initialized, give up */
//	return Default_pwmigpath ;
//    } else {
//	me = argv[0] ;
//    }
//
//    if ( strchr(me, '/') == 0 ) {
//	/* bare name -- have to find along PATH */
//	char *path ;
//	if ( (path = datapath ( "PATH", 0, me, "" )) == 0 ) {
//	    return Default_pwmigpath ; /* can't figure out path to program */
//	} else {
//	    strncpy ( PWMIGpath, path, FILENAME_MAX ) ;
//	    free ( path ) ;
//	}
//    } else {
//      strcpy ( PWMIGpath, me ) ;
//      slash = strrchr(PWMIGpath, '/' ) ;
//      *slash = 0 ;
//      if ( *me != '/' ) {
//	  char *old, *pwd ;
//	  /* relative path! */
//	  old = getcwd(0, FILENAME_MAX) ;
//	  if ( chdir(PWMIGpath) != 0 ) {
//	    return Default_pwmigpath ; /* die ( 1, "Can't cd to bin directory '%s'\n", PWMIGpath ) ; */
//	  }
//	  pwd = getcwd(0, FILENAME_MAX);
//	  strcpy ( PWMIGpath, pwd ) ;
//	  free ( pwd ) ;
//	  if ( chdir(old) != 0 ) {
//	    die ( 1, "getpwmig can't return to original directory '%s'\n", old ) ;
//	  }
//	  free ( old ) ;
//      }
//    }
//
//    n = strlen(PWMIGpath) ;
//    if (   PWMIGpath[n-4] != '/'
//	|| PWMIGpath[n-3] != 'b'
//	|| PWMIGpath[n-2] != 'i'
//	|| PWMIGpath[n-1] != 'n' ) {
//	return Default_pwmigpath ; /*
//	    die ( 0, "bad configuration: '%s' has no final bin component.\n",
//	    PWMIGpath ) ; */
//	}
//    PWMIGpath[n-4] = 0 ;
//
//    return PWMIGpath;
//}

char *datapath (const char *envname,const char *dirname,
        const char *filename,const char *suffix)
{
    char            local[FILENAME_MAX];

    local[0] = '\0';
    if (dirname != NULL && *dirname != '\0') {
	strcpy (local, dirname);
	strcat (local, "/");
    }
    strcat (local, filename);
    if (suffix != NULL && *suffix != '\0') {
	strcat (local, ".");
	strcat (local, suffix);
    }
    return datafile (envname, local);
}


char *datafile (const char *envname,const char *filename)
{
    char           *path;
    char            local[FILENAME_MAX];
    struct stat     statbuf;
    char           *basename,
                   *save;

    basename = strrchr (filename, '/');
    if (basename)
	basename++;
    /* This chain of conditionals strikes me as a tad
       dangerous, but known to work.   To work it depends on
       an assumption about how this is executed that is pretty
       subtle.   That is, it requires the set of conditionals
       evaluate left to right and as soon as we get a false it
       stops trying.  */
    if (envname != NULL
	    && *envname != '\0'
	    && (path = getenv (envname)) != NULL) {
/* Original code allowed a colon separated list of directories
   to search for a file.   I am disabling that functionality as
   it is not needed for this implementation and the code (below)
   requires an datascope tbl list container that I do not want
   to convert for pwmig */
	//if (strchr (path, ':') == NULL) {
	    if (stat (path, &statbuf) == 0) {
		if (S_ISDIR (statbuf.st_mode)) {
		    strcpy (local, path);
		    strcat (local, "/");
		    save = local + strlen (local);
		    strcat (local, filename);
		    if (access (local, R_OK) == 0)
			return strdup (local);
		    if (basename) {
			*save = '\0';
			strcat (local, basename);
			if (access (local, R_OK) == 0)
			    return strdup (local);
		    }
		} else if (access (path, R_OK) == 0)
		    return strdup (path);
	    }
            /* Top of : list search logic.   Commented out for
               now.  When sure this is ok should delete this. */
            /*
	} else {
	    localpath = strdup (path);
	    pathtbl = split (localpath, ':');
	    n = maxtbl (pathtbl);
	    for (i = 0; i < n; i++) {
		strcpy (local, gettbl (pathtbl, i));
		if (stat (local, &statbuf) == 0) {
		    if (S_ISDIR (statbuf.st_mode)) {
			strcat (local, "/");
			save = local + strlen (local);
			strcat (local, filename);
			if (access (local, R_OK) == 0) {
			    freetbl (pathtbl, 0);
			    free(localpath) ;
			    return strdup (local);
			}
			if (basename) {
			    *save = 0;
			    strcat (local, basename);
			    if (access (local, R_OK) == 0) {
				freetbl (pathtbl, 0);
				free(localpath) ;
				return strdup (local);
			    }
			}
		    } else {
			if (access (local, R_OK) == 0) {
			    freetbl (pathtbl, 0);
			    free(localpath) ;
			    return strdup (path);
			}
		    }
		}
	    }
	    freetbl (pathtbl, 0);
	    free (localpath);
	}
        */
    }

    path = getpwmig() ;

    strcpy (local, path);
    strcat (local, "/data/");
    strcat (local, filename);
    if (access (local, R_OK) == 0)
	return strdup (local);

    return NULL;
}

/* $Id: datafile.c,v 1.4 1998/02/16 21:26:13 danq Exp $ */
