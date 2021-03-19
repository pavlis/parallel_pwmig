#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "pwmig/dsap/stock.h"
//#include "dsap_regex.h"
#include <errno.h>
#include <dirent.h>
/* Most of this file is commented out to get rid of Tbl
   dependency.   Some may need to be converted if we
   find a hidden dependency */


/*
static Tbl     *patterns = NULL;
typedef struct Replace {
    char           *pattern,
                   *replacement;
    struct re_pattern_buffer *pattern_buffer;
}               Replace;


static void free_replace (void *vreplace)
{
    Replace        *replace;
    replace = (Replace *) vreplace ;
    free (replace->pattern);
    free (replace->replacement);
    re_free (replace->pattern_buffer);
}

void setmapcf (char *cfname)
{
    FILE           *f;
    Replace        *replace;

    char 	    aline[STRSZ] ;
    char            pattern[STRSZ],
                   *s;
    char            replacement[STRSZ];
    struct re_pattern_buffer *pattern_buffer = NULL;

    re_syntax_options = RE_SYNTAX_EGREP;

    if (patterns != NULL)
	freetbl (patterns, free_replace);
    patterns = newtbl (0L);
    if (cfname == 0 || *cfname == 0)
	cfname = datafile ("MAPPATH", "wd.cf");

    if ((f = fopen (cfname, "r")) == NULL)
			return;
    while ( fgets(aline, STRSZ, f) != NULL ) {
	// skip comments
	for ( s=aline ; *s == ' ' ; s++ ) ;
	if ( *s == '#' )
	    continue ;

	if (sscanf (aline, "%s %s\n", pattern, replacement) != 2)
	    // skip bad lines
	    continue ;

	allot (struct re_pattern_buffer *, pattern_buffer, 1);
	pattern_buffer->buffer = 0;
	pattern_buffer->allocated = 0;
	pattern_buffer->translate = 0;
	pattern_buffer->regs_allocated = REGS_UNALLOCATED;
	pattern_buffer->fastmap = 0;

	if ((s = re_compile_pattern (pattern,
				     strlen (pattern), pattern_buffer)) != 0)
	    die (0, "canonical name pattern '%s' from %s did not compile.\n",
		 pattern, cfname);

	allot (Replace *, replace, 1);
	replace->pattern = strdup (pattern);
	replace->replacement = strdup (replacement);
	replace->pattern_buffer = pattern_buffer;
	pushtbl (patterns, replace);

    }
    fclose (f);
}

// * mappath - Map path from "pwd" name space to "canonical" name space.
int mappath (char *spath,char *upath)
{
    long             i,
                    n;
    Replace        *replace;
    if (patterns == NULL)
	setmapcf ((char *)NULL);
    n = maxtbl (patterns);
    for (i = 0; i < n; i++) {
	replace = (Replace *) gettbl (patterns, i);
	if (patsub (spath, replace->pattern_buffer, replace->replacement, upath))
	    return 0;
    }
    strcpy ( upath, spath )  ;
    return 0 ;
}




// * dirbase -- convert path to dirname, basename

void dirbase (char *path,char *dir,char *base)
{
    char           *cp;
    strcpy (dir, path);
    if ( (cp = strrchr (dir, '/')) ) {
	*cp++ = 0;
	strcpy (base, cp);
    } else {
	strcpy (dir, ".");
	strcpy (base, path);
    }
}


// * abspath -- convert path to absolute pathname

int abspath (char *relp,char *absp)
{
    char            dir[FILENAME_MAX],
                    base[FILENAME_MAX],
                    cwd[FILENAME_MAX],
                    ncwd[FILENAME_MAX];
    struct stat     statbuf;

    *absp = 0;
    if (*relp == '/')
	mappath (relp, absp);
    else {
	if (stat (relp, &statbuf) == -1
		|| S_ISDIR (statbuf.st_mode) == 0) {
	    dirbase (relp, dir, base);
	    if (stat (dir, &statbuf) == -1
		    || S_ISDIR (statbuf.st_mode) == 0) {
		register_error (0, "directory does not exist: %s\n", dir);
		return -1;
	    }
	} else {
	    strcpy (dir, relp);
	    *base = 0;
	}
	getcwd (cwd, FILENAME_MAX);
	if (chdir (dir)) {
	    register_error (1, "Can't cd to directory %s\n", dir);
	    return -1;
	}
	getcwd (ncwd, FILENAME_MAX);
	chdir (cwd);
	mappath (ncwd, absp);
	if (*base != 0) {
	    strcat (absp, "/");
	    strcat (absp, base);
	}
    }
    return 0;
}

*/
int makedir (const char *dir)
{
    struct stat     statbuf;
    char           *s;

    if (stat (dir, &statbuf) == -1 && errno == ENOENT)
	mkdir (dir, 0775);
    if (stat (dir, &statbuf) == -1
	    && errno == ENOENT
	    && (s = strrchr (dir, '/')) != NULL) {
	char            parent[FILENAME_MAX];
	int             n;

	n = s - dir;
	strncpy (parent, dir, n);
	parent[n] = 0;

	makedir (parent);
	mkdir (dir, 0775);
    }
    if (stat (dir, &statbuf) == -1
	    || S_ISDIR (statbuf.st_mode) == 0
	    || access (dir, W_OK) != 0) {
        fprintf(stderr,"Can't create writable directory %s\n", dir);
	return -1;
    }
    return 0;
}

/* $Id: wd.c,v 1.2 1997/07/21 00:41:31 danq Exp $ */
