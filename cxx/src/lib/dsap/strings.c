
#include <string.h>
#include "stock.h"

int blank ( char *s )
{
	/* This is not Dan Quinlan's original code.  This procedure is supposed to 
	return nonzero if s is all blank characters. He did this a very obscure
	way so I just rewrote it.   There rest of the file retains the odd method
	as it was known to work*/
	size_t n=strlen(s);
	size_t i;
	for(i=0;i<n;++i)
		if(s[i]!=' ') return(0);
    return(1);
}

int whitespace ( char *s )
{
    while ( *s == ' ' || *s == '\t' || *s == '\n' ) s++;
    return *s == 0 ;
}

/* Remove this function to avoid Tbl */
/*
Tbl *split (char *s,char *c ) 
{
    Tbl *tbl ; 

    tbl = newtbl(5) ;

    while ( *s != '\0' )
	{
	while ( *s == *c ) s++ ;
	if ( *s != '\0' ) 
	    pushtbl ( tbl, s ) ; 
	while ( *s != *c && *s != '\0' ) 
	    s++ ; 
	if ( *s == *c ) *s++ = '\0' ; 
	}
    
    return tbl ; 
}
*/

/* copy source to dest, padding with blanks, no trailing zero */
void sncopy(char *dest,const char *source,int n)
{
    int             i = 0;

    while (i++ < n && *source != '\0')
	*dest++ = *source++;
    for (--i ; i < n; i++)
	*dest++ = ' ';
}

/* copy source to dest, add trailing zero */
void szcopy(char *dest,const char *source,int n)
{
    int             i = 0;

    for (i = 0; i < n; i++)
	*dest++ = *source++;
    *dest = '\0';
}

/* copy source to dest, add trailing zero, strip leading, trailing blanks */
void copystrip(char *dest,const char *source,int n)
{
    int i, start, last ; 

    for (start=0 ; start<n && source[start] == ' ' ; start++ ) 
	;
    for (last=n-1 ; last>=0 && source[last] == ' ' ; last-- ) 
	;
    last++ ;

    for(i=start ; i<last ; i++ ) 
	*(dest++) = source[i] ; 

    *dest = '\0';
}

void *memdup ( void *a, int n ) 
{
    char *b ;
    allot ( char *, b, n ) ; 
    memcpy(b,a,n) ; 
    return (void *) b ;
}

