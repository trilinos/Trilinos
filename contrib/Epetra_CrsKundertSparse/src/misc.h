/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#ifndef MISC_H
#define MISC_H

#define BSIZE_SP      512

#ifdef HAS_EXIT1
#  define EXIT_NORMAL 1
#  define EXIT_BAD    0
#else
#  define EXIT_NORMAL 0
#  define EXIT_BAD    1
#endif

#ifdef HAS_CTYPE
#  ifndef isalpha
#    include <ctype.h>
#  endif
#endif

#define eq(a,b)  (!strcmp((a), (b)))
#define eqc(a,b)  (cieq((a), (b)))
#define isalphanum(c)   (isalpha(c) || isdigit(c))
#define hexnum(c) ((((c) >= '0') && ((c) <= '9')) ? ((c) - '0') : ((((c) >= \
        'a') && ((c) <= 'f')) ? ((c) - 'a' + 10) : ((((c) >= 'A') && \
        ((c) <= 'F')) ? ((c) - 'A' + 10) : 0)))

#include "strext.h"

extern char *tmalloc();
extern char *trealloc();
extern void txfree();

#define tfree(x)	(txfree(x), x = 0)

#define	alloc(TYPE)	((TYPE *) tmalloc(sizeof(TYPE)))

extern char *copy();
extern char *gettok();
extern void appendc();
extern int scannum();
extern int prefix();
extern int ciprefix();
extern int cieq();
extern void strtolower();
extern int substring();
extern char *tilde_expand( );
extern void cp_printword();

extern char *datestring();
extern char *date_only();
extern char *time_only();
extern double seconds();

extern char *smktemp();

/* Externs from libc */

#ifdef HAS_STDLIB

#  ifndef _STDLIB_INCLUDED
#    define _STDLIB_INCLUDED
#    include <stdlib.h>
#  endif
#  ifndef HAS_BSDRAND
#    define random	rand
#    define srandom	srand
#  endif
#  ifdef HAS_DOSDIRS
#include <unistd.h>
#  endif

#else

#  ifdef HAS_BSDRAND
extern long random();
extern void srandom();
#  else
#    define random	rand
#    define srandom	srand
#  endif

extern void *calloc();
extern void *malloc();
extern void *realloc();
extern char *getenv();
extern int errno;
#if !defined(LINUX) && !defined(FREEBSD) && !defined(APPLEOSX)
extern char *sys_errlist[];
#endif
extern char *getenv();
extern char *getwd();
extern int rand();
extern void srand();
extern int atoi();
extern int kill();
extern int getpid();
extern void qsort();
#  ifdef notdef
extern void exit();
#  endif

#  ifdef HAS_GETCWD
char *getcwd(char *, size_t);
#  endif

#  ifdef HAS_CLEARERR
#    ifndef clearerr
extern void clearerr();
#    endif /* clearerr */
#  endif /* HAS_CLEARERR */

#if defined(LINUX) || defined(FREEBSD) || defined(NMcplant) || defined(DEC) || defined(APPLEOSX)
extern void bzero();
#else
#ifndef bzero
#ifndef dec
extern int bzero();
#  endif
# endif
#endif

#  ifndef bcopy
extern void bcopy();
#  endif

#  ifndef index
#    ifdef HAS_INDEX
extern char *rindex();
extern char *index();
#    else
#      ifdef HAS_STRCHR
/* For some strange reason these lines screw up the compile on linux:
extern char *strchr();
extern char *strrchr();
*/
#      else
#      endif
#    endif
#  endif

#endif	/* else STDLIB */

#ifndef HAS_INDEX
#  ifndef index
#    ifdef HAS_STRCHR
#      define	index	strchr
#      define	rindex	strrchr
#    endif
#  endif
#endif

#ifdef HAS_VPERROR
extern void perror();
#endif

#ifdef HAS_TIME_
#  ifdef HAS_BSDTIME
extern char *timezone();
#  endif
extern char *asctime();
extern struct tm *localtime();
#endif

#ifndef HAS_MEMAVL
#  ifdef HAS_RLIMIT_
extern char *sbrk();
#  endif
#endif

#define false 0
#define true 1

#ifdef HAS_DOSDIRS
typedef	int	*DIR;
struct direct {
	int	d_reclen;
	short	d_ino;
	short	d_namelen;
	char	d_name[20];
	};

#  ifdef __STDC__
extern DIR *opendir(char *);
extern struct direct *readdir(DIR *);
#  else
extern DIR *opendir( );
extern struct direct *readdir( );
#  endif

#endif

#endif /* MISC_H */
