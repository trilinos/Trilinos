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
Copyright 1991 Regents of the University of California.  All rights reserved.
**********/

/*
 *	SGI Irix 6.4
 */

#include "os_unix.h"

#define AVAIL_X11		/* If the X11 Window System can work    */
#define HAS_ATRIGH		/* acosh( ), asinh( ), atanh( )         */
#define HAS_BCOPY		/* bcopy( ), bzero( )			*/
#define HAS_BSDRANDOM		/* srandom( ) and random( )		*/
#define HAS_POSIXTTY		/* <termios.h>				*/
#define HAS_BSDDIRS		/* <sys/dir.h>				*/
#define HAS_BSDRUSAGE		/* getrusage( )				*/
#define HAS_BSDRLIMIT		/* getrlimit( )				*/
#define HAS_BSDSOCKETS		/* <net/inet.h>, socket( ), etc.	*/
#define HAS_DUP2
#define HAS_GETWD		/* getwd(buf)				*/
#define HAS_STRCHR              /* strchr( ) instead of index( )        */
#define HAS_LIMITS_H
#define HAS_FLOAT_H
#define DEV_bsim3
/* #define SIMPLE_INPUT */           /* Use fread rather than select etc., and allow
                                   interactive interface to operate from stdin */
#ifdef SHARED_MEM
#undef HAS_ISATTY
#endif
/*
*/
