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
 *	Ultrix
 */

#include "os_unix.h"

#define AVAIL_X11               /* If the X11 Window System can work    */
#define HAS_ATRIGH              /* acosh( ), asinh( ), atanh( )         */
#define HAS_BCOPY               /* bcopy( ), bzero( )                   */
#define HAS_BSDRANDOM           /* srandom( ) and random( )             */
#define HAS_POSIXTTY            /* <termios.h>                          */
#define HAS_BSDDIRS             /* <sys/dir.h>                          */
#define HAS_BSDRUSAGE           /* getrusage( )                         */
#define HAS_BSDRLIMIT           /* getrlimit( )                         */
#define HAS_BSDSOCKETS          /* <net/inet.h>, socket( ), etc.        */
#define HAS_DUP2
#define HAS_GETPW
#define HAS_STRCHR              /* strchr( ) instead of index( )        */
#define HAS_LIMITS_H
#define HAS_FLOAT_H
#ifdef SHARED_MEM
#undef HAS_ISATTY
#endif

