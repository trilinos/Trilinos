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

/*
 *  Generic UNIX
 */

#define HAS_STAT		/* stat( ) returns info on files (inode) */
#define HAS_QSORT		/* qsort( )				*/
#define HAS_FCNTL
#define HAS_ISATTY		/* isatty( )				*/
#define HAS_GETPW		/* getpwuid( ), etc.			*/
#define HAS_ENVIRON		/* getenv( )				*/
#define HAS_ACCESS		/* access( )				*/
#define HAS_CTYPE		/* <ctype.h>, iswhite( ), etc.		*/
#define HAS_LONGJUMP		/* setjmp( ), longjmp( )		*/
#define HAS_CHDIR		/* for tree filesystems, chdir( )	*/
#define     DIR_PATHSEP		"/"
#define     DIR_TERM		'/'
#define     DIR_CWD		"."
#define HAS_UNIX_SIGS		/* signal( ), kill( )			*/
#define HAS_UNLINK		/* unlink( ), for removing files	*/
#define HAS_GETPID		/* getpid( ) to identify processes	*/
#define HAS_WAIT		/* wait( ) wait for processes		*/
#define HAS_POPEN		/* popen( ), pipe through shell command	*/
#define HAS_SYSTEM		/* system( ), execute system command	*/
#define     SYSTEM_PLOT5LPR	"lpr -P%s -g %s"
#define     SYSTEM_PSLPR	"lpr -P%s %s"
#define     SYSTEM_MAIL		"Mail -s \"%s (%s) Bug Report\" %s"
				/* simulator, version, address		*/
#define HAS_VPERROR		/* perror( )				*/
#define HAS_ASCII		/* eigth bit of a character is not used	*/
#define HAS_CLEARERR		/* clearerr( ), should be in stdio	*/
#define HAS_LOCALTIME		/* localtime( ), asctime( ), time( )    */
#define HAS_GETCWD

#define TEMPFORMAT	"/tmp/%s%d"
