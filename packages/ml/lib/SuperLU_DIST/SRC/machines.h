/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * These macros define which machine will be used.
 */

#ifndef __SUPERLU_MACHINES /* allow multiple inclusions */
#define __SUPERLU_MACHINES

#define SGI	        0
#define ORIGIN	        1
#define DEC	        2
#define CRAY_T3E	3
#define SUN             4
#define PTHREAD         5
#define IBM             6

#ifdef _SGI
#define MACH SGI 
#endif

#ifdef _ORIGIN
#define MACH ORIGIN 
#endif

#ifdef _DEC
#define MACH DEC 
#endif

#ifdef _CRAY
#define MACH CRAY_T3E 
#endif

#ifdef _SOLARIS
#define MACH SUN 
#endif

#ifdef _PTHREAD
#define MACH PTHREAD
#endif

#if ( defined(_SP2) || defined(_SP) )
#define MACH IBM
#endif

#endif /* __SUPERLU_MACHINES */
