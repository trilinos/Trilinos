/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ML memory management functions                                       */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : April, 1998                                          */
/* ******************************************************************** */

#ifndef __MLMEM__
#define __MLMEM__

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#define MAX_MALLOC_LOG 1000

#ifdef __cplusplus
extern "C"
{
#endif

extern  int  ML_memory_alloc( void **, int, char * );
extern  int  ML_memory_free( void ** );
extern  int  ML_memory_check_var(void *);
extern  int  ML_memory_inquire();
extern  int  ML_memory_inquire_short( int );
extern  int  ML_memory_clean( char *, int );

#ifdef __cplusplus
}
#endif

#endif

