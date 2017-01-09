/* ========================================================================= */
/* === amesos_amd_global ========================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

#include <stdlib.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

#ifndef NULL
#define NULL 0
#endif

/* ========================================================================= */
/* === Default AMD memory manager ========================================== */
/* ========================================================================= */

/* The user can redefine these global pointers at run-time to change the memory
 * manager used by AMD.  AMD only uses malloc and free; realloc and calloc are
 * include for completeness, in case another package wants to use the same
 * memory manager as AMD.
 *
 * If compiling as a MATLAB mexFunction, the default memory manager is mxMalloc.
 * You can also compile AMD as a standard ANSI-C library and link a mexFunction
 * against it, and then redefine these pointers at run-time, in your
 * mexFunction.
 *
 * If -DNMALLOC is defined at compile-time, no memory manager is specified at
 * compile-time.  You must then define these functions at run-time, before
 * calling AMD, for AMD to work properly.
 */

#ifndef NMALLOC
#ifdef MATLAB_MEX_FILE
/* MATLAB mexFunction: */
void *(*amesos_amd_malloc) (size_t) = mxMalloc ;
void (*amesos_amd_free) (void *) = mxFree ;
void *(*amesos_amd_realloc) (void *, size_t) = mxRealloc ;
void *(*amesos_amd_calloc) (size_t, size_t) = mxCalloc ;
#else
/* standard ANSI-C: */
void *(*amesos_amd_malloc) (size_t) = malloc ;
void (*amesos_amd_free) (void *) = free ;
void *(*amesos_amd_realloc) (void *, size_t) = realloc ;
void *(*amesos_amd_calloc) (size_t, size_t) = calloc ;
#endif
#else
/* no memory manager defined at compile-time; you MUST define one at run-time */
void *(*amesos_amd_malloc) (size_t) = NULL ;
void (*amesos_amd_free) (void *) = NULL ;
void *(*amesos_amd_realloc) (void *, size_t) = NULL ;
void *(*amesos_amd_calloc) (size_t, size_t) = NULL ;
#endif

/* ========================================================================= */
/* === Default AMD printf routine ========================================== */
/* ========================================================================= */

/* The user can redefine this global pointer at run-time to change the printf
 * routine used by AMD.  If NULL, no printing occurs.  
 *
 * If -DNPRINT is defined at compile-time, stdio.h is not included.  Printing
 * can then be enabled at run-time by setting amesos_amd_printf to a non-NULL function.
 */

#ifndef NPRINT
#ifdef MATLAB_MEX_FILE
int (*amesos_amd_printf) (const char *, ...) = mexPrintf ;
#else
#include <stdio.h>
int (*amesos_amd_printf) (const char *, ...) = printf ;
#endif
#else
int (*amesos_amd_printf) (const char *, ...) = NULL ;
#endif
