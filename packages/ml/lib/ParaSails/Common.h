/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * Common.h header file - common definitions and parameters; also hides
 * HYPRE-specific definitions
 *
 *****************************************************************************/

#include <stdio.h>

#if 0 /* HYPRE */
#include "HYPRE_config.h"
#include "utilities.h"
#include "fortran.h"
#ifdef HYPRE_USING_ESSL
#define ESSL
#endif
#else /* not HYPRE */
#include "mpi.h"
#define hypre_F90_NAME(name) name##_
#endif

#ifndef _COMMON_H
#define _COMMON_H

#define PARASAILS_MAXLEN  300000 /* maximum nz in a pattern - can grow */
#define PARASAILS_NROWS   300000 /* maximum rows stored per proc - can grow */

#ifndef ABS
#define ABS(x) (((x)<0)?(-(x)):(x))
#endif
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define PARASAILS_EXIT              \
{                                   \
   fprintf(stderr, "Exiting...\n"); \
   fflush(NULL);                    \
   MPI_Abort(MPI_COMM_WORLD, -1);   \
}

#endif /* _COMMON_H */
