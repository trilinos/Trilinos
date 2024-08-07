// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include "par_const.h"
#include "zz_const.h"


/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME                             TYPE
----------------------------------------------------------------------
	Zoltan_Print_Stats			void

******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_Print_Stats (MPI_Comm communicator, int debug_proc, double x, char *msg)
{
/****************************************************************/
/* Print max, sum, and imbalance for a variable over all procs. */
/****************************************************************/
  double sum, max;
  int proc, num_proc;

  MPI_Comm_rank(communicator, &proc);
  MPI_Comm_size(communicator, &num_proc);

  MPI_Reduce((void *)&x, (void *)&sum, 1, MPI_DOUBLE, MPI_SUM, debug_proc, 
             communicator);

  MPI_Reduce((void *)&x, (void *)&max, 1, MPI_DOUBLE, MPI_MAX, debug_proc, 
             communicator);

  if (proc == debug_proc) {
    if (sum <= 0.0)
      printf("%s: Max: %g, Sum: %g, Imbal.: N.A.\n",
              msg, max, sum);
    else /* sum > 0.0 */
      printf("%s: Max: %g, Sum: %g, Imbal.: %g\n",
              msg, max, sum, max*(num_proc)/sum);
  }

}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
