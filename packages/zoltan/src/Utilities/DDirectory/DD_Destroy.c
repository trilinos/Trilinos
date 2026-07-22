// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include <stdio.h>
#include <stdlib.h>

#include "zoltan_dd_const.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif



/**********************  Zoltan_DD_Destroy()  *********************/

void Zoltan_DD_Destroy (Zoltan_DD_Directory **dd)
{
   char *yo = "ZOLTAN_DD_Destroy";

   /* input sanity check */
   if (dd == NULL || *dd == NULL) {
      ZOLTAN_PRINT_ERROR (0, yo, "Input argument dd is NULL");
      return;
   }
   if ((*dd)->debug_level > 4)
      ZOLTAN_TRACE_IN ((*dd)->my_proc, yo, NULL);

   ZOLTAN_FREE(&((*dd)->nodelist));
   ZOLTAN_FREE(&((*dd)->nodedata));

   /* execute user registered cleanup function, if needed */
   if ((*dd)->cleanup != NULL)
       (*dd)->cleanup((*dd)->hashdata);

   MPI_Comm_free (&((*dd)->comm));    /* free MPI Comm, ignore errors */

   if ((*dd)->debug_level > 4)
      ZOLTAN_TRACE_OUT ((*dd)->my_proc, yo, NULL);

   ZOLTAN_FREE (dd);                  /* free directory structure     */
   return;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
