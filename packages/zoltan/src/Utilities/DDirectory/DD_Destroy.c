/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"


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
