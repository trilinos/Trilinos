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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"

/****************************************************************************/
int Zoltan_LB_Box_Assign (
 ZZ *zz,
 double xlo,
 double ylo,
 double zlo,
 double xhi,
 double yhi,
 double zhi,
 int *procs,
 int *count)
{
  char *yo = "Zoltan_LB_Box_Assign";
  int tmp = 0;

  if (zz->LB.Box_Assign == NULL) {
    /* function not supported by current decomposition method */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Box_Assign not supported by chosen partitioning method.");
    return ZOLTAN_FATAL;  
  }

  if (zz->LB.PartDist != NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Non-uniform distribution of partitions over processors is specified; "
      "use Zoltan_LB_Box_PP_Assign.");
    return ZOLTAN_FATAL;
  }

  /* Call appropriate method.  Pass procs and count in partition arguments
   * for greater efficiency in LB.Box_Assign (Zoltan is partition-based.) */
  return zz->LB.Box_Assign(zz, xlo, ylo, zlo, xhi, yhi, zhi, NULL, &tmp, 
                           procs, count);
}

/****************************************************************************/
int Zoltan_LB_Box_PP_Assign (
 ZZ *zz,
 double xlo,
 double ylo,
 double zlo,
 double xhi,
 double yhi,
 double zhi,
 int *procs,
 int *proc_count,
 int *parts,
 int *part_count)
{
  char *yo = "Zoltan_LB_Box_PP_Assign";

  if (zz->LB.Box_Assign == NULL) {
    /* function not supported by current decomposition method */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Box_Assign not supported by chosen partitioning method.");
    return ZOLTAN_FATAL;  
  }

  /* Call appropriate method.  Pass procs and count in partition arguments
   * for greater efficiency in LB.Box_Assign (Zoltan is partition-based.) */
  return zz->LB.Box_Assign(zz, xlo, ylo, zlo, xhi, yhi, zhi, procs, proc_count,
                           parts, part_count);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
