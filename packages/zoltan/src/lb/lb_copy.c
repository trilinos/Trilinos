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

#define COPY_BUFFER(buf, type, num) \
  if (from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((num) * sizeof(type)); \
    if (!to->buf) { \
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory."); \
      Zoltan_LB_Free_Struct(to); \
      return ZOLTAN_MEMERR; \
    } \
    memcpy(to->buf, from->buf, (num) * sizeof(type)); \
  } \
  else { \
    to->buf = NULL; \
  }

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for copying LB and Migrate structures.
 *  The target of the copy should be a valid structure.
 *  These routines should be called only by Zoltan.
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_LB_Copy_Struct(ZZ *toZZ, ZZ const *fromZZ)
{
char *yo = "Zoltan_LB_Copy_Struct";
int proc = fromZZ->Proc;

  struct Zoltan_LB_Struct *to = &(toZZ->LB);
  struct Zoltan_LB_Struct const *from = &(fromZZ->LB);

  Zoltan_LB_Free_Struct(&(toZZ->LB));

  if (!from){
    return ZOLTAN_OK;
  }

  *to = *from;

  COPY_BUFFER(Part_Info, struct Zoltan_part_info, to->Part_Info_Max_Len);

  COPY_BUFFER(Remap, int, to->Num_Global_Parts);

  COPY_BUFFER(OldRemap, int, to->Num_Global_Parts);

  COPY_BUFFER(PartDist, int, to->Num_Global_Parts + 1);

  COPY_BUFFER(ProcDist, int, fromZZ->Num_Proc + 1);

  COPY_BUFFER(Imbalance_Tol, float, to->Imb_Tol_Len);

  if (from->Data_Structure) {
    to->Data_Structure = NULL;
    if (!from->Copy_Structure)
      {
      /* 
       * Some Zoltan codes don't save their Data_Structure after
       * partitioning.  However if they do, they must define a
       * copy function.
       */
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "A Copy function must be defined");
      return ZOLTAN_WARN;
      }
    from->Copy_Structure(toZZ, fromZZ);
  } 

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
