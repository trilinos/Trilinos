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
#include "lb_init_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_Migrate_Init(struct Zoltan_Migrate_Struct *mig)
{
  mig->Auto_Migrate = ZOLTAN_AUTO_MIGRATE_DEF;
  mig->Only_Proc_Changes = ZOLTAN_MIGRATE_ONLY_PROC_CHANGES_DEF;
  mig->Pre_Migrate_PP = NULL;
  mig->Mid_Migrate_PP = NULL;
  mig->Post_Migrate_PP = NULL;
  mig->Pre_Migrate = NULL;
  mig->Mid_Migrate = NULL;
  mig->Post_Migrate = NULL;
  mig->Pre_Migrate_PP_Fort = NULL;
  mig->Mid_Migrate_PP_Fort = NULL;
  mig->Post_Migrate_PP_Fort = NULL;
  mig->Pre_Migrate_Fort = NULL;
  mig->Mid_Migrate_Fort = NULL;
  mig->Post_Migrate_Fort = NULL;
  mig->Pre_Migrate_PP_Data = NULL;
  mig->Mid_Migrate_PP_Data = NULL;
  mig->Post_Migrate_PP_Data = NULL;
  mig->Pre_Migrate_Data = NULL;
  mig->Mid_Migrate_Data = NULL;
  mig->Post_Migrate_Data = NULL;
}

void Zoltan_LB_Init(struct Zoltan_LB_Struct *lb, int num_proc)
{
  int i;

  lb->Num_Global_Parts = num_proc;
  lb->Num_Global_Parts_Param = -1;
  lb->Num_Local_Parts_Param = -1;
  lb->Prev_Global_Parts_Param = -2;
  lb->Prev_Local_Parts_Param = -2;
  lb->Single_Proc_Per_Part = 1;
  lb->PartDist = NULL;
  lb->ProcDist = NULL;
  lb->Part_Info_Max_Len = 0;
  lb->Part_Info_Len = 0;
  lb->Part_Info = NULL;
  lb->Method = RCB;
  lb->LB_Fn = Zoltan_RCB;
  lb->Remap_Flag = 1;
  lb->Remap = NULL;
  lb->OldRemap = NULL;
  lb->Return_Lists = ZOLTAN_LB_RETURN_LISTS_DEF;
  lb->Uniform_Parts = 1;
  lb->Data_Structure = NULL;
  lb->Free_Structure = Zoltan_RCB_Free_Structure;
  lb->Copy_Structure = Zoltan_RCB_Copy_Structure;
  lb->Point_Assign = Zoltan_RB_Point_Assign;
  lb->Box_Assign = Zoltan_RB_Box_Assign;
  lb->Imb_Tol_Len = 10;
  lb->Imbalance_Tol = (float *)ZOLTAN_MALLOC((lb->Imb_Tol_Len)*sizeof(float));
  for (i=0; i<lb->Imb_Tol_Len; i++)
    lb->Imbalance_Tol[i] = ZOLTAN_LB_IMBALANCE_TOL_DEF;
  strcpy(lb->Approach, ZOLTAN_LB_APPROACH_DEF);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
