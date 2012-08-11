

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
