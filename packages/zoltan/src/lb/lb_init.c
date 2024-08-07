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
  lb->Serialize_Structure_Size = Zoltan_RCB_Serialize_Structure_Size;
  lb->Serialize_Structure = Zoltan_RCB_Serialize_Structure;
  lb->Deserialize_Structure = Zoltan_RCB_Deserialize_Structure;
  lb->Point_Assign = Zoltan_RB_Point_Assign;
  lb->Box_Assign = Zoltan_RB_Box_Assign;
  lb->Imb_Tol_Len = 10;
  lb->Imbalance_Tol = (float *)ZOLTAN_MALLOC((lb->Imb_Tol_Len)*sizeof(float));
  for (i=0; i<lb->Imb_Tol_Len; i++)
    lb->Imbalance_Tol[i] = ZOLTAN_LB_IMBALANCE_TOL_DEF;
  strcpy(lb->Approach, ZOLTAN_LB_APPROACH_DEF);
}

/*****************************************************************************/
/*****************************************************************************/

size_t Zoltan_LB_Serialize_Size(struct Zoltan_Struct const *zz) 
{
  size_t bufSize = 0;
  const struct Zoltan_LB_Struct *lb = &(zz->LB);

  /* Copy 12 integers from zz->LB */
  bufSize += 12 * sizeof(int);

  /* LB_Method_Name */
  bufSize += MAX_PARAM_STRING_LEN;

  /* Part_Info array */
  bufSize += lb->Part_Info_Len * sizeof(struct Zoltan_part_info);

  /* Imbalance_Tol array */
  bufSize += lb->Imb_Tol_Len * sizeof(float);

  /* lb->Approach */
  bufSize += MAX_PARAM_STRING_LEN;

  /* Remap array */
  bufSize += sizeof(int);
  if (lb->Remap != NULL)
    bufSize += lb->Num_Global_Parts * sizeof(int);

  /* Method specific data */
  if (lb->Serialize_Structure_Size != NULL) 
    bufSize += lb->Serialize_Structure_Size(zz);

  return bufSize;
}

/*****************************************************************************/
void Zoltan_LB_Serialize(struct Zoltan_Struct const *zz, char **buf)
{
  char *bufptr = *buf;
  const struct Zoltan_LB_Struct *lb = &(zz->LB);

  /* Copy 12 integers; if add more, update Zoltan_LB_Serialize_Size */
  int *intptr = (int *) bufptr;
  *intptr = lb->Num_Global_Parts; intptr++;
  *intptr = lb->Num_Global_Parts_Param; intptr++;
  *intptr = lb->Num_Local_Parts_Param; intptr++;
  *intptr = lb->Prev_Global_Parts_Param; intptr++;
  *intptr = lb->Prev_Local_Parts_Param; intptr++;
  *intptr = lb->Single_Proc_Per_Part; intptr++;
  *intptr = lb->Part_Info_Max_Len; intptr++;
  *intptr = lb->Part_Info_Len; intptr++;
  *intptr = lb->Remap_Flag; intptr++;
  *intptr = lb->Return_Lists; intptr++;
  *intptr = lb->Uniform_Parts; intptr++;
  *intptr = lb->Imb_Tol_Len; intptr++;
  bufptr = (char *) intptr;

  /* Copy LB_Method name */
  strcpy(bufptr, lb->Method_Name);
  bufptr += MAX_PARAM_STRING_LEN;
   
  /* Copy Part_Info */
  if (lb->Part_Info_Len) {
    size_t tmpSize = lb->Part_Info_Len * sizeof(struct Zoltan_part_info);
    memcpy(bufptr, (char *)(lb->Part_Info), tmpSize);
    bufptr += tmpSize;
  }

  /* Copy Imbalance_Tol */
  if (lb->Imb_Tol_Len) {
    size_t tmpSize = lb->Imb_Tol_Len * sizeof(float);
    memcpy(bufptr, (char *)(lb->Imbalance_Tol), tmpSize);
    bufptr += tmpSize;
  }

  /* Copy lb->Approach */
  strcpy(bufptr, lb->Approach);
  bufptr += MAX_PARAM_STRING_LEN;
 
  /* Copy Remap array; needed by Point_Assign  */
  if (lb->Remap != NULL) {
    int *intptr = (int *)bufptr;
    *intptr = 1;  // Sending Remap data
    bufptr += sizeof(int);
    int nbytes = lb->Num_Global_Parts * sizeof(int);
    memcpy(bufptr, lb->Remap, nbytes);
    bufptr += nbytes;
  }
  else {
    int *intptr = (int *)bufptr;
    *intptr = 0;  // Not sending Remap data
    bufptr += sizeof(int);
  }

  /* Serialize the method-specific load balancing data; advance bufptr */
  if (lb->Serialize_Structure != NULL)
    lb->Serialize_Structure(zz, &bufptr);

  *buf = bufptr;
}

/*****************************************************************************/
void Zoltan_LB_Deserialize(struct Zoltan_Struct *zz, char **buf)
{
  char *bufptr = *buf;
  struct Zoltan_LB_Struct *lb = &(zz->LB);
  int orig_Imb_Tol_Len = lb->Imb_Tol_Len;

  /* Copy 12 integers into zz->LB */
  int *intptr = (int *) bufptr;
  lb->Num_Global_Parts = *intptr; intptr++;
  lb->Num_Global_Parts_Param = *intptr; intptr++;
  lb->Num_Local_Parts_Param = *intptr; intptr++;
  lb->Prev_Global_Parts_Param = *intptr; intptr++;
  lb->Prev_Local_Parts_Param = *intptr; intptr++;
  lb->Single_Proc_Per_Part = *intptr; intptr++;
  lb->Part_Info_Max_Len = *intptr; intptr++;
  lb->Part_Info_Len = *intptr; intptr++;
  lb->Remap_Flag = *intptr; intptr++;
  lb->Return_Lists = *intptr; intptr++;
  lb->Uniform_Parts = *intptr; intptr++;
  lb->Imb_Tol_Len = *intptr; intptr++;
  bufptr = (char *) intptr;

  /* Reset the functions (Point_Assign, etc.) associated with the LB_Method */
  strcpy(lb->Method_Name, bufptr);
  bufptr += MAX_PARAM_STRING_LEN;
  Zoltan_Set_Param(zz, "LB_METHOD", lb->Method_Name);

  /* Copy Part_Info */
  if (lb->Part_Info_Len) {
    size_t tmpSize = lb->Part_Info_Len * sizeof(struct Zoltan_part_info);
    lb->Part_Info = (struct Zoltan_part_info *) ZOLTAN_MALLOC(tmpSize);
    memcpy((char *)(lb->Part_Info), bufptr, tmpSize);
    bufptr += tmpSize;
  }

  /* Copy Imbalance_Tol */
  if (lb->Imb_Tol_Len) {
    size_t tmpSize = lb->Imb_Tol_Len * sizeof(float);
    if (lb->Imb_Tol_Len > orig_Imb_Tol_Len) {
      if (lb->Imbalance_Tol != NULL) ZOLTAN_FREE(&(lb->Imbalance_Tol));
      lb->Imbalance_Tol = (float *) ZOLTAN_MALLOC(tmpSize);
    }
    memcpy((char *)(lb->Imbalance_Tol), bufptr, tmpSize);
    bufptr += tmpSize;
  }

  /* Copy lb->Approach */
  strcpy(lb->Approach, bufptr);
  bufptr += MAX_PARAM_STRING_LEN;
 
  /* Copy Remap array; needed by Point_Assign  */
  int nbytes = lb->Num_Global_Parts * sizeof(int);
  intptr = (int *) bufptr;
  bufptr += sizeof(int);
  if (*intptr) {  // Remap data was sent
    lb->Remap = (int *) ZOLTAN_MALLOC(nbytes);
    memcpy(lb->Remap, bufptr, nbytes);
    bufptr += nbytes;
  }
  else {
    lb->Remap = NULL;
  }

  /* Deserialize the method-specific load balancing data; advance bufptr */
  if (lb->Deserialize_Structure != NULL)
    lb->Deserialize_Structure(zz, &bufptr);

  *buf = bufptr;
}

/****************************************************************************/
/* Functions for serializations that are not yet implemented.
 * These are placeholders until someone requests them.
 */
size_t Zoltan_Serialize_Structure_Size_Not_Implemented(ZZ const *zz)
{
  return 0;
}

void Zoltan_Serialize_Structure_Not_Implemented(ZZ const *zz, char **buf)
{
  char msg[1024];
  sprintf(msg, "Zoltan_Serialize_Structure not implemented for method %s; "
               "no data will be copied.\n"
               "Contact Zoltan developers to request serialization of this "
               "method", zz->LB.Method_Name);
  ZOLTAN_PRINT_WARN(zz->Proc,
                    "Zoltan_Serialize_Structure_Not_Implemented", msg);
}

void Zoltan_Deserialize_Structure_Not_Implemented(ZZ *zz, char **buf)
{
  char msg[1024];
  sprintf(msg, "Zoltan_Deserialize_Structure not implemented for method %s; "
               "no data will be copied.\n"
               "Contact Zoltan developers to request serialization of this "
               "method", zz->LB.Method_Name);
  ZOLTAN_PRINT_WARN(zz->Proc, 
                   "Zoltan_Deserialize_Structure_Not_Implemented", msg);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
