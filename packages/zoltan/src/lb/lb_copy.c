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
