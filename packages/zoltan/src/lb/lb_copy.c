/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
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

#define COPY_FIELD(f) to->f = from->f;

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

int Zoltan_LB_Copy_Struct(ZZ *toZZ, ZZ *fromZZ)
{
char *yo = "Zoltan_LB_Copy_Struct";
int proc = fromZZ->Proc;

  struct Zoltan_LB_Struct *to = &(toZZ->LB);
  struct Zoltan_LB_Struct *from = &(fromZZ->LB);

  Zoltan_LB_Free_Struct(&(toZZ->LB));

  COPY_FIELD(Part_Info_Len);
  COPY_FIELD(Part_Info_Max_Len);
  COPY_FIELD(Num_Global_Parts);
  COPY_FIELD(Num_Global_Parts_Param);
  COPY_FIELD(Num_Local_Parts_Param);
  COPY_FIELD(Prev_Global_Parts_Param);
  COPY_FIELD(Prev_Local_Parts_Param);
  COPY_FIELD(Single_Proc_Per_Part);
  COPY_FIELD(Remap_Flag);
  COPY_FIELD(Return_Lists);
  COPY_FIELD(Uniform_Parts);
  COPY_FIELD(Method);
  COPY_FIELD(LB_Fn);
  COPY_FIELD(Imb_Tol_Len);
  COPY_FIELD(Free_Structure);
  COPY_FIELD(Copy_Structure);
  COPY_FIELD(Point_Assign);
  COPY_FIELD(Box_Assign);

  COPY_BUFFER(Part_Info, struct Zoltan_part_info, to->Part_Info_Max_Len);

  COPY_BUFFER(Remap, int, to->Num_Global_Parts);

  COPY_BUFFER(PartDist, int, to->Num_Global_Parts + 1);

  COPY_BUFFER(ProcDist, int, fromZZ->Num_Proc + 1);

  COPY_BUFFER(Imbalance_Tol, float, to->Imb_Tol_Len);

  if (from->Data_Structure) {
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

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void Zoltan_Migrate_Copy_Struct(struct Zoltan_Migrate_Struct *to, 
                                struct Zoltan_Migrate_Struct *from)
{
  COPY_FIELD(Auto_Migrate);
  COPY_FIELD(Only_Proc_Changes);
  COPY_FIELD(Pre_Migrate_PP );
  COPY_FIELD(Pre_Migrate_PP_Fort);
  COPY_FIELD(Pre_Migrate_PP_Data);
  COPY_FIELD(Mid_Migrate_PP);
  COPY_FIELD(Mid_Migrate_PP_Fort);
  COPY_FIELD(Mid_Migrate_PP_Data);
  COPY_FIELD(Post_Migrate_PP );
  COPY_FIELD(Post_Migrate_PP_Fort);
  COPY_FIELD(Post_Migrate_PP_Data);
  COPY_FIELD(Pre_Migrate);
  COPY_FIELD(Pre_Migrate_Fort);
  COPY_FIELD(Pre_Migrate_Data);
  COPY_FIELD(Mid_Migrate);
  COPY_FIELD(Mid_Migrate_Fort );
  COPY_FIELD(Mid_Migrate_Data);
  COPY_FIELD(Post_Migrate);
  COPY_FIELD(Post_Migrate_Fort);
  COPY_FIELD(Post_Migrate_Data);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
