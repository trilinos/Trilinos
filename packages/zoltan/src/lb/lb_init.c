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

#include "zz_const.h"
#include "lb_init_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_Migrate_Init(struct Zoltan_Migrate_Struct *mig)
{
  mig->Auto_Migrate = ZOLTAN_AUTO_MIGRATE_DEF;
  mig->Pre_Migrate = NULL;
  mig->Mid_Migrate = NULL;
  mig->Post_Migrate = NULL;
  mig->Pre_Migrate_Fort = NULL;
  mig->Mid_Migrate_Fort = NULL;
  mig->Post_Migrate_Fort = NULL;
}

void Zoltan_LB_Init(struct Zoltan_LB_Struct *lb)
{
  lb->Method = RCB;
  lb->LB_Fn = Zoltan_RCB;
  lb->Return_Lists = ZOLTAN_LB_RETURN_LISTS_DEF;
  lb->Imbalance_Tol = ZOLTAN_LB_IMBALANCE_TOL_DEF;
  lb->Data_Structure = NULL;
  lb->Free_Structure = Zoltan_RCB_Free_Structure;
  lb->Point_Assign = Zoltan_RB_Point_Assign;
  lb->Box_Assign = Zoltan_RB_Box_Assign;
}
