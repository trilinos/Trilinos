/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __LB_SET_FN_CONST_H
#define __LB_SET_FN_CONST_H

extern int Zoltan_Set_Pre_Migrate_Fn(ZZ *, ZOLTAN_PRE_MIGRATE_FN *, void *);
extern int Zoltan_Set_Mid_Migrate_Fn(ZZ *, ZOLTAN_MID_MIGRATE_FN *, void *);
extern int Zoltan_Set_Post_Migrate_Fn(ZZ *, ZOLTAN_POST_MIGRATE_FN *, void *);

#endif
