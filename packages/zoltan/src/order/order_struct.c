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
#include "order_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines to create, access, and destroy 
 *  Zoltan Ordering Structs (ZOS).
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Order_Create(ZOS **order_info, ZZ *zz)
{
  int ierr = ZOLTAN_OK; /* Error code to return */
  static char *yo = "Zoltan_Order_Create";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz == NULL){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Pointer to Zoltan struct is NULL");
    ierr = ZOLTAN_FATAL;
    return (ierr);
  }

  *order_info = (ZOS *) ZOLTAN_MALLOC (sizeof(ZOS));
  if (!(*order_info)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Not enough memory");
    ierr = ZOLTAN_MEMERR;
    return (ierr);
  }

  /* Initialize ordering struct */
  (*order_info)->zz = zz;
  (*order_info)->num_objects = 0;
  (*order_info)->gids = NULL;
  (*order_info)->lids = NULL;
  (*order_info)->method = NULL;
  (*order_info)->num_separators = 0;
  (*order_info)->sep_sizes = NULL;

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

int Zoltan_Order_Destroy(ZOS **order_info)
{
  int ierr = ZOLTAN_OK; /* Error code to return */
  /* static char *yo = "Zoltan_Order_Destroy"; */

  /* ZOLTAN_TRACE_ENTER(zz, yo); */

  if ((*order_info)->gids)      ZOLTAN_FREE(&((*order_info)->gids));
  if ((*order_info)->lids)      ZOLTAN_FREE(&((*order_info)->lids));
  if ((*order_info)->method)    ZOLTAN_FREE(&((*order_info)->method));
  if ((*order_info)->sep_sizes) ZOLTAN_FREE(&((*order_info)->sep_sizes));

  ZOLTAN_FREE(order_info);
  order_info = NULL;

  /* ZOLTAN_TRACE_EXIT(zz, yo); */
  return (ierr);
}
