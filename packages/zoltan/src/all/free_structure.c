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

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "lb_const.h"
#include "rcb_const.h"
#include "rib_const.h"
#include "oct_util_const.h"
#include "reftree_const.h"


void Zoltan_Free_Structure(
  ZZ *zz)				/* Zoltan structure */
{
/*
 * Free any persistent memory associated with a method.
 */

  switch (zz->Method) {

    case RCB:
      Zoltan_RCB_Free_Structure(zz);
      break;

    case RIB:
      Zoltan_RIB_Free_Structure(zz);
      break;

    case OCTPART:
      Zoltan_Oct_Free_Structure(zz);
      break;

    case REFTREE:
      Zoltan_Reftree_Free_Structure(zz);
      break;
/*
 * Add calls to additional method-specific free routines here.
 */

    default: {
    }
  }
}
