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
#include "irb_const.h"
#include "oct_util_const.h"
#include "reftree_const.h"


void LB_Free_Structure(
LB *lb)				/* load balance structure */
{
/*
 * Free any persistent memory associated with a method.
 */

  switch (lb->Method) {

    case RCB:
      LB_RCB_Free_Structure(lb);
      break;

    case IRB:
      LB_IRB_Free_Structure(lb);
      break;

    case OCTPART:
      LB_OCT_Free_Structure(lb);
      break;

    case REFTREE:
      LB_Reftree_Free_Structure(lb);
      break;
/*
 * Add calls to additional method-specific free routines here.
 */

    default: {
    }
  }
}
