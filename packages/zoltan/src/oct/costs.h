/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __COSTS_H
#define __COSTS_H

#include "costs_const.h"

extern void  LB_costs_init(OCT_Global_Info *OCT_info,pOctant octree);
extern float LB_costs_subtree_compute(OCT_Global_Info *OCT_info,pOctant octant, int *seq);
extern float LB_costs_weight(pOctant octant);

#endif
