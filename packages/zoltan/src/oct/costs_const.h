/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __COSTS_CONST_H
#define __COSTS_CONST_H

extern void  LB_costs_free(OCT_Global_Info *OCT_info,pOctant octree); 
extern float LB_costs_value(pOctant octant);
extern float LB_costs_global_compute(OCT_Global_Info *OCT_info);

#endif
