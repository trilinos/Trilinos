/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef __COSTS_CONST_H
#define __COSTS_CONST_H

#ifndef lint
static char *cvs_costsconsth_id = "$Id$";
#endif

extern void  LB_costs_free(pOctant octree); 
extern float LB_costs_value(pOctant octant);
extern float LB_costs_global_compute();

#endif
