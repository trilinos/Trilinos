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

extern void  costs_free(pOctant octree); 
extern float costs_value(pOctant octant);
extern float costs_global_compute();

#endif
