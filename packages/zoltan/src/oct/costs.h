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
#ifndef __COSTS_H
#define __COSTS_H

#ifndef lint
static char *cvs_costsh_id = "$Id$";
#endif

#include "costs_const.h"

extern void  costs_init(pOctant octree);
extern float costs_subtree_compute(pOctant octant, int *seq);
extern float costs_weight(pOctant octant);

#endif
