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
#ifndef __OCTUPDATE_H
#define __OCTUPDATE_H

#ifndef lint
static char *cvs_octantupdateh_id = "$Id$";
#endif

#include "octupdate_const.h"

void    LB_get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		      COORD min, COORD max, int wgtflag, float *c4);
int     LB_oct_fix(LB *lb, pRegion Region_array, int num_objs);
int     LB_oct_global_insert_object(LB *, pRegion Region_array, int num_objs);
pOctant LB_oct_global_find(COORD point);
pOctant LB_oct_findOctant(pOctant oct, COORD coord);
void    LB_oct_global_dref(void);
int     LB_oct_subtree_dref(pOctant oct);
void    LB_oct_terminal_coarsen(pOctant oct);
void    LB_oct_set_maxregions(int max);

#endif
