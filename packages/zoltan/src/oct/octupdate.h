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

void    get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		   COORD min, COORD max, float *c4);
int     oct_fix(LB *lb, pRegion Region_array, int num_objs);
int     oct_global_insert_object(pRegion Region_array, int num_objs);
pOctant oct_global_find(COORD point);
pOctant oct_findOctant(pOctant oct, COORD coord);
void    oct_global_dref(void);
int     oct_subtree_dref(pOctant oct);
void    oct_terminal_coarsen(pOctant oct);
void    oct_set_maxregions(int max);

#endif
