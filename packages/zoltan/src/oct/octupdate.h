/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __OCTUPDATE_H
#define __OCTUPDATE_H

#include "octupdate_const.h"

static void    LB_get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		      COORD min, COORD max, int wgtflag, float *c4);
static int     LB_oct_fix(LB *lb, pRegion Region_array, int num_objs);
static int     LB_oct_global_insert_object(LB *, pRegion Region_array, int num_objs);
static pOctant LB_oct_global_find(OCT_Global_Info *OCT_info,COORD point);
static pOctant LB_oct_findOctant(OCT_Global_Info *OCT_info,pOctant oct, COORD coord);
static void    LB_oct_global_dref(OCT_Global_Info *OCT_info);
static int     LB_oct_subtree_dref(OCT_Global_Info *OCT_info,pOctant oct);
static void    LB_oct_terminal_coarsen(OCT_Global_Info *OCT_info,pOctant oct);
static void    LB_oct_set_maxregions(int max);

#endif
