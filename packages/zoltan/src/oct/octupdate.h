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

static void    Zoltan_Oct_get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		      COORD min, COORD max, int wgtflag, float *c4);
static int     Zoltan_Oct_fix(LB *lb, pRegion Region_array, int num_objs);
static int     Zoltan_Oct_global_insert_object(LB *, pRegion Region_array, int num_objs);
static pOctant Zoltan_Oct_global_find(OCT_Global_Info *OCT_info,COORD point);
static pOctant Zoltan_Oct_findOctant(OCT_Global_Info *OCT_info,pOctant oct, COORD coord);
static void    Zoltan_Oct_global_dref(LB *, OCT_Global_Info *OCT_info);
static int     Zoltan_Oct_subtree_dref(LB *, OCT_Global_Info *OCT_info,pOctant oct);
static void    Zoltan_Oct_terminal_coarsen(LB *, OCT_Global_Info *OCT_info,pOctant oct);
static void    Zoltan_Oct_set_maxregions(int max);
static void    Zoltan_Oct_set_minregions(int min);
static void    Zoltan_Oct_global_clear(OCT_Global_Info * OCT_info);
#endif
