/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __OCTUPDATE_CONST_H
#define __OCTUPDATE_CONST_H

#include "octant_const.h"
#define MINOCTREGIONS 1              /* minimum number of regions per octant */

#ifdef LGG_MIGOCT
extern void    LB_oct_roots_in_order(pOctant **roots_ret,int *nroots_ret);
extern void    LB_oct_resetIdCount(int start_count);
extern int     LB_oct_nextId(void);
#endif /* LGG_MIGOCT */
extern pOctant LB_oct_findId(int i);
extern int     LB_oct_subtree_insert(LB *, pOctant oct, pRegion region);

extern void LB_oct_print_stats(LB *lb, double timetotal, double *timers, 
                           int *counters, float *c, int STATS_TYPE);
extern int LB_Set_Octpart_Param(char *name, char *val);
#endif
