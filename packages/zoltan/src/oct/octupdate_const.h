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

#ifndef __OCTUPDATE_CONST_H
#define __OCTUPDATE_CONST_H

#include "octant_const.h"
#define MINOCTREGIONS 1              /* minimum number of regions per octant */

extern void    LB_oct_gen_tree_from_input_data(LB *lb, int, int *c1, int *c2,
					       int *c3, float *c0);
#ifdef LGG_MIGOCT
extern void    LB_oct_roots_in_order(pOctant **roots_ret,int *nroots_ret);
extern void    LB_oct_resetIdCount(int start_count);
extern int     LB_oct_nextId(void);
#endif /* LGG_MIGOCT */
extern void    LB_oct_terminal_refine(LB *, pOctant oct,int count);
extern pOctant LB_oct_findId(int i);
extern pOctant LB_oct_global_insert(LB *, pRegion region);
extern int     LB_oct_subtree_insert(LB *, pOctant oct, pRegion region);

extern void LB_oct_print_stats(LB *lb, double timetotal, double *timers, 
                           int *counters, float *c, int STATS_TYPE);
extern int LB_Set_Octpart_Param(char *name, char *val);
#endif
