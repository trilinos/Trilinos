/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
#ifndef __OCTUPDATE_CONST_H
#define __OCTUPDATE_CONST_H

#include "octree_const.h"

#ifdef LGG_MIGOCT
extern void    LB_oct_resetIdCount(int start_count);
extern int     LB_oct_nextId(void);
#endif /* LGG_MIGOCT */
extern int     LB_oct_subtree_insert(LB *, pOctant oct, pRegion region);

extern void LB_oct_print_stats(LB *lb, double timetotal, double *timers, 
                           int *counters, float *c, int STATS_TYPE);
extern int LB_Set_Octpart_Param(char *name, char *val);
#endif
