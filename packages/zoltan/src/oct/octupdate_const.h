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
extern void    Zoltan_Oct_resetIdCount(int start_count);
extern int     Zoltan_Oct_nextId(void);
#endif /* LGG_MIGOCT */
extern int     Zoltan_Oct_subtree_insert(ZZ *, pOctant oct, pRegion region);

extern void Zoltan_Oct_print_stats(ZZ *zz, double timetotal, double *timers, 
                           int *counters, float *c, int STATS_TYPE);
extern int Zoltan_Oct_Set_Param(char *name, char *val);
#endif
