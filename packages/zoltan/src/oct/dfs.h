/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __DFS_H
#define __DFS_H

#include "dfs_const.h"

extern void LB_dfs_set_visit_criterion(int visit);
extern int  LB_dfs_SetIds(pOctant oct, int nprevoct);
extern void LB_visit_all_subtrees(OCT_Global_Info *OCT_info);
extern void LB_visit(OCT_Global_Info *OCT_info,pOctant octant);
extern void LB_visit_by_dist(OCT_Global_Info *OCT_info,pOctant octant, pOctant children[8]);
extern void LB_tag_subtree(OCT_Global_Info *OCT_info,pOctant octant, int partition);

#endif
