/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __DFS_CONST_H
#define __DFS_CONST_H

extern void LB_dfs_partition(LB *lb, int *counter, float *c1);
extern void LB_dfs_migrate(LB *lb, int *nsentags,
			   pRegion *import_tags, int *nrectags, 
			   float *c2, float *c3, int *counter3, int *counter4);


#endif
