/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __DFS_CONST_H
#define __DFS_CONST_H

extern void Zoltan_Oct_dfs_partition(ZZ *zz, int *counter, float *c1);
extern void Zoltan_Oct_dfs_migrate(ZZ *zz, int *nsentags,
			   pRegion *import_tags, int *nrectags, 
			   float *c2, float *c3, int *counter3, int *counter4);


#endif
