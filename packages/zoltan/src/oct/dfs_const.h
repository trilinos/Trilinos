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
#ifndef __DFS_CONST_H
#define __DFS_CONST_H

#ifndef lint
static char *cvs_dfsconsth_id = "$Id$";
#endif

extern void dfs_partition(int *counter, float *c1);
extern void dfs_migrate(pRegion *export_tags, int *nsentags,
			pRegion *import_tags, int *nrectags, 
			float *c2, float *c3, int *counter3, int *counter4);


#endif
