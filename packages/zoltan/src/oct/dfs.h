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
#ifndef __DFS_H
#define __DFS_H

#ifndef lint
static char *cvs_dfsh_id = "$Id$";
#endif

#include "dfs_const.h"

extern void LB_dfs_set_visit_criterion(int visit);
extern int  LB_dfs_SetIds(pOctant oct, int nprevoct);
extern void LB_visit_all_subtrees();
extern void LB_visit(pOctant octant);
extern void LB_visit_by_dist(pOctant octant, pOctant children[8]);
extern void LB_tag_subtree(pOctant octant, int partition);

/* global variables for dfs.c !! */
int partition;              /* Partition number we are working on */
float total;                /* Cost of all complete partitions so far */
float pcost;                /* Current partition cost */
float optcost;              /* Optimal partition cost */
float pmass;                /* octant volume for partition */
double pcoord[3];           /* Sum of octant position-volume products */

#endif
