#include "dfs_const.h"

extern void dfs_set_visit_criterion(int visit);
extern int  dfs_SetIds(pOctant oct, int nprevoct);
extern void visit_all_subtrees();
extern void visit(pOctant octant);
extern void visit_by_dist(pOctant octant, pOctant children[8]);
extern void tag_subtree(pOctant octant, int partition);

/* global variables for dfs.c !! */
int partition;              /* Partition number we are working on */
float total;                /* Cost of all complete partitions so far */
float pcost;                /* Current partition cost */
float optcost;              /* Optimal partition cost */
float pmass;                /* octant volume for partition */
double pcoord[3];           /* Sum of octant position-volume products */
