#include "costs_const.h"

void  costs_init(pOctant octree);
float costs_subtree_compute(pOctant octant, int *seq);
float costs_weight(pOctant octant);
