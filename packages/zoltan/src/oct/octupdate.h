#include "octupdate_const.h"

void    get_bounds(LB *lb, pRegion *ptr1, int *num_objs, 
		   double min[3], double max[3], int *c4);
void    initialize(LB *lb, pRegion *ret, int local_id, LB_ID global_id);
int     oct_fix(LB *lb, pRegion Region_array, int num_objs);
int     oct_global_insert_object(pRegion Region_array, int num_objs);
pOctant oct_global_find(double point[3]);
pOctant oct_findOctant(pOctant oct, double coord[3]);
void    oct_global_dref(void);
int     oct_subtree_dref(pOctant oct);
void    oct_terminal_coarsen(pOctant oct);
