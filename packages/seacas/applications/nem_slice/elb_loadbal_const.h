#ifndef _ELB_LOADBAL_CONST_H_
#define _ELB_LOADBAL_CONST_H_

#include "elb_const.h"

extern
int generate_loadbal(
  MACHINE_PTR   machine,
  PROB_INFO_PTR problem,
  MESH_INFO_PTR mesh,
  LB_INFO_PTR lb,
  SOLVE_INFO_PTR solve,
  GRAPH_INFO_PTR graph,
  WEIGHT_INFO_PTR weight,
  SPHERE_INFO_PTR sphere,
  int argc,
  char *argv[]
);

extern
int generate_maps(
  MACHINE_PTR machine,
  PROB_INFO_PTR problem,
  MESH_INFO_PTR mesh,
  LB_INFO_PTR lb,
  GRAPH_INFO_PTR graph
);

extern int extract_connected_lists( 
       int nrow, 
       const int* columns,
       const int* rows, 
       int* list,
       int** list_ptr 
);

#endif /* _ELB_LOADBAL_CONST_H_ */
