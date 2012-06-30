#ifndef _ELB_LOADBAL_CONST_H_
#define _ELB_LOADBAL_CONST_H_

#include "elb.h"

template <typename INT>
int generate_loadbal(
  Machine_Description*   machine,
  Problem_Description* problem,
  Mesh_Description<INT>* mesh,
  LB_Description<INT>* lb,
  Solver_Description* solve,
  Graph_Description<INT>* graph,
  Weight_Description<INT>* weight,
  Sphere_Info* sphere,
  int argc,
  char *argv[]
);


template <typename INT>
int generate_maps(
  Machine_Description* machine,
  Problem_Description* problem,
  Mesh_Description<INT>* mesh,
  LB_Description<INT>* lb,
  Graph_Description<INT>* graph
);
#endif /* _ELB_LOADBAL_CONST_H_ */
