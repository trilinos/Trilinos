/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef _ELB_LOADBAL_CONST_H_
#define _ELB_LOADBAL_CONST_H_

struct Machine_Description;
struct Problem_Description;
struct Solver_Description;
struct Sphere_Info;
template <typename INT> struct Graph_Description;
template <typename INT> struct LB_Description;
template <typename INT> struct Mesh_Description;
template <typename INT> struct Weight_Description;

template <typename INT>
int generate_loadbal(Machine_Description *machine, Problem_Description *problem,
                     Mesh_Description<INT> *mesh, LB_Description<INT> *lb,
                     Solver_Description *solve, Graph_Description<INT> *graph,
                     Weight_Description<INT> *weight, Sphere_Info *sphere, int argc, char *argv[]);

template <typename INT>
int generate_maps(Machine_Description *machine, Problem_Description *problem,
                  Mesh_Description<INT> *mesh, LB_Description<INT> *lb,
                  Graph_Description<INT> *graph);
#endif /* _ELB_LOADBAL_CONST_H_ */
