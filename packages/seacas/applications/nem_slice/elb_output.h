/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _ELB_OUTPUT_CONST_H_
#define _ELB_OUTPUT_CONST_H_
#include <string> // for string
struct Machine_Description;
struct Problem_Description;
struct Sphere_Info;
template <typename INT> struct LB_Description;
template <typename INT> struct Mesh_Description;

template <typename INT>
int write_nemesis(std::string &nemI_out_file, Machine_Description *machine,
                  Problem_Description *problem, Mesh_Description<INT> *mesh,
                  LB_Description<INT> *lb, Sphere_Info *sphere);

template <typename INT>
int write_vis(std::string &nemI_out_file, std::string &exoII_inp_file, Machine_Description *machine,
              Problem_Description *prob, Mesh_Description<INT> *mesh, LB_Description<INT> *lb);

#endif /* _ELB_OUTPUT_CONST_H_ */
