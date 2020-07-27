/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _ELB_GROUPS_CONST_H_
#define _ELB_GROUPS_CONST_H_

#include <cstddef> // for size_t
struct Machine_Description;
struct Problem_Description;
template <typename INT> struct Graph_Description;
template <typename INT> struct Mesh_Description;

/* Function prototypes */
template <typename INT>
int parse_groups(Mesh_Description<INT> *mesh, /* Mesh information structure */
                 Problem_Description *  prob);  /* Problem information */

template <typename INT>
int get_group_info(Machine_Description *machine, Problem_Description *prob,
                   Mesh_Description<INT> *mesh, Graph_Description<INT> *graph, int elem2grp[],
                   int nprocg[], int nelemg[], size_t *max_vtx, size_t *max_adj);

#endif /* _ELB_GROUPS_CONST_H_ */
