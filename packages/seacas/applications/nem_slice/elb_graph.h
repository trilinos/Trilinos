/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _ELB_GRAPH_CONST_H_
#define _ELB_GRAPH_CONST_H_

struct Problem_Description;
struct Sphere_Info;
template <typename INT> struct Graph_Description;
template <typename INT> struct Mesh_Description;
template <typename INT> struct Weight_Description;

template <typename INT>
int generate_graph(Problem_Description *problem,    /* Pointer to structure containing information
                                                     * about the type of decomposition */
                   Mesh_Description<INT> * mesh,    /* Pointer to structure containing the mesh */
                   Graph_Description<INT> *graph,   /* Pointer to structure in which to store
                                                     * the graph. */
                   Weight_Description<INT> *weight, /* Pointer to structure for graph weighting */
                   Sphere_Info *            sphere  /* Pointer to sphere adjustment structure */
);

#endif /* _ELB_GRAPH_CONST_H_ */
