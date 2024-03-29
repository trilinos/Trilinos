/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

struct Problem_Description;
struct Sphere_Info;
template <typename INT> struct Graph_Description;
template <typename INT> struct Mesh_Description;
struct Weight_Description;

template <typename INT>
int generate_graph(Problem_Description *problem,  /* Pointer to structure containing information
                                                   * about the type of decomposition */
                   Mesh_Description<INT>  *mesh,  /* Pointer to structure containing the mesh */
                   Graph_Description<INT> *graph, /* Pointer to structure in which to store
                                                   * the graph. */
                   Weight_Description *weight,    /* Pointer to structure for graph weighting */
                   Sphere_Info        *sphere     /* Pointer to sphere adjustment structure */
);
