/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_edata
#include <stdio.h>      // for NULL

struct refine_edata *
find_edge_mesh(int                  vertex,    /* vertex in comm_graph */
               int                  dim,       /* direction of edge from node */
               struct refine_edata *edata,     /* data structure for edge preferences */
               int *                mesh_dims, /* dimensions of mesh */
               int *                vtx2node   /* maps comm_graph vtxs to processors */
)
{
  struct refine_edata *eguy;      /* returned pointer to edge info */
  int                  dir;       /* higher or lower direction? */
  int                  my_node;   /* processor vertex assigned to */
  int                  my_loc[3]; /* location of my processor */
  int                  index = 0; /* computed index into edata */

  if (dim < 0) {
    dir = -1;
    dim = -(dim + 1);
  }
  else {
    dir = 1;
  }

  my_node   = vtx2node[vertex];
  my_loc[0] = my_node % mesh_dims[0];
  my_loc[1] = (my_node / mesh_dims[0]) % mesh_dims[1];
  my_loc[2] = my_node / (mesh_dims[0] * mesh_dims[1]);

  if ((my_loc[dim] == 0 && dir == -1) || (my_loc[dim] == mesh_dims[dim] - 1 && dir == 1)) {
    eguy = NULL;
  }

  else { /* Figure out where edge is in data structure. */
    /* Note: indexing must match with that in init_mesh_edata. */
    if (dir < 0) {
      --my_loc[dim];
    }

    if (dim == 0) { /* Edge in x-direction. */
      index = (mesh_dims[0] - 1) * mesh_dims[1] * my_loc[2] + (mesh_dims[0] - 1) * my_loc[1] +
              my_loc[0];
    }

    else if (dim == 1) { /* Edge in y-direction. */
      index = (mesh_dims[0] - 1) * mesh_dims[1] * mesh_dims[2] +
              mesh_dims[0] * (mesh_dims[1] - 1) * my_loc[2] + mesh_dims[0] * my_loc[1] + my_loc[0];
    }
    else if (dim == 2) { /* Edge in z-direction. */
      index = (mesh_dims[0] - 1) * mesh_dims[1] * mesh_dims[2] +
              mesh_dims[0] * (mesh_dims[1] - 1) * mesh_dims[2] +
              mesh_dims[0] * mesh_dims[1] * my_loc[2] + mesh_dims[0] * my_loc[1] + my_loc[0];
    }
    eguy = &(edata[index]);
  }

  return (eguy);
}
