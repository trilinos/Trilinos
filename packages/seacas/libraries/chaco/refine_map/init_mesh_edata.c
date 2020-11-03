/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_edata

/* Initialize the mapping of sets to endpoints of wires in the mesh. */
void init_mesh_edata(struct refine_edata *edata,       /* desire data for all edges */
                     int                  mesh_dims[3] /* dimensions of mesh */
)
{
  int wire;    /* loops through wires */
  int i, j, k; /* loop counters */

  wire = 0;
  /* First do all the x-axis wires. */
  for (k = 0; k < mesh_dims[2]; k++) {
    for (j = 0; j < mesh_dims[1]; j++) {
      for (i = 0; i < mesh_dims[0] - 1; i++) {
        edata[wire].node1 = i + mesh_dims[0] * (j + k * mesh_dims[1]);
        edata[wire].node2 = i + 1 + mesh_dims[0] * (j + k * mesh_dims[1]);
        edata[wire].dim   = 0;
        wire++;
      }
    }
  }

  /* Now do all the y-axis wires. */
  for (k = 0; k < mesh_dims[2]; k++) {
    for (j = 0; j < mesh_dims[1] - 1; j++) {
      for (i = 0; i < mesh_dims[0]; i++) {
        edata[wire].node1 = i + mesh_dims[0] * (j + k * mesh_dims[1]);
        edata[wire].node2 = i + mesh_dims[0] * (j + 1 + k * mesh_dims[1]);
        edata[wire].dim   = 1;
        wire++;
      }
    }
  }

  /* Finally, do all the z-axis wires. */
  for (k = 0; k < mesh_dims[2] - 1; k++) {
    for (j = 0; j < mesh_dims[1]; j++) {
      for (i = 0; i < mesh_dims[0]; i++) {
        edata[wire].node1 = i + mesh_dims[0] * (j + k * mesh_dims[1]);
        edata[wire].node2 = i + mesh_dims[0] * (j + (k + 1) * mesh_dims[1]);
        edata[wire].dim   = 2;
        wire++;
      }
    }
  }
}
