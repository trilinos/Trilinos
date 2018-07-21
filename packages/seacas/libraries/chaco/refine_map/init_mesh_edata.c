/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
