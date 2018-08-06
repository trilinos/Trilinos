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

#include "refine_map.h" // for refine_vdata
#include <stdio.h>      // for NULL

void update_mesh_vdata(int    old_loc, /* previous node for moved vertex in moved dimension */
                       int    new_loc, /* new node for moved vertex in moved dimension */
                       int    dim,     /* dimension that was changed */
                       double ewgt,    /* weight of edge */
                       struct refine_vdata *vdata,        /* array of vertex data */
                       int                  mesh_dims[3], /* size of processor mesh */
                       int                  neighbor,     /* vertex impacted by flip */
                       int *vtx2node /* mapping from comm_graph vtxs to processors */
)
{
  struct refine_vdata *vptr          = NULL; /* correct element in vdata */
  int                  offset        = 0;    /* index into vdata array */
  int                  my_loc        = 0;    /* my location in relevant dimension */
  int                  neighbor_node = 0;    /* processor neighbor assigned to */

  neighbor_node = vtx2node[neighbor];

  if (dim == 0) {
    offset = 0;
    my_loc = neighbor_node % mesh_dims[0];
  }
  else if (dim == 1) {
    offset = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
    my_loc = (neighbor_node / mesh_dims[0]) % mesh_dims[1];
  }
  else if (dim == 2) {
    offset = 2 * mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
    my_loc = neighbor_node / (mesh_dims[0] * mesh_dims[1]);
  }

  vptr = &vdata[offset + neighbor];

  /* If I'm far away from the flipped edge, I'm not effected. */
  if (!((old_loc < my_loc && new_loc < my_loc) || (old_loc > my_loc && new_loc > my_loc))) {
    if (old_loc < my_loc) { /* Old moves to the right to line up with me. */
      vptr->same += ewgt;
      vptr->below -= ewgt;
    }

    else if (old_loc > my_loc) { /* Old moves to the left to line up with me. */
      vptr->same += ewgt;
      vptr->above -= ewgt;
    }

    else if (new_loc < my_loc) { /* Old moves to the left to pass me. */
      vptr->same -= ewgt;
      vptr->below += ewgt;
    }

    else if (new_loc > my_loc) { /* Old moves to the right to pass me. */
      vptr->same -= ewgt;
      vptr->above += ewgt;
    }
  }
}
