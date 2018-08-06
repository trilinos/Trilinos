/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

#include "exodusII.h" // for MAX_LINE_LENGTH, ex_close, etc
#include "globals.h"
#include "nem_spread.h"
#include "rf_format.h"
#include "rf_io_const.h" // for Debug_Flag, ExoFile
#include <cstdio>        // for printf, fprintf, stderr
#include <cstdlib>       // for exit
#include <cstring>       // for memset

/* need to hang on to this to write it out to the proc 0 file */
char GeomTitle[MAX_LINE_LENGTH + 1];

template void NemSpread<double, int>::read_mesh_param();
template void NemSpread<float, int>::read_mesh_param();
template void NemSpread<double, int64_t>::read_mesh_param();
template void NemSpread<float, int64_t>::read_mesh_param();

template <typename T, typename INT> void NemSpread<T, INT>::read_mesh_param()

/*
 *       Function which reads the sizing parameters from the EXODUS II database,
 *       which contains the mesh.  This is used to cross-check against the
 *       load balancer file.
 *
 *      ------------------------------------------------------------------------
 *
 *       Functions called:
 *
 *       check_exodus_error -- function which handles the error code returned by
 *                             calls to EXODUS II API routines.
 *
 *      ------------------------------------------------------------------------
 */
{

  /* Local variables */

  const char *yo = "read_mesh_param";
  char *      exofile;
  int         exoid, error;
  int         cpu_ws;
  float       version;

  /**************************** execution begins *******************************/
  cpu_ws = sizeof(float);

  /*
   *         Generate the name of the parallel mesh file for this processor.
   *         Note: Default case is the scalar mesh file.
   */

  exofile = ExoFile;

  /* initialize the io word size on all of the processors */
  io_ws = 0;

  int mode = EX_READ | int64api;

  /* Open the EXODUS II mesh file */
  exoid = ex_open(exofile, mode, &cpu_ws, &io_ws, &version);
  if (exoid == -1) {
    fprintf(stderr, "%s: ERROR opening up the mesh exoII file, %s\n", yo, exofile);
    exit(-1);
  }

  /* Read the initialization parameters */
  memset(GeomTitle, '\0', MAX_LINE_LENGTH * sizeof(char));
  ex_init_params info{};
  info.title[0] = '\0';
  error         = ex_get_init_ext(exoid, &info);
  check_exodus_error(error, "ex_get_init");

  strncpy(GeomTitle, info.title, MAX_LINE_LENGTH);
  GeomTitle[MAX_LINE_LENGTH] = '\0';
  globals.Num_Dim            = info.num_dim;
  globals.Num_Node           = info.num_nodes;
  globals.Num_Elem           = info.num_elem;
  globals.Num_Elem_Blk       = info.num_elem_blk;
  globals.Num_Node_Set       = info.num_node_sets;
  globals.Num_Side_Set       = info.num_side_sets;

  printf("\nExodus file (%s)\n", exofile);
  printf("\tTitle of file: %s\n", GeomTitle);
  printf("\tDimensionality of problem = %d\n", globals.Num_Dim);
  printf("\tNumber of nodes           = " ST_ZU "\n", globals.Num_Node);
  printf("\tNumber of elements        = " ST_ZU "\n", globals.Num_Elem);
  printf("\tNumber of element blocks  = %d\n", globals.Num_Elem_Blk);
  printf("\tNumber of node sets       = %d\n", globals.Num_Node_Set);
  printf("\tNumber of side sets       = %d\n\n", globals.Num_Side_Set);

  /* Close the file */
  error = ex_close(exoid);
  check_exodus_error(error, "ex_close");

} /* END of routine read_mesh_params *****************************************/
