/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "copy_string_cpp.h"
#include "exodusII.h" // for MAX_LINE_LENGTH, ex_close, etc
#include "fmt/ostream.h"
#include "globals.h"
#include "nem_spread.h"
#include "rf_io_const.h" // for Debug_Flag, ExoFile
#include <cstdio>        // for stderr
#include <cstdlib>       // for exit
#include <cstring>       // for memset

/* need to hang on to this to write it out to the proc 0 file */
std::string GeomTitle;

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
  std::string exofile;
  int         exoid;
  int         error;
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
  exoid = ex_open(exofile.c_str(), mode, &cpu_ws, &io_ws, &version);
  if (exoid == -1) {
    fmt::print(stderr, "{}: ERROR opening up the mesh exoII file, {}\n", __func__, exofile.c_str());
    exit(-1);
  }

  /* Read the initialization parameters */
  ex_init_params info{};
  info.title[0] = '\0';
  error         = ex_get_init_ext(exoid, &info);
  check_exodus_error(error, "ex_get_init");

  GeomTitle            = info.title;
  globals.Num_Dim      = info.num_dim;
  globals.Num_Node     = info.num_nodes;
  globals.Num_Elem     = info.num_elem;
  globals.Num_Elem_Blk = info.num_elem_blk;
  globals.Num_Node_Set = info.num_node_sets;
  globals.Num_Side_Set = info.num_side_sets;

  fmt::print("\nExodus file ({})\n", exofile);
  fmt::print("\tTitle of file: '{}'\n", GeomTitle);
  fmt::print("\tDimensionality of problem = {:14n}\n", globals.Num_Dim);
  fmt::print("\tNumber of nodes           = {:14n}\n", globals.Num_Node);
  fmt::print("\tNumber of elements        = {:14n}\n", globals.Num_Elem);
  fmt::print("\tNumber of element blocks  = {:14n}\n", globals.Num_Elem_Blk);
  fmt::print("\tNumber of node sets       = {:14n}\n", globals.Num_Node_Set);
  fmt::print("\tNumber of side sets       = {:14n}\n\n", globals.Num_Side_Set);

  /* Close the file */
  error = ex_close(exoid);
  check_exodus_error(error, "ex_close");

} /* END of routine read_mesh_params *****************************************/
