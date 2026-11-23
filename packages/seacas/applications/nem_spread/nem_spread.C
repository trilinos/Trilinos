/*
 * Copyright(C) 1999-2021, 2023, 2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*--------------------------------------------------------------------------*/
/*    This program takes an Exodus II FEM database file and a parallel      */
/* load-balance file and creates a set of parallel Exodus II FEM files for  */
/* use by SALSA.                                                            */
/*--------------------------------------------------------------------------*/

#include "add_to_log.h" // for add_to_log
#include "exodusII.h"   // for ex_opts, ex_int64_status, etc
#include "fmt/ostream.h"
#include "nem_spread.h"     // for NemSpread, second, etc
#include "ps_pario.h"       // for PIO_Info
#include "ps_pario_const.h" // for Parallel_IO
#include "rf_allo.h"        // for safe_free
#include "rf_io.h"          // for ExoFile, Debug_Flag, etc
#include <cstdint>          // for int64_t
#include <cstdio>           // for stderr, etc
#include <cstdlib>          // for exit
#include <cstring>
#include <unistd.h> // for getopt, optarg, optind

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#include "XGetopt.h"
#else
#include "getopt.h" // for getopt
#endif

template <typename T, typename INT>
int read_pexoII_info(NemSpread<T, INT> &spreader, const char *filename);

int check_change_sets(int exoid, int &selected_change_set)
{
  // First, see if the file contains change sets...
  int num_change_sets = ex_inquire_int(exoid, EX_INQ_NUM_CHILD_GROUPS);

  if (num_change_sets == 0) {
    if (selected_change_set > 0) {
      fmt::print(stderr,
                 "ERROR: Reading from change set {} was specified, but the file contains no "
                 "change sets.\n",
                 selected_change_set);
      return 0;
    }
    return exoid;
  }

  // File contains change sets...
  if (selected_change_set == 0) {
    fmt::print(stderr,
               "\nWARNING: Exodus database contains {} change sets.\n         Setting to read from "
               "first change set since `-change_set #` option not specified.\n\n",
               num_change_sets);
    selected_change_set = 1;
    return exoid + selected_change_set;
  }

  if (selected_change_set > num_change_sets) {
    fmt::print(stderr,
               "ERROR: Change set {} was specified for reading, but mesh only contains {} change "
               "sets.\n",
               selected_change_set, num_change_sets);
    return 0;
  }
  // Contains change sets and selected change set is in range...
  fmt::print(stderr, "NOTE: Mesh data will be read from change set {} of {}\n", selected_change_set,
             num_change_sets);
  return exoid + selected_change_set;
}

template <typename T, typename INT>
int nem_spread(NemSpread<T, INT> &spreader, const char *salsa_cmd_file, int subcycles, int cycle)
{
  static char yo[] = "nem_spread";
  /* Local declarations. */
  double start_t;
  double end_t;

  fmt::print("Using {} byte integers and {} byte floating point values.\n", sizeof(INT), sizeof(T));

  /*
   * Read in the ASCII input file from the front end.
   * NOTE: In this case we read only the information needed to create the
   *       parallel exodus files.
   */
  fmt::print("Reading the command file, {}\n", salsa_cmd_file);
  if (read_pexoII_info(spreader, salsa_cmd_file) < 0) {
    fmt::print(stderr,
               "{} ERROR: Could not read in the the I/O command file"
               " \"{}\"!\n",
               yo, salsa_cmd_file);
    exit(1);
  }

  if (!spreader.check_inp()) {
    fmt::print(stderr, "{} ERROR: Error in user specified parameters.\n", yo);
    exit(1);
  }

  /* If debug is on the turn on netCDF/Exodus information as well */
  if (Debug_Flag > 0) {
    ex_opts(EX_VERBOSE | EX_DEBUG);
  }
  else {
    ex_opts(EX_VERBOSE);
  }

  /*
   * Read initial information from the mesh file.
   *  - This provides error checking against the load balance file
   *  - Broadcast the information
   */
  spreader.read_mesh_param();

  /*
   * Process the load balance information
   *  - Read info from the lb file with Proc 0
   *  - Distribute the information to all processors
   */
  start_t = second();
  spreader.load_lb_info();
  end_t = second() - start_t;
  fmt::print("\nLoad load balance information time: {} (sec.)\n\n", end_t);

  /* Process subcycle and cycle if specified. */
  if ((subcycles > 0 && cycle == -1) || (cycle != -1 && subcycles == 0)) {
    fmt::print(stderr, "ERROR: Only one of the -subcycle and -cycle options was specified.\n");
    fmt::print(stderr, "       Either both or neither are required.\n");
    exit(1);
  }

  if (subcycles > 0) {
    /* Determine number of processors per subcycle. */
    int part_count        = (spreader.Proc_Info[0] + subcycles - 1) / subcycles;
    int start_part        = part_count * cycle;
    spreader.Proc_Info[4] = start_part;
    spreader.Proc_Info[5] = part_count;
  }

  /*
   * Verify parameters in case spreading a subset of mesh...
   */
  if (spreader.Proc_Info[4] < 0) {
    spreader.Proc_Info[4] = 0;
  }
  if (spreader.Proc_Info[5] <= 0) {
    spreader.Proc_Info[5] = spreader.Proc_Info[0];
  }

  if (spreader.Proc_Info[4] + spreader.Proc_Info[5] > spreader.Proc_Info[0]) {
    spreader.Proc_Info[5] = spreader.Proc_Info[0] - spreader.Proc_Info[4];
  }

  if (spreader.Proc_Info[4] != 0 || spreader.Proc_Info[5] != spreader.Proc_Info[0]) {
    fmt::print(
        "\nSpreading subset of mesh.  Starting with processor {} and outputting {} processors.\n",
        spreader.Proc_Info[4], spreader.Proc_Info[5]);
  }

  /*
   * Get any restart parameter information
   *  - Read the parameters from the input ExodusII file
   */
  if (spreader.Restart_Info.Flag > 0) {
    fmt::print("Load exoII restart param info to each proc.\n\n");
    start_t = second();
    spreader.read_restart_params();
    end_t = second() - start_t;
    fmt::print("Load restart parameters time: {} (sec.)\n\n", end_t);
  }

  /*
   * Read the ExodusII mesh file and distribute information
   * contained in it to all processors
   *  - Each processor only gets the information that it needs
   *    to solve its piece of the mesh
   */
  fmt::print("Load exoII mesh info to each proc.\n\n");
  start_t = second();
  spreader.load_mesh();
  end_t = second() - start_t;

  fmt::print("Load mesh time: {} (sec.)\n\n", end_t);

  /*
   * Get any restart variable data
   *  - Read the restart data from the input ExodusII file, distribute
   *    it to the processors, and write it to the parallel files
   */
  if (spreader.Restart_Info.Flag > 0) {
    fmt::print("Load exoII restart data info to each proc.\n\n");
    start_t = second();
    spreader.read_restart_data();
    end_t = second() - start_t;
    fmt::print("Load restart data time: {} (sec.)\n\n", end_t);
  }

  fmt::print("Write of parallel exodus complete\n");

  safe_free(reinterpret_cast<void **>(&(PIO_Info.RDsk_List)));

  for (int i = 0; i < spreader.Proc_Info[0]; i++) {
    if (spreader.globals.Elem_Type != nullptr) {
      safe_free((void **)&spreader.globals.Elem_Type[i]);
    }
    safe_free((void **)&spreader.globals.Proc_SS_Ids[i]);
    safe_free((void **)&spreader.globals.Proc_SS_GEMap_List[i]);
    safe_free((void **)&spreader.globals.Proc_NS_Ids[i]);
    safe_free((void **)&spreader.globals.Proc_NS_GNMap_List[i]);
    safe_free((void **)&spreader.globals.Proc_Nodes_Per_Elem[i]);
    safe_free((void **)&spreader.globals.GElem_Blks[i]);
    spreader.globals.Proc_Global_Elem_Id_Map[i].clear();
    spreader.globals.Proc_Global_Node_Id_Map[i].clear();
  }
  safe_free((void **)&spreader.globals.Elem_Type);
  return 0;
}

template int nem_spread(NemSpread<float, int> &spreader, const char *salsa_cmd_file, int subcycles,
                        int cycle);
template int nem_spread(NemSpread<double, int> &spreader, const char *salsa_cmd_file, int subcycles,
                        int cycle);
template int nem_spread(NemSpread<float, int64_t> &spreader, const char *salsa_cmd_file,
                        int subcycles, int cycle);
template int nem_spread(NemSpread<double, int64_t> &spreader, const char *salsa_cmd_file,
                        int subcycles, int cycle);
