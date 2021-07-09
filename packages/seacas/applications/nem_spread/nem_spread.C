/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
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

#ifdef _MSC_VER
#include "XGetopt.h"
#else
#include "getopt.h" // for getopt
#endif

extern int read_mesh_file_name(const char *filename);

template <typename T, typename INT>
int read_pexoII_info(NemSpread<T, INT> &spreader, const char *filename);

template <typename T, typename INT>
int nem_spread(NemSpread<T, INT> &spreader, const char *salsa_cmd_file, int subcycles, int cycle);

int main(int argc, char *argv[])
{
  const char *salsa_cmd_file;
  int         c;

  double g_start_t    = second();
  bool   force_64_bit = false;
  int    start_proc   = 0;
  int    num_proc     = 0;
  int    subcycles    = 0;
  int    cycle        = -1;
  while ((c = getopt(argc, argv, "64Vhp:r:s:n:S:c:")) != -1) {
    switch (c) {
    case 'h':
      fmt::print(stderr, " usage:\n");
      fmt::print(stderr,
                 "\tnem_spread  [-s <start_proc>] [-n <num_proc>] [-S <subcycles> -c <cycle>] "
                 "[command_file]\n");
      fmt::print(stderr, "\t\tDecompose for processors <start_proc> to <start_proc>+<num_proc>\n");
      fmt::print(stderr, "\t\tDecompose for cycle <cycle> of <subcycle> groups\n");
      fmt::print(stderr, "\tnem_spread  [-V] [-h] (show version or usage info)\n");
      fmt::print(stderr, "\tnem_spread  [command file] [<-p Proc> <-r raid #>]\n");
      exit(1);
      break;
    case 'V':
      fmt::print("{} version {}\n", UTIL_NAME, VER_STR);
      exit(0);
      break;
    case 'p': /* Which proc to use? Also for compatibility */ break;
    case 'r': /* raid number.  Seems to be unused; left around for compatibility */ break;
    case 's': /* Start with processor <x> */ sscanf(optarg, "%d", &start_proc); break;
    case 'n': /* Number of processors to output files for */ sscanf(optarg, "%d", &num_proc); break;
    case '6':
    case '4':
      force_64_bit = true; /* Force storing output mesh using 64bit integers */
      break;
    case 'S': /* Number of subcycles to use (see below) */ sscanf(optarg, "%d", &subcycles); break;
    case 'c': /* Which cycle to spread (see below) */ sscanf(optarg, "%d", &cycle); break;
    }
  }

  if (optind >= argc) {
    salsa_cmd_file = "nem_spread.inp";
  }
  else {
    salsa_cmd_file = argv[optind];
  }

  fmt::print("{} version {}\n", UTIL_NAME, VER_STR);

  /* initialize some variables */
  ExoFile[0]      = '\0';
  Exo_Res_File[0] = '\0';
  Debug_Flag      = -1;

  Num_Nod_Var  = -1;
  Num_Elem_Var = -1;
  Num_Glob_Var = -1;
  Num_Nset_Var = -1;
  Num_Sset_Var = -1;

  PIO_Info.Dsk_List_Cnt        = -1;
  PIO_Info.Num_Dsk_Ctrlrs      = -1;
  PIO_Info.PDsk_Add_Fact       = -1;
  PIO_Info.Zeros               = -1;
  PIO_Info.NoSubdirectory      = 0;
  PIO_Info.Par_Dsk_Root[0]     = '\0';
  PIO_Info.Par_Dsk_SubDirec[0] = '\0';
  PIO_Info.Staged_Writes       = true;

  // Read the ASCII input file and get the name of the mesh file
  // so we can determine the floating point and integer word sizes
  // needed to instantiate the templates...
  if (read_mesh_file_name(salsa_cmd_file) < 0) {
    static char yo[] = "nem_spread";
    fmt::print(stderr,
               "{} ERROR: Could not read in the the I/O command file"
               " \"{}\"!\n",
               yo, salsa_cmd_file);
    exit(1);
  }

  // Open the mesh file and determine word sizes...
  int   io_ws  = 0;
  int   cpu_ws = sizeof(float);
  float version;

  int exoid = ex_open(ExoFile.c_str(), EX_READ, &cpu_ws, &io_ws, &version);

  // See if any 64-bit integers stored on database...
  int int64api = 0;
  int int64db  = ex_int64_status(exoid) & EX_ALL_INT64_DB;
  if (int64db != 0) {
    int64api     = EX_ALL_INT64_API;
    force_64_bit = true;
  }

  int status;
  if (io_ws == 4) {
    if (int64api != 0) {
      NemSpread<float, int64_t> spreader;
      spreader.io_ws        = io_ws;
      spreader.int64db      = int64db;
      spreader.int64api     = int64api;
      spreader.force64db    = force_64_bit;
      spreader.Proc_Info[4] = start_proc;
      spreader.Proc_Info[5] = num_proc;
      status                = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
    else {
      NemSpread<float, int> spreader;
      spreader.io_ws        = io_ws;
      spreader.int64db      = int64db;
      spreader.int64api     = int64api;
      spreader.force64db    = force_64_bit;
      spreader.Proc_Info[4] = start_proc;
      spreader.Proc_Info[5] = num_proc;
      status                = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
  }
  else {
    if (int64api != 0) {
      NemSpread<double, int64_t> spreader;
      spreader.io_ws        = io_ws;
      spreader.int64db      = int64db;
      spreader.int64api     = int64api;
      spreader.force64db    = force_64_bit;
      spreader.Proc_Info[4] = start_proc;
      spreader.Proc_Info[5] = num_proc;
      status                = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
    else {
      NemSpread<double, int> spreader;
      spreader.io_ws        = io_ws;
      spreader.int64db      = int64db;
      spreader.int64api     = int64api;
      spreader.force64db    = force_64_bit;
      spreader.Proc_Info[4] = start_proc;
      spreader.Proc_Info[5] = num_proc;
      status                = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
  }
  double g_end_t = second() - g_start_t;
  fmt::print("The average run time was: {:4f}s\n", g_end_t);

  ex_close(exoid);
  add_to_log(argv[0], g_end_t);
  return status;
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

  safe_free((void **)&(spreader.Proc_Ids));
  safe_free(reinterpret_cast<void **>(&(PIO_Info.RDsk_List)));

  for (int i = 0; i < spreader.Proc_Info[0]; i++) {
    safe_free((void **)&spreader.globals.GNodes[i]);
    if (spreader.globals.Elem_Type != nullptr) {
      safe_free((void **)&spreader.globals.Elem_Type[i]);
    }
    if (spreader.globals.Proc_Global_Node_Id_Map != nullptr) {
      safe_free((void **)&spreader.globals.Proc_Global_Node_Id_Map[i]);
    }
    safe_free((void **)&spreader.globals.Proc_SS_Ids[i]);
    safe_free((void **)&spreader.globals.Proc_SS_GEMap_List[i]);
    safe_free((void **)&spreader.globals.Proc_NS_Ids[i]);
    safe_free((void **)&spreader.globals.Proc_NS_GNMap_List[i]);
    safe_free((void **)&spreader.globals.Proc_Nodes_Per_Elem[i]);
    safe_free((void **)&spreader.globals.GElem_Blks[i]);
    safe_free((void **)&spreader.globals.Proc_Global_Elem_Id_Map[i]);
  }
  safe_free((void **)&spreader.globals.Elem_Type);
  safe_free((void **)&spreader.globals.GNodes);
  return 0;
}
