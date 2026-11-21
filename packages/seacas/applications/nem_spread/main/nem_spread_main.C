/*
 * Copyright(C) 1999-2025 National Technology & Engineering Solutions
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
#include "ps_pario_const.h" // for Parallel_IO
#include "rf_allo.h"        // for safe_free
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

extern int read_mesh_file_name(const char *filename);

template <typename T, typename INT>
int nem_spread(NemSpread<T, INT> &spreader, const char *salsa_cmd_file, int subcycles, int cycle);

int main(int argc, char *argv[])
{
  const char *salsa_cmd_file;
  int         c;

  double g_start_t           = second();
  bool   force_64_bit        = false;
  int    start_proc          = 0;
  int    num_proc            = 0;
  int    subcycles           = 0;
  int    cycle               = -1;
  int    selected_change_set = 0;

  while ((c = getopt(argc, argv, "64Vhp:r:s:n:S:c:C:")) != -1) {
    switch (c) {
    case 'h':
      fmt::print(stderr, " usage:\n");
      fmt::print(stderr,
                 "\tnem_spread  [-s <start_proc>] [-n <num_proc>] [-S <subcycles> -c <cycle>] "
                 "[-C change_set_#] [command_file]\n");
      fmt::print(stderr, "\t\tDecompose for processors <start_proc> to <start_proc>+<num_proc>\n");
      fmt::print(stderr, "\t\tDecompose for cycle <cycle> of <subcycle> groups\n");
      fmt::print(stderr, "\t\tRead from change_set `#` (1-based) if specified.\n");
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
    case 'C': /* change_set index <x> */ sscanf(optarg, "%d", &selected_change_set); break;
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

  static char yo[] = "nem_spread";
  // Read the ASCII input file and get the name of the mesh file
  // so we can determine the floating point and integer word sizes
  // needed to instantiate the templates...
  if (read_mesh_file_name(salsa_cmd_file) < 0) {
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
  if (exoid <= 0) {
    fmt::print(stderr, "{} ERROR: Could not open the mesh file '{}'!\n", yo, ExoFile);
    exit(1);
  }

  // Determine whether there are any change sets in file.
  // If there are, check whether user specified a specific
  // change set index and set to that one (if valid), or
  // if not specified, set to the first.
  exoid = check_change_sets(exoid, selected_change_set);
  if (exoid == 0) {
    exit(1);
  }

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
      spreader.selected_change_set = selected_change_set;
      spreader.io_ws               = io_ws;
      spreader.int64db             = int64db;
      spreader.int64api            = int64api;
      spreader.force64db           = force_64_bit;
      spreader.Proc_Info[4]        = start_proc;
      spreader.Proc_Info[5]        = num_proc;
      status                       = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
    else {
      NemSpread<float, int> spreader;
      spreader.selected_change_set = selected_change_set;
      spreader.io_ws               = io_ws;
      spreader.int64db             = int64db;
      spreader.int64api            = int64api;
      spreader.force64db           = force_64_bit;
      spreader.Proc_Info[4]        = start_proc;
      spreader.Proc_Info[5]        = num_proc;
      status                       = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
  }
  else {
    if (int64api != 0) {
      NemSpread<double, int64_t> spreader;
      spreader.selected_change_set = selected_change_set;
      spreader.io_ws               = io_ws;
      spreader.int64db             = int64db;
      spreader.int64api            = int64api;
      spreader.force64db           = force_64_bit;
      spreader.Proc_Info[4]        = start_proc;
      spreader.Proc_Info[5]        = num_proc;
      status                       = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
    else {
      NemSpread<double, int> spreader;
      spreader.selected_change_set = selected_change_set;
      spreader.io_ws               = io_ws;
      spreader.int64db             = int64db;
      spreader.int64api            = int64api;
      spreader.force64db           = force_64_bit;
      spreader.Proc_Info[4]        = start_proc;
      spreader.Proc_Info[5]        = num_proc;
      status                       = nem_spread(spreader, salsa_cmd_file, subcycles, cycle);
    }
  }
  double g_end_t = second() - g_start_t;
  fmt::print("The average run time was: {:4f}s\n", g_end_t);

  ex_close(exoid);
  add_to_log(argv[0], g_end_t);
  return status;
}
