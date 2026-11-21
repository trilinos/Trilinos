/*
 * Copyright(C) 1999-2021, 2023, 2024, 2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Functions contained in this file:
 *      main()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <cstddef> // for size_t
#include <cstdlib> // for free, exit, malloc
#include <cstring> // for strcmp
#include <iostream>
#include <stdexcept>

#include "add_to_log.h"  // for add_to_log
#include "elb.h"         // for LB_Description<INT>, get_time, etc
#include "elb_allo.h"    // for array_alloc
#include "elb_elem.h"    // for ElementType, ::NULL_EL
#include "elb_err.h"     // for error_report, Gen_Error, etc
#include "elb_exo.h"     // for init_weight_struct, etc
#include "elb_graph.h"   // for generate_graph
#include "elb_inp.h"     // for check_inp_specs, etc
#include "elb_loadbal.h" // for generate_loadbal, etc
#include "elb_nem_slice.h"
#include "elb_output.h" // for write_nemesis, write_vis
#include "fmt/ostream.h"

#ifdef USE_ZOLTAN
#include <mpi.h> // for MPI_Finalize, etc
#endif

int main(int argc, char *argv[])
{
  fmt::print(stderr, "Beginning nem_slice execution.\n");

  double start_time = get_time();

  /* Initialize to just reporting the error */
  error_lev = 1;

  /* check if the user just wants to know the version (or forcing 64-bit mode)*/
  bool int64com = false;
  bool int32com = false;
  int  int64db  = 0;

  for (int cnt = 0; cnt < argc; cnt++) {
    if (strcmp(argv[cnt], "-V") == 0) {
      fmt::print("{} version {}\n", UTIL_NAME, ELB_VERSION);
      exit(0);
    }
    if (strcmp(argv[cnt], "-64") == 0) {
      int64com = true;
    }

    if (strcmp(argv[cnt], "-32") == 0) {
      int32com = true;
    }
  }

  if (argc > 2) {
    /* Get the input mesh file so we can determine the integer size... */
    /* Should be the last item on the command line */
    char *mesh_file_name = argv[argc - 1];

    /* Open file and get the integer size... */
    int   cpu_ws = 0;
    int   io_ws  = 0;
    float vers   = 00.0;
    fmt::print(stderr, "Input Mesh File = '{}\n", mesh_file_name);
    int exoid = ex_open(mesh_file_name, EX_READ, &cpu_ws, &io_ws, &vers);
    if (exoid < 0) {
      std::string error("fatal: unable to open input ExodusII file ");
      error += mesh_file_name;
      Gen_Error(0, error);
      return 0;
    }

    int64db = ex_int64_status(exoid) & EX_ALL_INT64_DB;
    ex_close(exoid);

    ex_opts(EX_VERBOSE);
  }

  int status;
  if (int32com && (int64db != 0)) {
    fmt::print(stderr,
               "Forcing 32-bit integer mode for decomposition even though database is 64-bit.\n");
    NemSlice<int> nem_slice(argc, argv);
    status = nem_slice.run();
  }
  else if ((int64db != 0) || int64com) {
    fmt::print(stderr,
               "Using 64-bit integer mode for decomposition...\n"
               "NOTE: Only 'linear' and 'scattered' methods are supported for 64-bit models\n");
    NemSlice<int64_t> nem_slice(argc, argv);
    status = nem_slice.run();
  }
  else {
    fmt::print(stderr, "Using 32-bit integer mode for decomposition...\n");
    NemSlice<int> nem_slice(argc, argv);
    status = nem_slice.run();
  }

  /* Report any non-fatal errors that may have occurred */
  error_report();

  /* Get ending time */
  double end_time = get_time();
  fmt::print(stderr, "The entire load balance took {} seconds.\n", end_time - start_time);
  add_to_log(argv[0], end_time - start_time);
  return status;
}
