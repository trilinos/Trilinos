/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "exodusII.h" // for ex_close, etc
#include "fmt/ostream.h"
#include "nem_spread.h"     // for NemSpread, etc
#include "ps_pario_const.h" // for PIO_Info, etc
#include "rf_allo.h"        // for array_alloc, safe_free
#include "rf_io_const.h"    // for Exo_Res_File, ExoFile, etc
#include <cassert>          // for assert
#include <climits>          // for INT_MAX
#include <cstddef>          // for size_t
#include <cstdio>           // for stderr, nullptr, etc
#include <cstdlib>          // for exit, free, malloc
#include <string>
#include <unistd.h>
#include <vector> // for vector

namespace {
  int get_free_descriptor_count();

  template <typename INT>
  size_t find_gnode_inter(INT *intersect, size_t num_g_nodes, INT *glob_vec, size_t num_int_nodes,
                          size_t num_bor_nodes, size_t num_ext_nodes, INT *loc_vec);
} // namespace

/*****************************************************************************/
/*****************************************************************************/
template void NemSpread<double, int>::read_restart_params();
template void NemSpread<float, int>::read_restart_params();

template void NemSpread<double, int64_t>::read_restart_params();
template void NemSpread<float, int64_t>::read_restart_params();

template <typename T, typename INT> void NemSpread<T, INT>::read_restart_params()

/* Function which reads the restart variable parameters for the EXODUS II
 * database which contains the results information. Allocate necessary
 * memory on each processor.
 *
 *----------------------------------------------------------------------------
 *
 * Functions called:
 *
 * compare_mesh_param -- function which checks that parameters in
 *                       the restart EXODUS II file are the same as in
 *                       the mesh EXODUS II file
 * read_var_param -- function which reads the time indices, number
 *                   of variables, and their names from the restart file
 *
 *----------------------------------------------------------------------------
 */

{
  int   exoid;
  int   cpu_ws = 0;
  float vers;
  int   max_name_length = 0;

  /* Open the ExodusII file */
  cpu_ws   = io_ws;
  int mode = EX_READ | int64api;
  if ((exoid = ex_open(Exo_Res_File.c_str(), mode, &cpu_ws, &io_ws, &vers)) < 0) {
    fmt::print(stderr, "{}: Could not open file {} for restart info\n", __func__,
               Exo_Res_File.c_str());
    exit(1);
  }

  max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  ex_set_max_name_length(exoid, max_name_length);

  /*
   * Just do a rudimentary check to figure out if the mesh parameters
   * in the results file are the same as the mesh parameters in the
   * mesh file.
   */
  if (ExoFile != Exo_Res_File) {
    if (!compare_mesh_param(exoid)) {
      fmt::print(stderr,
                 "{}: Mesh parameters in mesh and result files"
                 " differ\n",
                 __func__);
      exit(1);
    }
  }

  /* get the time, and the variable names */
  if (read_var_param(exoid, max_name_length) < 0) {
    fmt::print(stderr, "{}: Error occurred while reading variable parameters\n", __func__);
    exit(1);
  }

  /* Close the ExodusII file */
  ex_close(exoid);
}

template void NemSpread<double, int>::read_restart_data();
template void NemSpread<float, int>::read_restart_data();

template void NemSpread<double, int64_t>::read_restart_data();
template void NemSpread<float, int64_t>::read_restart_data();

template <typename T, typename INT> void NemSpread<T, INT>::read_restart_data()

/* Function which reads the restart variable data from the EXODUS II
 * database which contains the results information. Then distribute
 * it to the processors, and write it to the parallel exodus files.
 *
 *----------------------------------------------------------------------------
 *
 * Functions called:
 *
 * read_vars -- function which reads the variable values from the restart
 *              file, and then distributes them to the processors
 *
 * write_var_timestep -- function which writes out the variables for a
 *                       to a parallel ExodusII file.
 *
 *----------------------------------------------------------------------------
 */

{
  /* need to get the element block ids and counts */
  std::vector<INT> eb_ids_global(globals.Num_Elem_Blk);
  std::vector<INT> eb_cnts_global(globals.Num_Elem_Blk);
  std::vector<INT> ss_ids_global(globals.Num_Side_Set);
  std::vector<INT> ss_cnts_global(globals.Num_Side_Set);
  std::vector<INT> ns_ids_global(globals.Num_Node_Set);
  std::vector<INT> ns_cnts_global(globals.Num_Node_Set);

  INT ***eb_map_ptr    = nullptr;
  INT ** eb_cnts_local = nullptr;
  int    exoid         = 0;
  int *  par_exoid     = nullptr;

  float       vers;
  std::string cTemp;

  /* computing precision should be the same as the database precision
   *
   * EXCEPTION: if the io_ws is smaller than the machine precision,
   * ie - database with io_ws == 4 on a Cray (sizeof(float) == 8),
   * then the cpu_ws must be the machine precision.
   */
  int cpu_ws;
  if (io_ws < static_cast<int>(sizeof(float))) {
    cpu_ws = sizeof(float);
  }
  else {
    cpu_ws = io_ws;
  }

  /* Open the ExodusII file */
  {
    cpu_ws   = io_ws;
    int mode = EX_READ | int64api;
    if ((exoid = ex_open(Exo_Res_File.c_str(), mode, &cpu_ws, &io_ws, &vers)) < 0) {
      fmt::print(stderr, "{}: Could not open file {} for restart info\n", __func__,
                 Exo_Res_File.c_str());
      exit(1);
    }
  }

  /* allocate space for the global variables */
  Restart_Info.Glob_Vals.resize(Restart_Info.NVar_Glob);

  if (Restart_Info.NVar_Elem > 0) {

    /* allocate storage space */
    Restart_Info.Elem_Vals.resize(Proc_Info[2]);

    /* now allocate storage for the values */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Elem *
                          (globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc]);
      Restart_Info.Elem_Vals[iproc].resize(array_size);
    }

    /*
     * at this point, I need to broadcast the global element block ids
     * and counts to the processors. I know that this is redundant data
     * since they will all receive this information in read_mesh, but
     * the variables which contain that information are static in
     * el_exoII_io.c, and cannot be used here. So, take a second and
     * broadcast all of this out.
     *
     * I want to do this here so that it is done only once no matter
     * how many time steps are retrieved
     */

    /* Get the Element Block IDs from the input file */
    if (ex_get_ids(exoid, EX_ELEM_BLOCK, eb_ids_global.data()) < 0) {
      fmt::print(stderr, "{}: unable to get element block IDs", __func__);
      exit(1);
    }

    /* Get the count of elements in each element block */
    for (int cnt = 0; cnt < globals.Num_Elem_Blk; cnt++) {
      char blk_name[MAX_STR_LENGTH];
      if (ex_get_block(exoid, EX_ELEM_BLOCK, eb_ids_global[cnt], blk_name, &(eb_cnts_global[cnt]),
                       nullptr, nullptr, nullptr, nullptr) < 0) {
        fmt::print(stderr, "{}: unable to get element count for block id {}", __func__,
                   (size_t)eb_ids_global[cnt]);
        exit(1);
      }
    }

    /*
     * in order to speed up finding matches in the global element
     * number map, set up an array of pointers to the start of
     * each element block's global element number map. That way
     * only entries for the current element block have to be searched
     */
    eb_map_ptr = (INT ***)array_alloc(__FILE__, __LINE__, 2, Proc_Info[2], globals.Num_Elem_Blk,
                                      sizeof(INT *));
    if (!eb_map_ptr) {
      fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
      exit(1);
    }
    eb_cnts_local =
        (INT **)array_alloc(__FILE__, __LINE__, 2, Proc_Info[2], globals.Num_Elem_Blk, sizeof(INT));
    if (!eb_cnts_local) {
      fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
      exit(1);
    }

    /*
     * for now, assume that element blocks have been
     * stored in the same order as the global blocks
     */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      int    ifound = 0;
      size_t offset = 0;
      int    ilocal;
      for (int cnt = 0; cnt < globals.Num_Elem_Blk; cnt++) {
        for (ilocal = ifound; ilocal < globals.Proc_Num_Elem_Blk[iproc]; ilocal++) {
          if (globals.Proc_Elem_Blk_Ids[iproc][ilocal] == eb_ids_global[cnt]) {
            break;
          }
        }

        if (ilocal < globals.Proc_Num_Elem_Blk[iproc]) {
          eb_map_ptr[iproc][cnt]    = &globals.GElems[iproc][offset];
          eb_cnts_local[iproc][cnt] = globals.Proc_Num_Elem_In_Blk[iproc][ilocal];
          offset += globals.Proc_Num_Elem_In_Blk[iproc][ilocal];
          ifound = ilocal; /* don't search the same part of the list over */
        }
        else {
          eb_map_ptr[iproc][cnt]    = nullptr;
          eb_cnts_local[iproc][cnt] = 0;
        }
      }
    }

  } /* End: "if (Restart_Info.NVar_Elem > 0 )" */

  if (Restart_Info.NVar_Node > 0) {
    /* allocate storage space */
    Restart_Info.Node_Vals.resize(Proc_Info[2]);

    /* now allocate storage for the values */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Node *
                          (globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                           globals.Num_External_Nodes[iproc]);
      Restart_Info.Node_Vals[iproc].resize(array_size);
    }
  }

  if (Restart_Info.NVar_Sset > 0) {

    /* allocate storage space */
    Restart_Info.Sset_Vals.resize(Proc_Info[2]);

    /* now allocate storage for the values */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Sset * globals.Proc_SS_Elem_List_Length[iproc];

      Restart_Info.Sset_Vals[iproc].resize(array_size);
    }

    /*
     * at this point, I need to broadcast the ids and counts to the
     * processors. I know that this is redundant data since they will
     * all receive this information in read_mesh, but the variables
     * which contain that information are static in el_exoII_io.c, and
     * cannot be used here. So, take a second and broadcast all of
     * this out.
     *
     * I want to do this here so that it is done only once no matter
     * how many time steps are retrieved
     */

    /* Get the Sideset IDs from the input file */
    if (ex_get_ids(exoid, EX_SIDE_SET, ss_ids_global.data()) < 0) {
      fmt::print(stderr, "{}: unable to get sideset IDs", __func__);
      exit(1);
    }

    /* Get the count of elements in each sideset */
    for (int cnt = 0; cnt < globals.Num_Side_Set; cnt++) {
      if (ex_get_set_param(exoid, EX_SIDE_SET, ss_ids_global[cnt], &(ss_cnts_global[cnt]),
                           nullptr) < 0) {
        fmt::print(stderr, "{}: unable to get element count for sideset id {}", __func__,
                   (size_t)ss_ids_global[cnt]);
        exit(1);
      }
    }
  } /* End: "if (Restart_Info.NVar_Sset > 0 )" */

  if (Restart_Info.NVar_Nset > 0) {

    /* allocate storage space */
    Restart_Info.Nset_Vals.resize(Proc_Info[2]);

    /* now allocate storage for the values */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      size_t array_size = Restart_Info.NVar_Nset * globals.Proc_NS_List_Length[iproc];
      Restart_Info.Nset_Vals[iproc].resize(array_size);
    }

    /*
     * at this point, I need to broadcast the ids and counts to the
     * processors. I know that this is redundant data since they will
     * all receive this information in read_mesh, but the variables
     * which contain that information are static in el_exoII_io.c, and
     * cannot be used here. So, take a second and broadcast all of
     * this out.
     *
     * I want to do this here so that it is done only once no matter
     * how many time steps are retrieved
     */

    /* Get the Nodeset IDs from the input file */
    if (ex_get_ids(exoid, EX_NODE_SET, ns_ids_global.data()) < 0) {
      fmt::print(stderr, "{}: unable to get nodeset IDs", __func__);
      exit(1);
    }

    /* Get the count of elements in each nodeset */
    for (int cnt = 0; cnt < globals.Num_Node_Set; cnt++) {
      if (ex_get_set_param(exoid, EX_NODE_SET, ns_ids_global[cnt], &(ns_cnts_global[cnt]),
                           nullptr) < 0) {
        fmt::print(stderr, "{}: unable to get element count for nodeset id {}", __func__,
                   (size_t)ns_ids_global[cnt]);
        exit(1);
      }
    }
  } /* End: "if (Restart_Info.NVar_Nset > 0 )" */

  /*
   * NOTE: A possible place to speed this up would be to
   * get the global node and element lists here, and broadcast
   * them out only once.
   */

  par_exoid = (int *)malloc(Proc_Info[2] * sizeof(int));
  if (par_exoid == nullptr) {
    fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
    exit(1);
  }

  /* See if any '/' in the name.  IF present, isolate the basename of the file */
  size_t found = Output_File_Base_Name.find_last_of('/');
  if (found != std::string::npos) {
    /* There is a path separator.  Get the portion after the
     * separator
     */
    cTemp = Output_File_Base_Name.substr(found + 1);
  }
  else {
    /* No separator; this is already just the basename... */
    cTemp = Output_File_Base_Name;
  }

  if (PIO_Info.Exo_Extension.empty()) {
    cTemp += ".par";
  }
  else {
    cTemp += PIO_Info.Exo_Extension;
  }

  int open_file_count = get_free_descriptor_count();
  if (open_file_count > Proc_Info[5]) {
    fmt::print("All output files opened simultaneously.\n");
    for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
      std::string Parallel_File_Name = gen_par_filename(cTemp, Proc_Ids[iproc], Proc_Info[0]);

      /* Open the parallel Exodus II file for writing */
      cpu_ws   = io_ws;
      int mode = EX_WRITE | int64api | int64db;
      if ((par_exoid[iproc] = ex_open(Parallel_File_Name.c_str(), mode, &cpu_ws, &io_ws, &vers)) <
          0) {
        fmt::print(stderr, "[{}] {} Could not open parallel Exodus II file: {}\n", iproc, __func__,
                   Parallel_File_Name.c_str());
        exit(1);
      }
    }
  }
  else {
    fmt::print("All output files opened one-at-a-time.\n");
  }

  /* Now loop over the number of time steps */
  for (int time_idx = 0; time_idx < Restart_Info.Num_Times; time_idx++) {

    double start_t = second();

    /* read and distribute the variables for this time step */
    if (read_vars(exoid, Restart_Info.Time_Idx[time_idx], eb_ids_global.data(),
                  eb_cnts_global.data(), eb_map_ptr, eb_cnts_local, ss_ids_global.data(),
                  ss_cnts_global.data(), ns_ids_global.data(), ns_cnts_global.data()) < 0) {
      fmt::print(stderr, "{}: Error occurred while reading variables\n", __func__);
      exit(1);
    }
    double end_t = second() - start_t;
    fmt::print("\tTime to read  vars for timestep {}: {} (sec.)\n", (time_idx + 1), end_t);

    start_t = second();
    for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {

      if (open_file_count < Proc_Info[5]) {
        std::string Parallel_File_Name = gen_par_filename(cTemp, Proc_Ids[iproc], Proc_Info[0]);

        /* Open the parallel Exodus II file for writing */
        cpu_ws   = io_ws;
        int mode = EX_WRITE | int64api | int64db;
        if ((par_exoid[iproc] = ex_open(Parallel_File_Name.c_str(), mode, &cpu_ws, &io_ws, &vers)) <
            0) {
          fmt::print(stderr, "[{}] {} Could not open parallel Exodus II file: {}\n", iproc,
                     __func__, Parallel_File_Name.c_str());
          exit(1);
        }
      }

      /*
       * Write out the variable data for the time steps in this
       * block to each parallel file.
       */
      write_var_timestep(par_exoid[iproc], iproc, (time_idx + 1), eb_ids_global.data(),
                         ss_ids_global.data(), ns_ids_global.data());

      if (iproc % 10 == 0 || iproc == Proc_Info[2] - 1) {
        fmt::print("{}", iproc);
      }
      else {
        fmt::print(".");
      }

      if (open_file_count < Proc_Info[5]) {
        if (ex_close(par_exoid[iproc]) == -1) {
          fmt::print(stderr, "[{}] {} Could not close the parallel Exodus II file.\n", iproc,
                     __func__);
          exit(1);
        }
      }
    } /* End "for (iproc=0; iproc <Proc_Info[2]; iproc++)" */

    end_t = second() - start_t;
    fmt::print("\n\tTime to write vars for timestep {}: {} (sec.)\n", (time_idx + 1), end_t);
  }
  if (Restart_Info.NVar_Elem > 0) {
    safe_free((void **)&eb_map_ptr);
    safe_free((void **)&eb_cnts_local);
  }

  /* Close the restart exodus II file */
  if (ex_close(exoid) == -1) {
    fmt::print(stderr, "{}Could not close the restart Exodus II file\n", __func__);
    exit(1);
  }

  if (open_file_count > Proc_Info[5]) {
    for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
      /* Close the parallel exodus II file */
      if (ex_close(par_exoid[iproc]) == -1) {
        fmt::print(stderr, "[{}] {} Could not close the parallel Exodus II file.\n", iproc,
                   __func__);
        exit(1);
      }
    }
  }
  free(par_exoid);
  par_exoid = nullptr;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_var_param(int exoid, int max_name_length)
{
  /* Get the number of time indices contained in the file */
  int ret_int = ex_inquire_int(exoid, EX_INQ_TIME);

  /* see if the user want to get all of the time indices */
  if (Restart_Info.Num_Times == -1) {

    Restart_Info.Num_Times = ret_int;

    if (ret_int > 0) {
      /* allocate array space */
      Restart_Info.Time_Idx.resize(Restart_Info.Num_Times);

      for (int cnt = 0; cnt < Restart_Info.Num_Times; cnt++) {
        Restart_Info.Time_Idx[cnt] = cnt + 1;
      }
    }
  }
  else {
    /* Check to see if the requested indices are valid */
    for (int cnt = 0; cnt < Restart_Info.Num_Times; cnt++) {

      /* if the user wants the last time, then set it */
      if (Restart_Info.Time_Idx[cnt] == 0) {
        Restart_Info.Time_Idx[cnt] = ret_int;
      }

      if (Restart_Info.Time_Idx[cnt] > ret_int) {
        fmt::print(stderr, "{}: Requested time index, {}, out of range.\n", __func__,
                   Restart_Info.Time_Idx[cnt]);
        fmt::print(stderr, "{}: Valid time indices in {} are from 1 to {}.\n", __func__,
                   Exo_Res_File.c_str(), ret_int);
        return -1;
      }
    }
  }

  /* if there are not any time steps, then return here without an error */
  if (Restart_Info.Num_Times == 0) {
    Restart_Info.Flag      = 0;
    Restart_Info.NVar_Glob = 0;
    Restart_Info.NVar_Node = 0;
    Restart_Info.NVar_Elem = 0;
    return 0;
  }

  /***************** Global Variables ********************/
  if (ex_get_variable_param(exoid, EX_GLOBAL, &(Restart_Info.NVar_Glob)) < 0) {
    fmt::print(stderr, "{}: Could not get global variable parameter from file\n", __func__);
    return -1;
  }

  /* allocate space for the global variable names */
  if (Restart_Info.NVar_Glob > 0) {
    Restart_Info.GV_Name = (char **)array_alloc(__FILE__, __LINE__, 2, Restart_Info.NVar_Glob,
                                                max_name_length + 1, sizeof(char));

    /* get the global variable names */
    if (ex_get_variable_names(exoid, EX_GLOBAL, Restart_Info.NVar_Glob, Restart_Info.GV_Name) < 0) {
      fmt::print(stderr, "{}: Could not get global variable names from file\n", __func__);
      return -1;
    }
  }

  /***************** Elemental Variables ********************/
  if (ex_get_variable_param(exoid, EX_ELEM_BLOCK, &(Restart_Info.NVar_Elem)) < 0) {
    fmt::print(stderr, "{}: Could not get elemental variable param from file\n", __func__);
    return -1;
  }

  /* allocate space for the elemental variable names */
  if (Restart_Info.NVar_Elem > 0) {
    Restart_Info.EV_Name = (char **)array_alloc(__FILE__, __LINE__, 2, Restart_Info.NVar_Elem,
                                                max_name_length + 1, sizeof(char));

    /* get the elemental variable names */
    if (ex_get_variable_names(exoid, EX_ELEM_BLOCK, Restart_Info.NVar_Elem, Restart_Info.EV_Name) <
        0) {
      fmt::print(stderr, "{}: Could not get elemental variable names from file\n", __func__);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GElem_TT.resize(globals.Num_Elem_Blk * Restart_Info.NVar_Elem);

    check_exodus_error(ex_get_truth_table(exoid, EX_ELEM_BLOCK, globals.Num_Elem_Blk,
                                          Restart_Info.NVar_Elem, Restart_Info.GElem_TT.data()),
                       "ex_get_truth_table");
  }

  /******************* Nodal Variables **********************/
  if (ex_get_variable_param(exoid, EX_NODAL, &(Restart_Info.NVar_Node)) < 0) {
    fmt::print(stderr, "{}: Could not get nodal variable param from file\n", __func__);
    return -1;
  }

  /* allocate space for the nodal variable names */
  if (Restart_Info.NVar_Node > 0) {
    Restart_Info.NV_Name = (char **)array_alloc(__FILE__, __LINE__, 2, Restart_Info.NVar_Node,
                                                max_name_length + 1, sizeof(char));

    /* get the nodal variable names */
    if (ex_get_variable_names(exoid, EX_NODAL, Restart_Info.NVar_Node, Restart_Info.NV_Name) < 0) {
      fmt::print(stderr, "{}: Could not get nodal variable names from file\n", __func__);
      return -1;
    }
  }

  /******************* Sideset Variables **********************/
  if (ex_get_variable_param(exoid, EX_SIDE_SET, &(Restart_Info.NVar_Sset)) < 0) {
    fmt::print(stderr, "{}: Could not get sideset variable param from file\n", __func__);
    return -1;
  }

  /* allocate space for the variable names */
  if (Restart_Info.NVar_Sset > 0) {
    Restart_Info.SSV_Name = (char **)array_alloc(__FILE__, __LINE__, 2, Restart_Info.NVar_Sset,
                                                 max_name_length + 1, sizeof(char));

    /* get the variable names */
    if (ex_get_variable_names(exoid, EX_SIDE_SET, Restart_Info.NVar_Sset, Restart_Info.SSV_Name) <
        0) {
      fmt::print(stderr, "{}: Could not get sideset variable names from file\n", __func__);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GSset_TT.resize(globals.Num_Side_Set * Restart_Info.NVar_Sset);

    check_exodus_error(ex_get_truth_table(exoid, EX_SIDE_SET, globals.Num_Side_Set,
                                          Restart_Info.NVar_Sset, Restart_Info.GSset_TT.data()),
                       "ex_get_truth_table");
  }

  /******************* Nodeset Variables **********************/
  if (ex_get_variable_param(exoid, EX_NODE_SET, &(Restart_Info.NVar_Nset)) < 0) {
    fmt::print(stderr, "{}: Could not get nodeset variable param from file\n", __func__);
    return -1;
  }

  /* allocate space for the variable names */
  if (Restart_Info.NVar_Nset > 0) {
    Restart_Info.NSV_Name = (char **)array_alloc(__FILE__, __LINE__, 2, Restart_Info.NVar_Nset,
                                                 max_name_length + 1, sizeof(char));

    /* get the variable names */
    if (ex_get_variable_names(exoid, EX_NODE_SET, Restart_Info.NVar_Nset, Restart_Info.NSV_Name) <
        0) {
      fmt::print(stderr, "{}: Could not get nodeset variable names from file\n", __func__);
      return -1;
    }

    /* and get the truth table */
    Restart_Info.GNset_TT.resize(globals.Num_Node_Set * Restart_Info.NVar_Nset);

    check_exodus_error(ex_get_truth_table(exoid, EX_NODE_SET, globals.Num_Node_Set,
                                          Restart_Info.NVar_Nset, Restart_Info.GNset_TT.data()),
                       "ex_get_var_tab");
  }

#ifdef DEBUG
  if (Debug_Flag >= 2) {
    fmt::print("\n\nRestart Parameters:\n");
    fmt::print("\tNumber of time indices: {}\n", Restart_Info.Num_Times);
    for (int cnt = 0; cnt < Restart_Info.Num_Times; cnt++)
      fmt::print("\t\tTime index: {}\n", Restart_Info.Time_Idx[cnt]);
    fmt::print("\tNumber of global variables: {}\n", Restart_Info.NVar_Glob);
    for (int cnt = 0; cnt < Restart_Info.NVar_Glob; cnt++)
      fmt::print("\t\tGlobal variable {}: {}\n", (cnt + 1), Restart_Info.GV_Name[cnt]);
    fmt::print("\tNumber of elental variables: {}\n", Restart_Info.NVar_Elem);
    for (int cnt = 0; cnt < Restart_Info.NVar_Elem; cnt++)
      fmt::print("\t\tElemental variable {}: {}\n", (cnt + 1), Restart_Info.EV_Name[cnt]);
    fmt::print("\tNumber of nodal variables: {}\n", Restart_Info.NVar_Node);
    for (int cnt = 0; cnt < Restart_Info.NVar_Node; cnt++)
      fmt::print("\t\tNodal variable {}: {}\n", (cnt + 1), Restart_Info.NV_Name[cnt]);
  }
#endif

  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_vars(int exoid, int index, INT *eb_ids, INT *eb_cnts, INT ***eb_map_ptr,
                                 INT **eb_cnts_local, INT *ss_ids, INT *ss_cnts, INT *ns_ids,
                                 INT *ns_cnts)
{
  /* first read the time */
  if (ex_get_time(exoid, index, &Restart_Info.Time) < 0) {
    fmt::print(stderr, "{}: ERROR, unable to get time for restart index {}!\n", __func__, index);
    return -1;
  }

  /***************** Global Variables ********************/
  /* allocate space for the global variables */
  if (Restart_Info.NVar_Glob > 0) {
    /* get the global variables */
    if (ex_get_var(exoid, index, EX_GLOBAL, 1, 1, Restart_Info.NVar_Glob,
                   Restart_Info.Glob_Vals.data()) < 0) {
      fmt::print(stderr, "{}: Could not get global variables from file\n", __func__);
      return -1;
    }
  }

  if (Restart_Info.NVar_Elem > 0) {
    fmt::print("Reading {} element variables...\n", Restart_Info.NVar_Elem);
    if (read_elem_vars(exoid, index, eb_ids, eb_cnts, eb_map_ptr, eb_cnts_local) < 0) {
      fmt::print(stderr, "{}: Error distributing elemental variables.\n", __func__);
      return -1;
    }
  }

  if (Restart_Info.NVar_Node > 0) {
    fmt::print("Reading {} nodal variables...\n", Restart_Info.NVar_Node);
    if (read_nodal_vars(exoid, index) < 0) {
      fmt::print(stderr, "{}: Error distributing nodal variables.\n", __func__);
      return -1;
    }
  }

  if (Restart_Info.NVar_Sset > 0) {
    fmt::print("Reading {} sideset variables...\n", Restart_Info.NVar_Sset);
    if (read_sset_vars(exoid, index, ss_ids, ss_cnts) < 0) {
      fmt::print(stderr, "{}: Error distributing sideset variables.\n", __func__);
      return -1;
    }
  }

  if (Restart_Info.NVar_Nset > 0) {
    fmt::print("Reading {} nodeset variables...\n", Restart_Info.NVar_Nset);
    if (read_nset_vars(exoid, index, ns_ids, ns_cnts) < 0) {
      fmt::print(stderr, "{}: Error distributing nodeset variables.\n", __func__);
      return -1;
    }
  }

  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_elem_vars(int exoid, int index, INT *eb_ids, INT *eb_cnts,
                                      INT ***eb_map_ptr, INT **eb_cnts_local)
{

  /* to speed up searches, keep track of element blocks offset on each proc */
  std::vector<INT> local_offset(Proc_Info[2]);

  /* loop over the number of element blocks */
  INT eb_offset = 0;
  for (int iblk = 0; iblk < globals.Num_Elem_Blk; iblk++) {
    read_elem_vars_1(exoid, index, eb_ids, eb_cnts, eb_map_ptr, eb_cnts_local, iblk, eb_offset,
                     local_offset.data());

    /* need to keep track of this for the element number map */
    eb_offset += eb_cnts[iblk];

    /* need to set up local offsets for next block */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      local_offset[iproc] += eb_cnts_local[iproc][iblk];
    }

  } /* End "for (iblk = 0; iblk < globals.Num_Elem_Blk; iblk++)" */

  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_elem_vars_1(int exoid, int index, INT *eb_ids, INT *eb_cnts,
                                        INT ***eb_map_ptr, INT **eb_cnts_local, int iblk,
                                        int eb_offset, INT *local_offset)
{
  /* Allocate memory for temporary storage */
  std::vector<T> vals(eb_cnts[iblk]);

  /* now loop over each variable */
  for (int ivar = 0; ivar < Restart_Info.NVar_Elem; ivar++) {

    /* check if this variable exists for this element block */
    if (Restart_Info.GElem_TT[iblk * Restart_Info.NVar_Elem + ivar]) {

      /*
       * Read in the specified element variable values and their associated
       * global FEM element numbers.
       */

      check_exodus_error(ex_get_var(exoid, index, EX_ELEM_BLOCK, (ivar + 1), eb_ids[iblk],
                                    eb_cnts[iblk], vals.data()),
                         "ex_get_var");

      /*
       * Find out which FEM elements belong on this processor and copy
       * them to the restart vector.
       */
      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

        /* check to see if this element block needs this variable */
        if (Restart_Info.GElem_TT[iblk * Restart_Info.NVar_Elem + ivar]) {

          /* calculate the offset for this variable */
          size_t var_offset =
              ivar * (globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc]);

          INT *  elem_map = eb_map_ptr[iproc][iblk];
          size_t num_elem = eb_cnts_local[iproc][iblk];

          for (size_t i1 = 0; i1 < num_elem; i1++) {
            size_t elem_loc = var_offset + i1 + local_offset[iproc];

            Restart_Info.Elem_Vals[iproc][elem_loc] = vals[elem_map[i1] - eb_offset];
          }
        }
      }
    } /* End "if (Restart_Info.GElem_TT[...])" */
  }
  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_sset_vars(int exoid, int index, INT *ss_ids, INT *ss_cnts)
{
  /* loop over the number of side sets */
  for (int iset = 0; iset < globals.Num_Side_Set; iset++) {
    read_sset_vars_1(exoid, index, ss_ids, ss_cnts, iset);
  }
  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_sset_vars_1(int exoid, int index, INT *ss_ids, INT *ss_cnts, int iset)
{
  /* Allocate memory for temporary storage */
  std::vector<T> vals(ss_cnts[iset]);

  /* now loop over each variable */
  for (int ivar = 0; ivar < Restart_Info.NVar_Sset; ivar++) {

    /* check if this variable exists for this set */
    if (Restart_Info.GSset_TT[iset * Restart_Info.NVar_Sset + ivar]) {

      /* Read in the specified variable values */
      check_exodus_error(ex_get_var(exoid, index, EX_SIDE_SET, (ivar + 1), ss_ids[iset],
                                    ss_cnts[iset], vals.data()),
                         "ex_get_var");

      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        size_t ss_offset  = 0;
        size_t var_offset = ivar * globals.Proc_SS_Elem_List_Length[iproc];
        for (int i = 0; i < globals.Proc_Num_Side_Sets[iproc]; i++) {
          if (globals.Proc_SS_Ids[iproc][i] == ss_ids[iset]) {

            size_t num_elem = globals.Proc_SS_Elem_Count[iproc][i];
            for (size_t i1 = 0; i1 < num_elem; i1++) {
              INT gelem_loc = globals.Proc_SS_GEMap_List[iproc][i1 + ss_offset];
              assert(gelem_loc < ss_cnts[iset]);
              Restart_Info.Sset_Vals[iproc][i1 + ss_offset + var_offset] = vals[gelem_loc];
            }
            break;
          }
          ss_offset += globals.Proc_SS_Elem_Count[iproc][i];
        }
      }
    }
  }
  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_nset_vars(int exoid, int index, INT *ns_ids, INT *ns_cnts)
{
  /* loop over the number of node sets */
  for (int iset = 0; iset < globals.Num_Node_Set; iset++) {
    read_nset_vars_1(exoid, index, ns_ids, ns_cnts, iset);
  }
  return 0;
}

template <typename T, typename INT>
int NemSpread<T, INT>::read_nset_vars_1(int exoid, int index, INT *ns_ids, INT *ns_cnts, int iset)
{
  /* Allocate memory for temporary storage */
  std::vector<T> vals(ns_cnts[iset]);

  /* now loop over each variable */
  for (int ivar = 0; ivar < Restart_Info.NVar_Nset; ivar++) {

    /* check if this variable exists for this set */
    if (Restart_Info.GNset_TT[iset * Restart_Info.NVar_Nset + ivar]) {

      /* Read in the specified variable values */
      check_exodus_error(ex_get_var(exoid, index, EX_NODE_SET, (ivar + 1), ns_ids[iset],
                                    ns_cnts[iset], vals.data()),
                         "ex_get_nset_var");

      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        size_t ns_offset  = 0;
        size_t var_offset = ivar * globals.Proc_NS_List_Length[iproc];
        for (int i = 0; i < globals.Proc_Num_Node_Sets[iproc]; i++) {
          if (globals.Proc_NS_Ids[iproc][i] == ns_ids[iset]) {

            size_t num_elem = globals.Proc_NS_Count[iproc][i];
            for (size_t i1 = 0; i1 < num_elem; i1++) {
              INT gelem_loc = globals.Proc_NS_GNMap_List[iproc][i1 + ns_offset];
              assert(gelem_loc < ns_cnts[iset]);
              Restart_Info.Nset_Vals[iproc][i1 + ns_offset + var_offset] = vals[gelem_loc];
            }
            break;
          }
          ns_offset += globals.Proc_NS_Count[iproc][i];
        }
      }
    }
  }
  return 0;
}

template <typename T, typename INT> int NemSpread<T, INT>::read_nodal_vars(int exoid, int index)
{
  /* Allocate memory for temporary storage */
  std::vector<T> vals(globals.Num_Node);

  /* Loop over each auxiliary variable */
  for (int var_num = 0; var_num < Restart_Info.NVar_Node; var_num++) {
    /*
     * Read in the specified nodal variable values and their associated
     * global FEM node numbers.
     */
    check_exodus_error(
        ex_get_var(exoid, index, EX_NODAL, (var_num + 1), 1, globals.Num_Node, vals.data()),
        "ex_get_var");

    /*
     * Find out which FEM nodes belong on this processor and copy
     * them to the restart vector.
     */
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

      /* calculate the offset for this variable */
      size_t loc_count = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                         globals.Num_External_Nodes[iproc];

      size_t var_offset = var_num * loc_count;

      for (size_t i2 = 0; i2 < loc_count; i2++) {
        size_t node_loc                         = var_offset + i2;
        Restart_Info.Node_Vals[iproc][node_loc] = vals[globals.GNodes[iproc][i2] - 1];
      }
    }

  } /* End "for (var_num = 0; var_num < Restart_Info.NVar_Node; var_num++)" */
  return 0;
}

template <typename T, typename INT> int NemSpread<T, INT>::compare_mesh_param(int exoid)
{
  int ret = 1;

  ex_init_params info{};
  info.title[0] = '\0';
  int error     = ex_get_init_ext(exoid, &info);
  check_exodus_error(error, "ex_get_init");

  /* now check that the parameters match those retrieved from the mesh file */
  if (info.num_dim != globals.Num_Dim) {
    ret = 0;
  }
  else if (static_cast<size_t>(info.num_nodes) != globals.Num_Node) {
    ret = 0;
  }
  else if (static_cast<size_t>(info.num_elem) != globals.Num_Elem) {
    ret = 0;
  }
  else if (info.num_elem_blk != globals.Num_Elem_Blk) {
    ret = 0;
  }
  else if (info.num_node_sets != globals.Num_Node_Set) {
    ret = 0;
  }
  else if (info.num_side_sets != globals.Num_Side_Set) {
    ret = 0;
  }

  return (ret);
}

/*****************************************************************************/
namespace {
  template <typename INT>
  size_t find_gnode_inter(INT *intersect, size_t num_g_nodes, INT *glob_vec, size_t num_int_nodes,
                          size_t num_bor_nodes, size_t num_ext_nodes, INT *loc_vec)

  /*
   * This function assumes that glob_vec is monotonic and that loc_vec is
   * monotonic for each of the internal, border and external node IDs it
   * contains.
   */
  {
    size_t count = 0;

    /* Initialize the intersect vector */
    for (size_t i1 = 0; i1 < num_g_nodes; i1++) {
      intersect[i1] = -1;
    }

    /* Check for the possibility of an intersection */
    size_t min_set1 = glob_vec[0];
    size_t max_set1 = glob_vec[num_g_nodes - 1];

    /* Search through the internal nodes */
    if (num_int_nodes > 0) {
      size_t min_set2 = loc_vec[0];
      size_t max_set2 = loc_vec[num_int_nodes - 1];

      if ((max_set2 >= min_set1) && (min_set2 <= max_set1)) {
        for (size_t i1 = 0, i2 = 0; i1 < num_g_nodes; i1++) {
          while ((i2 < (num_int_nodes - 1)) && (glob_vec[i1] > loc_vec[i2])) {
            i2++;
          }
          if (glob_vec[i1] == loc_vec[i2]) {
            intersect[i1] = i2;
            count++;
          }
        }
      }
    }

    /* Search through the border nodes */
    if (num_bor_nodes > 0) {
      size_t min_set2 = loc_vec[num_int_nodes];
      size_t max_set2 = loc_vec[num_int_nodes + num_bor_nodes - 1];

      size_t offset = num_int_nodes;

      if ((max_set2 >= min_set1) && (min_set2 <= max_set1)) {
        for (size_t i1 = 0, i2 = 0; i1 < num_g_nodes; i1++) {
          while ((i2 < (num_bor_nodes - 1)) && (glob_vec[i1] > loc_vec[offset + i2])) {
            i2++;
          }

          if (glob_vec[i1] == loc_vec[offset + i2]) {
            intersect[i1] = offset + i2;
            count++;
          }
        }
      }
    }

    /* Search through the external nodes */
    if (num_ext_nodes > 0) {
      size_t min_set2 = loc_vec[num_int_nodes + num_bor_nodes];
      size_t max_set2 = loc_vec[num_int_nodes + num_bor_nodes + num_ext_nodes - 1];

      size_t offset = num_int_nodes + num_bor_nodes;

      if ((max_set2 >= min_set1) && (min_set2 <= max_set1)) {
        for (size_t i1 = 0, i2 = 0; i1 < num_g_nodes; i1++) {
          while ((i2 < (num_ext_nodes - 1)) && (glob_vec[i1] > loc_vec[offset + i2])) {
            i2++;
          }

          if (glob_vec[i1] == loc_vec[offset + i2]) {
            intersect[i1] = offset + i2;
            count++;
          }
        }
      }
    }

#ifdef DEBUG
    assert(count < num_g_nodes &&
           "find_gnode_inter ERROR: Catastrophic error in node structure observed\n");
#endif

    return count;
  }

  int get_free_descriptor_count()
  {
/* Returns maximum number of files that one process can have open
 * at one time. (POSIX)
 */
#ifndef _MSC_VER
    int fdmax = sysconf(_SC_OPEN_MAX);
    if (fdmax == -1) {
      /* POSIX indication that there is no limit on open files... */
      fdmax = INT_MAX;
    }
#else
    int fdmax = _getmaxstdio();
#endif
    /* File descriptors are assigned in order (0,1,2,3,...) on a per-process
     * basis.
     *
     * Assume that we have stdin, stdout, stderr, and input exodus
     * file (4 total).
     */

    return fdmax - 4;

    /* Could iterate from 0..fdmax and check for the first EBADF (bad
     * file descriptor) error return from fcntl, but that takes too long
     * and may cause other problems.  There is a isastream(filedes) call
     * on Solaris that can be used for system-dependent code.
     *
     * Another possibility is to do an open and see which descriptor is
     * returned -- take that as 1 more than the current count of open files.
     */
  }
} // namespace
