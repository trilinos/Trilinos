/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "copy_string_cpp.h"
#include "el_check_monot.h" // for check_monot
#include "el_elm.h"         // for HEXSHELL, NN_SIDE, etc
#include "exodusII.h"       // for ex_inquire_int, etc
#include "fmt/ostream.h"
#include "nem_spread.h"        // for NemSpread, second, etc
#include "netcdf.h"            // for nc_set_fill, NC_NOFILL
#include "pe_common.h"         // for MAX_CHUNK_SIZE
#include "pe_str_util_const.h" // for string_to_lower
#include "ps_pario_const.h"    // for PIO_Time_Array, etc
#include "rf_allo.h"           // for safe_free, array_alloc
#include "rf_io_const.h"       // for Debug_Flag, etc
#include "rf_util.h"           // for print_line, my_sort
#include "sort_utils.h"        // for gds_qsort
#include "vector_data.h"
#include <cassert> // for assert
#include <climits> // for INT_MAX
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, free
#include <string>  // for string
#include <vector>  // for vector
template <typename T, typename INT> class Globals;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
namespace {
  template <typename T, typename INT>
  int check_sizes(int num_proc, std::vector<INT> &Elem_Blk_Ids, std::vector<INT> &Node_Set_Ids,
                  std::vector<INT> &Side_Set_Ids, Globals<T, INT> &globals)
  {
    /* Do we need 64-bit to write the entity (block/set) ids? */
    bool ids64bit = false;
    if (!Elem_Blk_Ids.empty()) {
      for (int i = 0; i < globals.Num_Elem_Blk && !ids64bit; i++) {
        if (Elem_Blk_Ids[i] >= INT_MAX) {
          ids64bit = true;
        }
      }
    }
    if (!Node_Set_Ids.empty()) {
      for (int i = 0; i < globals.Num_Node_Set && !ids64bit; i++) {
        if (Node_Set_Ids[i] >= INT_MAX) {
          ids64bit = true;
        }
      }
    }
    if (!Side_Set_Ids.empty()) {
      for (int i = 0; i < globals.Num_Side_Set && !ids64bit; i++) {
        if (Side_Set_Ids[i] >= INT_MAX) {
          ids64bit = true;
        }
      }
    }

    /* Do we need 64-bit to write the global node/element maps? */
    bool global_counts_64bit = false;
    if (globals.Num_Node >= INT_MAX) {
      global_counts_64bit = true;
    }
    if (globals.Num_Elem >= INT_MAX) {
      global_counts_64bit = true;
    }

    /* Does any individual mesh have > INT_MAX nodes/elements? */
    bool local_counts_64bit = false;
    for (int iproc = 0; iproc < num_proc && !local_counts_64bit; iproc++) {
      if (globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
              globals.Num_External_Nodes[iproc] >=
          INT_MAX) {
        local_counts_64bit = true;
      }
      if (globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc] >= INT_MAX) {
        local_counts_64bit = true;
      }
    }

    int mode = 0;
    if (ids64bit) {
      mode |= EX_IDS_INT64_DB;
      fmt::print("-- ID output requires 64-bit integers.\n");
    }

    if (global_counts_64bit) {
      mode |= EX_MAPS_INT64_DB;
      fmt::print("-- Global Id map output requires 64-bit integers.\n");
    }

    if (local_counts_64bit) {
      mode |= EX_BULK_INT64_DB;
      fmt::print("-- Bulk data output requires 64-bit integers.\n");
    }

    return mode;
  }

  void find_message_info(size_t iunit_size, size_t max_units, size_t *num_units_per_message,
                         size_t *num_messages, size_t *num_units_left_over)
  {
    // Function which determines information regarding the size and number of
    //    messages to be used by the functions which broadcast the mesh information.

    if (iunit_size > 0) {
      *num_units_per_message = MAX_CHUNK_SIZE / (2 * iunit_size);
    }
    else {
      fmt::print(stderr, "ERROR:  find_message_info called with unit_size = 0.\n");
      exit(1);
    }

    if (*num_units_per_message > max_units) {
      *num_units_per_message = max_units;
    }

    *num_messages        = max_units / (*num_units_per_message);
    *num_units_left_over = max_units - *num_messages * (*num_units_per_message);

    if (max_units % (*num_units_per_message) != 0) {
      (*num_messages)++;
    }
  }
} // namespace

template void NemSpread<double, int>::load_mesh();
template void NemSpread<float, int>::load_mesh();
template void NemSpread<double, int64_t>::load_mesh();
template void NemSpread<float, int64_t>::load_mesh();

template <typename T, typename INT> void NemSpread<T, INT>::load_mesh()
{

  /* Function which reads the EXODUS database which contains the mesh.  It
   * does this under the control of another EXODUS database file which
   * contains the load-balance information.  The load balance information must
   * have been already read.
   *
   *----------------------------------------------------------------------------
   *
   * Functions called:
   *
   * check_exodus_error -- function which handles the error code returned by
   *                        calls to EXODUS API routines.
   * construct_lb_filename --  function which appends the string '-#' where
   *                            '#' is the number of processors to the name of
   *                            the EXODUS filename.
   * read_lb          -- function which reads the load-balance information
   *                     from an EXODUS database for a given processor.
   * read_coord       -- function which reads the nodal coordinates information
   *                     from an EXODUS database for a given processor.
   * read_elem_blk_ids-- Function which read the element block ids and other
   *                     global information specific to the element block level.
   * read_elem_blk    -- Function which reads the element level information
   *                     for each element block in the mesh
   * read_node_sets   -- function which reads the node sets information from an
   *                     EXODUS database for a given processor.
   * read_side_sets   -- function which reads the side sets information from an
   *                     EXODUS database for a given processor.
   *
   *----------------------------------------------------------------------------
   */

  /* Local variables */

  int    mesh_exoid = 0;
  int    cpu_ws;
  float  version;
  double start_time = 0.0;

  /* Allocate some memory for each processor read by this processor */
  globals.Proc_Num_Elem_Blk  = (int *)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int));
  globals.Proc_Num_Node_Sets = (int *)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int));
  globals.Proc_Num_Side_Sets = (int *)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(int));

  globals.Proc_NS_List_Length =
      (INT *)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(INT));
  globals.Proc_SS_Elem_List_Length =
      (INT *)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(INT));

  /* Initialize */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    globals.Proc_Num_Elem_Blk[iproc]        = 0;
    globals.Proc_Num_Node_Sets[iproc]       = 0;
    globals.Proc_Num_Side_Sets[iproc]       = 0;
    globals.Proc_NS_List_Length[iproc]      = 0;
    globals.Proc_SS_Elem_List_Length[iproc] = 0;
  }

  globals.GElem_Blks = (INT **)array_alloc(__FILE__, __LINE__, 1, 24 * Proc_Info[2], sizeof(INT *));
  globals.Proc_Nodes_Per_Elem   = globals.GElem_Blks + Proc_Info[2];
  globals.Proc_Elem_Blk_Ids     = globals.Proc_Nodes_Per_Elem + Proc_Info[2];
  globals.Proc_Elem_Blk_Types   = globals.Proc_Elem_Blk_Ids + Proc_Info[2];
  globals.Proc_Num_Attr         = globals.Proc_Elem_Blk_Types + Proc_Info[2];
  globals.Proc_Num_Elem_In_Blk  = globals.Proc_Num_Attr + Proc_Info[2];
  globals.Proc_Connect_Ptr      = globals.Proc_Num_Elem_In_Blk + Proc_Info[2];
  globals.Proc_Elem_Connect     = globals.Proc_Connect_Ptr + Proc_Info[2];
  globals.Proc_NS_Ids           = globals.Proc_Elem_Connect + Proc_Info[2];
  globals.Proc_NS_Count         = globals.Proc_NS_Ids + Proc_Info[2];
  globals.Proc_NS_DF_Count      = globals.Proc_NS_Count + Proc_Info[2];
  globals.Proc_NS_Pointers      = globals.Proc_NS_DF_Count + Proc_Info[2];
  globals.Proc_NS_List          = globals.Proc_NS_Pointers + Proc_Info[2];
  globals.GNode_Sets            = globals.Proc_NS_List + Proc_Info[2];
  globals.Proc_NS_GNMap_List    = globals.GNode_Sets + Proc_Info[2];
  globals.Proc_SS_Ids           = globals.Proc_NS_GNMap_List + Proc_Info[2];
  globals.Proc_SS_Elem_Count    = globals.Proc_SS_Ids + Proc_Info[2];
  globals.Proc_SS_DF_Count      = globals.Proc_SS_Elem_Count + Proc_Info[2];
  globals.Proc_SS_Elem_Pointers = globals.Proc_SS_DF_Count + Proc_Info[2];
  globals.Proc_SS_Elem_List     = globals.Proc_SS_Elem_Pointers + Proc_Info[2];
  globals.Proc_SS_Side_List     = globals.Proc_SS_Elem_List + Proc_Info[2];
  globals.Proc_SS_DF_Pointers   = globals.Proc_SS_Side_List + Proc_Info[2];
  globals.GSide_Sets            = globals.Proc_SS_DF_Pointers + Proc_Info[2];
  globals.Proc_SS_GEMap_List    = globals.GSide_Sets + Proc_Info[2];

  globals.Proc_Global_Elem_Id_Map.resize(Proc_Info[2]);
  globals.Proc_Global_Node_Id_Map.resize(Proc_Info[2]);

  /* Initialize */
  for (int iproc = 0; iproc < 24 * Proc_Info[2]; iproc++) {
    globals.GElem_Blks[iproc] = nullptr;
  }

  globals.Coor = (T ***)array_alloc(__FILE__, __LINE__, 1, Proc_Info[2], sizeof(T **));

  /* Initialize */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    globals.Coor[iproc] = nullptr;
  }

  globals.Proc_Elem_Attr.resize(Proc_Info[2]);
  globals.Proc_NS_Dist_Fact.resize(Proc_Info[2]);
  globals.Proc_SS_Dist_Fact.resize(Proc_Info[2]);

  /* Check for a problem which has too many processors for a given mesh */

  if (globals.Num_Node / Proc_Info[0] < 1) {
    fmt::print(stderr,
               "[{}]: ERROR: Problem divided among too many "
               "processors.\n",
               __func__);
    exit(1);
  }
  else if (globals.Num_Elem / Proc_Info[0] < 1) {
    fmt::print(stderr,
               "[{}]: ERROR: Problem divided among too many "
               "processors.\n",
               __func__);
    exit(1);
  }

  /* Open the EXODUS mesh file */

  /* computing precision should be the same as the database precision
   *
   * EXCEPTION: if the io_ws is smaller than the machine precision,
   * ie - database with io_ws == 4 on a Cray (sizeof(float) == 8),
   * then the cpu_ws must be the machine precision.
   */
  {
    cpu_ws     = io_ws;
    int mode   = EX_READ | int64api;
    mesh_exoid = ex_open(ExoFile.c_str(), mode, &cpu_ws, &io_ws, &version);
    if (mesh_exoid < 0) {
      fmt::print(stderr, "{}Exodus returned error opening mesh file, {}\n", __func__, ExoFile);
      exit(1);
    }
  }

  auto max_name_length = ex_inquire_int(mesh_exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  ex_set_max_name_length(mesh_exoid, max_name_length);

  globals.Num_QA_Recs = ex_inquire_int(mesh_exoid, EX_INQ_QA);

  /* Add a QA record for the spreader */
  globals.Num_QA_Recs++;

  int length_qa = 4 * (globals.Num_QA_Recs);
  globals.QA_Record =
      reinterpret_cast<char **>(array_alloc(__FILE__, __LINE__, 1, length_qa, sizeof(char *)));
  for (int i1 = 0, index = 0; i1 < globals.Num_QA_Recs; i1++) {
    for (int i2 = 0; i2 < 4; i2++) {
      globals.QA_Record[index++] = reinterpret_cast<char *>(
          array_alloc(__FILE__, __LINE__, 1, MAX_STR_LENGTH + 1, sizeof(char)));
    }
  }

  if (ex_get_qa(mesh_exoid, (char *(*)[4])globals.QA_Record) < 0) {
    fmt::print(stderr, "[{}]: ERROR, could not get QA record(s)\n", __func__);
    exit(1);
  }

  /* Read in the information records */
  globals.Num_Info_Recs = ex_inquire_int(mesh_exoid, EX_INQ_INFO);

  if (globals.Num_Info_Recs > 0) {
    /* Allocate the Information records */
    globals.Info_Record = (char **)array_alloc(__FILE__, __LINE__, 2, globals.Num_Info_Recs,
                                               MAX_LINE_LENGTH + 1, sizeof(char));
    if (!globals.Info_Record) {
      fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
      exit(1);
    }

    if (ex_get_info(mesh_exoid, globals.Info_Record) < 0) {
      fmt::print(stderr, "[{}]: ERROR, could not get Info record(s)\n", __func__);
      exit(1);
    }
  }

  /* Read in the assembly information */
  {
    ex_init_params ex_info;
    ex_get_init_ext(mesh_exoid, &ex_info);
    globals.Num_Assemblies = ex_info.num_assembly;
  }

  if (globals.Num_Assemblies > 0) {
    globals.Assemblies.resize(globals.Num_Assemblies);
    for (auto &assembly : globals.Assemblies) {
      assembly.name        = nullptr;
      assembly.entity_list = nullptr;
    }
    ex_get_assemblies(mesh_exoid, Data(globals.Assemblies));

    for (auto &assembly : globals.Assemblies) {
      assembly.entity_list = new int64_t[assembly.entity_count];
    }

    // Now get the assembly entity lists...
    ex_get_assemblies(mesh_exoid, Data(globals.Assemblies));
  }

  /* Read in the coordinate frame information */
  globals.Num_Coordinate_Frames = ex_inquire_int(mesh_exoid, EX_INQ_COORD_FRAMES);

  if (globals.Num_Coordinate_Frames > 0) {

    /* Allocate the Coordinate Frame records */
    globals.Coordinate_Frame_Ids =
        (INT *)array_alloc(__FILE__, __LINE__, 1, globals.Num_Coordinate_Frames, sizeof(INT));
    if (!globals.Coordinate_Frame_Ids) {
      fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
      exit(1);
    }

    globals.Coordinate_Frame_Coordinates =
        (T *)array_alloc(__FILE__, __LINE__, 1, 9 * globals.Num_Coordinate_Frames, sizeof(T));
    if (!globals.Coordinate_Frame_Coordinates) {
      fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
      exit(1);
    }

    globals.Coordinate_Frame_Tags =
        (char *)array_alloc(__FILE__, __LINE__, 1, globals.Num_Coordinate_Frames, sizeof(char));
    if (!globals.Coordinate_Frame_Tags) {
      fmt::print(stderr, "[{}]: ERROR, insufficient memory!\n", __func__);
      exit(1);
    }

    int num_frames = 0;
    if (ex_get_coordinate_frames(mesh_exoid, &num_frames, globals.Coordinate_Frame_Ids,
                                 globals.Coordinate_Frame_Coordinates,
                                 globals.Coordinate_Frame_Tags) < 0) {
      fmt::print(stderr, "[{}]: ERROR, could not get Coordinate Frame record(s)\n", __func__);
      exit(1);
    }
    if (num_frames != globals.Num_Coordinate_Frames) {
      fmt::print(stderr, "[{}]: ERROR, frame count inconsistency\n", __func__);
      exit(1);
    }
  }

  /*
   * Allocate Temporary Arrays that are only used in el_exoII_io.c (Note: Calls
   * to array_alloc are minimized by using pointer arithmetic. Resulting memory
   * structure may be assumed during later broadcast routines)
   */

  if (globals.Num_Elem_Blk <= 0) {
    fmt::print(stderr, "ERROR, globals.Num_Elem_Blk = {}\n", globals.Num_Elem_Blk);
    exit(1);
  }

  Num_Elem_In_Blk.resize(globals.Num_Elem_Blk);
  Num_Nodes_Per_Elem.resize(globals.Num_Elem_Blk);
  Num_Attr_Per_Elem.resize(globals.Num_Elem_Blk);
  Elem_Blk_Ids.resize(globals.Num_Elem_Blk);

  Elem_Blk_Types = (char **)array_alloc(__FILE__, __LINE__, 2, globals.Num_Elem_Blk,
                                        MAX_STR_LENGTH + 1, sizeof(char));
  Elem_Blk_Names = (char **)array_alloc(__FILE__, __LINE__, 2, globals.Num_Elem_Blk,
                                        max_name_length + 1, sizeof(char));
  Elem_Blk_Attr_Names =
      (char ***)array_alloc(__FILE__, __LINE__, 1, globals.Num_Elem_Blk, sizeof(char **));

  Node_Set_Ids.resize(globals.Num_Node_Set);
  std::vector<INT> num_nodes_in_node_set(globals.Num_Node_Set);
  std::vector<INT> num_df_in_nsets(globals.Num_Node_Set);

  if (globals.Num_Node_Set > 0) {
    Node_Set_Names = (char **)array_alloc(__FILE__, __LINE__, 2, globals.Num_Node_Set,
                                          max_name_length + 1, sizeof(char));
  }

  Side_Set_Ids.resize(globals.Num_Side_Set);
  std::vector<INT> num_elem_in_ssets(globals.Num_Side_Set);
  std::vector<INT> num_df_in_ssets(globals.Num_Side_Set);

  if (globals.Num_Side_Set > 0) {
    Side_Set_Names = (char **)array_alloc(__FILE__, __LINE__, 2, globals.Num_Side_Set,
                                          max_name_length + 1, sizeof(char));
  }

  /*
   * Process the Global Element Block IDs and associated information related to
   * element blocks (i.e., read them, broadcast them, and check them against
   * the input file).
   */
  start_time = second();

  read_elem_blk_ids(mesh_exoid, max_name_length);

  fmt::print("\tTime to read element block IDs:            {:.4f}\n", second() - start_time);

  /*
   * Process the Node Set IDs and associated information related to node sets
   * (i.e., read them, broadcast them, and check them against the input file).
   */

  if (globals.Num_Node_Set > 0) {
    start_time = second();
    read_node_set_ids(mesh_exoid, num_nodes_in_node_set, num_df_in_nsets, max_name_length);

    fmt::print("\tTime to read node set IDs:                 {:.4f}\n", second() - start_time);
  }

  /*
   * Process the Side Set IDs and associated information related to side sets
   * (i.e., read them, broadcast them, and check them against the input file).
   */
  if (globals.Num_Side_Set > 0) {
    start_time = second();
    read_side_set_ids(mesh_exoid, num_elem_in_ssets, num_df_in_ssets, max_name_length);
    fmt::print("\tTime to read side set IDs:                 {:.4f}\n", second() - start_time);
  }

  /*
   * Process the element block information.  Find out which element blocks have
   * elements needed by which processors.  Set-up and fill in vectors which
   * have length equal to the number of element blocks on the processor,
   * globals.Num_Elem_Blk.
   */
  start_time = second();

  extract_elem_blk();

  fmt::print("\tTime to extract element block information: {:.4f}\n", second() - start_time);

  /*
   * Read the mesh information from the exodus file and broadcast it to the
   * processors.  NOTE: this will only read as much information as will fit in
   * the communication buffer.  It will then send the information to all the
   * processors so they can extract the information they need. Then, the next
   * block of mesh data will be read and sent, etc., until the entire mesh has
   * been read and distributed to the processors.  Note, that the load balance
   * information has already been distributed by this point.  Therefore, each
   * node already nodes which global nodes it owns, which elements it has
   * global nodes in that it owns, and which ghost nodes it needs information
   * about.
   *
   *        The following information is defined as the mesh information:
   *
   *            1 - Coordinate values of global nodes
   *            2 - Element connectivity graph
   *            3 - Attributes of the elements
   *            4 - Composition of Node Sets
   *            5 - Composition of Side Sets
   */

  /*      Read the coordinates, broadcast, and sort */

  start_time = second();

  read_coord(mesh_exoid, max_name_length);

  fmt::print("\tTime to read nodal coordinates:            {:.4f}\n", second() - start_time);

  /*
   * Read the element connectivity and attributes information.  Broadcast it
   * and extract only the information that the current processor needs.
   */
  start_time = second();
  read_elem_blk(mesh_exoid);
  fmt::print("\tTime to read element blocks:               {:.4f}\n", second() - start_time);

  /* Read the node sets and sort */
  if (globals.Num_Node_Set > 0) {
    start_time = second();
    read_node_sets(mesh_exoid, num_nodes_in_node_set, num_df_in_nsets);
    fmt::print("\tTime to read node sets:                    {:.4f}\n", second() - start_time);
  }

  /* Assign the element types. */
  if (globals.Num_Side_Set > 0) {
    start_time = second();
    create_elem_types(); /* globals.Elem_Type only used in sideset read */
    fmt::print("\tTime to categorize element types:          {:.4f}\n", second() - start_time);

    /* Read the side sets and sort */
    start_time = second();
    read_side_sets(mesh_exoid, num_elem_in_ssets, num_df_in_ssets);
    fmt::print("\tTime to read side sets:                    {:.4f}\n", second() - start_time);
  }

  /* Close the EXODUS  mesh file */
  check_exodus_error(ex_close(mesh_exoid), "ex_close");

  /************* Output the information to the parallel disks *****************/
  /*==========================================================================*/
  /*==========================================================================*/

  /* Generate the processor to disk map */
  gen_disk_map(&PIO_Info, Proc_Info, 0, 1);

  /* Generate the parallel exodus file name */
  /* See if any '/' in the name.  IF present, isolate the basename of the file */
  std::string cTemp;

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

  /* Check sizes to see if need to store using 64-bit on the database */
  int db_mode = check_sizes(Proc_Info[2], Elem_Blk_Ids, Node_Set_Ids, Side_Set_Ids, globals);
  if (force64db) {
    fmt::print("-- Command-line option forcing output to use all 64-bit integers.\n");
    db_mode |= EX_ALL_INT64_DB;
  }

  for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
    std::string Parallel_File_Name = gen_par_filename(cTemp, Proc_Ids[iproc], Proc_Info[0]);

    /* Create the parallel Exodus file for writing */
    if (Debug_Flag >= 7) {
      fmt::print("{}Parallel mesh file name is {}\n", __func__, Parallel_File_Name);
    }
    else {
      if (iproc % 10 == 0 || iproc == Proc_Info[2] - 1) {
        fmt::print("{}", iproc);
      }
      else {
        fmt::print(".");
      }
    }

    int mode = EX_CLOBBER;
    mode |= int64api;
    mode |= db_mode;
    if ((mesh_exoid = ex_create(Parallel_File_Name.c_str(), mode, &cpu_ws, &io_ws)) == -1) {

      fmt::print(stderr, "[{}] Could not create parallel Exodus file:\n\t{}\n", __func__,
                 Parallel_File_Name);
      exit(1);
    }

    ex_set_max_name_length(mesh_exoid, max_name_length);

    /* Set fill mode off... */
    {
      int old_fill;
      nc_set_fill(mesh_exoid, NC_NOFILL, &old_fill);
    }

    if (Debug_Flag >= 7) {
      fmt::print("{}Parallel mesh file id is {}\n", __func__, mesh_exoid);
    }

    /* Write out a parallel mesh file local to each processor */
    write_parExo_data(mesh_exoid, max_name_length, iproc, num_nodes_in_node_set, num_elem_in_ssets,
                      Num_Elem_In_Blk);

    /* Close the parallel exodus file */
    if (ex_close(mesh_exoid) == -1) {
      fmt::print(stderr, "{}Could not close the parallel Exodus file\n", __func__);
      exit(1);
    }
  } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */

  if (Debug_Flag >= 4) {
    fmt::print("\n\n\t\tTIMING TABLE FOR PROCESSORS\n");
    fmt::print("===========================================================\n");
    fmt::print("[0]: {:6.2f}\n", PIO_Time_Array[0]);
    for (int i1 = 1; i1 < 25; i1++) {
      fmt::print("\t\t[{}]: {:6.2f}\n", i1, PIO_Time_Array[i1]);
    }
    fmt::print("\nOutput rate: {:6.2f} kB/s\n", PIO_Time_Array[25]);
  }
  else if (Debug_Flag >= 1) {
    fmt::print("\n\nOutput rate: {:6.2f} kB/s\n", PIO_Time_Array[25]);
  }

  /*--------------------------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/

  /* Free QA records */
  for (int i1 = 0; i1 < 4 * (globals.Num_QA_Recs); i1++) {
    safe_free((void **)&(globals.QA_Record[i1]));
  }

  safe_free((void **)&globals.QA_Record);

  for (int i1 = 0; i1 < globals.Num_Elem_Blk; i1++) {
    if (Num_Attr_Per_Elem[i1] > 0) {
      safe_free((void **)&(Elem_Blk_Attr_Names[i1]));
    }
  }
  safe_free((void **)&Elem_Blk_Attr_Names);

  /* Free Coordinate Frames */
  if (globals.Num_Coordinate_Frames > 0) {
    safe_free((void **)&globals.Coordinate_Frame_Ids);
    safe_free((void **)&globals.Coordinate_Frame_Coordinates);
    safe_free((void **)&globals.Coordinate_Frame_Tags);
  }

  // Clear out allocated memory for assembly entity_list...
  for (auto &assembly : globals.Assemblies) {
    free(assembly.name);
    delete[] assembly.entity_list;
    assembly.name        = nullptr;
    assembly.entity_list = nullptr;
  }

  /* done with the Coordinate names */
  for (int i1 = 0; i1 < globals.Num_Dim; i1++) {
    safe_free((void **)&(Coord_Name[i1]));
  }

  /* Free some local arrays */
  Num_Elem_In_Blk.clear();
  Node_Set_Ids.clear();
  Side_Set_Ids.clear();

  safe_free((void **)&Elem_Blk_Types);
  safe_free((void **)&Node_Set_Names);
  safe_free((void **)&Side_Set_Names);
  safe_free((void **)&Elem_Blk_Names);

  safe_free((void **)&Restart_Info.NV_Name);
  safe_free((void **)&Restart_Info.EV_Name);
  safe_free((void **)&Restart_Info.GV_Name);
  safe_free((void **)&Restart_Info.NSV_Name);
  safe_free((void **)&Restart_Info.SSV_Name);

  /*
   * free up some other memory so that there is more
   * memory available for reading restart variables
   */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    safe_free((void **)&(globals.Proc_Connect_Ptr[iproc]));
    safe_free((void **)&(globals.Proc_Elem_Connect[iproc]));

    safe_free((void **)&(globals.Coor[iproc]));
    globals.Proc_Elem_Attr[iproc].clear();
    globals.Proc_NS_Dist_Fact[iproc].clear();
    globals.Proc_SS_Dist_Fact[iproc].clear();
  }
  globals.N_Comm_Map.clear();
  globals.E_Comm_Map.clear();
  globals.Proc_Elem_Attr.clear();
  globals.Proc_NS_Dist_Fact.clear();
  globals.Proc_SS_Dist_Fact.clear();
  fmt::print("\n");

} /* END of routine load_mesh () *********************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_elem_blk_ids(int mesh_exoid, int max_name_length)
{

  /* This function reads part of the element block info from the EXODUS file.
   * It reads all information having a length equal to the number of elements
   * blocks, specifically.
   *
   * The function then broadcasts this information to all processors.
   *
   * The function also calls check_matrl_blks, which checks the element block ids
   * in the input file and in the mesh file for consistency (and prints out a
   * table for a large enough value of Debug_Flag).
   */

  /* Get the Element Block IDs from the input file */
  check_exodus_error(ex_get_ids(mesh_exoid, EX_ELEM_BLOCK, Data(Elem_Blk_Ids)), "ex_get_ids");

  check_exodus_error(ex_get_names(mesh_exoid, EX_ELEM_BLOCK, Elem_Blk_Names), "ex_get_names");

  /*
   *     Get from the input file:
   *         Number of Elements               in each element block
   *         Number of nodes per element      in each element block
   *         Number of attributes per element in each element block
   *         The element type for elements    in each element block
   */
  for (int i = 0; i < globals.Num_Elem_Blk; i++) {

    check_exodus_error(ex_get_block(mesh_exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[i], Elem_Blk_Types[i],
                                    &Num_Elem_In_Blk[i], &Num_Nodes_Per_Elem[i], nullptr, nullptr,
                                    &Num_Attr_Per_Elem[i]),
                       "ex_get_elem_block");

    /* Convert element block types to lower case here */
    string_to_lower(Elem_Blk_Types[i], '\0');

    /* Allocate space for attribute names (if any), read and store. */
    if (Num_Attr_Per_Elem[i] > 0) {
      Elem_Blk_Attr_Names[i] = (char **)array_alloc(__FILE__, __LINE__, 2, Num_Attr_Per_Elem[i],
                                                    max_name_length + 1, sizeof(char));
      check_exodus_error(
          ex_get_attr_names(mesh_exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[i], Elem_Blk_Attr_Names[i]),
          "ex_get_attr_names");
    }
    else {
      Elem_Blk_Attr_Names[i] = nullptr;
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_node_set_ids(int mesh_exoid, std::vector<INT> &num_nodes_in_node_set,
                                          std::vector<INT> &num_df_in_nsets,
                                          int /*max_name_length*/)

/*
 * This function reads part of the node set info from the EXODUS file.  It
 * reads all information having a length equal to the number of node sets,
 * specifically.
 *
 * The function then broadcasts this information to all processors.
 *
 */

{
  if (globals.Num_Node_Set > 0) {
    int error = ex_get_ids(mesh_exoid, EX_NODE_SET, Data(Node_Set_Ids));
    check_exodus_error(error, "ex_get_node_set_ids");

    error = ex_get_names(mesh_exoid, EX_NODE_SET, Node_Set_Names);
    check_exodus_error(error, "ex_get_node_set_ids");

    for (int i = 0; i < globals.Num_Node_Set; i++) {
      check_exodus_error(ex_get_set_param(mesh_exoid, EX_NODE_SET, Node_Set_Ids[i],
                                          &num_nodes_in_node_set[i], &num_df_in_nsets[i]),
                         "ex_get_set_param");
    }
  }

  /* Output debug info */
  if ((Debug_Flag > 1)) {
    fmt::print("\n\n");
    print_line("=", 79);
    fmt::print("\tTABLE OF NODE SET ID\'s\n\n");
    fmt::print("Node_Set_Num   ID  globals.Num_Nodes\n");
    print_line("-", 79);

    if (globals.Num_Node_Set > 0) {
      for (int i = 0; i < globals.Num_Node_Set; i++) {
        fmt::print("{:6d}{:11d}{:12d}\n", i, Node_Set_Ids[i], num_nodes_in_node_set[i]);
      }
    }

    else {
      fmt::print("\tNO NODE SETS ARE DEFINED IN THE MESH FILE\n");
    }
    print_line("=", 79);
    fmt::print("\n");
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_side_set_ids(int mesh_exoid, std::vector<INT> &num_elem_in_ssets,
                                          std::vector<INT> &num_df_in_ssets,
                                          int /*max_name_length*/)

/* This function reads part of the side set info from the EXODUS file.  It
 * reads all information having a length equal to the number of side sets,
 * specifically.
 *
 * The function then broadcasts this information to all processors.
 *
 */

{
  if (globals.Num_Side_Set > 0) {
    int error = ex_get_ids(mesh_exoid, EX_SIDE_SET, Data(Side_Set_Ids));
    check_exodus_error(error, "ex_get_side_set_ids");

    error = ex_get_names(mesh_exoid, EX_SIDE_SET, Side_Set_Names);
    check_exodus_error(error, "ex_get_side_set_ids");

    for (int i = 0; i < globals.Num_Side_Set; i++) {
      check_exodus_error(ex_get_set_param(mesh_exoid, EX_SIDE_SET, Side_Set_Ids[i],
                                          &num_elem_in_ssets[i], &num_df_in_ssets[i]),
                         "ex_get_set_param");
    }
  }

  /* Output debug information */
  if (Debug_Flag > 1) {
    fmt::print("\n\n");
    print_line("=", 79);
    fmt::print("\tTABLE OF SIDE SET ID\'s\n\n");
    fmt::print("Side_Set_Num   ID   Number Elements\n");
    print_line("-", 79);

    if (globals.Num_Side_Set > 0) {
      for (int i = 0; i < globals.Num_Side_Set; i++) {
        fmt::print("{:6d}{:11d}  {:12}\n", i, Side_Set_Ids[i],
                   fmt::group_digits(num_elem_in_ssets[i]));
      }
    }

    else {
      fmt::print("\tNO SIDE SETS ARE DEFINED IN THE MESH FILE\n");
    }
    print_line("=", 79);
    fmt::print("\n");
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_coord(int exoid, int max_name_length)
{

  /* Function which reads the nodal coordinates information from an
   * EXODUS database for a given processor.
   */

  /*
   * Calculate the size of the coordinate space for each processors coordinate matrix.
   */
  for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
    size_t itotal_nodes = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                          globals.Num_External_Nodes[iproc];

    /* Allocate permanent storage for the coordinates */
    if (itotal_nodes > 0) {
      globals.Coor[iproc] =
          (T **)array_alloc(__FILE__, __LINE__, 2, globals.Num_Dim, itotal_nodes, sizeof(T));
    }
    else {
      globals.Coor[iproc] = nullptr;
    }
  }

  /* Allocate temporary space to hold 1 dimensions worth of coordinates... */
  std::vector<T> coord(globals.Num_Node);

  /* Read in the coordinates and broadcast to the processors */
  for (int idim = 0; idim < globals.Num_Dim; idim++) {
    switch (idim) {
    case 0:
      check_exodus_error(ex_get_coord(exoid, Data(coord), nullptr, nullptr), "ex_get_coord");
      break;
    case 1:
      check_exodus_error(ex_get_coord(exoid, nullptr, Data(coord), nullptr), "ex_get_coord");
      break;
    case 2:
      check_exodus_error(ex_get_coord(exoid, nullptr, nullptr, Data(coord)), "ex_get_coord");
      break;
    }

    for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
      size_t itotal_nodes = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                            globals.Num_External_Nodes[iproc];

      for (size_t j = 0; j < itotal_nodes; j++) {
        size_t inode                 = globals.GNodes[iproc][j];
        globals.Coor[iproc][idim][j] = coord[inode];
      }
    }
  }
  coord.clear();

  for (int i = 0; i < globals.Num_Dim; i++) {
    Coord_Name[i] = reinterpret_cast<char *>(
        array_alloc(__FILE__, __LINE__, 1, max_name_length + 1, sizeof(char)));
  }

  /* Get the coordinate names */
  if (ex_get_coord_names(exoid, Coord_Name) < 0) {
    fmt::print(stderr, "ERROR:Unable to obtain coordinate names\n");
    exit(1);
  }

  /* Handle global node ids... */
  std::vector<INT> global_node_ids(globals.Num_Node);

  check_exodus_error(ex_get_id_map(exoid, EX_NODE_MAP, Data(global_node_ids)), "ex_get_id_map");
  /*
   * Check whether map is sequential (1..globals.Num_Node). If it is, then it
   * provides no information and we don't need to store it in the
   * output databases.
   */
  {
    bool sequential = true;
    for (size_t i = 0; i < globals.Num_Node; i++) {
      if ((size_t)global_node_ids[i] != i + 1) {
        sequential = false;
        break;
      }
    }

    // Check that map is valid 1 <= global_node_id[*]
    // If not, output a warning and disable the map.
    for (size_t i = 0; i < globals.Num_Node; i++) {
      if (global_node_ids[i] <= 0) {
        fmt::print(stderr,
                   "---------------------------------------------------------------------\n"
                   "ERROR: Local node {} has a global id of {} which is invalid.\n"
                   "       All global ids must be greater than 0. The map will be ignored.\n"
                   "---------------------------------------------------------------------\n",
                   fmt::group_digits(i + 1), fmt::group_digits(global_node_ids[i]));
        sequential = true; // Map is invalid, ignore it.
        break;
      }
    }

    if (!sequential) {
      for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {

        size_t itotal_nodes = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                              globals.Num_External_Nodes[iproc];

        globals.Proc_Global_Node_Id_Map[iproc].resize(itotal_nodes);
        extract_global_node_ids(global_node_ids, globals.Num_Node, iproc);
      }
    }
  }
} /* END of routine read_coord */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT> void NemSpread<T, INT>::extract_elem_blk()

/* Function which calculates the element block information for the current
 * processor, given the global element block information.  If an error is
 * tripped, the program will set a flag, and abort at the end of the routine.
 *
 * The following global variables are allocated and calculated in this routine:
 *
 *           globals.Num_Elem_Blk = Number of element blocks on the current
 *                               processor
 *           globals.GElem_Blks [globals.Num_Elem_Blk]
 *                             = Map from the local element block number
 *                               to the global element block number.
 *           globals.Proc_Nodes_Per_Elem [globals.Num_Elem_Blk]
 *                             = Number of nodes per element for each
 *                               block on the current processor
 *           globals.Proc_Elem_Blk_Ids [globals.Num_Elem_Blk]
 *                             = Element block id's for the processor's
 *                               element blocks
 *           globals.Proc_Elem_Blk_Types [globals.Num_Elem_Blk]
 *                             = Element block types for the processor's
 *                               element blocks
 *                               (this is a unique integer number)
 *           globals.Proc_Num_Attr [globals.Num_Elem_Blk]
 *                             = Number of attributes for each block on the
 *                               current processor
 *           globals.Proc_Num_Elem_In_Blk [globals.Num_Elem_Blk]
 *                             = Number of elements in the processor's
 *                               element blocks
 *
 *    ------------------------------------------------------------------------
 *
 *       Functions called:
 *
 *       find_elem_block -- function which finds the element block which
 *                          owns each element in this processor.  In addition,
 *                          it determines a map from the current element to an
 *                          element block.
 *
 *    ------------------------------------------------------------------------
 */

{
  std::vector<INT> proc_elem_blk{};

  /* Element blocks local to the current processor */

  /**************************** execution begins ******************************/

  /*
   * Allocate temporary array for a mapping between the local element number
   * and the corresponding element block id for the element block that the
   * element belongs to.
   */

  for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {

    proc_elem_blk.resize(globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc]);

    /* Find out which element block each element in this processor belongs to.
     * Fill this information into the temporary vector, proc_elem_blk.  Also,
     * calculate:
     *
     *              globals.Num_Elem_Blk = Number of element blocks defined on
     *                                  the current processor.
     *              globals.GElem_Blks        = Map from the local element block number
     *                                  to the global element block number.
     */
    find_elem_block(proc_elem_blk, iproc, Proc_Ids[iproc]);

    /* Allocate Permanent integer arrays, which have lengths equal to the
     * number of element blocks defined on the processor.  This is done with a
     * single malloc call.  Initialize this space to zero.
     */
    if (globals.Num_Elem_Blk > 0) {
      globals.Proc_Nodes_Per_Elem[iproc] = (INT *)array_alloc(
          __FILE__, __LINE__, 1, (4 * globals.Num_Elem_Blk + globals.Proc_Num_Elem_Blk[iproc]),
          sizeof(INT));
      globals.Proc_Elem_Blk_Ids[iproc] = globals.Proc_Nodes_Per_Elem[iproc] + globals.Num_Elem_Blk;
      globals.Proc_Elem_Blk_Types[iproc] = globals.Proc_Elem_Blk_Ids[iproc] + globals.Num_Elem_Blk;
      globals.Proc_Num_Attr[iproc] =
          globals.Proc_Elem_Blk_Types[iproc] + globals.Proc_Num_Elem_Blk[iproc];
      globals.Proc_Num_Elem_In_Blk[iproc] = globals.Proc_Num_Attr[iproc] + globals.Num_Elem_Blk;
      /* Initialize */
      for (int i = 0; i < (4 * globals.Num_Elem_Blk + globals.Proc_Num_Elem_Blk[iproc]); i++) {
        globals.Proc_Nodes_Per_Elem[iproc][i] = 0;
      }
    }
    else {
      fmt::print(stderr, "ERROR globals.Num_Elem_Blk = {}\n", globals.Num_Elem_Blk);
      exit(1);
    }

    /*
     * Fill in the local processor element block arrays from a lookup from the
     * global element block arrays.  Note that the global element block arrays
     * will be freed at the end of the file.
     */
    for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {

      size_t iglobal_blk                    = globals.GElem_Blks[iproc][i];
      globals.Proc_Nodes_Per_Elem[iproc][i] = Num_Nodes_Per_Elem[iglobal_blk];
      globals.Proc_Elem_Blk_Ids[iproc][i]   = Elem_Blk_Ids[iglobal_blk];
      globals.Proc_Num_Attr[iproc][i]       = Num_Attr_Per_Elem[iglobal_blk];

      /* Determine the element type integer id for this current element block.
       * The element type integer ID is a unique identifier (unique for all
       * element types that the program knows about).  It is determined from
       * the character string name, Elem_Blk_Types, and the number of nodes in
       * the element.
       */
      globals.Proc_Elem_Blk_Types[iproc][i] = get_type(
          Elem_Blk_Types[iglobal_blk], globals.Proc_Nodes_Per_Elem[iproc][i], globals.Num_Dim);
    }

    /*
     * Determine the number of elements in each of the element blocks defined
     * on the local processor.
     */

    for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {
      for (INT j = 0; j < globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];
           j++) {
        if (proc_elem_blk[j] == globals.Proc_Elem_Blk_Ids[iproc][i]) {
          (globals.Proc_Num_Elem_In_Blk[iproc][i])++;
        }
      }
    }

    /* Sort globals.GElems so that each element block is monotonic */
    size_t j = 0;
    for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {
      gds_qsort((Data(globals.GElems[iproc])) + j, globals.Proc_Num_Elem_In_Blk[iproc][i]);
      j += globals.Proc_Num_Elem_In_Blk[iproc][i];
    }
  } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */

  if (Debug_Flag >= 5) {
    for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {

      /* Printout the Element Block Information defined on each Processor */
      print_line("=", 79);
      fmt::print("\t\tLocal Element Block information for Proc = {}\n", Proc_Ids[iproc]);
      fmt::print("\t\tNumber of Elem blocks on processor = {}\n", globals.Proc_Num_Elem_Blk[iproc]);
      fmt::print("{}{}\n", "Local_Block_Num  Global_Block_Num  Block_ID Nodes_Per_Elem ",
                 "Num_Attributes  Elem_Blk_Type  globals.Proc_Num_Elem_In_Blk "
                 "Glb_Elm_In_Blk");
      print_line("-", 79);
      for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {
        fmt::print("{:4d}\t\t{:5}\t{:8}\t{:8}\t{:8}\t{:8}\t{:8}\t{:8}\n", i,
                   globals.GElem_Blks[iproc][i],
                   fmt::group_digits(globals.Proc_Elem_Blk_Ids[iproc][i]),
                   fmt::group_digits(globals.Proc_Nodes_Per_Elem[iproc][i]),
                   fmt::group_digits(globals.Proc_Num_Attr[iproc][i]),
                   globals.Proc_Elem_Blk_Types[iproc][i],
                   fmt::group_digits(globals.Proc_Num_Elem_In_Blk[iproc][i]),
                   fmt::group_digits(Num_Elem_In_Blk[globals.GElem_Blks[iproc][i]]));
      }
      print_line("=", 79);
    }
  }

} /* END extract_elem_blk() */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT> void NemSpread<T, INT>::read_elem_blk(int exoid)

/* Function which reads the element block information from an EXODUS
 * database.  If an error is tripped, the program will set a flag, and abort at
 * the end of the routine.
 *
 *    ------------------------------------------------------------------------
 *
 *       extract_elem_connect
 *                        - Function which sorts through the element block
 *                          connectivity information and extracts the
 *                          information required by each processor.
 *
 *    ------------------------------------------------------------------------
 */

{

  /* Local variables */

  T *elem_attr = nullptr;
#ifdef DEBUG
  size_t ielem_count;
#endif
  INT   *elem_blk             = nullptr, ipos;
  size_t num_elem_per_message = 0;
  size_t num_attr_per_message = 0;
  size_t num_attr_left_over   = 0;
  size_t num_elem_messages    = 0;
  size_t num_attr_messages    = 0;
  size_t num_elem_left_over   = 0;

  size_t istart_elem;
  size_t iend_elem;
  size_t istart_attr;
  size_t iend_attr;
  int    local_ielem_blk;

  /**************************** execution begins ******************************/

  /*
   * Allocate memory for the vector of global element types in the entire
   * mesh. Used in sideset manipulations
   */
  GM_Elem_Types = (int *)array_alloc(__FILE__, __LINE__, 1, globals.Num_Elem, sizeof(int));

  std::vector<INT> global_ids(globals.Num_Elem);

  for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {

    /*
     * Build globals.Proc_Connect_Ptr, a vector of pointers to the start of each
     * element's connectivity list in globals.Proc_Elem_Connect.  The last entry in
     * globals.Proc_Connect_Ptr is the length of globals.Proc_Elem_Connect.
     */
    globals.Proc_Connect_Ptr[iproc] = (INT *)array_alloc(
        __FILE__, __LINE__, 1,
        globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc] + 1, sizeof(INT));

    size_t iptr_count      = 0;
    size_t iconnect_length = 0;
    for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {
      for (INT j = 0; j < globals.Proc_Num_Elem_In_Blk[iproc][i]; j++) {
        globals.Proc_Connect_Ptr[iproc][iptr_count++] = iconnect_length;
        iconnect_length += globals.Proc_Nodes_Per_Elem[iproc][i];
      }
    }
    globals.Proc_Connect_Ptr[iproc][iptr_count] = iconnect_length;

    /*
     * Allocate the processor's connectivity list vector
     * - This is a global vector
     */
    globals.Proc_Elem_Connect[iproc] =
        (INT *)array_alloc(__FILE__, __LINE__, 1, iconnect_length, sizeof(INT));

#ifdef DEBUG
    for (size_t i = 0; i < iconnect_length; i++)
      globals.Proc_Elem_Connect[iproc][i] = -1111111;
#endif

    /*
     * Allocate the processor's attribute list vector, if its length
     * is nonzero
     */
    size_t iattr_length = 0;
    for (int i = 0; i < globals.Proc_Num_Elem_Blk[iproc]; i++) {
      iattr_length += globals.Proc_Num_Attr[iproc][i] * globals.Proc_Num_Elem_In_Blk[iproc][i];
    }
    if (iattr_length > 0) {
      globals.Proc_Elem_Attr[iproc].resize(iattr_length);
    }
  } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc)" */

  size_t icount = 0;

  /*
   * READ THE ELEMENT CONNECTIVITY AND ATTRIBUTE INFORMATION
   * FROM THE EXODUS FILE FOR EACH ELEMENT BLOCK
   */
  for (int ielem_blk = 0; ielem_blk < globals.Num_Elem_Blk; ielem_blk++) {
    if (Num_Elem_In_Blk[ielem_blk] > 0) {

      /*
       * Calculate the size of a single element's connectivity list and a
       * single element's attribute list
       */
      size_t iconnect_size = Num_Nodes_Per_Elem[ielem_blk] * sizeof(INT);
      size_t iattr_size    = Num_Attr_Per_Elem[ielem_blk] * sizeof(T);

      /*
       * Determine how many element connectivity list and how many element
       * attributes can be send in single message
       */
      find_message_info(iconnect_size, Num_Elem_In_Blk[ielem_blk], &num_elem_per_message,
                        &num_elem_messages, &num_elem_left_over);

      if (iattr_size > 0) {
        find_message_info(iattr_size, Num_Elem_In_Blk[ielem_blk], &num_attr_per_message,
                          &num_attr_messages, &num_attr_left_over);
      }
      else {
        num_attr_messages = 0;
      }

      if (Debug_Flag > 1) {
        fmt::print("\n\nMessage summary for Element Block number {}, ", ielem_blk);
        fmt::print("having a block id of {}:\n", Elem_Blk_Ids[ielem_blk]);
        fmt::print("\tNumber of messages needed for the element connectivity "
                   "vector = {}\n",
                   num_elem_messages);
        fmt::print("\tNumber of elements per message = {}\n", num_elem_per_message);
        fmt::print("\tNumber of nodes per element = {}\n", Num_Nodes_Per_Elem[ielem_blk]);
        fmt::print("\tLength of each message = {} bytes\n",
                   (Num_Nodes_Per_Elem[ielem_blk] * num_elem_per_message * sizeof(INT)));
        if (num_attr_messages > 0) {
          fmt::print("\tNumber of attribute messages: {}\n\tNumber "
                     "of attributes per message: {}\n\n",
                     num_attr_messages, num_attr_per_message);
        }
      }

      /*
       * Allocate the arrays for reading from the ExodusII file and
       * broadcasting
       */
      elem_blk = (INT *)array_alloc(
          __FILE__, __LINE__, 1, Num_Nodes_Per_Elem[ielem_blk] * num_elem_per_message, sizeof(INT));
      if (num_attr_messages > 0) {
        elem_attr = (T *)array_alloc(
            __FILE__, __LINE__, 1, Num_Attr_Per_Elem[ielem_blk] * num_attr_per_message, sizeof(T));
      }
      else {
        elem_attr = nullptr;
      }

      /*
       * Read in the element connectivity list for the current element block,
       * ielem_blk, from Proc 0, and then, broadcast it to all of the processors
       */
      for (size_t i = 0; i < num_elem_messages; i++) {

        if (Debug_Flag >= 2) {
          fmt::print("\telem block message: {} of {}\n", i + 1, num_elem_messages);
        }

        /* Initialize the element connectivity list to a value of -1.0 */
        for (size_t j      = 0; j < Num_Nodes_Per_Elem[ielem_blk] * num_elem_per_message;
             elem_blk[j++] = -1) {
          ;
        }

        size_t num_to_get = 0;
        istart_elem       = i * num_elem_per_message;
        if (num_elem_left_over == 0 || i < num_elem_messages - 1) {
          num_to_get = num_elem_per_message;
        }
        else {
          num_to_get = num_elem_left_over;
        }
        iend_elem = istart_elem + num_to_get;

        int el_type;
        check_exodus_error(ex_get_partial_conn(exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[ielem_blk],
                                               (istart_elem + 1), num_to_get, elem_blk, nullptr,
                                               nullptr),
                           "ex_get_partial_conn");
        if (Debug_Flag >= 2) {
          fmt::print("\t\tread connectivity\n");
        }

        el_type =
            get_type(Elem_Blk_Types[ielem_blk], Num_Nodes_Per_Elem[ielem_blk], globals.Num_Dim);
        for (size_t ielem = 0; ielem < num_to_get; ielem++) {
          GM_Elem_Types[icount++] = el_type;
        }

        if (Debug_Flag >= 2) {
          fmt::print("\t\tgot element types\n");
        }

        /* PRINT OUT THE ELEMENT CONNECTIVITY TABLE IF IN DEBUGGING MODE */
        if (Debug_Flag >= 6) {
          fmt::print("\n\n\n");
          print_line("=", 79);
          fmt::print("Printout of Element connectivity list obtained from "
                     "Exodus file:\n");
          fmt::print("\tGlobal element block number = {}\n", ielem_blk);
          fmt::print("\tElement ID number     = {}\n", Elem_Blk_Ids[ielem_blk]);
          fmt::print("\tMessage number        = {}\n", i);
          print_line("-", 79);
          ipos = 0;
          for (size_t j = 0; j < num_to_get; j++) {
            fmt::print("\t elem: {}, nodes:", j);
            for (int k = 0; k < Num_Nodes_Per_Elem[ielem_blk]; k++) {
              fmt::print(" {}", elem_blk[ipos++]);
            }
            fmt::print("\n");
          }
          print_line("=", 79);
        }

        /*
         * On each processor, extract the element connectivity lists that the
         * processor needs
         */
        for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
          extract_elem_connect(elem_blk, ielem_blk, istart_elem, iend_elem, &local_ielem_blk,
                               iproc);
        }
        if (Debug_Flag >= 2) {
          fmt::print("\t\textract connectivity\n");
        }

      } /* End "for(i=0; i < num_elem_messages; i++)" */

      /* Read in the element attribute lists and broadcast to the processors */
      for (size_t i = 0; i < num_attr_messages; i++) {

        if (Debug_Flag >= 2) {
          fmt::print("\tattribute message: {} of {}\n", i + 1, num_attr_messages);
        }

        /* Initialize */
        for (size_t j = 0; j < Num_Attr_Per_Elem[ielem_blk] * num_attr_per_message; j++) {
          elem_attr[j] = 0.0;
        }

        istart_attr = i * num_attr_per_message;

        if (num_attr_left_over == 0 || i < (num_attr_messages - 1)) {
          iend_attr = istart_attr + num_attr_per_message;
        }
        else {
          iend_attr = istart_attr + num_attr_left_over;
        }

        if (num_attr_left_over == 0 || i < (num_attr_messages - 1)) {
          check_exodus_error(ex_get_partial_attr(exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[ielem_blk],
                                                 (istart_attr + 1), num_attr_per_message,
                                                 elem_attr),
                             "ex_get_partial_attr");
        }
        else {
          check_exodus_error(ex_get_partial_attr(exoid, EX_ELEM_BLOCK, Elem_Blk_Ids[ielem_blk],
                                                 (istart_attr + 1), num_attr_left_over, elem_attr),
                             "ex_get_partial_attr");
        }

        for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {

          if (Debug_Flag > 6) {
            fmt::print("\t\tExtract attributes for processor {}\n", Proc_Ids[iproc]);
          }
          /*
           * On each processor, extract the element attributes that the
           * processor needs
           */
          extract_elem_attr(elem_attr, ielem_blk, istart_attr, iend_attr,
                            Num_Attr_Per_Elem[ielem_blk], iproc);
        }
      }

      /*       Free vectors */
      safe_free((void **)&elem_blk);
      safe_free((void **)&elem_attr);

    } /* END "if (globals.Proc_Num_Elem_In_Blk[ielem_blk] > 0)" */

  } /* End "for(ielem_blk=0; ielem_blk < globals.Num_Elem_Blk; ielem_blk++)" */

  /* Handle global element ids... */
  check_exodus_error(ex_get_id_map(exoid, EX_ELEM_MAP, Data(global_ids)), "ex_get_id_map");

  /*
   * Check whether map is sequential (1..globals.Num_Elem). If it is, then it
   * provides no information and we don't need to store it in the
   * output databases.
   */
  {
    bool sequential = true;
    for (size_t i = 0; i < globals.Num_Elem; i++) {
      if ((size_t)global_ids[i] != i + 1) {
        sequential = false;
        break;
      }
    }

    // Check that map is valid 1 <= global_node_id[*]
    // If not, output a warning and disable the map.
    for (size_t i = 0; i < globals.Num_Elem; i++) {
      if (global_ids[i] <= 0) {
        fmt::print(stderr,
                   "---------------------------------------------------------------------\n"
                   "ERROR: Local element {} has a global id of {} which is invalid.\n"
                   "       All global ids must be greater than 0. The map will be ignored.\n"
                   "---------------------------------------------------------------------\n",
                   fmt::group_digits(i + 1), fmt::group_digits(global_ids[i]));
        sequential = true; // Map is invalid, ignore it.
        break;
      }
    }

    if (!sequential) {
      for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
        globals.Proc_Global_Elem_Id_Map[iproc].resize(globals.Num_Internal_Elems[iproc] +
                                                      globals.Num_Border_Elems[iproc]);
        extract_global_element_ids(global_ids, globals.Num_Elem, iproc);
      }
    }
    else {
      /* Should be nullptr already, but make it more clear */
      for (int iproc = Proc_Info[4]; iproc < Proc_Info[4] + Proc_Info[5]; iproc++) {
        globals.Proc_Global_Elem_Id_Map[iproc].clear();
      }
    }
  }
} /* read_elem_blk */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::extract_global_element_ids(const std::vector<INT> &global_ids,
                                                   size_t /*Num_Elem*/, int iproc)
{
  /*
   * globals.Num_Elem -- number of elements in serial mesh.
   * global_ids[] (size globals.Num_Elem) = global id corresponding to index
   *                                 of element in serial file
   * globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc] = number of element on this
   processor
   * globals.GElems[iproc][i] = index of element in serial file corresponding
   to i'th element on this processor. Both are zero-based?
  */
  size_t num_local_element = globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];
  for (size_t i = 0; i < num_local_element; i++) {
    globals.Proc_Global_Elem_Id_Map[iproc][i] = global_ids[globals.GElems[iproc][i]];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::extract_global_node_ids(const std::vector<INT> &global_ids,
                                                size_t /*Num_Node*/, int iproc)
{
  /*
   * Num_Elem -- number of elements in serial mesh.
   * global_ids[] (size Num_Elem) = global id corresponding to index
   *                                 of element in serial file
   * Num_Internal_Elems[iproc] + Num_Border_Elems[iproc] = number of element on this processor
   * GElems[iproc][i] = index of element in serial file corresponding
   * to i'th element on this processor. Both are zero-based?
   */
  size_t num_local_node = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                          globals.Num_External_Nodes[iproc];
  for (size_t i = 0; i < num_local_node; i++) {
    globals.Proc_Global_Node_Id_Map[iproc][i] = global_ids[globals.GNodes[iproc][i]];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
size_t NemSpread<T, INT>::extract_elem_connect(INT elem_blk[], int icurrent_elem_blk,
                                               size_t istart_elem, size_t iend_elem,
                                               int *local_ielem_blk, int iproc)

/*
 * Function which extracts the element connectivity list from the current
 * element connectivity message.  The results are put into the global
 * vector, globals.Proc_Elem_Connect.
 *
 *   The function returns the number of elements in the current message
 *   that are defined on the current processor.
 */

{
  size_t count_hits = 0;
  bool   found      = false;

  /* Match the Element Block Id of the global block with the local block number
   * The end result is the current value of ielem_blk, Store this value in the
   * output variable, local_ielem_blk.
   */
  for (int ielem_blk = 0; ielem_blk < globals.Proc_Num_Elem_Blk[iproc]; ielem_blk++) {
    if (globals.Proc_Elem_Blk_Ids[iproc][ielem_blk] == Elem_Blk_Ids[icurrent_elem_blk]) {
      *local_ielem_blk = ielem_blk;
      found            = true;

      /* Some elements in the current element block passed in from
         read_elem_blk may be on this processor */

      /* Calculate the number of elements in the current message */
      size_t num_elem = iend_elem - istart_elem;

      /* Calculate iglobal_offset - The sum of the elements in all the global
       *                            element blocks preceding this one
       */
      size_t iglobal_offset = 0;
      for (int i = 0; i < globals.GElem_Blks[iproc][ielem_blk]; i++) {
        iglobal_offset += Num_Elem_In_Blk[i];
      }

      /* Calculate iproc_offset - The sum of the elements in all the
       *                          processor's
       *                          element blocks preceding this one.
       * Calculate iproc_start  - The starting position in the connectivity
       *                          vector for this element block,
       *                          globals.Proc_Elem_Connect[]
       */
      size_t iproc_offset = 0;
      size_t iproc_start  = 0;
      for (int i = 0; i < ielem_blk; i++) {
        iproc_offset += globals.Proc_Num_Elem_In_Blk[iproc][i];
        iproc_start +=
            globals.Proc_Num_Elem_In_Blk[iproc][i] * globals.Proc_Nodes_Per_Elem[iproc][i];
      }

      INT iglobal_begin = iglobal_offset + istart_elem;
      INT iglobal_end   = iglobal_begin + num_elem;

      size_t ibegin = iproc_offset;
      size_t iend   = iproc_offset + globals.Proc_Num_Elem_In_Blk[iproc][ielem_blk];

      /* OPTIMIZE: Check that max globals.GElems is >min < max...*/
      /* Note that globals.GElems is sorted, so we can check that the
       * ranges overlap and skip the loop if they don't */
      if (globals.GElems[iproc][ibegin] < iglobal_end &&
          globals.GElems[iproc][iend - 1] >= iglobal_begin) {

        for (size_t i = ibegin; i < iend; i++) {
          INT iglobal_elem = globals.GElems[iproc][i];

          /* globals.GElems is sorted, so if we are outside the iglobal range,
           * break out of this loop...
           */
          if (iglobal_elem > iglobal_end) {
            break;
          }

          if (iglobal_begin <= iglobal_elem && iglobal_elem < iglobal_end) {
            count_hits++;

            size_t ipos =
                iproc_start + (i - iproc_offset) * globals.Proc_Nodes_Per_Elem[iproc][ielem_blk];

            size_t iglobal_pos =
                (iglobal_elem - iglobal_begin) * globals.Proc_Nodes_Per_Elem[iproc][ielem_blk];

            assert(iglobal_pos < num_elem * globals.Proc_Nodes_Per_Elem[iproc][ielem_blk]);
            /*
             * Store the connectivity information for the current element into
             * globals.Proc_Elem_Connect [].  At the same time, Decrement by one the
             * value of all nodes in the connectivity matrix.
             */
            for (int j = 0; j < globals.Proc_Nodes_Per_Elem[iproc][ielem_blk]; j++) {
              globals.Proc_Elem_Connect[iproc][ipos++] = elem_blk[iglobal_pos++] - 1;
            }
          }
        }
      }
    } /* END if (globals.Proc_Elem_Blk_Ids[ielem_blk]== ... */
  }

  if (!found) {
    *local_ielem_blk = -1;
  }

  return count_hits;

} /* extract_elem_connect */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::extract_elem_attr(T *elem_attr, int icurrent_elem_blk, size_t istart_elem,
                                          size_t iend_elem, int natt_p_elem, int iproc)
{
  /*
   * Function which extracts the element connectivity list from the current
   * element connectivity message.  The results are put into the global
   * vector, globals.Proc_Elem_Connect.  This function is run on all processors,
   * include Proc = 0
   *
   *       Author:          Scott Hutchinson (1421)
   *
   *      ------------------------------------------------------------------
   *
   *      input
   *    -----------
   *       elem_attr  = Vector containing the element attributes for this message.
   *                    This message consists of the either the whole or
   *                    part of the element attribute array for a single
   *                    element block.
   *
   *      istart_elem = This is the starting value of the element number for
   *                    the element attribute message.  i.e., the first
   *                    message in an element block will have
   *                    istart_elem = 0.  The second message for that
   *                    element block will have the value:
   *                         istart_elem = num_elem_per_message.
   *      iend_elem   = This is the value + 1 of the last element number for
   *                    the element attribute message.  i.e., the first
   *                    message in an element block will have
   *                          iend_elem = num_elem_per_message
   */

  /* Match the Element Block Id of the global block with the local block number */
  for (int ielem_blk = 0; ielem_blk < globals.Proc_Num_Elem_Blk[iproc]; ielem_blk++) {
    if (globals.Proc_Elem_Blk_Ids[iproc][ielem_blk] == Elem_Blk_Ids[icurrent_elem_blk]) {

      /* Some elements in the current element block passed in from
         read_elem_blk may be on this processor */

      /* Calculate the number of elements in the current message */
      size_t num_elem = iend_elem - istart_elem;

      /* Calculate iglobal_offset - The sum of the elements in all the global
       *                            element blocks preceding this one
       */
      size_t iglobal_offset = 0;
      for (int i = 0; i < globals.GElem_Blks[iproc][ielem_blk]; i++) {
        iglobal_offset += Num_Elem_In_Blk[i];
      }

      /* Calculate iproc_offset - The sum of the elements in all the
       *                          processor's
       *                          element blocks preceding this one.
       * Calculate iproc_start  - The starting position in the connectivity
       *                          vector for this element block,
       *                          globals.Proc_Elem_Connect[]
       */
      size_t iproc_offset = 0;
      size_t iproc_start  = 0;
      for (int i = 0; i < ielem_blk; i++) {
        iproc_offset += globals.Proc_Num_Elem_In_Blk[iproc][i];
        iproc_start += globals.Proc_Num_Elem_In_Blk[iproc][i] * globals.Proc_Num_Attr[iproc][i];
      }

      INT iglobal_begin = iglobal_offset + istart_elem;
      INT iglobal_end   = iglobal_begin + num_elem;

      size_t ibegin = iproc_offset;
      size_t iend   = iproc_offset + globals.Proc_Num_Elem_In_Blk[iproc][ielem_blk];

      /* OPTIMIZE: Check that max globals.GElems is >min < max...*/
      /* Note that globals.GElems is sorted, so we can check that the
       * ranges overlap and skip the loop if they don't */
      if (globals.GElems[iproc][ibegin] < iglobal_end &&
          globals.GElems[iproc][iend - 1] >= iglobal_begin) {

        for (size_t i = ibegin; i < iend; i++) {
          INT iglobal_elem = globals.GElems[iproc][i];

          /* globals.GElems is sorted, so if we are outside the iglobal range,
           * break out of this loop...
           */
          if (iglobal_elem > iglobal_end) {
            break;
          }

          if (iglobal_begin <= iglobal_elem && iglobal_elem < iglobal_end) {
            size_t ipos        = iproc_start + (i - iproc_offset) * natt_p_elem;
            size_t iglobal_pos = (iglobal_elem - iglobal_begin) * natt_p_elem;

            assert(iglobal_pos < num_elem * natt_p_elem);
            /*
             * Store the connectivity information for the current element into
             * globals.Proc_Elem_Connect [].  At the same time, Decrement by one the
             * value of all nodes in the connectivity matrix.
             */
            for (int j = 0; j < natt_p_elem; j++) {
              globals.Proc_Elem_Attr[iproc][ipos++] = elem_attr[iglobal_pos++];
            }
          }
        }
      }
    } /* END if (globals.Proc_Elem_Blk_Ids[ielem_blk]== ... */
  }
} /* extract_elem_connect */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::find_elem_block(std::vector<INT> &proc_elem_blk, int iproc,
                                        int /*proc_for*/)

{
  /* Function which finds the element block which owns each element on the
   * current processor.  In addition, a map from the local element block number
   * to the global element block number, *globals.GElem_Blks, is created.
   * globals.Num_Elem_Blk, the number of element blocks on the current processor,
   * is calculated This function is called by every processor.
   *
   * Author(s):          Scott Hutchinson (9221)
   *                     Gary Hennigan (9221)
   *
   *       Output:
   *
   *         globals.Num_Elem_Blk  = Number of element blocks on the processor
   *                              (Global int variable)
   *
   *        *proc_elem_blk      = Vector of element block ids for the local
   *                              elements defined on the processor.
   *                              (Local int vector of length
   *                              globals.Num_Internal_Elems + globals.Num_Border_Elems)
   *
   *        *globals.GElem_Blks         = Map from the local element block number to
   *                              the global element block number.
   *                              (Global int vector of length globals.Num_Elem_Blk)
   */
  if (globals.Num_Elem_Blk == 0) {
    return;
  }

  /* Boolean vector of length globals.Num_Elem_Blk If the i'th
     element block exists on the current processor, the i'th entry
     is set to TRUE  */
  std::vector<bool> elem_in_blk(globals.Num_Elem_Blk);

  /* Vector of integer offsets into the vector globals.GElems.
     It has a length of globals.Num_Elem_Blk+1.  The i'th entry
     points to the beginning of of the element map information
     for the first element in the i'th element block.  The
     (i+1)th entry points to the last element for the i'th block
     in the vector globals.GElem. */
  std::vector<INT> elem_blk_point(globals.Num_Elem_Blk + 1);

  /*
   * Construct a vector of pointers to the beginning of element blocks in an
   * element connectivity list (vector)
   */
  elem_blk_point[0] = 0;
  for (int i = 0; i < globals.Num_Elem_Blk; i++) {
    elem_blk_point[i + 1] = elem_blk_point[i] + Num_Elem_In_Blk[i];
  }

  /*
   * Find out which elements belonging to the current processor are in which
   * local element block.  Store this in *proc_elem_blk.
   */

  /* Internal Elements */
  if (check_monot(&globals.GElems[iproc][0], globals.Num_Internal_Elems[iproc])) {
    size_t tmp_cnt = globals.Num_Internal_Elems[iproc];
    int    j       = 0;
    size_t i       = 0;
    while (i < tmp_cnt && j < globals.Num_Elem_Blk) {
      while (i < tmp_cnt && globals.GElems[iproc][i] < elem_blk_point[j + 1]) {
        assert(globals.GElems[iproc][i] >= elem_blk_point[j]);
        proc_elem_blk[i++] = j;
        elem_in_blk[j]     = true;
      }
      j++;
    }
  }
  else {
    for (INT i = 0; i < globals.Num_Internal_Elems[iproc]; i++) {
      bool found = false;
      for (int j = 0; j < globals.Num_Elem_Blk && !found; j++) {
        if (globals.GElems[iproc][i] < elem_blk_point[j + 1] &&
            globals.GElems[iproc][i] >= elem_blk_point[j]) {
          proc_elem_blk[i] = j;
          elem_in_blk[j]   = true;
          found            = true;
        }
      }
      if (!found) {
        fmt::print(stderr, "find_elem_block: Error!:\n");
        fmt::print(stderr,
                   "\tElement {} not found in any element "
                   "block.\n",
                   i);
        exit(1);
      }
    }
  }

  /* Border Elements */
  if (check_monot(&globals.GElems[iproc][globals.Num_Internal_Elems[iproc]],
                  globals.Num_Border_Elems[iproc])) {
    INT tmp_cnt = globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];
    int j       = 0;
    INT i       = globals.Num_Internal_Elems[iproc];
    while (i < tmp_cnt && j < globals.Num_Elem_Blk) {
      while (i < tmp_cnt && globals.GElems[iproc][i] < elem_blk_point[j + 1]) {
        assert(globals.GElems[iproc][i] >= elem_blk_point[j]);
        proc_elem_blk[i++] = j;
        elem_in_blk[j]     = true;
      }
      j++;
    }
  }
  else {
    INT tmp_cnt = globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];
    for (INT i = globals.Num_Internal_Elems[iproc]; i < tmp_cnt; i++) {
      bool found = false;
      for (int j = 0; j < globals.Num_Elem_Blk && !found; j++) {
        if (globals.GElems[iproc][i] < elem_blk_point[j + 1] &&
            globals.GElems[iproc][i] >= elem_blk_point[j]) {
          proc_elem_blk[i] = j;
          elem_in_blk[j]   = true;
          found            = true;
        }
      }
      if (!found) {
        fmt::print(stderr, "find_elem_block: Error!:\n");
        fmt::print(stderr,
                   "\tElement {} not found in any element "
                   "block.\n",
                   i);
        exit(1);
      }
    }
  }

  /*
   * Reorder globals.GElems based on the element block information. This is necessary
   * to insure that internal and border elements that are in the same element
   * block are sequentially contained in globals.GElems, which is assumed in some
   * later operations.
   */

  my_sort(globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc], Data(proc_elem_blk),
          Data(globals.GElems[iproc]));

  /* Now change proc_elem_blk to be a list of global element block IDs */
  for (INT i = 0; i < globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc]; i++) {
    size_t j         = proc_elem_blk[i];
    proc_elem_blk[i] = Elem_Blk_Ids[j];
  }

  /* Count the number of element blocks defined on this processor */
  globals.Proc_Num_Elem_Blk[iproc] = 0;
  for (int i = 0; i < globals.Num_Elem_Blk; i++) {
    if (elem_in_blk[i]) {
      globals.Proc_Num_Elem_Blk[iproc]++;
    }
  }

  /*
   * Create a map from the current processor's element block number to the
   * 'global' element block number, globals.GElem_Blks.  This is a permanent Global
   * Map.
   */
  globals.GElem_Blks[iproc] =
      (INT *)array_alloc(__FILE__, __LINE__, 1, globals.Num_Elem_Blk, sizeof(INT));

  for (int i = 0, icount = 0; i < globals.Num_Elem_Blk; i++) {
    if (elem_in_blk[i]) {
      globals.GElem_Blks[iproc][icount++] = i;
    }
  }
} /* END of routine find_elem_block() ****************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_node_sets(int exoid, std::vector<INT> &num_nodes_in_node_set,
                                       std::vector<INT> &num_df_in_nsets)

/* Function which reads the node sets information from an EXODUS database.
 * It read information in chunks.  Then, it broadcasts the chunk to all
 * processors.  Each processor, then searches for the chunk for what it is
 * responsible for.  This function is called by all processors.
 *
 *     Global Variables which are set by this routine
 *     -----------------------------------------------
 *
 *       globals.Proc_Num_Node_Sets = Number of node sets on the current processor
 *
 *       globals.GNode_Sets [globals.Proc_Num_Node_Sets]
 *                          = Mapping between the local node set number and the
 *                            global node set number.
 *
 *       globals.S_Ids [globals.Proc_Num_Node_Sets]
 *                          = Set IDs for the node sets defined on the current
 *                            processor
 *
 *       globals.S_Count [globals.Proc_Num_Node_Sets]
 *                          = Number of nodes in the node sets defined on
 *                            the current processor
 *
 *       globals.S_Pointers [globals.Proc_Num_Node_Sets]
 *                          = Pointers into the node set list for the node sets
 *                            defined on the current processor.
 *
 *       globals.NS_List_Length = Length of the local node set list.
 *
 *       globals.NS_List [globals.NS_List_Length]
 *                          = Concatenated node set list for the current
 *                            processor.
 */

{
  size_t num_messages;
  size_t num_left_over;
  size_t num_node_per_message;

  /* Allocate arrays */
  std::vector<INT>  list_length(Proc_Info[2]);
  std::vector<INT>  proc_num_ns(Proc_Info[2]);
  std::vector<INT>  ns_cntr(Proc_Info[2]);
  std::vector<bool> ns_on_proc(Proc_Info[2]);

  /* pointers into the concatenated node set list which locate the start of the node sets     */
  std::vector<std::vector<INT>> proc_list_pointer(Proc_Info[2]);

  /* mapping from the processor node set numbers to the global node set numbers              */
  std::vector<std::vector<INT>> ns_proc2glob(Proc_Info[2]);

  /* the number of nodes in each processor's node set                                         */
  std::vector<std::vector<INT>> proc_list_count(Proc_Info[2]);

  /* the number of df's in each processor's node set (usually either zero or proc_list_count)*/
  std::vector<std::vector<INT>> proc_list_df_count(Proc_Info[2]);

  /* the node number lists for the current node set for the current processor               */
  std::vector<std::vector<std::vector<INT>>> proc_list(Proc_Info[2]);

  /*  Vector of pointers to nodeset global nodeeset node map lists.
      There is one pointer for each nodeset     */
  std::vector<std::vector<std::vector<INT>>> proc_gmap_list(Proc_Info[2]);

  /* the dist. factors for the current processor */
  std::vector<std::vector<std::vector<T>>> proc_list_df(Proc_Info[2]);

  /* Initialize */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    ns_proc2glob[iproc].resize(globals.Num_Node_Set);
    proc_list_pointer[iproc].resize(globals.Num_Node_Set);
    proc_list_count[iproc].resize(globals.Num_Node_Set);
    proc_list_df_count[iproc].resize(globals.Num_Node_Set);

    proc_list[iproc].resize(globals.Num_Node_Set);
    proc_gmap_list[iproc].resize(globals.Num_Node_Set);
    proc_list_df[iproc].resize(globals.Num_Node_Set);
  }

  std::vector<INT> node_proc_index(globals.Num_Node + 1);
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    size_t size = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc];
    for (size_t i = 0; i < size; i++) {
      size_t ind = globals.GNodes[iproc][i];
      node_proc_index[ind]++;
    }
  }

  /* Now convert the node_proc_index counts to an index... */
  size_t index = 0;
  for (size_t i = 0; i < globals.Num_Node; i++) {
    size_t count       = node_proc_index[i];
    node_proc_index[i] = index;
    index += count;
  }
  node_proc_index[globals.Num_Node] = index;

  /* Create a map from 'global node id' to processor */
  std::vector<INT> gmap(index, -1);

  /* Fill it in....
   * Note that since nodes can be on multiple processors, we can't
   * normally do a simple one-to-one map of node to its processor.
   *
   * The values in gmap[] from node_proc_index[node] to node_proc_index[node+1]
   * contains the listing of the processors that node 'node' are on.
   */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    size_t size = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc];
    for (size_t i = 0; i < size; i++) {
      size_t node = globals.GNodes[iproc][i];
      size_t beg  = node_proc_index[node];
      size_t end  = node_proc_index[node + 1];
      for (size_t j = beg; j < end; j++) {
        if (gmap[j] == -1) {
          gmap[j] = iproc;
          break;
        }
      }
    }
  }

  /*-------------------------------------------------------------------------*/
  /*                    LOOP OVER THE NODE SETS                              */
  /*-------------------------------------------------------------------------*/

  for (int i = 0; i < globals.Num_Node_Set; i++) {
    if (num_nodes_in_node_set[i] > 0) {

      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        proc_num_ns[iproc] = 0;
        ns_cntr[iproc]     = 0;
        ns_on_proc[iproc]  = false;
      }

      size_t iss_size = sizeof(INT);
      if (num_df_in_nsets[i] > 0) {
        iss_size += sizeof(T);
      }

      find_message_info(iss_size, num_nodes_in_node_set[i], &num_node_per_message, &num_messages,
                        &num_left_over);

      if (Debug_Flag > 1) {
        fmt::print("\nMessage summary for Node Set number {}, with an ID of {}:\n", i,
                   Node_Set_Ids[i]);
        fmt::print("\tNumber of messages need for node set = {}\n", num_messages);
        fmt::print("\tNumber of node IDs and dist. factors per message = {}\n",
                   num_node_per_message);
        fmt::print("\tLength of each message = {}\n", num_node_per_message * iss_size);
      }

      /* pointers into 'node_set' for the nodes that are common to
         both 'node_set' and the internal and border nodes for the
         current processor  */

      std::vector<INT> node_set(num_node_per_message);
      std::vector<T>   node_set_df;
      if (num_df_in_nsets[i] > 0) {
        node_set_df.resize(num_node_per_message);
      }

      /*---------------------------------------------------------------------*/
      /*                  LOOP OVER THE MESSAGES                             */
      /*---------------------------------------------------------------------*/

      for (size_t imess = 0; imess < num_messages; imess++) {

        std::vector<INT> proc_ns_node_count(Proc_Info[2]); // Count nodes on this processor.

        size_t istart_ns = imess * num_node_per_message;
        if (num_left_over != 0 && imess == num_messages - 1) {
          num_node_per_message = num_left_over;
        }

        /* Read in the part of the node set that will fit in the message */
        check_exodus_error(ex_get_partial_set(exoid, EX_NODE_SET, Node_Set_Ids[i], (istart_ns + 1),
                                              num_node_per_message, Data(node_set), nullptr),
                           "ex_get_partial_set");

        if (num_df_in_nsets[i] > 0) {
          check_exodus_error(ex_get_partial_set_dist_fact(exoid, EX_NODE_SET, Node_Set_Ids[i],
                                                          (istart_ns + 1), num_node_per_message,
                                                          Data(node_set_df)),
                             "ex_get_partial_node_set_df");
        }

        /* Renumber nodes to start at node '0' instead of node '1' */
        for (size_t j = 0; j < num_node_per_message; j++) {
          node_set[j]--;
        }

        for (size_t j = 0; j < num_node_per_message; j++) {
          size_t node = node_set[j];
          size_t beg  = node_proc_index[node];
          size_t end  = node_proc_index[node + 1];
          for (size_t k = beg; k < end; k++) {
            proc_ns_node_count[gmap[k]]++;
          }
        }

        std::vector<std::vector<INT>> proc_ns_pointer(Proc_Info[2]);
        for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
          proc_ns_pointer[iproc].reserve(proc_ns_node_count[iproc]);
        }

        for (size_t j = 0; j < num_node_per_message; j++) {
          size_t node = node_set[j];
          size_t beg  = node_proc_index[node];
          size_t end  = node_proc_index[node + 1];
          for (size_t k = beg; k < end; k++) {
            proc_ns_pointer[gmap[k]].push_back(j);
          }
        }

        /* Loop over the number of processors being handled */
        for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
          /*
           * If the message node set and either the Internal node or the Border
           * node sets intersect, then add this intersection to the list
           */
          if (proc_ns_node_count[iproc] > 0) {
            ns_on_proc[iproc] = true;
            proc_num_ns[iproc] += proc_ns_node_count[iproc];

            /* Allocate and store node information in a temporary vector */
            int nset = globals.Proc_Num_Node_Sets[iproc];
            proc_list[iproc][nset].resize(proc_num_ns[iproc]);
            proc_gmap_list[iproc][nset].resize(proc_num_ns[iproc]);
            if (num_df_in_nsets[i] > 0) {
              proc_list_df[iproc][nset].resize(proc_num_ns[iproc]);
            }

            for (INT j = 0; j < proc_ns_node_count[iproc]; j++) {

              proc_list[iproc][nset][ns_cntr[iproc]]      = node_set[proc_ns_pointer[iproc][j]];
              proc_gmap_list[iproc][nset][ns_cntr[iproc]] = proc_ns_pointer[iproc][j] + istart_ns;

              if (num_df_in_nsets[i] > 0) {
                proc_list_df[iproc][nset][ns_cntr[iproc]] = node_set_df[proc_ns_pointer[iproc][j]];
              }

              ns_cntr[iproc]++;
            }
          } /* End "if(intersection)" */
        } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */
      } /* End "for(imess=0; imess < num_messages; imess++)" */

      /*
       * If any part of this node-set is on the processor, update the various
       * pointers, lengths, etc.
       */
      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        if (ns_on_proc[iproc]) {
          ns_proc2glob[iproc][globals.Proc_Num_Node_Sets[iproc]] = i;

          proc_list_pointer[iproc][globals.Proc_Num_Node_Sets[iproc]] = list_length[iproc];

          proc_list_count[iproc][globals.Proc_Num_Node_Sets[iproc]] = proc_num_ns[iproc];

          if (num_df_in_nsets[i] > 0) {
            proc_list_df_count[iproc][globals.Proc_Num_Node_Sets[iproc]] = proc_num_ns[iproc];
          }
          else {
            proc_list_df_count[iproc][globals.Proc_Num_Node_Sets[iproc]] = 0;
          }

          (globals.Proc_Num_Node_Sets[iproc])++;

          list_length[iproc] += proc_num_ns[iproc];
        }
      }
    } /* END "if (num_nodes_in_node_set[i] > 0)" */
  } /* END "for (i = 0; i < globals.Num_Node_Set; i++)" */

  /*--------------------------------------------------------------------------*/
  /*             WRITE PERMANENT ARRAYS FOR NODE SET INFO                     */
  /*--------------------------------------------------------------------------*/

  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

    /* Allocate Permanent Arrays for node sets in one long malloc */
    if (globals.Num_Node_Set > 0) {

      /*
       * Note that the alloc is done based on globals.Num_Node_Set, rather than
       * globals.Proc_Num_Node_Sets[] due to the fact that nullptr entities are
       * stored on processors not having a particular node set.
       */
      globals.Proc_NS_Ids[iproc] = (INT *)array_alloc(
          __FILE__, __LINE__, 1,
          (3 * globals.Num_Node_Set + 2 * globals.Proc_Num_Node_Sets[iproc] + list_length[iproc]),
          sizeof(INT));

      globals.Proc_NS_Count[iproc]    = globals.Proc_NS_Ids[iproc] + globals.Num_Node_Set;
      globals.Proc_NS_DF_Count[iproc] = globals.Proc_NS_Count[iproc] + globals.Num_Node_Set;
      globals.Proc_NS_Pointers[iproc] = globals.Proc_NS_DF_Count[iproc] + globals.Num_Node_Set;
      globals.GNode_Sets[iproc] =
          globals.Proc_NS_Pointers[iproc] + globals.Proc_Num_Node_Sets[iproc];
      globals.Proc_NS_List[iproc] = globals.GNode_Sets[iproc] + globals.Proc_Num_Node_Sets[iproc];

      globals.Proc_NS_GNMap_List[iproc] =
          (INT *)array_alloc(__FILE__, __LINE__, 1, list_length[iproc], sizeof(INT));

      if (list_length[iproc] > 0) {
        globals.Proc_NS_Dist_Fact[iproc].resize(list_length[iproc]);
      }
    }
    else {
      globals.Proc_NS_Ids[iproc]      = nullptr;
      globals.Proc_NS_Count[iproc]    = nullptr;
      globals.Proc_NS_Pointers[iproc] = nullptr;
      globals.Proc_NS_List[iproc]     = nullptr;
    }

    /*
     * Fill in the permanent node set arrays which have length,
     * globals.Proc_Num_Node_Sets, the total number of node sets defined on the
     * current processor.
     */
    globals.Proc_NS_List_Length[iproc] = 0;
    for (int i = 0; i < globals.Proc_Num_Node_Sets[iproc]; i++) {
      globals.GNode_Sets[iproc][i]       = ns_proc2glob[iproc][i];
      globals.Proc_NS_Ids[iproc][i]      = Node_Set_Ids[ns_proc2glob[iproc][i]];
      globals.Proc_NS_Count[iproc][i]    = proc_list_count[iproc][i];
      globals.Proc_NS_Pointers[iproc][i] = proc_list_pointer[iproc][i];
      globals.Proc_NS_List_Length[iproc] += globals.Proc_NS_Count[iproc][i];
      globals.Proc_NS_DF_Count[iproc][i] = proc_list_df_count[iproc][i];
    }

    /* Construct the concatenated node list */
    for (int i = 0; i < globals.Proc_Num_Node_Sets[iproc]; i++) {
      for (INT j = 0; j < globals.Proc_NS_Count[iproc][i]; j++) {
        globals.Proc_NS_List[iproc][globals.Proc_NS_Pointers[iproc][i] + j] =
            proc_list[iproc][i][j];
        globals.Proc_NS_GNMap_List[iproc][globals.Proc_NS_Pointers[iproc][i] + j] =
            proc_gmap_list[iproc][i][j];

        if (globals.Proc_NS_DF_Count[iproc][i] > 0) {
          globals.Proc_NS_Dist_Fact[iproc][globals.Proc_NS_Pointers[iproc][i] + j] =
              proc_list_df[iproc][i][j];
        }
      }
    }
  } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */

} /* END of routine read_node_sets () ****************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_side_sets(int exoid, std::vector<INT> &num_elem_in_ssets,
                                       std::vector<INT> &num_df_in_ssets)
{
  /*
   * Function which reads the side sets information from an EXODUS database
   * for a given processor. It then broadcasts all information to every
   * processor. Each processor extracts and keeps only the information that it
   * needs to run its local problem.
   *
   *     Global Variables which are set by this routine
   *     -----------------------------------------------
   *
   *       globals.Proc_Num_Side_Sets = Number of side sets on the current processor
   *
   *       globals.GSide_Sets [globals.Proc_Num_Side_Sets]
   *                          = Mapping between the local side set number and the
   *                            global side set number.
   *
   *       globals.S_Ids [globals.Proc_Num_Side_Sets]
   *                          = Set IDs for the side sets defined on the current
   *                            processor
   *
   *       globals.Proc_SS_Elem_Count [globals.Proc_Num_Side_Sets]
   *                          = Number of elements in the side set for the
   *                            current processor
   *
   *       globals.S_Side_Count [globals.Proc_Num_Side_Sets]
   *                          = Number of nodes in the node set list for the side
   *                            set for the current processor.
   *
   *       globals.Proc_SS_Elem_Pointers [globals.Proc_Num_Side_Sets]
   *                          = Side sets pointer record for elements for
   *                            the side sets in a given processor
   *
   *       globals.Proc_SS_Elem_List_Length
   *                          = Length of the element list length for side sets
   *                            on the current processor
   *
   *       globals.Proc_SS_Elem_List [globals.Proc_SS_Elem_List_Length]
   *                          = Concatenated vector of elements that comprise
   *                            the side set definitions on the current
   *                            processor.
   *
   *       globals.S_Side_List [globals.Proc_SS_Elem_List_Length]
   *                          = Concatenated vector of sides that comprise
   *                            the side set definitions on the current
   *                            processor.
   */

  /*  Vector of pointers to side-set element lists. There is one pointer for each side-set     */
  std::vector<std::vector<std::vector<INT>>> proc_elem_list(Proc_Info[2]);

  /*  Vector of pointers to side-set side lists. There is one pointer for each side-set     */
  std::vector<std::vector<std::vector<INT>>> proc_side_list(Proc_Info[2]);

  /*  Vector of pointers to side-set global sideset element map lists. There is one pointer for each
   * side-set     */
  std::vector<std::vector<std::vector<INT>>> proc_gmap_list(Proc_Info[2]);

  /*  Vector of pointers to side-set df pointers There is one pointer for each side-set     */
  std::vector<std::vector<std::vector<INT>>> proc_df_ptr(Proc_Info[2]);
  std::vector<std::vector<std::vector<INT>>> proc_df_indx(Proc_Info[2]);
  std::vector<std::vector<std::vector<INT>>> proc_df_indx_cnt(Proc_Info[2]);

  std::vector<std::vector<std::vector<T>>> proc_ss_df(Proc_Info[2]);

  std::vector<INT>  ss_elem_cntr(Proc_Info[2]);
  std::vector<bool> ss_on_proc(Proc_Info[2]); /* Flag to indicate that the ss is on the proc */
  std::vector<INT>  elem_list_length(Proc_Info[2]); /* length of the element side-set list for the
                                                       current processor                           */
  std::vector<INT> proc_num_sides_ss(Proc_Info[2]); /* Number of sides in the current side set that
                                                       exist for the current processor */
  std::vector<INT> ntotal(Proc_Info[2]);
  std::vector<INT> ss_num(Proc_Info[2]);

  /* Allocate temporary arrays */
  /* Array of counts of the dist. factors.        */
  std::vector<std::vector<INT>> proc_ss_df_cnt(Proc_Info[2]);

  /* mapping from the processor side-set numbers to the global side-set numbers              */
  std::vector<std::vector<INT>> ss_proc2glob(Proc_Info[2]);

  /* pointers into the concatenated element list which locate the start of the element lists */
  std::vector<std::vector<INT>> proc_elem_list_ptr(Proc_Info[2]);

  /* the number of elements in each processor's element set */
  std::vector<std::vector<INT>> proc_elem_list_cnt(Proc_Info[2]);

  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    proc_elem_list[iproc].resize(globals.Num_Side_Set);
    proc_side_list[iproc].resize(globals.Num_Side_Set);
    proc_df_ptr[iproc].resize(globals.Num_Side_Set);
    proc_df_indx[iproc].resize(globals.Num_Side_Set);
    proc_df_indx_cnt[iproc].resize(globals.Num_Side_Set);
    proc_gmap_list[iproc].resize(globals.Num_Side_Set);

    proc_ss_df[iproc].resize(globals.Num_Side_Set);

    ss_proc2glob[iproc].resize(globals.Num_Side_Set);
    proc_elem_list_ptr[iproc].resize(globals.Num_Side_Set);
    proc_elem_list_cnt[iproc].resize(globals.Num_Side_Set);
    proc_ss_df_cnt[iproc].resize(globals.Num_Side_Set);
  }

  /* Create a map from 'global element id' to processor */
  std::vector<INT> gmap(globals.Num_Elem, -1);

  /* Fill it in.... */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    size_t size = globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];
    for (size_t i = 0; i < size; i++) {
      size_t ind = globals.GElems[iproc][i];
      gmap[ind]  = iproc;
    }
  }

  /*-------------------------------------------------------------------------*/
  /*                    LOOP OVER THE SIDE SETS                              */
  /*-------------------------------------------------------------------------*/
  for (int i = 0; i < globals.Num_Side_Set; i++) {

    if (num_elem_in_ssets[i] > 0) {

      size_t ilast      = 0;
      int    ilast_side = 0;

      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        ss_on_proc[iproc]        = false;
        proc_num_sides_ss[iproc] = 0;
        ss_elem_cntr[iproc]      = 0;
      }

      /* One element ID + One side ID */
      size_t iss_size = 3 * sizeof(INT);

      size_t num_messages;
      size_t num_left_over;
      size_t num_elem_per_message;
      find_message_info(iss_size, num_elem_in_ssets[i], &num_elem_per_message, &num_messages,
                        &num_left_over);

      if (Debug_Flag >= 2) {
        fmt::print("Message summary for Side Set number {}, with an ID of {}:\n", i,
                   Side_Set_Ids[i]);
        fmt::print("\tNumber of messages needed for element and side list = {}\n", num_messages);
        fmt::print("\tNumber of element and side IDs per message = {}\n", num_elem_per_message);
        fmt::print("\tLength of each message = {}\n", iss_size * num_elem_per_message);
      }

      /* Allocate temporary storage for the current message for the
       * current side set on each proc
       */
      std::vector<INT> ss_elem_list(num_elem_per_message); /* side-set element list */
      std::vector<INT> ss_side_list(num_elem_per_message); /* side-set node list */
      std::vector<INT> ss_df_ptr(num_elem_per_message);    /* pointer into the dist. factor list */

      /*---------------------------------------------------------------------*/
      /*                  LOOP OVER THE MESSAGES                             */
      /*---------------------------------------------------------------------*/
      for (size_t imess = 0; imess < num_messages; imess++) {

        if (Debug_Flag >= 2) {
          fmt::print("\tside set message: {} of {}\n", imess + 1, num_messages);
        }

        size_t istart_ss = imess * num_elem_per_message;

        if (num_left_over != 0 && imess == num_messages - 1) {
          num_elem_per_message = num_left_over;
        }

        /* Read in the part of the side set that will fit in the message. */

        check_exodus_error(ex_get_partial_set(exoid, EX_SIDE_SET, Side_Set_Ids[i], (istart_ss + 1),
                                              num_elem_per_message, Data(ss_elem_list),
                                              Data(ss_side_list)),
                           "ex_get_partial_set");

        /* Fill in the distribution factor pointer vector */
        if (imess == 0) {
          ilast      = 0;
          ilast_side = 0;
        }

        ss_df_ptr[0] = ilast + ilast_side;

        for (size_t j = 1; j < num_elem_per_message; j++) {

          /*
           * Kluge for "special" sidesest which has a single distribution
           * factor for each element instead of 1 per face node..
           */
          if (num_elem_in_ssets[i] == num_df_in_ssets[i]) {
            ilast_side = 1;
          }
          else {

            ilast_side =
                elem_info(NN_SIDE, GM_Elem_Types[(ss_elem_list[j - 1]) - 1], ss_side_list[j - 1]);

            /*
             * kludge for HEXSHELL's
             * where distribution factors are concerned, only count
             * 4 nodes on the 6 node faces (they are just like HEX's)
             */
            if (GM_Elem_Types[(ss_elem_list[j - 1]) - 1] == HEXSHELL) {
              ilast_side = 4;
            }
          }

          ss_df_ptr[j] = ss_df_ptr[j - 1] + ilast_side;
          ilast        = ss_df_ptr[j];
        }

        /* Renumber elements to start at '0' instead of '1' */
        for (size_t j = 0; j < num_elem_per_message; ss_elem_list[j++]--) {
          ;
        }

#if 1
        /* This method may be slightly faster, but uses lots of push backs
           that may be bad in the long wrong.  Keep both sets of code in case
           this turns out to be bad...  The other branch is also faster than
           the older code, so it isn't bad to turn it on instead of this... */

        /* Map the ss_elem_list elements to the processor that they are on */
        /* Note that an element may be in ss_elem_list multiple times */
        for (size_t j = 0; j < num_elem_per_message; j++) {
          int iproc = gmap[ss_elem_list[j]];
          int sset  = globals.Proc_Num_Side_Sets[iproc];
          proc_elem_list[iproc][sset].push_back(ss_elem_list[j]);
          proc_side_list[iproc][sset].push_back(ss_side_list[j]);
          proc_gmap_list[iproc][sset].push_back(j + istart_ss);
          proc_df_ptr[iproc][sset].push_back(ss_df_ptr[j]);
        }

        for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
          if (!proc_elem_list[iproc][globals.Proc_Num_Side_Sets[iproc]].empty()) {
            ss_on_proc[iproc] = true;
            proc_num_sides_ss[iproc] =
                proc_elem_list[iproc][globals.Proc_Num_Side_Sets[iproc]].size();
          }
        } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */
#else
        /* Map the ss_elem_list elements to the processor that they are on */
        /* Note that an element may be in ss_elem_list multiple times */

        /* pointers into the global element side-set list for elements
           which are common to the processor and the side-set  */
        std::vector<std::vector<INT>> proc_es_pointer(Proc_Info[2]);

        for (size_t j = 0; j < num_elem_per_message; j++) {
          int iproc = gmap[ss_elem_list[j]];
          int sset  = globals.Proc_Num_Side_Sets[iproc];
          proc_es_pointer[iproc].push_back(j);
        }

        for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

          /*
           * Find the intersection between the elements in the side set and the
           * elements owned by the current processor.
           */

          if (!proc_es_pointer[iproc].empty()) {
            proc_num_sides_ss[iproc] += proc_es_pointer[iproc].size();
            ss_on_proc[iproc] = true;

            /*
             * Store the information for the current side set in temporary arrays.
             *
             * This part of the routine only gets executed if the side set is
             * defined on the current processor.
             */

            /* Allocate and store element information in a temporary vector */
            int sset = globals.Proc_Num_Side_Sets[iproc];
            proc_elem_list[iproc][sset].resize(proc_num_sides_ss[iproc]);
            proc_side_list[iproc][sset].resize(proc_num_sides_ss[iproc]);
            proc_gmap_list[iproc][sset].resize(proc_num_sides_ss[iproc]);
            proc_df_ptr[iproc][sset].resize(proc_num_sides_ss[iproc]);

            /* Transfer the element, side number, and df pointers  into "permanent" storage... */
            for (size_t j = 0; j < proc_es_pointer[iproc].size(); j++) {
              proc_elem_list[iproc][sset][ss_elem_cntr[iproc]] =
                  ss_elem_list[proc_es_pointer[iproc][j]];
              proc_side_list[iproc][sset][ss_elem_cntr[iproc]] =
                  ss_side_list[proc_es_pointer[iproc][j]];
              proc_gmap_list[iproc][sset][ss_elem_cntr[iproc]] =
                  proc_es_pointer[iproc][j] + istart_ss;
              proc_df_ptr[iproc][sset][(ss_elem_cntr[iproc])++] =
                  ss_df_ptr[proc_es_pointer[iproc][j]];
            }
          } /* End "if(ipos_elem > 0)" */
        } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */
#endif
      } /* End "for(imess=0; imess < num_messages; imess++)" */

      /* Entire sideset has been read at this point */
      /*
       * If any part of this side-set is on the processor, update the various
       * pointers, lengths, etc.
       */

      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        if (ss_on_proc[iproc]) {
          ss_proc2glob[iproc][globals.Proc_Num_Side_Sets[iproc]] = i;

          proc_elem_list_ptr[iproc][globals.Proc_Num_Side_Sets[iproc]] = elem_list_length[iproc];

          proc_elem_list_cnt[iproc][globals.Proc_Num_Side_Sets[iproc]] = proc_num_sides_ss[iproc];

          (globals.Proc_Num_Side_Sets[iproc])++;
          elem_list_length[iproc] += proc_num_sides_ss[iproc];
        }
      }

      /* Process any distribution factors in the side set */
      if (num_df_in_ssets[i] > 0) {

        iss_size = sizeof(T);
        find_message_info(iss_size, num_df_in_ssets[i], &num_elem_per_message, &num_messages,
                          &num_left_over);

        if (Debug_Flag >= 4) {
          fmt::print("Message summary for Side Set number {}, with ID of {}:\n", i,
                     Side_Set_Ids[i]);
          fmt::print("\tNumber of messages needed for distribution "
                     "factors = {}\n",
                     num_messages);
          fmt::print("\tNumber of dist. factors in each message = {}\n", num_elem_per_message);
          fmt::print("\tLength of each message = {}\n", (num_elem_per_message * sizeof(T)));
        }

        std::vector<T> ss_dist_fact(num_elem_per_message); /* side-set distribution factors */

        /* set up an array to help speed up the searches below */
        std::vector<std::vector<INT>> proc_df_map(Proc_Info[2]);

        /* Loop over the messages needed for the distribution factors */
        for (size_t imess = 0; imess < num_messages; imess++) {

          size_t istart_ss = imess * num_elem_per_message;
          if (num_left_over != 0 && imess == (num_messages - 1)) {
            num_elem_per_message = num_left_over;
          }

          /* Read in the part of the side set df's that will fit in the msg. */
          check_exodus_error(ex_get_partial_set_dist_fact(exoid, EX_SIDE_SET, Side_Set_Ids[i],
                                                          (istart_ss + 1), num_elem_per_message,
                                                          Data(ss_dist_fact)),
                             "ex_get_partial_set_dist_fact");

          /*
           * At this point a processor has the list of global element IDs
           * that it owns that are contained in the side set. It also has
           * a vector of pointers into the distribution factor vector telling
           * where the distribution factors for a particular element in
           * the side set begin.
           */

          /*
           * First figure out how many df's there are for this side set
           * on this processor.
           */
          for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

            if (imess == 0) {
              ntotal[iproc] = 0;
            }

            if (imess == 0 && ss_on_proc[iproc]) {

              ss_num[iproc] = globals.Proc_Num_Side_Sets[iproc];
              ntotal[iproc] = 0;

              size_t indx     = 0;
              int    ntotal_s = 0; // Nodes per element side...

              size_t cnt = proc_elem_list_cnt[iproc][ss_num[iproc] - 1];
              proc_df_map[iproc].resize(cnt);
              proc_df_indx[iproc][ss_num[iproc] - 1].resize(cnt);
              proc_df_indx_cnt[iproc][ss_num[iproc] - 1].resize(cnt);

              for (INT j = 0; j < proc_elem_list_cnt[iproc][ss_num[iproc] - 1]; j++) {

                size_t loc_elem = proc_elem_list[iproc][ss_num[iproc] - 1][j];
                size_t loc_side = proc_side_list[iproc][ss_num[iproc] - 1][j];

                /*
                 * Kluge for "special" sidesest which has a single distribution
                 * factor for each element instead of 1 per face node..
                 */
                if (num_elem_in_ssets[i] == num_df_in_ssets[i]) {
                  ntotal_s = 1;
                }
                else {

                  ntotal_s = elem_info(NN_SIDE, GM_Elem_Types[loc_elem], loc_side);

                  /*
                   * kludge for HEXSHELL's
                   * where distribution factors are concerned, only count
                   * 4 nodes on the 6 node faces (they are just like HEX's)
                   */
                  if (GM_Elem_Types[loc_elem] == HEXSHELL) {
                    ntotal_s = 4;
                  }
                }

                ntotal[iproc] += ntotal_s;

                proc_df_indx[iproc][ss_num[iproc] - 1][j]     = indx;
                proc_df_indx_cnt[iproc][ss_num[iproc] - 1][j] = 0;

                indx += ntotal_s;

                /* and set up map for sort in a minute */
                proc_df_map[iproc][j] = j;
              }

              /*
               * now sort the proc_df_ptr array so that it is monotonic,
               * and can be searched more easily
               */
              my_sort(proc_elem_list_cnt[iproc][ss_num[iproc] - 1],
                      &proc_df_ptr[iproc][ss_num[iproc] - 1][0], &proc_df_map[iproc][0]);
            }

            /* Allocate memory for the dfs */
            if (ntotal[iproc] > 0) {

              if (imess == 0) {
                proc_ss_df[iproc][ss_num[iproc] - 1].resize(ntotal[iproc]);
              }

              INT k       = 0;
              INT ss_indx = ss_num[iproc] - 1;

              for (size_t j = 0; j < num_elem_per_message; j++) {

                INT indx = istart_ss + j;

                /*
                 * Find out if this index fits in the list of this
                 * processors indices.
                 */

                /* 'istart_ss' <= indx <= 'istart_ss+num_elem_per_message' */
                /* 'indx' increments by 1 each time through loop */

                /* 'proc_df_ptr[iproc][ss_indx]' is the list
                   and it is sorted....
                   ... So, don't need a binary search, can just increment
                   as we proceed through the loops...
                */
                if (indx >= proc_df_ptr[iproc][ss_indx][0]) {

                  while (k < proc_elem_list_cnt[iproc][ss_indx] - 1 &&
                         indx >= proc_df_ptr[iproc][ss_indx][k + 1]) {
                    k++;
                  }

                  /* now get the proper location from the map */
                  INT ipos_elem = proc_df_map[iproc][k];

                  INT idx_begin = proc_df_indx[iproc][ss_indx][ipos_elem];
                  INT idx_end;
                  if ((ipos_elem + 1) < proc_elem_list_cnt[iproc][ss_indx]) {
                    idx_end = proc_df_indx[iproc][ss_indx][ipos_elem + 1];
                  }
                  else {
                    idx_end = ntotal[iproc];
                  }

                  INT idx_diff = idx_end - idx_begin;

                  /*
                   * check to see if indx really fits in this spot
                   *
                   * need to make sure that indx > the k-th value since
                   * bin_search_min will return 0 for anything less than
                   * the first entry in the list
                   */
                  if (indx >= proc_df_ptr[iproc][ss_indx][k] &&
                      indx < proc_df_ptr[iproc][ss_indx][k] + idx_diff) {

                    size_t df_loc = proc_df_indx[iproc][ss_indx][ipos_elem] +
                                    proc_df_indx_cnt[iproc][ss_indx][ipos_elem];

                    proc_ss_df[iproc][ss_indx][df_loc] = ss_dist_fact[indx - istart_ss];

                    proc_ss_df_cnt[iproc][ss_indx]++;

                    proc_df_indx_cnt[iproc][ss_indx][ipos_elem]++;
                  }
                }
              } /* End "for(j=0; j < num_elem_per_message; j++)" */
            } /* End "if(ntotal[iproc] > 0)" */
          } /* End "for(iproc=0; iproc <Proc_Info[2]; iproc++)" */
        } /* End "for(imess=0; imess < num_messages; imess++)" */
      } /* End "if(num_df_in_ssets[i] > 0)" */
    } /* END "if (num_elem_in_ssets[i] > 0)" */
  } /* END "for (i = 0; i < globals.Num_Side_Set; i++)" */

  /*-------------------------------------------------------------------------*/
  /*---------------- Store Structures back into a packed form ---------------*/
  /*-------------------------------------------------------------------------*/

  /*
   * Allocate storage for permanent side set integer info vectors using one
   * long malloc statement.
   */
  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    if (globals.Num_Side_Set > 0) {

      /*
       * Note that the alloc is done based on globals.Num_Side_Set, rather than
       * globals.Proc_Num_Side_Sets[] due to the fact that nullptr entities are
       * stored on processors not having a particular side set.
       */
      globals.Proc_SS_Ids[iproc] =
          (INT *)array_alloc(__FILE__, __LINE__, 1,
                             (3 * globals.Num_Side_Set + 3 * globals.Proc_Num_Side_Sets[iproc] +
                              2 * elem_list_length[iproc] + 1),
                             sizeof(INT));
      globals.Proc_SS_Elem_Count[iproc] = globals.Proc_SS_Ids[iproc] + globals.Num_Side_Set;
      globals.Proc_SS_Elem_Pointers[iproc] =
          globals.Proc_SS_Elem_Count[iproc] + globals.Num_Side_Set;
      globals.Proc_SS_DF_Count[iproc] =
          globals.Proc_SS_Elem_Pointers[iproc] + globals.Proc_Num_Side_Sets[iproc];
      globals.GSide_Sets[iproc] = globals.Proc_SS_DF_Count[iproc] + globals.Num_Side_Set;
      globals.Proc_SS_Elem_List[iproc] =
          globals.GSide_Sets[iproc] + globals.Proc_Num_Side_Sets[iproc];
      globals.Proc_SS_Side_List[iproc] = globals.Proc_SS_Elem_List[iproc] + elem_list_length[iproc];
      globals.Proc_SS_DF_Pointers[iproc] =
          globals.Proc_SS_Side_List[iproc] + elem_list_length[iproc];

      globals.Proc_SS_GEMap_List[iproc] =
          (INT *)array_alloc(__FILE__, __LINE__, 1, elem_list_length[iproc], sizeof(INT));
    }

    /* Construct the side sets global array information */
    globals.Proc_SS_Elem_List_Length[iproc] = 0;
    for (int i = 0; i < globals.Proc_Num_Side_Sets[iproc]; i++) {
      globals.GSide_Sets[iproc][i]            = ss_proc2glob[iproc][i];
      globals.Proc_SS_Ids[iproc][i]           = Side_Set_Ids[ss_proc2glob[iproc][i]];
      globals.Proc_SS_Elem_Count[iproc][i]    = proc_elem_list_cnt[iproc][i];
      globals.Proc_SS_Elem_Pointers[iproc][i] = proc_elem_list_ptr[iproc][i];
      globals.Proc_SS_Elem_List_Length[iproc] += globals.Proc_SS_Elem_Count[iproc][i];
      globals.Proc_SS_DF_Count[iproc][i] = proc_ss_df_cnt[iproc][i];

      if (i == 0) {
        globals.Proc_SS_DF_Pointers[iproc][i] = 0;
      }
      else {
        globals.Proc_SS_DF_Pointers[iproc][i] =
            globals.Proc_SS_DF_Pointers[iproc][i - 1] + proc_ss_df_cnt[iproc][i - 1];
      }
    }

    if (globals.Proc_Num_Side_Sets[iproc] > 0) {
      int i = globals.Proc_Num_Side_Sets[iproc];
      globals.Proc_SS_DF_Pointers[iproc][i] =
          globals.Proc_SS_DF_Pointers[iproc][i - 1] + proc_ss_df_cnt[iproc][i - 1];

      if (globals.Proc_SS_DF_Pointers[iproc][i] > 0) {
        globals.Proc_SS_Dist_Fact[iproc].resize(globals.Proc_SS_DF_Pointers[iproc][i]);
      }
    }

    /* Construct the concatenated element and side list */
    for (int i = 0; i < globals.Proc_Num_Side_Sets[iproc]; i++) {
      for (INT j = 0; j < globals.Proc_SS_Elem_Count[iproc][i]; j++) {
        globals.Proc_SS_Elem_List[iproc][globals.Proc_SS_Elem_Pointers[iproc][i] + j] =
            proc_elem_list[iproc][i][j];
        globals.Proc_SS_Side_List[iproc][globals.Proc_SS_Elem_Pointers[iproc][i] + j] =
            proc_side_list[iproc][i][j];
        globals.Proc_SS_GEMap_List[iproc][globals.Proc_SS_Elem_Pointers[iproc][i] + j] =
            proc_gmap_list[iproc][i][j];
      }
      for (INT j = 0; j < globals.Proc_SS_DF_Count[iproc][i]; j++) {
        globals.Proc_SS_Dist_Fact[iproc][globals.Proc_SS_DF_Pointers[iproc][i] + j] =
            proc_ss_df[iproc][i][j];
      }
    }
  } /* End "for(iproc)" */
} /*      END of routine read_side_sets                                      */
