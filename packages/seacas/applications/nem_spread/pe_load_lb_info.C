/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "el_check_monot.h" // for check_monot
#include "exodusII.h"       // for ex_inquire, ex_opts, etc
#include "fmt/ostream.h"
#include "globals.h"     // for ELEM_COMM_MAP, etc
#include "nem_spread.h"  // for NemSpread, etc
#include "pe_common.h"   // for PEX_MAX
#include "rf_allo.h"     // for array_alloc, safe_free
#include "rf_io_const.h" // for Debug_Flag, Exo_LB_File
#include "rf_util.h"     // for print_line
#include "sort_utils.h"  // for gds_qsort
#include "vector_data.h"
#include <array>
#include <cassert>
#include <cstddef> // for size_t
#include <cstdio>  // for stderr, etc
#include <cstdlib> // for exit
#include <string>

char **qa_record_ptr, **inf_record_ptr;
int    num_inf_rec = 0, num_qa_rec = 0, length_qa = 0;

/****************************************************************************
 ****************************************************************************
 ****************************************************************************
 * This function reads load-balance information for salsa on a parallel
 * computer.
 *
 * Author: Gary L. Hennigan (1421)
 *---------------------------------------------------------------------------
 * Revision History.
 *
 *   Gary Hennigan:
 *      08 November 1993 - Modified scalar version to write information out
 *                         to the parallel disk array once processed.
 ****************************************************************************
 ****************************************************************************
 ****************************************************************************/
template void NemSpread<float, int>::load_lb_info();
template void NemSpread<double, int>::load_lb_info();
template void NemSpread<float, int64_t>::load_lb_info();
template void NemSpread<double, int64_t>::load_lb_info();

template <typename T, typename INT> void NemSpread<T, INT>::load_lb_info()

/*
 *       load_lb_info:
 *
 *            This function reads and distributes the load balance information.
 *    This information is defined as the following scalars, which are specific
 *    to each processor:
 *
 *        Num_Internal_Nodes
 *        Num_Border_Nodes
 *        Num_External_Nodes
 *        globals.Num_Internal_Elems
 *        globals.Num_Border_Elems
 *
 *    and the following vectors, also specific to each processor:
 *
 *        globals.GNodes [Num_Internal_Nodes+Num_Border_Nodes+Num_External_Nodes]
 *        globals.GElems [globals.Num_Internal_Elems+globals.Num_Border_Elems]
 *
 *    For the 1 processor case, this routine fills in appropriate values
 *    for these quantities.
 */
{

  int   lb_exoid      = 0;
  INT   cmap_max_size = 0;
  char  Title[MAX_LINE_LENGTH + 1];
  float version;

  int cpu_ws = 0;

  /******************************** START EXECUTION ****************************/
  if (Debug_Flag != 0) {
    fmt::print("\nStart to read in and distribute the load balance info\n");
  }

  /* Open the Load Balance exoII file for reading */

  fmt::print("EXODUS II load-balance file: {}\n", Exo_LB_File);
  cpu_ws     = io_ws;
  int mode   = EX_READ | int64api;
  int iio_ws = 0; // Don't interfere with exodus files; this is the nemesis file.
  if ((lb_exoid = ex_open(Exo_LB_File.c_str(), mode, &cpu_ws, &iio_ws, &version)) == -1) {
    fmt::print(stderr, "[{}] ERROR: Could not open lb file, {}\n", __func__, Exo_LB_File);
    exit(1);
  }

  /* Read information about the processor configuration */
  read_proc_init(lb_exoid, Proc_Info, Proc_Ids);

  /* Allocate space for the counts */
  globals.Num_Internal_Nodes.resize(Proc_Info[2]);
  globals.Num_Border_Nodes.resize(Proc_Info[2]);
  globals.Num_External_Nodes.resize(Proc_Info[2]);
  globals.Num_Internal_Elems.resize(Proc_Info[2]);
  globals.Num_Border_Elems.resize(Proc_Info[2]);
  globals.Num_N_Comm_Maps.resize(Proc_Info[2]);
  globals.Num_E_Comm_Maps.resize(Proc_Info[2]);

  /* Allocate space for each processor entity */
  globals.GNodes.resize(Proc_Info[2]);
  globals.GElems.resize(Proc_Info[2]);
  globals.Elem_Map.resize(Proc_Info[2]);

  /* Allocate contiguous space for the pointer vectors on all processors */
  std::vector<INT> Int_Space(1);
  std::vector<INT> Int_Node_Num(Proc_Info[0]);
  std::vector<INT> Bor_Node_Num(Proc_Info[0]);
  std::vector<INT> Ext_Node_Num(Proc_Info[0]);
  std::vector<INT> Int_Elem_Num(Proc_Info[0]);
  std::vector<INT> Bor_Elem_Num(Proc_Info[0]);
  std::vector<INT> Node_Comm_Num(Proc_Info[0]);
  std::vector<INT> Elem_Comm_Num(Proc_Info[0]);

  /* Read the initial information contained in the load balance file */
  read_lb_init(lb_exoid, Int_Space, Int_Node_Num, Bor_Node_Num, Ext_Node_Num, Int_Elem_Num,
               Bor_Elem_Num, Node_Comm_Num, Elem_Comm_Num, Title);

  /* Allocate memory for the communication map arrays */
  globals.N_Comm_Map.resize(Proc_Info[2]);
  globals.E_Comm_Map.resize(Proc_Info[2]);

  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {

    /*
     * Error check:
     *  Currently a maximum of one nodal communication map and one
     * elemental communication map is supported.
     */
    if (globals.Num_N_Comm_Maps[iproc] > 1 || globals.Num_E_Comm_Maps[iproc] > 1) {
      fmt::print(stderr,
                 "[{}] ERROR. Only 1 nodal and elemental comm map "
                 "is supported\n",
                 __func__);
      exit(1);
    }
    else {
      /* Always allocate at least one and initialize the counts to 0 */
      globals.N_Comm_Map[iproc].node_cnt = 0;
      globals.E_Comm_Map[iproc].elem_cnt = 0;
    }
  } /* End "for (int iproc=0; iproc <Proc_Info[2]; iproc++)" */

  /* Set up each processor for the communication map parameters */
  read_cmap_params(lb_exoid, globals.E_Comm_Map, globals.N_Comm_Map, &cmap_max_size);

  /*
   * loop through the processors, one at a time, to read
   * their load balance information
   *
   * NOTE: From here on there are no provisions for multiple nodal
   *       or elemental communication maps.
   */

  for (int iproc = 0; iproc < Proc_Info[0]; iproc++) {

    /* Get the node map for processor "iproc" */
    size_t itotal_nodes = globals.Num_Internal_Nodes[iproc] + globals.Num_Border_Nodes[iproc] +
                          globals.Num_External_Nodes[iproc];
    globals.GNodes[iproc].resize(itotal_nodes);
    if (ex_get_processor_node_maps(
            lb_exoid, &globals.GNodes[iproc][0], &globals.GNodes[iproc][Int_Node_Num[iproc]],
            &globals.GNodes[iproc][Int_Node_Num[iproc] + Bor_Node_Num[iproc]], iproc) < 0) {
      fmt::print(stderr, "[{}] ERROR, failed to get node map for Proc {}!\n", __func__, iproc);
      exit(1);
    }

    /* Get the element map for processor number "iproc" */
    size_t itotal_elems = globals.Num_Internal_Elems[iproc] + globals.Num_Border_Elems[iproc];
    globals.GElems[iproc].resize(itotal_elems);
    globals.Elem_Map[iproc].resize(itotal_elems);

    if (ex_get_processor_elem_maps(lb_exoid, &globals.GElems[iproc][0],
                                   &globals.GElems[iproc][Int_Elem_Num[iproc]], iproc) < 0) {
      fmt::print(stderr, "[{}] ERROR, failed to get element map for Proc {}!\n", __func__, iproc);
      exit(1);
    }
    globals.Elem_Map[iproc] = globals.GElems[iproc];

    if (Node_Comm_Num[iproc] > 0) {
      globals.N_Comm_Map[iproc].node_ids.resize(globals.N_Comm_Map[iproc].node_cnt);
      globals.N_Comm_Map[iproc].proc_ids.resize(globals.N_Comm_Map[iproc].node_cnt);

      if (ex_get_node_cmap(lb_exoid, globals.N_Comm_Map[iproc].map_id,
                           Data(globals.N_Comm_Map[iproc].node_ids),
                           Data(globals.N_Comm_Map[iproc].proc_ids), iproc) < 0) {
        /*
         * If there are disconnected mesh pieces, then it is
         * possible that there is no communication between the
         * pieces and there will be no communication maps.  Normally
         * this is a problem, so output a warning, but don't abort.
         */
        fmt::print(stderr, "[{}] WARNING. Failed to get nodal comm map for Proc {}!\n", __func__,
                   iproc);
      }
    }

    if (Elem_Comm_Num[iproc] > 0) {
      globals.E_Comm_Map[iproc].elem_ids.resize(globals.E_Comm_Map[iproc].elem_cnt);
      globals.E_Comm_Map[iproc].side_ids.resize(globals.E_Comm_Map[iproc].elem_cnt);
      globals.E_Comm_Map[iproc].proc_ids.resize(globals.E_Comm_Map[iproc].elem_cnt);

      if (ex_get_elem_cmap(lb_exoid, globals.E_Comm_Map[iproc].map_id,
                           Data(globals.E_Comm_Map[iproc].elem_ids),
                           Data(globals.E_Comm_Map[iproc].side_ids),
                           Data(globals.E_Comm_Map[iproc].proc_ids), iproc) < 0) {
        fmt::print(stderr, "[{}] ERROR. Failed to get elemental comm map for Proc {}!\n", __func__,
                   iproc);
        exit(1);
      }
    }

    /*
     * Communicate load balance information to the correct processor
     *   - if iproc = Proc_Ids[*] then process the data instead.
     */
    assert(Proc_Ids[iproc] == iproc);

    /*
     * Sort the local element numbers in ascending global element numbers.
     * This means that globals.GElems will be monotonic.
     */
    gds_qsort(Data(globals.GElems[iproc]), globals.Num_Internal_Elems[iproc]);
    gds_qsort(Data(globals.Elem_Map[iproc]), globals.Num_Internal_Elems[iproc]);

    /* Check that globals.GNodes is monotonic, from i = 0 to Num_Internal_Nodes */
#ifdef DEBUG
    assert(check_monot(&globals.GNodes[iproc], globals.Num_Internal_Nodes[iproc]));

    /*
     * Check that globals.GNodes is monotonic, from i = Num_Internal_Nodes to
     *    (Num_Internal_Nodes + Num_Border_Nodes)
     */
    assert(check_monot(&(globals.GNodes[iproc][globals.Num_Internal_Nodes[iproc]]),
                       globals.Num_Border_Nodes[iproc]));
#endif
  }

  /* Close the load balance file - we are finished with it */
  if (ex_close(lb_exoid) == -1) {
    fmt::print(stderr, "[{}] ERROR: Error in closing load balance file\n", __func__);
    exit(1);
  }

  /************************* Cleanup and Printout Phase ***********************/

  if (num_qa_rec > 0) {
    for (int i = 0; i < length_qa; i++) {
      safe_free(reinterpret_cast<void **>(&(qa_record_ptr[i])));
    }
    safe_free(reinterpret_cast<void **>(&qa_record_ptr));
  }

  if (num_inf_rec > 0) {
    for (int i = 0; i < num_inf_rec; i++) {
      safe_free(reinterpret_cast<void **>(&(inf_record_ptr[i])));
    }
    safe_free(reinterpret_cast<void **>(&inf_record_ptr));
  }

  for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
    if (globals.Num_Internal_Nodes[iproc] == 0 && globals.Num_Border_Nodes[iproc] == 0 &&
        globals.Num_External_Nodes[iproc] == 0) {
      fmt::print(stderr, "\n[{}] WARNING, Processor {} has no nodes!\n", __func__, iproc);
    }
    if (globals.Num_Internal_Elems[iproc] == 0 && globals.Num_Border_Elems[iproc] == 0) {
      fmt::print(stderr, "\n[{}] WARNING, Processor {} has no elements!\n", __func__, iproc);
    }
  }

  /*========================================================================*/

  if (Debug_Flag != 0) {
    fmt::print("\nFinished distributing load balance info\n");
  }

  /* Output Detailed timing information for the program */
  /*
   * Print out a Large table of Load Balance Information if the debug_flag
   * setting is large enough
   */

  if (Debug_Flag >= 7) {
    fmt::print("\n\n");
    print_line("=", 79);
    for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
      fmt::print("\n\t***For Processor {}***\n", Proc_Ids[iproc]);
      fmt::print("\tInternal nodes owned by the current processor\n\t");
      for (INT i = 0; i < globals.Num_Internal_Nodes[iproc]; i++) {
        fmt::print(" {}", globals.GNodes[iproc][i]);
      }

      fmt::print("\n");

      fmt::print("\tBorder nodes owned by the current processor\n\t");
      for (INT i = 0; i < globals.Num_Border_Nodes[iproc]; i++) {
        fmt::print(" {}", globals.GNodes[iproc][i + globals.Num_Internal_Nodes[iproc]]);
      }

      fmt::print("\n");

      if (globals.Num_External_Nodes[iproc] > 0) {
        fmt::print("\tExternal nodes needed by the current processor\n\t");
        for (INT i = 0; i < globals.Num_External_Nodes[iproc]; i++) {
          fmt::print(" {}", globals.GNodes[iproc][i + globals.Num_Internal_Nodes[iproc] +
                                                  globals.Num_Border_Nodes[iproc]]);
        }

        fmt::print("\n");
      }

      fmt::print("\tInternal elements owned by the current processor\n\t");
      for (INT i = 0; i < globals.Num_Internal_Elems[iproc]; i++) {
        fmt::print(" {}", globals.GElems[iproc][i]);
      }

      fmt::print("\n");

      if (globals.Num_Border_Elems[iproc] > 0) {
        fmt::print("\tBorder elements owned by the current processor\n\t");
        for (INT i = 0; i < globals.Num_Border_Elems[iproc]; i++) {
          fmt::print(" {}", globals.GElems[iproc][i + globals.Num_Internal_Elems[iproc]]);
        }

        fmt::print("\n");
      }

      if (globals.Num_N_Comm_Maps[iproc] > 0) {
        fmt::print("\tNodal Comm Map for the current processor\n");
        fmt::print("\t\tnode IDs:");
        for (size_t i = 0; i < globals.N_Comm_Map[iproc].node_cnt; i++) {
          fmt::print(" {}", globals.N_Comm_Map[iproc].node_ids[i]);
        }
        fmt::print("\n\t\tproc IDs:");
        for (size_t i = 0; i < globals.N_Comm_Map[iproc].node_cnt; i++) {
          fmt::print(" {}", globals.N_Comm_Map[iproc].proc_ids[i]);
        }

        fmt::print("\n");
      }

      if (globals.Num_E_Comm_Maps[iproc] > 0) {
        fmt::print("\tElemental Comm Map for the current processor\n");
        fmt::print("\t\telement IDs:");
        for (size_t i = 0; i < globals.E_Comm_Map[iproc].elem_cnt; i++) {
          fmt::print(" {}", globals.E_Comm_Map[iproc].elem_ids[i]);
        }
        fmt::print("\n\t\tside IDs:");
        for (size_t i = 0; i < globals.E_Comm_Map[iproc].elem_cnt; i++) {
          fmt::print(" {}", globals.E_Comm_Map[iproc].side_ids[i]);
        }
        fmt::print("\n\t\tproc IDs:");
        for (size_t i = 0; i < globals.E_Comm_Map[iproc].elem_cnt; i++) {
          fmt::print(" {}", globals.E_Comm_Map[iproc].proc_ids[i]);
        }

        fmt::print("\n");
      }
    }
    fmt::print("\n");
    print_line("=", 79);
  }

} /* END of routine load_lb_info () ******************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_proc_init(int lb_exoid, std::array<int, 6> &proc_info,
                                       std::vector<int> &proc_ids)

/*----------------------------------------------------------------------------
 *  read_proc_init:
 *
 *      This function reads information about the processor configuration
 * which the load balance was generated for and makes assignments for each
 * processor.
 *----------------------------------------------------------------------------
 * Variable Glossary (after return):
 *
 *              proc_info[0] = # procs, from load balance file
 *              proc_info[1] = # procs for, from load balance file
 *              proc_info[2] = # procs this processor is responsible for
 *              proc_info[3] = # of extra procs
 *
 */
{
  char ftype[2];
  if (ex_get_init_info(lb_exoid, &proc_info[0], &proc_info[1], ftype) < 0) {
    fmt::print(stderr, "[{}] ERROR, could not get init info!\n", __func__);
    exit(1);
  }

  /* Calculate which processor is responsible for what */
  proc_info[2] = proc_info[0];
  proc_ids.resize(proc_info[2]);

  for (int i1 = 0; i1 < proc_info[2]; i1++) {
    proc_ids[i1] = i1;
  }
} /* End of read_proc_init() *************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_lb_init(int lb_exoid, std::vector<INT> &Int_Space,
                                     std::vector<INT> &Int_Node_Num, std::vector<INT> &Bor_Node_Num,
                                     std::vector<INT> &Ext_Node_Num, std::vector<INT> &Int_Elem_Num,
                                     std::vector<INT> &Bor_Elem_Num,
                                     std::vector<INT> &Node_Comm_Num,
                                     std::vector<INT> &Elem_Comm_Num, char * /*Title*/)

/*
 *    read_lb_init:
 *
 *         This function reads the initial information contained in the
 * load balance file
 *
 *
 */
{
  /*********************BEGIN EXECUTABLE STATEMENTS****************************/

  /* Read Set-up information from the load-balance file on Proc 0 */

  /*
   * If debugging is not on go ahead and report errors from init. This
   * will show version mismatch information by default.
   */
  int old_opt = 0;
  if (Debug_Flag == 0) {
    old_opt = ex_opts(EX_VERBOSE);
  }

  /* Read the title of the LB File and about the size of the mesh */
  INT num_nodes;
  INT num_elem;
  INT num_elem_blk;
  INT num_node_sets;
  INT num_side_sets;
  int error = ex_get_init_global(lb_exoid, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                                 &num_side_sets);

  check_exodus_error(error, "ex_get_init");

  if (Debug_Flag == 0) {
    ex_opts(old_opt);
  }

#ifdef DEBUG
  if (Debug_Flag >= 2) {
    fmt::print("---------------------------------------------------------\n"
               "\t\tLoad balance file global information\n"
               "---------------------------------------------------------\n"
               "\tNumber of nodes:          {}\n"
               "\tNumber of elements:       {}\n"
               "\tNumber of element blocks: {}\n"
               "---------------------------------------------------------\n",
               fmt::group_digits(num_nodes), fmt::group_digits(num_elem),
               fmt::group_digits(num_elem_blk));
  }
#endif

  /* Cross-check the load balance file info against the mesh info */
  if (((size_t)num_nodes != globals.Num_Node) || ((size_t)num_elem != globals.Num_Elem) ||
      (num_elem_blk != globals.Num_Elem_Blk)) {
    fmt::print(stderr, "[{}] ERROR: Problem dimensions in the LB File don't match with those \
in mesh file",
               __func__);
    exit(1);
  }

  /* Read the QA Records */
  num_qa_rec = ex_inquire_int(lb_exoid, EX_INQ_QA);
  if (num_qa_rec > 0) {
    length_qa = 4 * num_qa_rec;
    qa_record_ptr =
        reinterpret_cast<char **>(array_alloc(__FILE__, __LINE__, 1, length_qa, sizeof(char *)));
    for (int i = 0; i < length_qa; i++) {
      qa_record_ptr[i] = reinterpret_cast<char *>(
          array_alloc(__FILE__, __LINE__, 1, (MAX_STR_LENGTH + 1), sizeof(char)));
    }
    error = ex_get_qa(lb_exoid, reinterpret_cast<char *(*)[4]>(&qa_record_ptr[0]));
    check_exodus_error(error, "ex_get_qa");
  }

  /* Read the Info Records */
  num_inf_rec = ex_inquire_int(lb_exoid, EX_INQ_INFO);
  if (num_inf_rec > 0) {
    inf_record_ptr =
        reinterpret_cast<char **>(array_alloc(__FILE__, __LINE__, 1, num_inf_rec, sizeof(char *)));
    for (int i = 0; i < num_inf_rec; i++) {
      inf_record_ptr[i] = reinterpret_cast<char *>(
          array_alloc(__FILE__, __LINE__, 1, (MAX_LINE_LENGTH + 2), sizeof(char)));
    }
    error = ex_get_info(lb_exoid, inf_record_ptr);
    check_exodus_error(error, "ex_get_info");
  }

  Int_Space[0] = 0;

  for (int i = 0; i < Proc_Info[0]; i++) {

    if (ex_get_loadbal_param(lb_exoid, &Int_Node_Num[i], &Bor_Node_Num[i], &Ext_Node_Num[i],
                             &Int_Elem_Num[i], &Bor_Elem_Num[i], &Node_Comm_Num[i],
                             &Elem_Comm_Num[i], i) < 0) {
      fmt::print(stderr, "[{}] ERROR, could not get load balance params!\n", __func__);
      exit(1);
    }

#ifdef DEBUG
    if (Debug_Flag >= 5) {
      if (i == 0) {
        fmt::print("--------------------------------------------------------\n"
                   "\t\tLoad balance parameters as read by Processor 0\n"
                   "--------------------------------------------------------\n");
      }
      fmt::print("Read on processor 0 for processor {}\n"
                 "\tNumber internal nodes:         {}\n"
                 "\tNumber border nodes:           {}\n"
                 "\tNumber external nodes:         {}\n"
                 "\tNumber internal elements:      {}\n"
                 "\tNumber border elements:        {}\n"
                 "\tNumber of nodal comm maps:     {}\n"
                 "\tNumber of elemental comm maps: {}\n"
                 "--------------------------------------------------------\n",
                 fmt::group_digits(i), fmt::group_digits(Int_Node_Num[i]),
                 fmt::group_digits(Bor_Node_Num[i]), fmt::group_digits(Ext_Node_Num[i]),
                 fmt::group_digits(Int_Elem_Num[i]), fmt::group_digits(Bor_Elem_Num[i]),
                 fmt::group_digits(Node_Comm_Num[i]), fmt::group_digits(Elem_Comm_Num[i]));
    }
#endif /* DEBUG */

    Int_Space[0] = PEX_MAX(Int_Space[0], Int_Node_Num[i] + Bor_Node_Num[i] + Ext_Node_Num[i] +
                                             Int_Elem_Num[i] + Bor_Elem_Num[i]);
  }

  /* Each processor extracts the information that it needs */
  for (int i = 0; i < Proc_Info[2]; i++) {
    globals.Num_Internal_Nodes[i] = Int_Node_Num[i];
    globals.Num_Border_Nodes[i]   = Bor_Node_Num[i];
    globals.Num_External_Nodes[i] = Ext_Node_Num[i];
    globals.Num_Internal_Elems[i] = Int_Elem_Num[i];
    globals.Num_Border_Elems[i]   = Bor_Elem_Num[i];
    globals.Num_N_Comm_Maps[i]    = Node_Comm_Num[i];
    globals.Num_E_Comm_Maps[i]    = Elem_Comm_Num[i];
  }

  /*
   * Print Out a Summary of the Load Balance Information, as distributed
   * across the processors, for the appropriate value of Debug_Flag
   */
  if (Debug_Flag >= 3) {
    print_line("=", 79);
    fmt::print("\n\t\tTABLE OF LOAD BALANCE STATISTICS\n\n");
    fmt::print("{}{}\n", "globals. Int_Nodes Bor_Nodes Ext_Nodes",
               " Int_Elems Bor_Elems N_Comm_Maps E_Comm_Maps");
    print_line("-", 79);
    fmt::print("\n");
    for (int i = 0; i < Proc_Info[2]; i++) {
      fmt::print("{:6d}  {:6d}  {:6d}   {:6d}    {:6d}    {:6d}     {:6d}     {:6d}\n", Proc_Ids[i],
                 globals.Num_Internal_Nodes[i], globals.Num_Border_Nodes[i],
                 globals.Num_External_Nodes[i], globals.Num_Internal_Elems[i],
                 globals.Num_Border_Elems[i], globals.Num_N_Comm_Maps[i],
                 globals.Num_E_Comm_Maps[i]);
    }
    print_line("=", 79);
    fmt::print("\n\n");
  }

} /* END of read_lb_init () **************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template <typename T, typename INT>
void NemSpread<T, INT>::read_cmap_params(int lb_exoid, std::vector<ELEM_COMM_MAP<INT>> &E_Comm_Map,
                                         std::vector<NODE_COMM_MAP<INT>> &N_Comm_Map,
                                         INT                             *cmap_max_size)
{
  for (int iproc = 0; iproc < Proc_Info[0]; iproc++) {

    std::array<INT, 1> node_cm_ids{0};
    std::array<INT, 1> node_cm_cnts{0};
    std::array<INT, 1> elem_cm_ids{0};
    std::array<INT, 1> elem_cm_cnts{0};

    /* Read the communication map IDs for processor "iproc" */
    if (ex_get_cmap_params(lb_exoid, Data(node_cm_ids), Data(node_cm_cnts), Data(elem_cm_ids),
                           Data(elem_cm_cnts), iproc) < 0) {
      fmt::print(stderr, "[{}] ERROR, unable to read communication map params\n", __func__);
      exit(1);
    }

    N_Comm_Map[iproc].map_id   = node_cm_ids[0];
    N_Comm_Map[iproc].node_cnt = node_cm_cnts[0];
    E_Comm_Map[iproc].map_id   = elem_cm_ids[0];
    E_Comm_Map[iproc].elem_cnt = elem_cm_cnts[0];

    INT psum       = 2 * node_cm_cnts[0] + 3 * elem_cm_cnts[0];
    *cmap_max_size = std::max(*cmap_max_size, psum);
  }

  if (Debug_Flag >= 4) {
    print_line("=", 79);
    fmt::print("\t\tCOMMUNICATION MAP INFORMATION\n");
    fmt::print("\t\t   largest cmap = {} integers\n", *cmap_max_size);
    print_line("=", 79);

    int bprnt_n = 0;
    int bprnt_e = 0;

    for (int i1 = 0; i1 < Proc_Info[2]; i1++) {
      if (globals.Num_N_Comm_Maps[i1] > 0) {
        bprnt_n = 1;
      }
      if (globals.Num_E_Comm_Maps[i1] > 0) {
        bprnt_e = 1;
      }
    }

    if (bprnt_n > 0) {
      fmt::print("\tFor Proc\tNode Map ID\tNode Count\n");
      fmt::print("\t------------------------------------------------\n");
      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        if (globals.Num_N_Comm_Maps[iproc] > 0) {
          fmt::print("\t     {}\t\t    {}\t\t    {}\n", Proc_Ids[iproc], N_Comm_Map[iproc].map_id,
                     N_Comm_Map[iproc].node_cnt);
        }
      }
    }

    if (bprnt_e > 0) {
      fmt::print("\tFor Proc\tElem Map ID\tElem Count\n");
      fmt::print("\t------------------------------------------------\n");
      for (int iproc = 0; iproc < Proc_Info[2]; iproc++) {
        if (globals.Num_E_Comm_Maps[iproc] > 0) {
          fmt::print("\t     {}\t\t    {}\t\t    {}\n", Proc_Ids[iproc], E_Comm_Map[iproc].map_id,
                     E_Comm_Map[iproc].elem_cnt);
        }
      }
    }

    print_line("=", 79);
  }

} /* End of read_cmap_params() ***********************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
