/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "elb.h"      // for Machine_Description, etc
#include "elb_allo.h" // for array_alloc
#include "elb_elem.h" // for NNODES, get_elem_info
#include "elb_err.h"  // for Gen_Error, error_lev
#include "elb_output.h"
#include "elb_util.h" // for gds_qsort, qsort2, in_list, etc
#include "fmt/chrono.h"
#include "fmt/ostream.h"
#include "scopeguard.h"
#include <copy_string_cpp.h>
#include <cstddef>    // for size_t, nullptr
#include <cstdlib>    // for free, malloc, realloc
#include <cstring>    // for strcat, strlen, etc
#include <ctime>      // for asctime, localtime, time, etc
#include <exodusII.h> // for ex_close, ex_opts, etc
#include <iostream>
#include <vector> // for vector

namespace {
  std::string remove_extension(const std::string &filename)
  {
    // Strip off the extension
    size_t ind = filename.find_last_of('.', filename.size());
    if (ind != std::string::npos) {
      return filename.substr(0, ind);
    }
    return filename;
  }
} // namespace

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs a load balance file using the ExodusII and NemesisI
 * API.
 *****************************************************************************/
template int write_nemesis(std::string &nemI_out_file, Machine_Description *machine,
                           Problem_Description *problem, Mesh_Description<int> *mesh,
                           LB_Description<int> *lb, Sphere_Info *sphere);
template int write_nemesis(std::string &nemI_out_file, Machine_Description *machine,
                           Problem_Description *problem, Mesh_Description<int64_t> *mesh,
                           LB_Description<int64_t> *lb, Sphere_Info *sphere);

template <typename INT>
int write_nemesis(std::string &nemI_out_file, Machine_Description *machine,
                  Problem_Description *problem, Mesh_Description<INT> *mesh,
                  LB_Description<INT> *lb, Sphere_Info *sphere)
{
  int         exoid;
  std::string method1{}, method2{};
  std::string title;

  int cpu_ws = sizeof(float);
  int io_ws  = sizeof(float);

  fmt::print("Outputting load balance to file {}\n", nemI_out_file.c_str());

  /* Create the load balance file */
  /* Attempt to create a netcdf4-format file; if it fails, then assume
     that the netcdf library does not support that mode and fall back
     to classic netcdf3 format.  If that fails, issue an error and
     return failure.
  */
  int mode3 = EX_CLOBBER;
  int mode4 = mode3 | EX_NETCDF4 | EX_NOCLASSIC | problem->int64db | problem->int64api;

  ex_opts(EX_DEFAULT); // Eliminate misleading error if the first ex_create fails, but the second
                       // succeeds.
  if ((exoid = ex_create(nemI_out_file.c_str(), mode4, &cpu_ws, &io_ws)) < 0) {
    /* If int64api or int64db non-zero, then netcdf-4 format is required, so
       fail now...
    */
    if (problem->int64db | problem->int64api) {
      Gen_Error(0, "fatal: failed to create Nemesis netcdf-4 file");
      return 0;
    }
    if ((exoid = ex_create(nemI_out_file.c_str(), mode3, &cpu_ws, &io_ws)) < 0) {
      Gen_Error(0, "fatal: failed to create Nemesis file");
      return 0;
    }
  }
  ON_BLOCK_EXIT(ex_close, exoid);

  /* Set the error reporting value */
  if (error_lev > 1) {
    ex_opts(EX_VERBOSE | EX_DEBUG);
  }
  else {
    ex_opts(EX_VERBOSE);
  }

  /* Enable compression (if netcdf-4) */
  ex_set_option(exoid, EX_OPT_COMPRESSION_LEVEL, 1);
  ex_set_option(exoid, EX_OPT_COMPRESSION_SHUFFLE, 1);

  /* Create the title */
  if (problem->type == NODAL) {
    method1 = "nodal";
  }
  else {
    method1 = "elemental";
  }

  title = fmt::format("nem_slice {} load balance file", method1.c_str());

  method1 = "method1: ";
  method2 = "method2: ";

  switch (lb->type) {
  case MULTIKL:
    method1 += "Multilevel-KL decomposition";
    method2 += "With Kernighan-Lin refinement";
    break;
  case SPECTRAL: method1 += "Spectral decomposition"; break;
  case INERTIAL: method1 += "Inertial decomposition"; break;
  case ZPINCH: method1 += "ZPINCH decomposition"; break;
  case BRICK: method1 += "BRICK decomposition"; break;
  case ZOLTAN_RCB: method1 += "RCB decomposition"; break;
  case ZOLTAN_RIB: method1 += "RIB decomposition"; break;
  case ZOLTAN_HSFC: method1 += "HSFC decomposition"; break;
  case LINEAR: method1 += "Linear decomposition"; break;
  case RANDOM: method1 += "Random decomposition"; break;
  case SCATTERED: method1 += "Scattered decomposition"; break;
  }

  if (lb->refine == KL_REFINE && lb->type != MULTIKL) {
    method2 += "with Kernighan-Lin refinement";
  }
  else if (lb->type != MULTIKL) {
    method2 += "no refinement";
  }

  switch (lb->num_sects) {
  case 1: method1 += " via bisection"; break;
  case 2: method1 += " via quadrasection"; break;
  case 3: method1 += " via octasection"; break;
  }

  /* Do some sorting */
  for (int proc = 0; proc < machine->num_procs; proc++) {

    /* Sort node maps */
    gds_qsort(lb->int_nodes[proc].data(), lb->int_nodes[proc].size());
    if (problem->type == NODAL) {
      sort2(lb->ext_nodes[proc].size(), lb->ext_nodes[proc].data(), lb->ext_procs[proc].data());
    }

    /* Sort element maps */
    gds_qsort(lb->int_elems[proc].data(), lb->int_elems[proc].size());
  }

  /* Output the info records */
  char *info[3];
  info[0] = const_cast<char *>(title.c_str());
  info[1] = const_cast<char *>(method1.c_str());
  info[2] = const_cast<char *>(method2.c_str());

  if (ex_put_info(exoid, 3, info) < 0) {
    Gen_Error(0, "warning: output of info records failed");
  }

  /* Generate a QA record for the utility */
  time_t      time_val = time(nullptr);
  auto *      lt       = std::localtime(&time_val);
  std::string time     = fmt::format("{:%H:%M:%S}", *lt);
  std::string date     = fmt::format("{:%Y/%m/%d}", *lt);

  char qa_date[15];
  char qa_time[10];
  char qa_name[MAX_STR_LENGTH];
  char qa_vers[10];

  copy_string(qa_time, time);
  copy_string(qa_date, date);
  copy_string(qa_name, UTIL_NAME);
  copy_string(qa_vers, ELB_VERSION);

  char **lqa_record = reinterpret_cast<char **>(array_alloc(1, 4, sizeof(char *)));
  for (int i2 = 0; i2 < 4; i2++) {
    lqa_record[i2] = reinterpret_cast<char *>(array_alloc(1, MAX_STR_LENGTH + 1, sizeof(char)));
  }

  copy_string(lqa_record[0], qa_name, MAX_STR_LENGTH + 1);
  copy_string(lqa_record[1], qa_vers, MAX_STR_LENGTH + 1);
  copy_string(lqa_record[2], qa_date, MAX_STR_LENGTH + 1);
  copy_string(lqa_record[3], qa_time, MAX_STR_LENGTH + 1);

  fmt::print("QA Record:\n");
  for (int i2 = 0; i2 < 4; i2++) {
    fmt::print("\t{}\n", lqa_record[i2]);
  }

  if (ex_put_qa(exoid, 1, reinterpret_cast<char *(*)[4]>(&lqa_record[0])) < 0) {
    Gen_Error(0, "fatal: unable to output QA records");
    return 0;
  }

  /* free up memory */
  for (int i2 = 0; i2 < 4; i2++) {
    free(lqa_record[i2]);
  }

  free(lqa_record);

  /* Output the the initial Nemesis global information */
  if (ex_put_init_global(exoid, mesh->num_nodes, mesh->num_elems, mesh->num_el_blks, 0, 0) < 0) {
    Gen_Error(0, "fatal: failed to output initial Nemesis parameters");
    return 0;
  }

  ex_put_eb_info_global(exoid, mesh->eb_ids.data(), mesh->eb_cnts.data());

  /* Set up dummy arrays for output */
  std::vector<INT> num_nmap_cnts(machine->num_procs);
  std::vector<INT> num_emap_cnts(machine->num_procs);

  if (problem->type == NODAL) {
    /* need to check and make sure that there really are comm maps */
    for (int cnt = 0; cnt < machine->num_procs; cnt++) {
      if (!lb->bor_nodes[cnt].empty()) {
        num_nmap_cnts[cnt] = 1;
      }
    }
  }
  else { /* Elemental load balance */
    if (problem->num_vertices > sphere->num) {
      /* need to check and make sure that there really are comm maps */
      for (int cnt = 0; cnt < machine->num_procs; cnt++) {
        if (!lb->bor_nodes[cnt].empty()) {
          num_nmap_cnts[cnt] = 1;
        }
      }
      for (int cnt = 0; cnt < machine->num_procs; cnt++) {
        if (!lb->bor_elems[cnt].empty()) {
          num_emap_cnts[cnt] = 1;
        }
      }
    }
  }

  if (ex_put_init_info(exoid, machine->num_procs, machine->num_procs, (char *)"s") < 0) {
    Gen_Error(0, "fatal: unable to output init info");
    return 0;
  }

  // Need to create 5 arrays with the sizes of lb->int_nodes[i].size()...
  {
    std::vector<INT> ins(machine->num_procs);
    std::vector<INT> bns(machine->num_procs);
    std::vector<INT> ens(machine->num_procs);
    std::vector<INT> ies(machine->num_procs);
    std::vector<INT> bes(machine->num_procs);

    for (int iproc = 0; iproc < machine->num_procs; iproc++) {
      ins[iproc] = lb->int_nodes[iproc].size();
      bns[iproc] = lb->bor_nodes[iproc].size();
      ens[iproc] = lb->ext_nodes[iproc].size();
      ies[iproc] = lb->int_elems[iproc].size();
      bes[iproc] = lb->bor_elems[iproc].size();
    }

    if (ex_put_loadbal_param_cc(exoid, ins.data(), bns.data(), ens.data(), ies.data(), bes.data(),
                                num_nmap_cnts.data(), num_emap_cnts.data()) < 0) {
      Gen_Error(0, "fatal: unable to output load-balance parameters");
      return 0;
    }
  }

  if (problem->type == NODAL) /* Nodal load balance output */
  {
    /* Set up for the concatenated communication map parameters */
    std::vector<INT> node_proc_ptr(machine->num_procs + 1);
    std::vector<INT> node_cmap_ids_cc(machine->num_procs);
    std::vector<INT> node_cmap_cnts_cc(machine->num_procs);

    node_proc_ptr[0] = 0;
    for (int proc = 0; proc < machine->num_procs; proc++) {
      node_proc_ptr[proc + 1] = node_proc_ptr[proc] + 1;
      node_cmap_cnts_cc[proc] = lb->ext_nodes[proc].size();
      node_cmap_ids_cc[proc]  = 1;
    }

    /* Output the communication map parameters */
    if (ex_put_cmap_params_cc(exoid, node_cmap_ids_cc.data(), node_cmap_cnts_cc.data(),
                              node_proc_ptr.data(), nullptr, nullptr, nullptr) < 0) {
      Gen_Error(0, "fatal: unable to output communication map parameters");
      return 0;
    }

    /* Output the node and element maps */
    for (int proc = 0; proc < machine->num_procs; proc++) {
      /* Output the nodal map */
      if (ex_put_processor_node_maps(exoid, lb->int_nodes[proc].data(), lb->bor_nodes[proc].data(),
                                     lb->ext_nodes[proc].data(), proc) < 0) {
        Gen_Error(0, "fatal: failed to output node map");
        return 0;
      }

      /* Output the elemental map */
      if (ex_put_processor_elem_maps(exoid, lb->int_elems[proc].data(), nullptr, proc) < 0) {
        Gen_Error(0, "fatal: failed to output element map");
        return 0;
      }

      /*
       * Reorder the nodal communication maps so that they are ordered
       * by processor and then by global ID.
       */

      /* This is a 2-key sort */
      qsort2(lb->ext_procs[proc].data(), lb->ext_nodes[proc].data(), lb->ext_nodes[proc].size());

      /* Output the nodal communication map */
      if (ex_put_node_cmap(exoid, 1, lb->ext_nodes[proc].data(), lb->ext_procs[proc].data(), proc) <
          0) {
        Gen_Error(0, "fatal: failed to output nodal communication map");
        return 0;
      }

    } /* End "for(proc=0; proc < machine->num_procs; proc++)" */
  }
  else if (problem->type == ELEMENTAL) /* Elemental load balance output */
  {
    std::vector<INT> node_proc_ptr(machine->num_procs + 1);
    std::vector<INT> node_cmap_ids_cc(machine->num_procs);
    std::vector<INT> node_cmap_cnts_cc(machine->num_procs);

    node_proc_ptr[0] = 0;
    for (int proc = 0; proc < machine->num_procs; proc++) {
      node_proc_ptr[proc + 1] = node_proc_ptr[proc] + 1;

      node_cmap_cnts_cc[proc] = 0;
      for (size_t cnt = 0; cnt < lb->bor_nodes[proc].size(); cnt++) {
        node_cmap_cnts_cc[proc] += lb->born_procs[proc][cnt].size();
      }

      node_cmap_ids_cc[proc] = 1;
    }

    std::vector<INT> elem_proc_ptr(machine->num_procs + 1);
    std::vector<INT> elem_cmap_ids_cc(machine->num_procs);
    std::vector<INT> elem_cmap_cnts_cc(machine->num_procs);

    elem_proc_ptr[0] = 0;
    for (int proc = 0; proc < machine->num_procs; proc++) {
      elem_proc_ptr[proc + 1] = elem_proc_ptr[proc] + 1;
      elem_cmap_cnts_cc[proc] = lb->e_cmap_elems[proc].size();
      elem_cmap_ids_cc[proc]  = 1;
    }

    /* Output the communication map parameters */
    if (ex_put_cmap_params_cc(exoid, node_cmap_ids_cc.data(), node_cmap_cnts_cc.data(),
                              node_proc_ptr.data(), elem_cmap_ids_cc.data(),
                              elem_cmap_cnts_cc.data(), elem_proc_ptr.data()) < 0) {
      Gen_Error(0, "fatal: unable to output communication map parameters");
      return 0;
    }

    /* Output the node and element maps */
    for (int proc = 0; proc < machine->num_procs; proc++) {
      /* Output the nodal map */
      if (ex_put_processor_node_maps(exoid, lb->int_nodes[proc].data(), lb->bor_nodes[proc].data(),
                                     nullptr, proc) < 0) {
        Gen_Error(0, "fatal: failed to output node map");
        return 0;
      }

      /* Output the elemental map */
      if (ex_put_processor_elem_maps(exoid, lb->int_elems[proc].data(), lb->bor_elems[proc].data(),
                                     proc) < 0) {
        Gen_Error(0, "fatal: failed to output element map");
        return 0;
      }

      /*
       * Build a nodal communication map from the list of border nodes
       * and their associated processors and side IDs.
       */
      size_t nsize = 0;
      for (size_t cnt = 0; cnt < lb->bor_nodes[proc].size(); cnt++) {
        nsize += lb->born_procs[proc][cnt].size();
      }

      if (nsize > 0) {
        std::vector<INT> n_cmap_nodes(nsize);
        std::vector<INT> n_cmap_procs(nsize);

        size_t cnt3 = 0;
        for (size_t cnt = 0; cnt < lb->bor_nodes[proc].size(); cnt++) {
          for (size_t cnt2 = 0; cnt2 < lb->born_procs[proc][cnt].size(); cnt2++) {
            n_cmap_nodes[cnt3]   = lb->bor_nodes[proc][cnt];
            n_cmap_procs[cnt3++] = lb->born_procs[proc][cnt][cnt2];
          }
        }

        /*
         * Reorder the nodal communication maps so that they are ordered
         * by processor and then by global ID.
         */
        /* This is a 2-key sort */
        qsort2(n_cmap_procs.data(), n_cmap_nodes.data(), cnt3);

        /* Output the nodal communication map */
        if (ex_put_node_cmap(exoid, 1, n_cmap_nodes.data(), n_cmap_procs.data(), proc) < 0) {
          Gen_Error(0, "fatal: unable to output nodal communication map");
          return 0;
        }
      } /* End "if (nsize > 0)" */

      /* Output the elemental communication map */
      if (!lb->e_cmap_elems[proc].empty()) {
        if (ex_put_elem_cmap(exoid, 1, lb->e_cmap_elems[proc].data(), lb->e_cmap_sides[proc].data(),
                             lb->e_cmap_procs[proc].data(), proc) < 0) {
          Gen_Error(0, "fatal: unable to output elemental communication map");
          return 0;
        }
      }

    } /* End "for(proc=0; proc < machine->num_procs; proc++)" */
  }
  return 1;
} /*------------------------End write_nemesis()------------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function outputs an ExodusII file for the purposes of visualizing
 * the load balance.
 *****************************************************************************/
template int write_vis(std::string &nemI_out_file, std::string &exoII_inp_file,
                       Machine_Description *machine, Problem_Description *prob,
                       Mesh_Description<int> *mesh, LB_Description<int> *lb);
template int write_vis(std::string &nemI_out_file, std::string &exoII_inp_file,
                       Machine_Description *machine, Problem_Description *prob,
                       Mesh_Description<int64_t> *mesh, LB_Description<int64_t> *lb);

template <typename INT>
int write_vis(std::string &nemI_out_file, std::string &exoII_inp_file, Machine_Description *machine,
              Problem_Description *prob, Mesh_Description<INT> *mesh, LB_Description<INT> *lb)
{
  int exid_vis;
  int exid_inp;

  std::string title;
  const char *coord_names[] = {"X", "Y", "Z"};

  /*-----------------------------Execution Begins------------------------------*/

  /* Generate the file name for the visualization file */
  std::string vis_file_name = remove_extension(nemI_out_file);
  vis_file_name += "-vis.exoII";

  /* Generate the title for the file */
  title = fmt::format("{} {}  load balance visualization file", UTIL_NAME, ELB_VERSION);

  /*
   * If the vis technique is to be by element block then calculate the
   * number of element blocks.
   */
  int vis_nelem_blks;
  if (prob->type == ELEMENTAL) {
    vis_nelem_blks = machine->num_procs;
  }
  else {
    vis_nelem_blks = machine->num_procs + 1;
  }

  /* Create the ExodusII file */
  fmt::print("Outputting load balance visualization file {}\n", vis_file_name.c_str());
  int cpu_ws = 0;
  int io_ws  = 0;
  int mode   = EX_CLOBBER;
  if (prob->int64db | prob->int64api) {
    mode |= EX_NETCDF4 | EX_NOCLASSIC | prob->int64db | prob->int64api;
  }
  if ((exid_vis = ex_create(vis_file_name.c_str(), mode, &cpu_ws, &io_ws)) < 0) {
    Gen_Error(0, "fatal: unable to create visualization output file");
    return 0;
  }
  ON_BLOCK_EXIT(ex_close, exid_vis);

  /*
   * Open the original input ExodusII file, read the values for the
   * element blocks and output them to the visualization file.
   */
  int   icpu_ws = 0;
  int   iio_ws  = 0;
  float vers    = 0.0;
  mode          = EX_READ | prob->int64api;
  if ((exid_inp = ex_open(exoII_inp_file.c_str(), mode, &icpu_ws, &iio_ws, &vers)) < 0) {
    Gen_Error(0, "fatal: unable to open input ExodusII file");
    return 0;
  }
  ON_BLOCK_EXIT(ex_close, exid_inp);

  std::vector<std::string> elem_type(mesh->num_el_blks);
  std::vector<INT>         el_blk_ids(mesh->num_el_blks);
  std::vector<INT>         el_cnt_blk(mesh->num_el_blks);
  std::vector<INT>         node_pel_blk(mesh->num_el_blks);
  std::vector<INT>         nattr_el_blk(mesh->num_el_blks);

  if (ex_get_ids(exid_inp, EX_ELEM_BLOCK, el_blk_ids.data()) < 0) {
    Gen_Error(0, "fatal: unable to get element block IDs");
    return 0;
  }

  int acc_vis = ELB_TRUE; // Output a different element block per processor
  if (prob->vis_out == 2) {
    acc_vis = ELB_FALSE; // Output a nodal/element variable showing processor
  }

  size_t nsize = 0;

  /*
   * Find out if the mesh consists of mixed elements. If not then
   * element blocks will be used to visualize the partitioning. Otherwise
   * nodal/element results will be used.
   */
  for (size_t ecnt = 0; ecnt < mesh->num_el_blks; ecnt++) {
    char type[MAX_STR_LENGTH + 1];
    if (ex_get_block(exid_inp, EX_ELEM_BLOCK, el_blk_ids[ecnt], type, &el_cnt_blk[ecnt],
                     &node_pel_blk[ecnt], nullptr, nullptr, &nattr_el_blk[ecnt]) < 0) {
      Gen_Error(0, "fatal: unable to get element block parameters");
      return 0;
    }
    elem_type[ecnt] = type;
    nsize += el_cnt_blk[ecnt] * node_pel_blk[ecnt];

    if (elem_type[0] == elem_type[ecnt]) {
      if (node_pel_blk[0] != node_pel_blk[ecnt]) {
        acc_vis = ELB_FALSE;
      }
    }
    else {
      acc_vis = ELB_FALSE;
    }
  }

  if (acc_vis == ELB_TRUE) {
    /* Output the initial information */
    if (ex_put_init(exid_vis, title.c_str(), mesh->num_dims, mesh->num_nodes, mesh->num_elems,
                    vis_nelem_blks, 0, 0) < 0) {
      Gen_Error(0, "fatal: unable to output initial params to vis file");
      return 0;
    }

    /* Output the nodal coordinates */
    float *xptr = nullptr;
    float *yptr = nullptr;
    float *zptr = nullptr;
    switch (mesh->num_dims) {
    case 3: zptr = (mesh->coords) + 2 * mesh->num_nodes; FALL_THROUGH;
    case 2: yptr = (mesh->coords) + mesh->num_nodes; FALL_THROUGH;
    case 1: xptr = mesh->coords;
    }
    if (ex_put_coord(exid_vis, xptr, yptr, zptr) < 0) {
      Gen_Error(0, "fatal: unable to output coords to vis file");
      return 0;
    }

    if (ex_put_coord_names(exid_vis, (char **)coord_names) < 0) {
      Gen_Error(0, "fatal: unable to output coordinate names");
      return 0;
    }

    std::vector<INT> elem_block(mesh->num_elems);
    std::vector<INT> elem_map(mesh->num_elems);
    std::vector<INT> tmp_connect(nsize);
    for (size_t ecnt = 0; ecnt < mesh->num_elems; ecnt++) {
      elem_map[ecnt] = ecnt + 1;
      if (prob->type == ELEMENTAL) {
        elem_block[ecnt] = lb->vertex2proc[ecnt];
      }
      else {
        int proc         = lb->vertex2proc[mesh->connect[ecnt][0]];
        int nnodes       = get_elem_info(NNODES, mesh->elem_type[ecnt]);
        elem_block[ecnt] = proc;
        for (int ncnt = 1; ncnt < nnodes; ncnt++) {
          if (lb->vertex2proc[mesh->connect[ecnt][ncnt]] != proc) {
            elem_block[ecnt] = machine->num_procs;
            break;
          }
        }
      }
    }

    int              ccnt = 0;
    std::vector<INT> vis_el_blk_ptr(vis_nelem_blks + 1);
    for (INT bcnt = 0; bcnt < vis_nelem_blks; bcnt++) {
      vis_el_blk_ptr[bcnt] = ccnt;
      int    pos           = 0;
      int    old_pos       = 0;
      INT *  el_ptr        = elem_block.data();
      size_t ecnt          = mesh->num_elems;
      while (pos != -1) {
        pos = in_list(bcnt, ecnt, el_ptr);
        if (pos != -1) {
          old_pos += pos + 1;
          ecnt       = mesh->num_elems - old_pos;
          el_ptr     = elem_block.data() + old_pos;
          int nnodes = get_elem_info(NNODES, mesh->elem_type[old_pos - 1]);
          for (int ncnt = 0; ncnt < nnodes; ncnt++) {
            tmp_connect[ccnt++] = mesh->connect[old_pos - 1][ncnt] + 1;
          }
        }
      }
    }
    vis_el_blk_ptr[vis_nelem_blks] = ccnt;

    /* Output the element map */
    if (ex_put_map(exid_vis, elem_map.data()) < 0) {
      Gen_Error(0, "fatal: unable to output element number map");
      return 0;
    }

    /* Output the visualization element blocks */
    for (int bcnt = 0; bcnt < vis_nelem_blks; bcnt++) {
      /*
       * Note this assumes all the blocks contain the same type
       * element.
       */
      int ecnt = (vis_el_blk_ptr[bcnt + 1] - vis_el_blk_ptr[bcnt]) / node_pel_blk[0];
      if (ex_put_block(exid_vis, EX_ELEM_BLOCK, bcnt + 1, elem_type[0].c_str(), ecnt,
                       node_pel_blk[0], 0, 0, 0) < 0) {
        Gen_Error(0, "fatal: unable to output element block params");
        return 0;
      }

      /* Output the connectivity */
      if (ex_put_conn(exid_vis, EX_ELEM_BLOCK, bcnt + 1, &tmp_connect[vis_el_blk_ptr[bcnt]],
                      nullptr, nullptr) < 0) {
        Gen_Error(0, "fatal: unable to output element connectivity");
        return 0;
      }
    }
  }

  else { /* For nodal/element results visualization of the partitioning. */
    // Copy the mesh portion to the vis file.
    ex_copy(exid_inp, exid_vis);

    /* Set up the file for nodal/element results */
    float time_val = 0.0;
    if (ex_put_time(exid_vis, 1, &time_val) < 0) {
      Gen_Error(0, "fatal: unable to output time to vis file");
      return 0;
    }

    const char *var_names[] = {"proc"};
    if (prob->type == NODAL) {
      /* Allocate memory for the nodal values */
      std::vector<float> proc_vals(mesh->num_nodes);

      if (ex_put_variable_param(exid_vis, EX_NODAL, 1) < 0) {
        Gen_Error(0, "fatal: unable to output var params to vis file");
        return 0;
      }

      if (ex_put_variable_names(exid_vis, EX_NODAL, 1, (char **)var_names) < 0) {
        Gen_Error(0, "fatal: unable to output variable name");
        return 0;
      }

      /* Do some problem specific assignment */
      for (size_t ncnt = 0; ncnt < mesh->num_nodes; ncnt++) {
        proc_vals[ncnt] = lb->vertex2proc[ncnt];
      }

      for (int pcnt = 0; pcnt < machine->num_procs; pcnt++) {
        for (auto &elem : lb->bor_nodes[pcnt]) {
          proc_vals[elem] = machine->num_procs + 1;
        }
      }

      /* Output the nodal variables */
      if (ex_put_var(exid_vis, 1, EX_NODAL, 1, 1, mesh->num_nodes, proc_vals.data()) < 0) {
        Gen_Error(0, "fatal: unable to output nodal variables");
        return 0;
      }
    }
    else if (prob->type == ELEMENTAL) {
      /* Allocate memory for the element values */
      std::vector<float> proc_vals(mesh->num_elems);

      if (ex_put_variable_param(exid_vis, EX_ELEM_BLOCK, 1) < 0) {
        Gen_Error(0, "fatal: unable to output var params to vis file");
        return 0;
      }

      if (ex_put_variable_names(exid_vis, EX_ELEM_BLOCK, 1, (char **)var_names) < 0) {
        Gen_Error(0, "fatal: unable to output variable name");
        return 0;
      }

      /* Do some problem specific assignment */
      for (int proc = 0; proc < machine->num_procs; proc++) {
        for (size_t e = 0; e < lb->int_elems[proc].size(); e++) {
          size_t ecnt     = lb->int_elems[proc][e];
          proc_vals[ecnt] = proc;
        }

        for (size_t e = 0; e < lb->bor_elems[proc].size(); e++) {
          size_t ecnt     = lb->bor_elems[proc][e];
          proc_vals[ecnt] = proc;
        }
      }

      /* Output the element variables */
      size_t offset = 0;
      for (size_t i = 0; i < mesh->num_el_blks; i++) {
        if (ex_put_var(exid_vis, 1, EX_ELEM_BLOCK, 1, el_blk_ids[i], el_cnt_blk[i],
                       &proc_vals[offset]) < 0) {
          Gen_Error(0, "fatal: unable to output nodal variables");
          return 0;
        }
        offset += el_cnt_blk[i];
      }
    }
  }
  return 1;
} /*---------------------------End write_vis()-------------------------------*/
