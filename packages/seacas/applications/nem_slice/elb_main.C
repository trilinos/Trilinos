/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Functions contained in this file:
 *      main()
 *      print_input()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <cstddef> // for size_t
#include <cstdlib> // for free, exit, malloc
#include <cstring> // for strcmp
#include <iostream>
#include <stdexcept>

#include "add_to_log.h"  // for add_to_log
#include "elb.h"         // for LB_Description<INT>, get_time, etc
#include "elb_allo.h"    // for array_alloc
#include "elb_elem.h"    // for E_Type, ::NULL_EL
#include "elb_err.h"     // for error_report, Gen_Error, etc
#include "elb_exo.h"     // for init_weight_struct, etc
#include "elb_graph.h"   // for generate_graph
#include "elb_inp.h"     // for check_inp_specs, etc
#include "elb_loadbal.h" // for generate_loadbal, etc
#include "elb_output.h"  // for write_nemesis, write_vis
#include "fmt/ostream.h"

#ifdef USE_ZOLTAN
#include <mpi.h> // for MPI_Finalize, etc
#endif

#ifdef SGI10K
#include <sys/resource.h>
#endif

namespace {
  template <typename INT>
  void print_input(Machine_Description * /*machine*/, LB_Description<INT> * /*lb*/,
                   Problem_Description * /*prob*/, Solver_Description * /*solver*/,
                   Weight_Description<INT> * /*weight*/);
} // namespace

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Routine which reads in a EXODUS II mesh database and load balances the
 * mesh by either a nodal or element assignment.
 *----------------------------------------------------------------------------
 * Functions called:
 *      cmd_line_arg_parse()
 *      read_cmd_file()
 *      check_inp_specs()
 *      print_input()
 *      read_exo_weights()
 *      read_mesh_params()
 *      read_mesh()
 *      generate_gaph()
 *      generate_loadbal()
 *      write_vis()
 *      generate_maps()
 *      write_nemesis()
 *----------------------------------------------------------------------------
 * Variable Index:
 *      machine       - structure containing information about the machine
 *                      for which the load balance is to be generated.
 *      lb            - structure containing information about what type
 *                      of load balance is to be performed.
 *      problem       - structure containing information about the problem
 *                      type. Currently whether the problem is a nodal or
 *                      elemental based decomposition.
 *      solver        - structure containing various parameters for the
 *                      eigensolver used by Chaco.
 *      weight        - structure used to store parameters about how to
 *                      weight the resulting graph before it is fed into
 *                      Chaco.
 *      graph         - structure containing the graph of the problem.
 *      mesh          - structure containing a description of the FEM mesh.
 *      exoII_inp_file - name of the ExodusII file containing the problem
 *                       geometry.
 *      ascii_inp_file - name of the input command file.
 *      nemI_out_file  - name of the output NemesisI file.
 *
 *****************************************************************************/
template <typename INT> int internal_main(int argc, char *argv[], INT /* dummy */);

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
    status = internal_main(argc, argv, int(0));
  }
  else if ((int64db != 0) || int64com) {
    fmt::print(stderr,
               "Using 64-bit integer mode for decomposition...\n"
               "NOTE: Only 'linear' and 'scattered' methods are supported for 64-bit models\n");

    status = internal_main(argc, argv, int64_t(0));
  }
  else {
    fmt::print(stderr, "Using 32-bit integer mode for decomposition...\n");
    status = internal_main(argc, argv, int(0));
  }

  /* Report any non-fatal errors that may have occurred */
  error_report();

  /* Get ending time */
  double end_time = get_time();
  fmt::print(stderr, "The entire load balance took {} seconds.\n", end_time - start_time);
  add_to_log(argv[0], end_time - start_time);
  return status;
}

template <typename INT> int internal_main(int argc, char *argv[], INT /* dummy */)
{
  /* Local variables */
  std::string exoII_inp_file;
  std::string ascii_inp_file;
  std::string nemI_out_file;

  Machine_Description     machine;
  LB_Description<INT>     lb;
  Problem_Description     problem;
  Solver_Description      solver;
  Weight_Description<INT> weight;
  Mesh_Description<INT>   mesh;
  Sphere_Info             sphere;
  Graph_Description<INT>  graph;

/*-----------------------------Execution Begins------------------------------*/
#ifdef USE_ZOLTAN
  MPI_Init(&argc, &argv);
#endif

  mesh.title[0] = '\0';

  if (sizeof(INT) == 8) {
    problem.int64api = EX_ALL_INT64_API;
  }

  /* Parse the command line */
  if (!cmd_line_arg_parse(argc, argv, exoII_inp_file, ascii_inp_file, nemI_out_file, &machine, &lb,
                          &problem, &solver, &weight)) {
    fmt::print(stderr, "error parsing command line\n");
    error_report();
    exit(1);
  }

  /*
   *If the information is to be taken from an ASCII input file then
   * read that file.
   */
  if (!ascii_inp_file.empty()) {
    if (!read_cmd_file(ascii_inp_file, exoII_inp_file, nemI_out_file, &machine, &lb, &problem,
                       &solver, &weight)) {
      fmt::print(stderr, "error parsing command file\n");
      error_report();
      exit(1);
    }
  }

  /* make sure that this type is set */
  if (weight.type < 0) {
    weight.type = NO_WEIGHT;
  }

  /*
   * Perform at least some rudimentary error checks on the user
   * specified input.
   */
  if (!check_inp_specs(exoII_inp_file, nemI_out_file, &machine, &lb, &problem, &solver, &weight)) {
    fmt::print(stderr, "Error in user specified parameters\n");
    error_report();
    exit(1);
  }

  /* Output the parameters for the run to the screen */
  print_input(&machine, &lb, &problem, &solver, &weight);

  /* Read in mesh parameters */
  double time1 = get_time();
  if (!read_mesh_params(exoII_inp_file, &problem, &mesh, &sphere)) {
    fmt::print(stderr, "Error reading mesh parameters\n");
    error_report();
    exit(1);
  }
  double time2 = get_time();
  fmt::print("Time to read mesh parameters: {}s\n", time2 - time1);

  /* Check for error conditions */
  if ((mesh.num_nodes) / (machine.num_procs) < 1) {
    Gen_Error(0, "fatal: problem divided among too many processors");
    error_report();
    exit(1);
  }
  else if ((problem.type == ELEMENTAL) && ((mesh.num_elems) / (machine.num_procs) < 1)) {
    Gen_Error(0, "fatal: problem divided among too many processors");
    error_report();
    exit(1);
  }

  /* If vertex weighting is turned on, prepare weight struct */
  if ((weight.type & READ_EXO) || (weight.type & EL_BLK)) {
    time1 = get_time();
    if (!init_weight_struct(&problem, &mesh, &weight)) {
      fmt::print(stderr, "Error during initialization of weight struct\n");
      error_report();
      exit(1);
    }
    time2 = get_time();
    fmt::print("Time to initialize weight struct: {}s\n", time2 - time1);
  }

  /* If desired, read in the weighting factors from the ExodusII file */
  if (weight.type & READ_EXO) {
    time1 = get_time();
    if (!read_exo_weights(&problem, &weight)) {
      fmt::print(stderr, "Error during read of ExodusII weights\n");
      error_report();
      exit(1);
    }
    time2 = get_time();
    fmt::print("Time to read ExodusII weights: {}s\n", time2 - time1);
  }

  /* Initialize various parameters */
  if (lb.type == INERTIAL || lb.type == ZPINCH || lb.type == BRICK || lb.type == ZOLTAN_RCB ||
      lb.type == ZOLTAN_RIB || lb.type == ZOLTAN_HSFC || problem.vis_out == 1 ||
      problem.vis_out == 2) {
    problem.read_coords = ELB_TRUE;
  }
  else {
    problem.read_coords = ELB_FALSE;
  }

  if (lb.type != SPECTRAL) {
    problem.coarse_flag = ELB_FALSE;
  }
  else {
    problem.coarse_flag = ELB_TRUE;
  }

  if (lb.refine == KL_REFINE) {
    problem.alloc_graph = ELB_TRUE;
  }
  else if (lb.type == SPECTRAL) {
    problem.alloc_graph = ELB_TRUE;
  }
  else {
    problem.alloc_graph = ELB_FALSE;
  }

  /* if fix_columns is on, we need the face adjacency graph. So if
   * nothing else has asked for the full adjacency graph, ask for the
   * face adjacency graph. If something else did ask for the adjacency
   * graph, we don't know if its full or face adjacency only, so leave
   * the option as is */

  if (problem.fix_columns) {
    if (problem.alloc_graph == ELB_FALSE) {
      problem.alloc_graph = ELB_TRUE;
      problem.face_adj    = ELB_TRUE;
    }
  }

  /* Allocate necessary memory */
  if (problem.type == NODAL) {
    problem.num_vertices = mesh.num_nodes;
  }
  else if (problem.type == ELEMENTAL) {
    problem.num_vertices = (mesh.num_elems - sphere.num);
  }

  if (problem.read_coords == ELB_TRUE) {
    size_t mem_req = (size_t)(mesh.num_dims) * (mesh.num_nodes) * sizeof(float);
    mesh.coords    = reinterpret_cast<float *>(malloc(mem_req));
    if (!(mesh.coords)) {
      Gen_Error(0, "fatal: insufficient memory for coordinates");
      error_report();
      exit(1);
    }
  }
  else {
    mesh.coords = nullptr;
  }

  mesh.elem_type = static_cast<E_Type *>(array_alloc(1, mesh.num_elems, sizeof(E_Type)));
  mesh.connect = static_cast<INT **>(array_alloc(2, mesh.num_elems, mesh.max_np_elem, sizeof(INT)));
  if (!(mesh.elem_type) || !(mesh.connect)) {
    Gen_Error(0, "fatal: insufficient memory");
    error_report();
    exit(1);
  }

  /* Read the mesh */
  time1 = get_time();
  if (!read_mesh(exoII_inp_file, &problem, &mesh, &weight)) {
    Gen_Error(0, "fatal: could not read mesh");
    error_report();
    exit(1);
  }
  time2 = get_time();
  fmt::print("Time to read mesh: {}s\n", time2 - time1);

  /* free unneeded memory */
  vec_free(weight.ow);
  vec_free(weight.elemblk);
  vec_free(weight.elemblk_wgt);

  /* Generate the graph for the mesh */
  time1 = get_time();
  if (!generate_graph(&problem, &mesh, &graph, &weight, &sphere)) {
    Gen_Error(0, "fatal: could not generate graph");
    error_report();
    exit(1);
  }
  time2 = get_time();
  fmt::print("Time to generate graph: {}s\n", time2 - time1);

  /* Generate load balance */
  try {
    time1 = get_time();
    if (!generate_loadbal(&machine, &problem, &mesh, &lb, &solver, &graph, &weight, &sphere, argc,
                          argv)) {
      Gen_Error(0, "fatal: could not generate load balance");
      error_report();
      exit(1);
    }
    time2 = get_time();
    fmt::print("Time to generate load balance: {}s\n", time2 - time1);
  }
  catch (const std::exception &e) {
    fmt::print(stderr, "NEM_SLICE: Exception in generate_loadbal: {}\n", e.what());
  }

  /* free up memory */
  if (sphere.adjust) {
    free(sphere.adjust);
  }

#ifdef PRINT_VERT
  for (size_t cnt = 0; cnt < problem.num_vertices; cnt++)
    fmt::print("element = {}, proc = {}\n", cnt, lb.vertex2proc[cnt]);
#endif

  /*
   * NOTE: in Chaco, if FREE_GRAPH is set to 1, the following arrays
   * are freed: graph.start
   *            graph.adj
   *            weight.vertices
   *            weight.edges
   *
   * Need to take into account special case where the mesh contains
   * only spheres. In this case Chaco is not called, and the above
   * arrays are not freed
   */

  if (sphere.num >= mesh.num_elems) {
    vec_free(graph.start);
    vec_free(graph.adj);
    vec_free(weight.vertices);
    vec_free(weight.edges);
  }

  /* Generate the load balance maps */
  time1 = get_time();
  if (!generate_maps(&machine, &problem, &mesh, &lb, &graph)) {
    Gen_Error(0, "fatal: could not generate load-balance maps");
    error_report();
    exit(1);
  }
  time2 = get_time();
  fmt::print("Time to generate load-balance maps: {}s\n", time2 - time1);

  /* Output the visualization file */
  if (problem.vis_out == 1 || problem.vis_out == 2) {
    time1 = get_time();
    if (!write_vis(nemI_out_file, exoII_inp_file, &machine, &problem, &mesh, &lb)) {
      Gen_Error(0, "warning: unable to output visualization file");
    }
    time2 = get_time();
    fmt::print("Time to output visualization file: {}s\n", time2 - time1);
  }

  /* Free up memory no longer needed */
  if (problem.read_coords == ELB_TRUE) {
    free(mesh.coords);
  }

  free(mesh.elem_type);
  free(mesh.connect);

  if (!graph.sur_elem.empty()) {
    for (size_t cnt = 0; cnt < mesh.num_nodes; cnt++) {
      vec_free(graph.sur_elem[cnt]);
    }
    vec_free(graph.sur_elem);
  }

  /* Output a Nemesis load balance file */
  time1 = get_time();
  if (!write_nemesis(nemI_out_file, &machine, &problem, &mesh, &lb, &sphere)) {
    Gen_Error(0, "fatal: could not output Nemesis file");
    error_report();
    exit(1);
  }
  time2 = get_time() - time1;
  fmt::print("Time to write Nemesis file: {}s\n", time2);

  /* Free up unused memory for leak checking */
  for (int cnt = 0; cnt < machine.num_procs; cnt++) {
    vec_free(lb.int_nodes[cnt]);
    vec_free(lb.bor_nodes[cnt]);
    if (problem.type == NODAL) {
      vec_free(lb.ext_nodes[cnt]);
      vec_free(lb.ext_procs[cnt]);
    }

    vec_free(lb.int_elems[cnt]);
    if (problem.type == ELEMENTAL) {
      vec_free(lb.bor_elems[cnt]);
    }
  }

  vec_free(lb.int_nodes);

  if (problem.type == NODAL) {
    vec_free(lb.ext_nodes);
    vec_free(lb.ext_procs);
  }

  vec_free(lb.int_elems);

  if (problem.type == ELEMENTAL) {
    vec_free(lb.bor_elems);
    for (int cnt = 0; cnt < machine.num_procs; cnt++) {
      for (size_t cnt1 = 0; cnt1 < lb.bor_nodes[cnt].size(); cnt1++) {
        vec_free(lb.born_procs[cnt][cnt1]);
      }

      if (!lb.born_procs[cnt].empty()) {
        vec_free(lb.born_procs[cnt]);
      }
      vec_free(lb.ext_procs[cnt]);
      vec_free(lb.e_cmap_elems[cnt]);
      vec_free(lb.e_cmap_sides[cnt]);
      vec_free(lb.e_cmap_procs[cnt]);
      vec_free(lb.e_cmap_neigh[cnt]);
    }
    vec_free(lb.e_cmap_elems);
    vec_free(lb.born_procs);
    vec_free(lb.ext_procs);
  }

  vec_free(lb.bor_nodes);
  free(lb.vertex2proc);

#ifdef USE_ZOLTAN
  MPI_Finalize();
#endif

  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function prints out parameters as read from the command line and/or
 * the input ASCII load-balance file.
 *---------------------------------------------------------------------------*/
namespace {
  template <typename INT>
  void print_input(Machine_Description *machine, LB_Description<INT> *lb, Problem_Description *prob,
                   Solver_Description *solver, Weight_Description<INT> *weight)
  {
    fmt::print("{} version {}\n", UTIL_NAME, ELB_VERSION);

    fmt::print("Performing ");
    switch (prob->type) {
    case NODAL: fmt::print("a nodal "); break;

    case ELEMENTAL: fmt::print("an elemental "); break;
    }

    fmt::print("load balance with the following parameters...\n");

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                         Machine_Description PARAMETERS                                */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    fmt::print("Machine Description\n");
    if (machine->num_boxes > 1) {
      fmt::print("\tCluster of {} boxes\n", machine->num_boxes);
      fmt::print("\tarchitecture of each box: ");
    }
    else {
      fmt::print("\tarchitecture: ");
    }
    switch (machine->type) {
    case HCUBE: fmt::print("hypercube\n"); break;

    case MESH: fmt::print("mesh\n"); break;
    }
    if (machine->num_boxes > 1) {
      fmt::print("\tdimension(s) of each box: ");
    }
    else {
      fmt::print("\tdimension(s): ");
    }
    switch (machine->type) {
    case HCUBE: fmt::print("{}\n", machine->dim[0]); break;

    case MESH:
      for (int cnt = 0; cnt < (machine->num_dims) - 1; cnt++) {
        fmt::print("{}x", machine->dim[cnt]);
      }

      fmt::print("{}\n", machine->dim[(machine->num_dims) - 1]);
      break;
    }
    fmt::print("\ttotal number of processors: {}\n", machine->num_procs);

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                       LOAD BALANCE PARAMETERS                             */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    fmt::print("Load Balance Parameters\n");
    switch (lb->type) {
    case MULTIKL:
      fmt::print("\ttype: multilevel\n");
      fmt::print("\tnumber of sections: {}\n", lb->num_sects);
      break;

    case SPECTRAL:
      fmt::print("\ttype: spectral\n");
      fmt::print("\tnumber of sections: {}\n", lb->num_sects);
      break;

    case INERTIAL: fmt::print("\ttype: inertial\n"); break;

    case ZPINCH: fmt::print("\ttype: zpinch\n"); break;

    case BRICK: fmt::print("\ttype: brick\n"); break;

    case ZOLTAN_RCB: fmt::print("\ttype: rcb\n"); break;

    case ZOLTAN_RIB: fmt::print("\ttype: rib\n"); break;

    case ZOLTAN_HSFC: fmt::print("\ttype: hsfc\n"); break;

    case LINEAR: fmt::print("\ttype: linear\n"); break;

    case RANDOM: fmt::print("\ttype: random\n"); break;

    case SCATTERED: fmt::print("\ttype: scattered\n"); break;

    case INFILE:
      fmt::print("\ttype: input from file\n");
      fmt::print("\tfile name: {}\n", lb->file.c_str());
      break;
    }
    fmt::print("\trefinement: ");
    switch (lb->refine) {
    case KL_REFINE: fmt::print("Kernighan-Lin\n"); break;

    case NO_REFINE: fmt::print("none\n"); break;
    }
    if (lb->cnctd_dom) {
      fmt::print("\tConnected Domain enforced\n");
    }
    if (lb->outfile) {
      fmt::print("\toutput assignment vector\n");
      fmt::print("\tfile name: {}\n", lb->file.c_str());
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                         EIGENSOLVER PARAMETERS                            */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    if (lb->type == MULTIKL || lb->type == SPECTRAL) {
      fmt::print("Eigensolver Parameters\n");
      fmt::print("\teignsolver tolerance: {}\n", solver->tolerance);
      if (solver->rqi_flag == USE_RQI) {
        fmt::print("\tusing RQI/Symmlq eigensolver\n");
        fmt::print("\tnumber of vertices to coarsen down to: {}\n", solver->vmax);
      }
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                          WEIGHTING PARAMETERS                             */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    fmt::print("Weighting Parameters\n");
    if (weight->type == NO_WEIGHT) {
      fmt::print("\tno weighting\n");
    }
    else if (weight->type & READ_EXO) {
      fmt::print("\tweights from: ExodusII file\n");
      fmt::print("\ttime index: {}\n", weight->exo_tindx);
      fmt::print("\tvariable index: {}\n", weight->exo_vindx);
      fmt::print("\tvariable name: {}\n", weight->exo_varname.c_str());
    }
    else if (weight->type & EL_BLK) {
      fmt::print("\tElement Block weights specified\n");
      for (size_t cnt = 0; cnt < weight->elemblk.size(); cnt++) {
        fmt::print("\telement block: {}, weight: {}\n", weight->elemblk[cnt],
                   weight->elemblk_wgt[cnt]);
      }
    }
    else if (weight->type & EDGE_WGT) {
      fmt::print("\tedge weights turned on\n");
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                          WEIGHTING PARAMETERS                             */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    fmt::print("Miscellaneous Options\n");
    if (prob->face_adj == 1) {
      fmt::print("\tusing face definition of adjacency\n");
    }
    if (prob->global_mech == 1) {
      fmt::print("\tlooking for global mechanisms\n");
    }
    if (prob->local_mech == 1) {
      fmt::print("\tlooking for local mechanisms\n");
    }
    if (prob->find_cnt_domains == 1) {
      fmt::print("\tidentifying the number of disconnected element blocks on a subdomain\n");
    }
    if (prob->mech_add_procs == 1) {
      fmt::print("\tincreasing the number of processors if mechanisms are found\n");
    }
    if (prob->dsd_add_procs == 1) {
      fmt::print("\tincreasing the number of processors if disconnected sudomains are found\n");
    }
    if (prob->no_sph == 1) {
      fmt::print("\tSPHERES are being treated as concentrated mass - connectivity exists\n");
    }
    if (prob->skip_checks == 1) {
      fmt::print("\tWARNING: side id error checks turned off\n");
    }
    if (prob->groups != nullptr) {
      fmt::print("\telement block groups defined\n");
      fmt::print("\tgroup string: \"{}\"\n", prob->groups);
    }
  }
} // namespace
