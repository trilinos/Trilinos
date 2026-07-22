// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>   // must appear before stdio or iostream

#ifdef TFLOP
#include <iostream.h>
#include <ostream.h>
#else
#include <iostream>
#include <ostream>
using namespace std;
#endif

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_loadbal_const.h"
#include "dr_output_const.h"
#include "dr_err_const.h"
#include "dr_elem_util_const.h"
#include "dr_dd.h"
#include "dr_compress_const.h"

#include "zoltan_cpp.h"

int Debug_Driver = 1;
int Number_Iterations = 1;
int Driver_Action = 1;	/* Flag indicating coloring, load-balancing or ordering. */
int Debug_Chaco_Input = 0;
int Chaco_In_Assign_Inv = 0;
struct Test_Flags Test;
struct Output_Flags Output;

double Total_Partition_Time = 0.0;  /* Total over Number_Iterations */

static int read_mesh(int, int, PROB_INFO_PTR, PARIO_INFO_PTR, MESH_INFO_PTR);
static void initialize_mesh(MESH_INFO_PTR);
static void print_input_info(ostream &os, int Num_Proc, PROB_INFO_PTR prob, 
PARIO_INFO_PTR pio, float zoltan_version);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int main(int argc, char *argv[])
{
/* Local declarations. */
  const char  *cmd_file;
  char   cmesg[256]; /* for error messages */
  int    error, gerror;
  int    print_output = 1;

/***************************** BEGIN EXECUTION ******************************/

  /* Initialize MPI */

  // We must use the C bindings to MPI because the C++ bindings are
  // are not available or not complete on some of our platforms.
 
  MPI_Init(&argc, &argv);

  /* get some machine information */
  int Proc = 0, Num_Proc = 0;
  MPI_Comm_rank(zoltan_get_global_comm(), &Proc);
  MPI_Comm_size(zoltan_get_global_comm(), &Num_Proc);

  /* Initialize flags */
  Test.DDirectory = 0;
  Test.Local_Parts = 0;
  Test.Drops = 0;
  Test.RCB_Box = 0;
  Test.Multi_Callbacks = 0;
  Test.Gen_Files = 0;
  Test.Fixed_Objects = 0;
  Test.Null_Lists = NO_NULL_LISTS;

  Output.Text = 1;
  Output.Gnuplot = 0;
  Output.Nemesis = 0;
  Output.Plot_Partition = 0;
  Output.Mesh_Info_File = 0;

  /* Interpret the command line */
  switch(argc)
  {
  case 1:
    cmd_file = "zdrive.inp";
    break;

  case 2:
    cmd_file = argv[1];
    break;

  default:
    cerr << "MAIN: ERROR in command line " ;

    if (Proc == 0)
    {
      cerr << "usage:" << endl;
      cerr << "\t" << DRIVER_NAME << " [command file]";
    }
    cerr << endl;

    exit(1);
    break;
  }

  /* Initialize Zoltan 
   *  Not part of C++ interface at this time. (C++ wishlist)
   */

  float version;
  
  if ((error = Zoltan_Initialize(argc, argv, &version)) != ZOLTAN_OK) {
    sprintf(cmesg, "fatal: Zoltan_Initialize returned error code, %d", error);
    Gen_Error(0, cmesg);
    error_report(Proc);
    MPI_Finalize();
    return 1;
  }

  /*
   *  Create a Zoltan structure.  
   *  No exception handling at this time. (C++ wishlist)
   *  We must dynamically create the object so that we can delete it
   *  before MPI_Finalize().  (If zz is created on the stack, it will
   *  be deleted atexit, after MPI_Finalize().  Zoltan_Destroy calls
   *  MPI functions.)
   */

  Zoltan *zz = NULL;
  zz = new Zoltan();

  /* initialize some variables */
  MESH_INFO  mesh;
  initialize_mesh(&mesh);

  PARIO_INFO pio_info;
  pio_info.dsk_list_cnt		= -1;
  pio_info.num_dsk_ctrlrs	= -1;
  pio_info.pdsk_add_fact	= -1;
  pio_info.zeros		= -1;
  pio_info.file_type		= -1;
  pio_info.init_dist_type	= -1;
  pio_info.init_size		= ZOLTAN_ID_INVALID;
  pio_info.init_dim 		= -1;
  pio_info.init_vwgt_dim 	= -1;
  pio_info.init_dist_pins       = -1;
  pio_info.pdsk_root[0]		= '\0';
  pio_info.pdsk_subdir[0]	= '\0';
  pio_info.pexo_fname[0]	= '\0';
  pio_info.file_comp            = STANDARD;

  PROB_INFO  prob;
  prob.method[0]		= '\0';
  prob.num_params		= 0;
  prob.params			= NULL;

  /* Read in the ascii input file */
  error = 0;
  if (Proc == 0) {
    cout << "\n\nReading the command file, " << cmd_file << endl;
    if (!read_cmd_file(cmd_file, &prob, &pio_info, NULL)) {
      sprintf(cmesg,"fatal: Could not read in the command file"
              " \"%s\"!\n", cmd_file);
      Gen_Error(0, cmesg);
      error_report(Proc);
      print_output = 0;
      error = 1;
    }

    if (!check_inp(&prob, &pio_info)) {
      Gen_Error(0, "fatal: Error in user specified parameters.\n");
      error_report(Proc);
      print_output = 0;
      error = 1;
    }

    print_input_info(cout, Num_Proc, &prob, &pio_info, version);
  }

  MPI_Allreduce(&error, &gerror, 1, MPI_INT, MPI_MAX, zoltan_get_global_comm());

  if (gerror) goto End;

  /* broadcast the command info to all of the processor */
  brdcst_cmd_info(Proc, &prob, &pio_info, &mesh);

  zz->Set_Param("DEBUG_MEMORY", "1");
  print_output = Output.Text;

  if (!setup_zoltan(*zz, Proc, &prob, &mesh, &pio_info)) {
    Gen_Error(0, "fatal: Error returned from setup_zoltan\n");
    error_report(Proc);
    print_output = 0;
    goto End;
  }

  srand(Proc);

  /* Loop over read and balance for a number of iterations */
  /* (Useful for testing REUSE parameters in Zoltan.) */
  for (int iteration = 1; iteration <= Number_Iterations; iteration++) {

    /*
     * now read in the mesh and element information.
     * This is the only function call to do this. Upon return,
     * the mesh struct and the elements array should be filled.
     */
    if (iteration == 1) {
      if (!read_mesh(Proc, Num_Proc, &prob, &pio_info, &mesh)) {
        Gen_Error(0, "fatal: Error returned from read_mesh\n");
        error_report(Proc);
        print_output = 0;
        goto End;
      }
      /* 
       *  Create a Zoltan DD for tracking elements during repartitioning.
       */

      if (mesh.data_type == ZOLTAN_HYPERGRAPH && !build_elem_dd(&mesh)) {
        Gen_Error(0, "fatal: Error returned from build_elem_dd\n");
        error_report(Proc);
        print_output = 0;
        goto End;
      }
    }


#ifdef KDDKDD_COOL_TEST
/* KDD Cool test of changing number of partitions  */
    sprintf(cmesg, "%d", Num_Proc * iteration);
    zz->Set_Param("NUM_GLOBAL_PARTS", cmesg);
#endif

    /*
     * Produce files to verify input.
     */
    if (iteration == 1)
      if (Debug_Driver > 2) {
        if (!output_results(cmd_file,"in",Proc,Num_Proc,&prob,&pio_info,&mesh)){
          Gen_Error(0, "fatal: Error returned from output_results\n");
          error_report(Proc);
        }
        if (Output.Gnuplot)
          if (!output_gnu(cmd_file,"in",Proc,Num_Proc,&prob,&pio_info,&mesh)) {
            Gen_Error(0, "warning: Error returned from output_gnu\n");
            error_report(Proc);
          }
      }

    /*
     * now run Zoltan to get a new load balance and perform
     * the migration
     */
    if (!run_zoltan(*zz, Proc, &prob, &mesh, &pio_info)) {
      Gen_Error(0, "fatal: Error returned from run_zoltan\n");
      error_report(Proc);
      print_output = 0;
      goto End;
    }

    /* Reset the mesh data structure for next iteration. */
    if (iteration < Number_Iterations) {
      float twiddle = 0.01;
      char str[32];
      /* Perturb coordinates of mesh */
      if (mesh.data_type == ZOLTAN_GRAPH)
        for (int i = 0; i < mesh.num_elems; i++) {
          for (int j = 0; j < mesh.num_dims; j++) {
            /* tmp = ((float) rand())/RAND_MAX; *//* Equiv. to sjplimp's test */
            float tmp = (float) (i % 10) / 10.;
            mesh.elements[i].coord[0][j] += twiddle * (2.0*tmp-1.0);
            mesh.elements[i].avg_coord[j] = mesh.elements[i].coord[0][j];
          }
        }
      /* change the ParMETIS Seed */
      sprintf(str, "%d", iteration%10000);
#ifdef ZOLTAN_PARMETIS      
      zz->Set_Param("PARMETIS_SEED", str);
#endif
    }

  } /* End of loop over read and balance */

  if (Proc == 0) {
    cout << "FILE " << cmd_file << ":  Total:  " ;
    cout << Total_Partition_Time << " seconds in Partitioning" << endl;

    cout << "FILE " << cmd_file << ":  Average:  ";
    cout << Total_Partition_Time/Number_Iterations << " seconds per Iteration";
    cout << endl;
  }

End:
  
  if (mesh.data_type == ZOLTAN_HYPERGRAPH)
  {
    destroy_elem_dd();
  }

  delete zz;

  Zoltan_Memory_Stats();

  /*
   * output the results
   */
  if (print_output) {
    if (!output_results(cmd_file,"out",Proc,Num_Proc,&prob,&pio_info,&mesh)) {
      Gen_Error(0, "fatal: Error returned from output_results\n");
      error_report(Proc);
    }

    if (Output.Gnuplot) {
      if (!output_gnu(cmd_file,"out",Proc,Num_Proc,&prob,&pio_info,&mesh)) {
        Gen_Error(0, "warning: Error returned from output_gnu\n");
        error_report(Proc);
      }
    }
  }

  free_mesh_arrays(&mesh);
  if (prob.params != NULL) free(prob.params);

  MPI_Finalize();

  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function determines which input file type is being used,
 * and calls the appropriate read function. If a new type of input
 * file is added to the driver, then a section needs to be added for
 * it here.
 *---------------------------------------------------------------------------*/
static int read_mesh(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh)
{
/* local declarations */
/*-----------------------------Execution Begins------------------------------*/
  if (pio_info->file_type == CHACO_FILE) {
    if (!read_chaco_file(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_chaco_mesh\n");
        return 0;
    }
  }
  else if (pio_info->file_type == NEMESIS_FILE) {
    if (!read_exoII_file(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_exoII_mesh\n");
        return 0;
    }
  }
  else if (pio_info->file_type == HYPERGRAPH_FILE) {
    if (!read_hypergraph_file(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_hypergraph_file\n");
        return 0;
    }
  }
  else if (pio_info->file_type == MATRIXMARKET_FILE) {
    if (!read_mm_file(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_mm_file\n");
        return 0;
    }
  }
  else if (pio_info->file_type == MATRIXMARKET_PLUS_FILE) {
    if (!read_mtxplus_file(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_mtxplus_file\n");
        
    }
  }
  else if (pio_info->file_type == NO_FILE_POINTS) {
    if (!create_random_input(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from create_random_input\n");
        return 0;
    }
  }
  else {
    Gen_Error(0, "fatal: Invalid file type.\n");
    return 0;
  }
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
static void print_input_info(ostream &os, int Num_Proc, PROB_INFO_PTR prob, 
PARIO_INFO_PTR pio, float zoltan_version)
{
int i;

  os << "Input values: " << endl;
  os << "  Zoltan version " << zoltan_version << endl;
  os << DRIVER_NAME << " version " <<  VER_STR << endl;
  os << "  Total number of Processors = " << Num_Proc << endl;

  os << endl << "  Performing load balance using " << prob->method << endl;
  os << "\tParameters:\n" ;
  for (i = 0; i < prob->num_params; i++)
    os << "\t\t" << prob->params[i].Name << " " << prob->params[i].Val << endl;

  if ((pio->init_dist_procs > 0) && (pio->init_dist_procs != Num_Proc)){
    os << "\n  Distribute input objects to only " << pio->init_dist_procs;
    os << " processes initially.\n";
  }
  if (pio->init_dist_type >= 0){
    os << "\n  Initially distribute input objects" ;
    switch (pio->init_dist_type){
      case INITIAL_FILE:
      case INITIAL_OWNER:
        os << " according to assignments in file.";
        break;
      case INITIAL_LINEAR:
        os << " in linear fashion (first n/p to process 0, etc).";
        break;
      case INITIAL_CYCLIC:
        os << " in cyclic (round robin) fashion.";
        break;
    }
    os << endl;
  }
  if (pio->init_dist_pins >= 0){
    os << "\n  Distribute pins";
    switch (pio->init_dist_pins){
      case INITIAL_FILE:
        os << " according to assignments in file.";
        break;
      case INITIAL_LINEAR:
        os << " in linear fashion (first n/p to process 0, etc).";
        break;
      case INITIAL_CYCLIC:
        os << " in cyclic (round robin) fashion.";
        break;
      case INITIAL_ROW:
        os << " so each process gets full rows.";
        break;
      case INITIAL_COL:
        os << " so each process gets full columns.";
        break;
      case INITIAL_ZERO:
        os << " all to process zero.";
        break;
    }
    os << endl;
  }


  os << "##########################################################" << endl;
}
/*****************************************************************************/
/*****************************************************************************/
static void initialize_mesh(MESH_INFO_PTR mesh)
{
/* Initializes mesh variables */

  mesh->dd = NULL;
  mesh->data_type = MESH;
  mesh->num_elems = mesh->num_nodes
                  = mesh->num_dims
                  = mesh->num_el_blks
                  = mesh->num_node_sets
                  = mesh->num_side_sets
                  = mesh->necmap
                  = mesh->elem_array_len
                  = mesh->gnhedges
                  = mesh->nhedges
                  = mesh->vwgt_dim
                  = mesh->ewgt_dim
                  = mesh->hewgt_dim
                  = mesh->blank_count
                  = mesh->global_blank_count
                  = 0;
  mesh->eb_names       = NULL;
  mesh->eb_ids         = NULL;
  mesh->eb_cnts        = NULL;
  mesh->eb_nnodes      = NULL;
  mesh->eb_nattrs      = NULL;
  mesh->ecmap_id       = NULL;
  mesh->ecmap_cnt      = NULL;
  mesh->ecmap_elemids  = NULL;
  mesh->ecmap_sideids  = NULL;
  mesh->ecmap_neighids = NULL;
  mesh->elements       = NULL;
  mesh->blank          = NULL;
  mesh->format         = ZOLTAN_COMPRESSED_EDGE;
  mesh->hgid           = NULL;
  mesh->hindex         = NULL;
  mesh->hvertex        = NULL;
  mesh->hvertex_proc   = NULL;
  mesh->heWgtId        = NULL;
  mesh->hewgts         = NULL;
}
