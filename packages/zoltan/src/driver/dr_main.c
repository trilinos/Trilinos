/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <mpi.h>
/*--------------------------------------------------------------------------*/
/* Purpose: Driver for dynamic load-balance library, ZOLTAN.                */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
/* Supported Environment(s):    Intel Paragon                               */
/*                              Intel Teraflop                              */
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*                                                                          */
/*    30 March 1999:    Date of creation                                    */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_loadbal_const.h"
#include "dr_output_const.h"
#include "dr_err_const.h"
#include "dr_elem_util_const.h"
#include "dr_dd.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int Debug_Driver = 1;
int Number_Iterations = 1;
int Driver_Action = 1;	/* Flag indicating load-balancing or ordering. */
int Debug_Chaco_Input = 0;
int Chaco_In_Assign_Inv = 0;
struct Test_Flags Test;
struct Output_Flags Output;

double Total_Partition_Time = 0.0;  /* Total over Number_Iterations */

static int read_mesh(int, int, PROB_INFO_PTR, PARIO_INFO_PTR, MESH_INFO_PTR);
static void print_input_info(FILE *, int, PROB_INFO_PTR, PARIO_INFO_PTR, float );
static void initialize_mesh(MESH_INFO_PTR);
#ifdef DEBUG_READ_MESH
static void print_mesh(MESH_INFO_PTR m, int *tp, int *the, int *tv);
#endif

#include <unistd.h>

#ifdef VAMPIR
#include <VT.h>
#endif
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int main(int argc, char *argv[])
{
/* Local declarations. */
  struct Zoltan_Struct *zz = NULL;

  char  *cmd_file;
  char   cmesg[256]; /* for error messages */

  float  version;

  int    Proc, Num_Proc;
  int    iteration;
  int    error, gerror;
  int    print_output = 1;

  MESH_INFO  mesh;             /* mesh information struct */
  PARIO_INFO pio_info;
  PROB_INFO  prob;

/***************************** BEGIN EXECUTION ******************************/

  /* initialize MPI */
  MPI_Init(&argc, &argv);

#ifdef VAMPIR
  VT_initialize(&argc, &argv);
#endif

  /* get some machine information */
  MPI_Comm_rank(MPI_COMM_WORLD, &Proc);
  MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);

  /* Initialize flags */
  Test.DDirectory = 0;
  Test.Local_Partitions = 0;
  Test.Fixed_Objects = 0;
  Test.Drops = 0;
  Test.RCB_Box = 0;
  Test.Multi_Callbacks = 0;
  Test.Gen_Files = 0;
  Test.Null_Lists = NONE;
  Test.Dynamic_Hgraph = .0;

  Output.Text = 1;
  Output.Gnuplot = 0;
  Output.Nemesis = 0;
  Output.Plot_Partitions = 0;
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
    fprintf(stderr, "MAIN: ERROR in command line,");
    if(Proc == 0)
    {
      fprintf(stderr, " usage:\n");
      fprintf(stderr, "\t%s [command file]", DRIVER_NAME);
    }
    exit(1);
    break;
  }

  /* initialize Zoltan */
  if ((error = Zoltan_Initialize(argc, argv, &version)) != ZOLTAN_OK) {
    sprintf(cmesg, "fatal: Zoltan_Initialize returned error code, %d", error);
    Gen_Error(0, cmesg);
    error_report(Proc);
    print_output = 0;
    goto End;
  }

  /* initialize some variables */
  initialize_mesh(&mesh);

  pio_info.dsk_list_cnt		= -1;
  pio_info.num_dsk_ctrlrs	= -1;
  pio_info.pdsk_add_fact	= -1;
  pio_info.zeros		= -1;
  pio_info.file_type		= -1;
  pio_info.init_dist_type	= -1;
  pio_info.init_size		= -1;
  pio_info.init_dim 		= -1;
  pio_info.init_vwgt_dim 	= -1;
  pio_info.init_dist_pins       = -1;
  pio_info.pdsk_root[0]		= '\0';
  pio_info.pdsk_subdir[0]	= '\0';
  pio_info.pexo_fname[0]	= '\0';

  prob.method[0]		= '\0';
  prob.num_params		= 0;
  prob.params			= NULL;

  /* Read in the ascii input file */
  error = gerror = 0;
  if (Proc == 0) {
    printf("\n\nReading the command file, %s\n", cmd_file);
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

    print_input_info(stdout, Num_Proc, &prob, &pio_info, version);
  }

  MPI_Allreduce(&error, &gerror, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (gerror) goto End;

  /* broadcast the command info to all of the processor */
  brdcst_cmd_info(Proc, &prob, &pio_info, &mesh);

  Zoltan_Set_Param(NULL, "DEBUG_MEMORY", "1");
  print_output = Output.Text;

  /*
   *  Create a Zoltan structure.
   */
  if ((zz = Zoltan_Create(MPI_COMM_WORLD)) == NULL) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_Create()\n");
    return 0;
  }

  if (!setup_zoltan(zz, Proc, &prob, &mesh)) {
    Gen_Error(0, "fatal: Error returned from setup_zoltan\n");
    error_report(Proc);
    print_output = 0;
    goto End;
  }

  srand(Proc);

  /* Loop over read and balance for a number of iterations */
  /* (Useful for testing REUSE parameters in Zoltan.) */
  for (iteration = 1; iteration <= Number_Iterations; iteration++) {

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

      if (mesh.data_type == HYPERGRAPH && !build_elem_dd(&mesh)) {
        Gen_Error(0, "fatal: Error returned from build_elem_dd\n");
        error_report(Proc);
        print_output = 0;
        goto End;
      }
    }


#ifdef KDDKDD_COOL_TEST
/* KDD Cool test of changing number of partitions  */
    sprintf(cmesg, "%d", Num_Proc * iteration);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTITIONS", cmesg);
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
    if (!run_zoltan(zz, Proc, &prob, &mesh, &pio_info)) {
      Gen_Error(0, "fatal: Error returned from run_zoltan\n");
      error_report(Proc);
      print_output = 0;
      goto End;
    }

    /* Reset the mesh data structure for next iteration. */
    if (iteration < Number_Iterations) {
      int i, j;
      float tmp;
      float twiddle = 0.01;
      char str[4];
      /* Perturb coordinates of mesh */
      if (mesh.data_type == ZOLTAN_GRAPH){
        for (i = 0; i < mesh.num_elems; i++) {
          for (j = 0; j < mesh.num_dims; j++) {
            /* tmp = ((float) rand())/RAND_MAX; *//* Equiv. to sjplimp's test */
            tmp = (float) (i % 10) / 10.;
            mesh.elements[i].coord[0][j] += twiddle * (2.0*tmp-1.0);
            mesh.elements[i].avg_coord[j] = mesh.elements[i].coord[0][j];
          }
        }
      }
      else if (mesh.data_type == HYPERGRAPH) {
        printf("Debug: dynamic hgraph pertubation=%f\n", Test.Dynamic_Hgraph);
        /* Randomly perturb pins. */
        for (i = 0; i < mesh.hindex[mesh.nhedges]; i++) {
          tmp = ((float) rand())/RAND_MAX; /* random in [0,1) */
          if (tmp < Test.Dynamic_Hgraph){
            mesh.hvertex[i] += tmp * mesh.num_elems;
            if (mesh.hvertex[i] > mesh.num_elems)
              mesh.hvertex[i] -= mesh.num_elems;
          }
        }
      }
      /* change the ParMETIS Seed */
      sprintf(str, "%d", iteration);
      Zoltan_Set_Param(zz, "PARMETIS_SEED", str);
    }

  } /* End of loop over read and balance */

  if (Proc == 0) {
    printf("FILE %s:  Total:    %e seconds in Partitioning\n", 
           cmd_file, Total_Partition_Time);
    printf("FILE %s:  Average:  %e seconds per Iteration\n", 
           cmd_file, Total_Partition_Time/Number_Iterations);
  }

End:
  Zoltan_Destroy(&zz);
  if (mesh.dd) Zoltan_DD_Destroy(&(mesh.dd));

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
  
#ifdef VAMPIR
  VT_finalize();
#endif

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
#ifdef DEBUG_READ_MESH
  int i, tpins, the, tverts, pins, he, verts;
#endif
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
#ifdef ZOLTAN_HG
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
        return 0;
    }
  }
#endif
  else if (pio_info->file_type == NO_FILE) {
    if (!create_random_input(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from create_random_input\n");
        return 0;
    }
  }
  else {
    Gen_Error(0, "fatal: Invalid file type.\n");
    return 0;
  }

#ifdef DEBUG_READ_MESH
  for (i=0; i<Num_Proc; i++){
    if (i == Proc){
      printf("Process %d:\n",i);
      print_mesh(mesh, &pins, &he, &verts);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Reduce(&pins, &tpins, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&he, &the, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&verts, &tverts, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (Proc == 0){
    if (mesh->format == ZOLTAN_COMPRESSED_EDGE){
      printf("Total pins %d, total vertices %d, total rows %d\n",
                 tpins, tverts, the);
    }
    else{
      printf("Total pins %d, total vertices %d, total columns %d\n",
                 tpins, tverts, the);
    }
  }
#endif
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
static void print_input_info(FILE *fp, int Num_Proc, PROB_INFO_PTR prob, 
PARIO_INFO_PTR pio, float zoltan_version)
{
int i;

  fprintf(fp, "Input values:\n");
  fprintf(fp, "  Zoltan version %g\n",zoltan_version);
  fprintf(fp, "  %s version %s\n", DRIVER_NAME, VER_STR);
  fprintf(fp, "  Total number of Processors = %d\n", Num_Proc);

  fprintf(fp, "\n  Performing load balance using %s.\n", prob->method);
  fprintf(fp, "\tParameters:\n");
  for (i = 0; i < prob->num_params; i++)
    fprintf(fp, "\t\t%s %s\n", prob->params[i].Name, prob->params[i].Val);

  if ((pio->init_dist_procs > 0) && (pio->init_dist_procs != Num_Proc)){
    fprintf(fp, "\n  Distribute input objects to only %d processes initially.\n", 
             pio->init_dist_procs);
  }
  if (pio->init_dist_type >= 0){
    fprintf(fp, "\n  Initially distribute input objects");
    switch (pio->init_dist_type){
      case INITIAL_FILE:
      case INITIAL_OWNER:
        fprintf(fp," according to assignments in file.");
        break;
      case INITIAL_LINEAR:
        fprintf(fp," in linear fashion (first n/p to process 0, etc).");
        break;
      case INITIAL_CYCLIC:
        fprintf(fp," in cyclic (round robin) fashion.");
        break;
    }
    fprintf(fp, "\n");
  }
  if (pio->init_dist_pins >= 0){
    fprintf(fp, "\n  Distribute pins");
    switch (pio->init_dist_pins){
      case INITIAL_FILE:
        fprintf(fp," according to assignments in file.");
        break;
      case INITIAL_LINEAR:
        fprintf(fp," in linear fashion (first n/p to process 0, etc).");
        break;
      case INITIAL_CYCLIC:
        fprintf(fp," in cyclic (round robin) fashion.");
        break;
      case INITIAL_ROW:
        fprintf(fp," so each process gets full rows.");
        break;
      case INITIAL_COL:
        fprintf(fp," so each process gets full columns.");
        break;
      case INITIAL_ZERO:
        fprintf(fp," all to process zero.");
        break;
    }
    fprintf(fp, "\n");
  }


  fprintf(fp, "##########################################################\n");
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
  mesh->format         = ZOLTAN_COMPRESSED_EDGE;
  mesh->hgid           = NULL;
  mesh->hindex         = NULL;
  mesh->hvertex        = NULL;
  mesh->hvertex_proc   = NULL;
  mesh->heWgtId        = NULL;
  mesh->hewgts         = NULL;
}

#ifdef DEBUG_READ_MESH
static void print_mesh(MESH_INFO_PTR m, int *tp, int *the, int *tv)
{
  int i, j;
  printf("Global number of hyperedges %d\n",m->gnhedges);
  if (m->format == ZOLTAN_COMPRESSED_EDGE){
    printf("Pins: %d edges\n",m->nhedges);
  }
  else{
    printf("Pins: %d vertices\n",m->nhedges);
  }
  for (i=0; i<m->nhedges; i++){
    printf("  %d: ", m->hgid[i]);
    for (j=m->hindex[i]; j<m->hindex[i+1]; j++){
      printf("%d ", m->hvertex[j]);
    }
    printf("\n");
  }
  
  printf("Total pins: %d\n", m->hindex[m->nhedges]);

  printf("%d vertices: ", m->num_elems);
  for (i=0; i<m->num_elems; i++){
    printf("%d ", m->elements[i].globalID);
  }
  printf("\n");
  fflush(stdout);

  *tp = m->hindex[m->nhedges];
  *the = m->nhedges;
  *tv = m->num_elems;
}
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
