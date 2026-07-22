// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>
/*--------------------------------------------------------------------------*/
/* Purpose: Driver for dynamic load-balance library, ZOLTAN.                */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
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
#include <signal.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_loadbal_const.h"
#include "dr_util_const.h"
#include "dr_output_const.h"
#include "dr_err_const.h"
#include "dr_elem_util_const.h"
#include "dr_dd.h"
#include "dr_compress_const.h"

/* #define IGNORE_FIRST_ITERATION_STATS */
/* #define RANDOM_DIST */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern void Zoltan_write_linux_meminfo(int append, char *msg, int committedOnly);
int Debug_Driver = 1;
int Number_Iterations = 1;
int Driver_Action = 1;	/* Flag indicating coloring, load-balancing or ordering. */
int Debug_Chaco_Input = 0;
int Chaco_In_Assign_Inv = 0;
struct Test_Flags Test;
struct Output_Flags Output;

double Total_Partition_Time = 0.0;  /* Total over Number_Iterations */
int CITESEER[200]; /* Hack for Citeseer Vtx_Inc experiments */

static int read_mesh(int, int, PROB_INFO_PTR, PARIO_INFO_PTR, MESH_INFO_PTR);
static void print_input_info(FILE *, int, PROB_INFO_PTR, PARIO_INFO_PTR, float );
static void initialize_mesh(MESH_INFO_PTR, int proc);
static void remove_random_vertices(MESH_INFO_PTR mesh, int iter, float factor);
#ifdef DEBUG_READ_MESH
static void print_mesh(int proc, MESH_INFO_PTR m, int *tp, int *the, int *tv);
#endif

int my_rank=-1;

void meminfo_signal_handler(int sig)
{
  char msg[128];

  sprintf(msg,"(%d) Received signal %d\n",my_rank,sig);

  /* Signal handler for Linux that helps us to understand */
  /* whether failure was due to insufficient memory.      */

  signal(SIGINT, SIG_IGN);
  signal(SIGTERM, SIG_IGN);
  signal(SIGABRT, SIG_IGN);
  signal(SIGSEGV, SIG_IGN);
  signal(SIGFPE, SIG_IGN);

  Zoltan_write_linux_meminfo(0, msg,0);

  exit(sig);
}


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
  MPI_Comm_rank(zoltan_get_global_comm(), &Proc);
  MPI_Comm_size(zoltan_get_global_comm(), &Num_Proc);

  my_rank = Proc;

#ifdef HOST_LINUX
  signal(SIGSEGV, meminfo_signal_handler);
  signal(SIGINT, meminfo_signal_handler);
  signal(SIGTERM, meminfo_signal_handler);
  signal(SIGABRT, meminfo_signal_handler);
  signal(SIGFPE, meminfo_signal_handler);
#endif

#ifdef ZOLTAN_PURIFY
  printf("%d of %d ZDRIVE LAUNCH pid = %d file = %s\n", 
         Proc, Num_Proc, getpid(), argv[1]);
#endif

  /* Initialize flags */
  Test.DDirectory = 0;
  Test.Local_Parts = 0;
  Test.Fixed_Objects = 0;
  Test.Drops = 0;
  Test.RCB_Box = 0;
  Test.Multi_Callbacks = 0;
  Test.Graph_Callbacks = 1;
  Test.Hypergraph_Callbacks = 1;
  Test.Gen_Files = 0;
  Test.Null_Lists = NO_NULL_LISTS;
  Test.Dynamic_Weights = .0;
  Test.Dynamic_Graph = .0;
  Test.Vtx_Inc = 0;

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
  initialize_mesh(&mesh, Proc);

  pio_info.dsk_list_cnt		= -1;
  pio_info.file_comp            = STANDARD;
  pio_info.num_dsk_ctrlrs	= -1;
  pio_info.pdsk_add_fact	= -1;
  pio_info.zeros		= -1;
  pio_info.file_type		= -1;
  pio_info.chunk_reader         = 0;
  pio_info.init_dist_type	= -1;
  pio_info.init_size		= ZOLTAN_ID_INVALID;
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

  MPI_Allreduce(&error, &gerror, 1, MPI_INT, MPI_MAX, zoltan_get_global_comm());
  if (gerror) goto End;

  /* broadcast the command info to all of the processor */
  brdcst_cmd_info(Proc, &prob, &pio_info, &mesh);

  Zoltan_Set_Param(NULL, "DEBUG_MEMORY", "1");
  print_output = Output.Text;

  /*
   *  Create a Zoltan structure.
   */
  if ((zz = Zoltan_Create(zoltan_get_global_comm())) == NULL) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_Create()\n");
    return 0;
  }

  if (!setup_zoltan(zz, Proc, &prob, &mesh, &pio_info)) {
    Gen_Error(0, "fatal: Error returned from setup_zoltan\n");
    error_report(Proc);
    print_output = 0;
    goto End;
  }

  /* srand(Proc); Different seeds on different procs. */
  srand(1);  /* Same seed everywhere. */

  if (Test.Dynamic_Weights){
    /* Set obj weight dim to 1; can be overridden by user parameter */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
  }

  /* Loop over read and balance for a number of iterations */
  /* (Useful for testing REUSE parameters in Zoltan.) */
  for (iteration = 1; iteration <= Number_Iterations; iteration++) {

    if (Proc == 0) {
      printf("Starting iteration %d\n", iteration); 
      fflush(stdout);
    }

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
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", cmesg);
#endif

    /*
     * Produce files to verify input.
     */
    if (iteration == 1) {
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
      if (Test.Vtx_Inc<0){
        /* Read Citeseer data from file */
        FILE *fp;
        int i=0;
        if (Proc==0){
          fp = fopen("months.txt", "r");
          if (!fp)
            printf("ERROR: Couldn't open file months.txt\n");
          while (fscanf(fp, "%d", &CITESEER[i])==1){
            ++i;
          }
          fclose(fp);
        }
        MPI_Bcast (CITESEER, 200, MPI_INT, 0, zoltan_get_global_comm());
      }
    }

    if (Test.Dynamic_Graph > 0.0){
      if (mesh.data_type == ZOLTAN_GRAPH) {
        remove_random_vertices(&mesh, iteration, Test.Dynamic_Graph); 
      }
      else{
        Gen_Error(0, "fatal: \"test dynamic graph\" only works on graphs, not hypergraphs\n");
        error_report(Proc);
        print_output = 0;
        goto End;
      }
    }

    if (Test.Vtx_Inc){
      if (mesh.data_type == ZOLTAN_HYPERGRAPH ) {
        if (Test.Vtx_Inc>0)
          mesh.visible_nvtx += Test.Vtx_Inc; /* Increment uniformly */
        else
          mesh.visible_nvtx = CITESEER[iteration-1]; /* Citeseer document matrix. */
      }
      else{
        Gen_Error(0, "fatal: \"vertex increment\" only works on hypergraphs\n");
        error_report(Proc);
        print_output = 0;
        goto End;
      }
    }

    /*
     * now run Zoltan to get a new load balance and perform
     * the migration
     */
  
#ifdef IGNORE_FIRST_ITERATION_STATS
if (iteration == 1) {
  /* Exercise partitioner once on Tbird because first run is slow. */
  /* Lee Ann suspects Tbird is loading shared libraries. */
  struct Zoltan_Struct *zzcopy;
  zzcopy = Zoltan_Copy(zz);
  /* Don't do any migration or accumulate any stats. */
  if (Proc == 0) printf("%d KDDKDD IGNORING FIRST ITERATION STATS\n", Proc);
  Zoltan_Set_Param(zzcopy, "RETURN_LISTS", "NONE");
  Zoltan_Set_Param(zzcopy, "FINAL_OUTPUT", "0");
  Zoltan_Set_Param(zzcopy, "USE_TIMERS", "0");
  if (!run_zoltan(zzcopy, Proc, &prob, &mesh, &pio_info)) {
    Gen_Error(0, "fatal: Error returned from run_zoltan\n");
    error_report(Proc);
    print_output = 0;
    goto End;
  }
  Zoltan_Destroy(&zzcopy);
}
#endif /* IGNORE_FIRST_ITERATION_STATS */
#ifdef RANDOM_DIST
 if (iteration % 2 == 0) {
   char LB_METHOD[1024];

  if (Proc == 0) printf("%d CCCC Randomizing the input\n", Proc);
   strcpy(LB_METHOD, prob.method);
   strcpy(prob.method, "RANDOM");
   Zoltan_Set_Param(zz, "LB_METHOD", "RANDOM");
   Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    if (!run_zoltan(zz, Proc, &prob, &mesh, &pio_info)) {
      Gen_Error(0, "fatal: Error returned from run_zoltan\n");
      error_report(Proc);
      print_output = 0;
      goto End;
    }
   Zoltan_Set_Param(zz, "RETURN_LISTS", "NONE");
   Zoltan_Set_Param(zz, "LB_METHOD", LB_METHOD);
   strcpy(prob.method, LB_METHOD);
  if (Proc == 0) printf("%d CCCC Randomizing the input -- END\n", Proc);
 }
#endif /* RANDOM_DIST */
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
      char str[32];
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
        /* Increase weights in some parts */
        if (Test.Dynamic_Weights){
          /* Randomly pick 10% of parts to "refine" */
          /* Note:  Assumes at least 10 parts!  */
          /* Increase vertex weight, and also edge weights? TODO */
          j = (int) ((10.0*rand())/RAND_MAX + .5);
          for (i = 0; i < mesh.num_elems; i++) {
            if ((mesh.elements[i].my_part%10) == j){
                mesh.elements[i].cpu_wgt[0] = Test.Dynamic_Weights*(1+rand()%5);
            }
          }
        }
      }
      /* change the ParMETIS Seed */
      sprintf(str, "%d", iteration);
#ifdef ZOLTAN_PARMETIS      
      Zoltan_Set_Param(zz, "PARMETIS_SEED", str);
#endif
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
    mm_cleanup(mesh);
  }
  else if (pio_info->file_type == MATRIXMARKET_PLUS_FILE) {
    if (!read_mtxplus_file(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_mtxplus_file\n");
        return 0;
    }
    /* KDDKDD 3/26/10:  
     * Eventually, we should do cleanup here to address bug 3346.
     * but doing so will change the answer files.
     * mm_cleanup(mesh);
     */
  }
  else if (pio_info->file_type == NO_FILE_POINTS) {
    if (!create_random_input(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from create_random_input\n");
        return 0;
    }
  }
  else if (pio_info->file_type == NO_FILE_TRIANGLES) {
    if (!create_random_triangles(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from create_random_triangles\n");
        return 0;
    }
  }
  else if (pio_info->file_type == NO_FILE_GRAPH) {
    if (!create_a_graph(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from create_a_graph\n");
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
      print_mesh(Proc, mesh, &pins, &he, &verts);
    }
    MPI_Barrier(zoltan_get_global_comm());
  }
  MPI_Reduce(&pins, &tpins, 1, MPI_INT, MPI_SUM, 0, zoltan_get_global_comm());
  MPI_Reduce(&he, &the, 1, MPI_INT, MPI_SUM, 0, zoltan_get_global_comm());
  MPI_Reduce(&verts, &tverts, 1, MPI_INT, MPI_SUM, 0, zoltan_get_global_comm());
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
int idtypesize;
char *idtypename;

  fprintf(fp, "Input values:\n");
  fprintf(fp, "  Zoltan version %g\n",zoltan_version);
  fprintf(fp, "  %s version %s\n", DRIVER_NAME, VER_STR);
  fprintf(fp, "  Total number of Processors = %d\n", Num_Proc);
  idtypesize = Zoltan_get_global_id_type(&idtypename);
  fprintf(fp, "  ZOLTAN_ID_TYPE size %d name %s\n", idtypesize, idtypename);

  fprintf(fp, "\n  Performing load balance using %s.\n", prob->method);
  fprintf(fp, "\tParameters:\n");
  for (i = 0; i < prob->num_params; i++)
    fprintf(fp, "\t\t%s %s\n", prob->params[i].Name, prob->params[i].Val);

  if ((pio->init_dist_procs > 0) && (pio->init_dist_procs != Num_Proc)){
    fprintf(fp, "\n  Distribute input objects to only %d processes initially.\n", 
             pio->init_dist_procs);
  }
  if (pio->chunk_reader > 0){
    fprintf(fp, "\n  Initially read input file in chunks (due to file size)\n");
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
      case INITIAL_NO_DIST:
        fprintf(fp," not at all.  They are created as distributed.");
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
static void initialize_mesh(MESH_INFO_PTR mesh, int proc)
{
/* Initializes mesh variables */

  mesh->dd = NULL;
  mesh->data_type = MESH;
  mesh->gnhedges = mesh->global_blank_count = 0;
  mesh->num_elems = mesh->num_nodes
                  = mesh->num_dims
                  = mesh->num_el_blks
                  = mesh->num_node_sets
                  = mesh->num_side_sets
                  = mesh->necmap
                  = mesh->elem_array_len
                  = mesh->nhedges
                  = mesh->vwgt_dim
                  = mesh->ewgt_dim
                  = mesh->hewgt_dim
                  = mesh->blank_count
                  = 0;
  mesh->eb_names       = NULL;
  mesh->eb_etypes      = NULL;
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
  mesh->proc           = proc;
  mesh->visible_nvtx   = 0;
}

static void remove_random_vertices(MESH_INFO_PTR mesh, int iteration, 
                                   float blank_factor) 
{
int i, j, blankmine = (mesh->proc % 2) == (iteration % 2);
ZOLTAN_ID_TYPE tmp, total_vertices;
ELEM_INFO *elem;

  for (i=0; i < mesh->num_elems; i++){
    safe_free((void **)(void *)&(mesh->elements[i].adj_blank));
  }

  safe_free((void **)(void *)&mesh->blank);
  mesh->blank_count = 0;
  mesh->global_blank_count = 0;

  /*
   * Mark some portion of vertices as blanked.  The graph callbacks
   * will not report blanked vertices.  
   */

  if ((blank_factor <= 0.0) || (blank_factor >= 1.0)){
    return;
  }

  mesh->blank = (int *)calloc(sizeof(int) , mesh->elem_array_len);
  if (mesh->elem_array_len && !mesh->blank){
    Gen_Error(0, "memory allocation in remove_random_vertices");
    error_report(mesh->proc);
    return;
  }
 
  for (i=0; i < mesh->num_elems; i++){
  
    /* Each vertex (element) has probability given by
     * by blank_factor of being blanked.  The blanked vertices should
     * vary somewhat in each iteration.
     */
    elem = mesh->elements + i;

    srand((long int)(elem->globalID * iteration));

    if (blankmine && (((double)rand()/RAND_MAX) <= (double)blank_factor)){
      mesh->blank[i] = 1.0;
      mesh->blank_count++;
    }

    /* adj_blank marks each off processor adjacent vertex as
     * either blanked (1) or not blanked (0).
     */
    elem->adj_blank = (int *)calloc(sizeof(int), elem->nadj);

    if (elem->nadj){
      if (!elem->adj_blank){
        Gen_Error(0, "memory allocation in remove_random_vertices");
        error_report(mesh->proc);
        return;
      }
      for (j=0; j<elem->nadj; j++){
        if (elem->adj_proc[j] != mesh->proc){
          srand((long int)(elem->adj[j] * iteration));
          if (((elem->adj_proc[j] % 2) == (iteration % 2) ) && (((double)rand()/RAND_MAX) <= (double)blank_factor)){
            elem->adj_blank[j] = 1.0;
          }
        }
      }
    }
  }

  MPI_Allreduce(&mesh->blank_count, &mesh->global_blank_count, 1, ZOLTAN_ID_MPI_TYPE, MPI_SUM, zoltan_get_global_comm());

  tmp = (ZOLTAN_ID_TYPE)mesh->num_elems;

  MPI_Reduce(&tmp, &total_vertices, 1, ZOLTAN_ID_MPI_TYPE, MPI_SUM, 0, zoltan_get_global_comm());

  if (mesh->proc == 0){
    printf("Dynamic graph factor %0.4f, " ZOLTAN_ID_SPEC " vertices, " ZOLTAN_ID_SPEC " blanked (%0.2f%%)\n",
            blank_factor, total_vertices, mesh->global_blank_count,
            ((double)mesh->global_blank_count*100.0/total_vertices));
  }
  fflush(stdout);
  if (Debug_Driver > 1) {
    MPI_Barrier(zoltan_get_global_comm());
    if (mesh->num_elems){
      printf("Proc %d: %d vertices, %d blanked (%0.2f%%)\n",
            mesh->proc, mesh->num_elems,  mesh->blank_count,
            ((float)mesh->blank_count*100.0/mesh->num_elems));
    }
    else{
      printf("Proc %d: 0 vertices\n", mesh->proc);
    }
    fflush(stdout);
    MPI_Barrier(zoltan_get_global_comm());
  }
}

#ifdef DEBUG_READ_MESH
static void print_mesh(int proc, MESH_INFO_PTR m, int *tp, int *the, int *tv)
{
  int i, j, ii;
  ZOLTAN_ID_TYPE globalID, adj;
  ELEM_INFO_PTR el;
  printf("Global number of hyperedges " ZOLTAN_ID_SPEC "\n",m->gnhedges);
  if (m->format == ZOLTAN_COMPRESSED_EDGE){
    printf("Pins: %d edges\n",m->nhedges);
  }
  else{
    printf("Pins: %d vertices\n",m->nhedges);
  }
  for (i=0; i<m->nhedges; i++){
    printf("  " ZOLTAN_ID_SPEC ": ", m->hgid[i]);
    for (j=m->hindex[i],ii=0; j<m->hindex[i+1]; j++,ii++){
      if (ii && (ii%15==0)) printf("\n       ");
      printf(ZOLTAN_ID_SPEC " ", m->hvertex[j]);
    }
    printf("\n");
  }
  
  printf("Total pins: %d\n", m->hindex[m->nhedges]);
  printf("%d vertices: ", m->num_elems);

  el = m->elements;
  
  for (i=0; i<m->num_elems; i++){
    printf(ZOLTAN_ID_SPEC " (%d adj: ", el->globalID, el->nadj);
    for (j=0; j<el->nadj; j++){
      adj = el->adj[j];
      if (el->adj_proc[j] == proc){
        globalID = m->elements[adj].globalID;
      }
      else{
        globalID = adj;
      } 
      if (j && (j%15==0)) printf("\n       ");
      printf(ZOLTAN_ID_SPEC " ",globalID);
    }
    printf(")\n");

    el++;
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
