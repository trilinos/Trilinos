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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


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

#include <mpi.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_loadbal_const.h"
#include "dr_output_const.h"
#include "dr_err_const.h"
#include "dr_elem_util_const.h"

int Debug_Driver = 1;
int Test_DDirectory = 1;       /* Flag indicating whether to exercise 
                                  DDirectory */
int Test_Multi_Callbacks = 0;  /* Flag indicating whether to use 
                                  list-based (MULTI) callbacks */
int Test_Null_Import_Lists = 0;/* Flag indicating whether to test passing
                                  NULL import lists to Help_Migrate */
int Gnuplot_Output = 0;
int Print_Mesh_Info_File = 0;
int Number_Iterations = 1;
int Driver_Action = 1;	/* Flag indicating load-balancing or ordering. */
int Debug_Chaco_Input = 0;
int Chaco_In_Assign_Inv = 0;

static int read_mesh(int, int, PROB_INFO_PTR, PARIO_INFO_PTR, MESH_INFO_PTR);
static void print_input_info(FILE *fp, int Num_Proc, PROB_INFO_PTR prob);
static void initialize_mesh(MESH_INFO_PTR);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int main(int argc, char *argv[])
{
/* Local declarations. */
  struct Zoltan_Struct *zz;

  char  *cmd_file;
  char   cmesg[256]; /* for error messages */

  float  version;

  int    Proc, Num_Proc;
  int    iteration;
  int    error;

  MESH_INFO  mesh;             /* mesh information struct */
  PARIO_INFO pio_info;
  PROB_INFO  prob;

/***************************** BEGIN EXECUTION ******************************/

  /* initialize MPI */
  MPI_Init(&argc, &argv);

  /* get some machine information */
  MPI_Comm_rank(MPI_COMM_WORLD, &Proc);
  MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);

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
    exit(1);
  }

  /* initialize some variables */
  initialize_mesh(&mesh);

  pio_info.dsk_list_cnt		= -1;
  pio_info.num_dsk_ctrlrs	= -1;
  pio_info.pdsk_add_fact	= -1;
  pio_info.zeros		= -1;
  pio_info.file_type		= -1;
  pio_info.init_dist_type	= -1;
  pio_info.pdsk_root[0]		= '\0';
  pio_info.pdsk_subdir[0]	= '\0';
  pio_info.pexo_fname[0]	= '\0';

  prob.method[0]		= '\0';
  prob.num_params		= 0;
  prob.params			= NULL;

  /* Read in the ascii input file */
  if(Proc == 0) {
    printf("\n\nReading the command file, %s\n", cmd_file);
    if(!read_cmd_file(cmd_file, &prob, &pio_info))
    {
      sprintf(cmesg,"fatal: Could not read in the command file"
              " \"%s\"!\n", cmd_file);
      Gen_Error(0, cmesg);
      error_report(Proc);
      exit(1);
    }

    if (!check_inp(&prob, &pio_info)) {
      Gen_Error(0, "fatal: Error in user specified parameters.\n");
      error_report(Proc);
      exit(1);
    }

    print_input_info(stdout, Num_Proc, &prob);
  }

  /* broadcast the command info to all of the processor */
  brdcst_cmd_info(Proc, &prob, &pio_info);

  Zoltan_Set_Param(NULL, "DEBUG_MEMORY", "1");

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
    exit(1);
  }

  /* Loop over read and balance for a number of iterations */
  /* (Useful for testing REUSE parameters in Zoltan.) */
  for (iteration = 1; iteration <= Number_Iterations; iteration++) {

    /*
     * now read in the mesh and element information.
     * This is the only function call to do this. Upon return,
     * the mesh struct and the elements array should be filled.
     */
    if (!read_mesh(Proc, Num_Proc, &prob, &pio_info, &mesh)) {
      Gen_Error(0, "fatal: Error returned from read_mesh\n");
      error_report(Proc);
      exit(1);
    }

/* KDD Cool test of changing number of partitions 
    sprintf(cmesg, "%d", Num_Proc * iteration);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTITIONS", cmesg);
*/

    /*
     * Produce files to verify input.
     */
    if (Debug_Driver > 2) {
      if (!output_results(cmd_file,"in",Proc,Num_Proc,&prob,&pio_info,&mesh)) {
        Gen_Error(0, "fatal: Error returned from output_results\n");
        error_report(Proc);
        exit(1);
      }
      if (Gnuplot_Output)
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
      exit(1);
    }

    /* Reset the mesh data structure for next iteration. */
    if (iteration < Number_Iterations) {
      free_mesh_arrays(&mesh);
      initialize_mesh(&mesh);
    }

  } /* End of loop over read and balance */

  Zoltan_Destroy(&zz);

  Zoltan_Memory_Stats();

  /*
   * output the results
   */
  if (!output_results(cmd_file,"out",Proc,Num_Proc,&prob,&pio_info,&mesh)) {
    Gen_Error(0, "fatal: Error returned from output_results\n");
    error_report(Proc);
    exit(1);
  }

  if (Gnuplot_Output) {
    if (!output_gnu(cmd_file,"out",Proc,Num_Proc,&prob,&pio_info,&mesh)) {
      Gen_Error(0, "warning: Error returned from output_gnu\n");
      error_report(Proc);
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
    if (!read_chaco_mesh(Proc, Num_Proc, prob, pio_info, mesh)) {
        Gen_Error(0, "fatal: Error returned from read_chaco_mesh\n");
        return 0;
    }
  }
  else if (pio_info->file_type == NEMESIS_FILE) {
    if (!read_exoII_mesh(Proc, Num_Proc, prob, pio_info, mesh)) {
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
#endif
  else {
    Gen_Error(0, "fatal: Invalid file type.\n");
    return 0;
  }
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
static void print_input_info(FILE *fp, int Num_Proc, PROB_INFO_PTR prob)
{
int i;

  fprintf(fp, "Input values:\n");
  fprintf(fp, "  %s version %s\n", DRIVER_NAME, VER_STR);
  fprintf(fp, "  Total number of Processors = %d\n\n", Num_Proc);

  fprintf(fp, "\n  Performing load balance using %s.\n", prob->method);
  fprintf(fp, "\tParameters:\n");
  for (i = 0; i < prob->num_params; i++)
    fprintf(fp, "\t\t%s %s\n", prob->params[i][0], prob->params[i][1]);

  fprintf(fp, "##########################################################\n");
}
/*****************************************************************************/
/*****************************************************************************/
static void initialize_mesh(MESH_INFO_PTR mesh)
{
/* Initializes mesh variables */

  mesh->data_type = MESH;
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
  mesh->hindex         = mesh->hvertex
                       = NULL;
  mesh->hewgts         = NULL;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
