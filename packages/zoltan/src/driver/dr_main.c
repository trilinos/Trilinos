
/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_dr_main = "$Id$";
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

#include <mpi.h>

#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_err_const.h"

/* global mesh information struct variable */
MESH_INFO Mesh;

extern int run_zoltan(int, PROB_INFO_PTR, ELEM_INFO_PTR *);
extern int output_results(int, int, PARIO_INFO_PTR, ELEM_INFO *);

int read_mesh(int, int, PROB_INFO_PTR, PARIO_INFO_PTR, ELEM_INFO_PTR *);

int main(int argc, char *argv[])
{
/* Local declarations. */
  char  *cmd_file;
  char   cmesg[256]; /* for error messages */

  float  version;

  int    Proc, Num_Proc;
  int    error, i;

  PARIO_INFO    pio_info;
  ELEM_INFO_PTR elements;
  PROB_INFO     prob;

  /* Local function prototype */
  void   print_input(int, PROB_INFO_PTR);
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
      fprintf(stderr, "\t%s [command file]", UTIL_NAME);
    }
    exit(1);
    break;
  }

  /* initialize Zoltan */
  if ((error = LB_Initialize(argc, argv, &version)) != DLB_OK) {
    sprintf(cmesg, "fatal: LB_Initialize returned error code, %d", error);
    Gen_Error(0, cmesg);
    error_report(Proc);
    exit(1);
  }

  /* initialize some variables */
  Mesh.eb_names			= NULL;
  Mesh.eb_ids			= NULL;
  Mesh.eb_cnts			= NULL;
  Mesh.eb_nnodes		= NULL;
  Mesh.eb_nattrs		= NULL;

  pio_info.dsk_list_cnt		= -1;
  pio_info.num_dsk_ctrlrs	= -1;
  pio_info.pdsk_add_fact	= -1;
  pio_info.zeros		= -1;
  pio_info.file_type		= -1;
  pio_info.pdsk_root[0]		= '\0';
  pio_info.pdsk_subdir[0]	= '\0';
  pio_info.pexo_fname[0]	= '\0';

  prob.method[0]		= '\0';
  if ((error = LB_Initialize_Params_Array(prob.params)) != DLB_OK) {
    Gen_Error(0, "fatal: error returned from LB_Initialize_Params_Array");
    error_report(Proc);
    exit(1);
  }
  prob.tol			= -1.0;
  prob.read_coord		= 0;
  prob.gen_graph		= 0;

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

    print_input(Num_Proc, &prob);
  }

  /* broadcast the command info to all of the processor */
  brdcst_cmd_info(Proc, &prob, &pio_info);

  /*
   * now read in the mesh and element information.
   * This is the only function call to do this. Upon return,
   * the mesh struct and the elements array should be filled.
   */
  if (!read_mesh(Proc, Num_Proc, &prob, &pio_info, &elements)) {
      Gen_Error(0, "fatal: Error returned from read_mesh\n");
      error_report(Proc);
      exit(1);
  }

  /*
   * now run zoltan to get a new load balance and perform
   * the migration
   */
  if (!run_zoltan(Proc, &prob, &elements)) {
      Gen_Error(0, "fatal: Error returned from run_zoltan\n");
      error_report(Proc);
      exit(1);
  }

  /*
   * output the results
   */
  if (!output_results(Proc, Num_Proc, &pio_info, elements)) {
      Gen_Error(0, "fatal: Error returned from output_results\n");
      error_report(Proc);
      exit(1);
  }

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
int read_mesh(int Proc,
              int Num_Proc,
              PROB_INFO_PTR prob,
              PARIO_INFO_PTR pio_info,
              ELEM_INFO **elements)
{
/* local declarations */
/*-----------------------------Execution Begins------------------------------*/
  if (pio_info->file_type == CHACO_FILE) {
    if (!read_chaco_mesh(Proc, Num_Proc, prob, pio_info, elements)) {
        Gen_Error(0, "fatal: Error returned from read_chaco_mesh\n");
        return 0;
    }
  }
#ifdef NEMESIS_IO
  else if (pio_info->file_type == NEMESIS_FILE) {
    if (!read_exoII_mesh(Proc, Num_Proc, prob, pio_info, elements)) {
        Gen_Error(0, "fatal: Error returned from read_exoII_mesh\n");
        return 0;
    }
  }
#endif
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function prints out parameters as read from the command file.
 *---------------------------------------------------------------------------*/
void print_input(int Num_Proc, PROB_INFO_PTR prob)
{
/* local declarations */
  int i;
/*-----------------------------Execution Begins------------------------------*/

  printf("%s version %s\n", UTIL_NAME, VER_STR);
  printf("Total number of Processors = %d\n\n", Num_Proc);

  printf("\nPerforming load balance using %s.\n", prob->method);
  printf("\tTolerance: %lf\n", prob->tol);
  printf("\tParameters:\n");
  for (i = 0; i < LB_PARAMS_MAX_SIZE; i++)
    printf("\t\t%d %lf\n", i, prob->params[i]);

  if (prob->gen_graph)
    printf("\tGenerating graph.\n");
  else
    printf("\tNot generating graph.\n");

  if (prob->read_coord)
    printf("\tReading coordinates.\n");
  else
    printf("\tNot reading coordinates.\n");
}
