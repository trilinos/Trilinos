/*****************************************************************************/
/* Copyright 1999, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Sample driver for AZTEC/ML package. The software is tested by reading     */ 
/* input to drive ML from a file (ml_inputfile). Matrices can be created     */
/* by using Aztec's AZ_matrix_capture capability and then running            */
/* ml/util/az_capt2read.c (see comments for ways using matlab to create      */
/* matrices).                                                                */ 
/*****************************************************************************/
/* Here is a sample input file:
#
#  Test input file ML 
#       
# Notes: 
#   1) Captilization should not matter
#   2) Lines starting with # are comments
#   3) Parallel Partitioning File is not used
#      when only 1 processor is present.
#   4) comments in [] indicate available options.
#   5) The matrix must be stored in a file '.data' according
#      to Aztec's AZ_read_msr() format.
#   6) Including the file 'rhsfile' will cause this 
#      data (stored in Aztec's AZ_read_msr() format)
#      to be used as righthand side.
#   7) Including the file 'initguessfile' will cause this 
#      data (stored in Aztec's AZ_read_msr() format)
#      to be used as righthand side.
#   8) rigid body mode information can be input by
#      keeping files 'rigid_body_mode%d' (stored
#      in Aztec's AZ_read_msr() format) where %d
#      starts from 0 and increases. Each file
#      should contain 1 rigid body mode.
#   9) The utility ml/util/az_capt2read.c (see comments)
#      can be used to convert matlab/Aztec type data into 
#      AZ_read_msr() format.
#
-----------------------------------------------
      General Problem Specifications
-----------------------------------------------
Number of DOF per node       = 1
Parallel Partitioning File   = myfile  
Output Frequency             = 2       
Tolerance                    = 1.0e-11

-----------------------------------------------
      Solution Specifications
-----------------------------------------------
Max Number of Levels         = 4
Type of Smoother             = SymGaussSeidel 
#                              [Parasails, GaussSeidel, SymGaussSeidel, 
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis]
Smoother steps per level     = 7
Coarse grid solver           = SuperLU
#                              [Parasails, GaussSeidel, SymGaussSeidel, 
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis,
#                               SuperLU]
Coarse Grid iterations       = 1
Outer Iteration              = Cg
#                              [Cg, Bicgstab, Tfqmr, Gmres] 

-----------------------------------------------
      Aggregation Specifications
-----------------------------------------------
Type of Aggregation          = Mis
#                              [Mis, Uncoupled, Coupled]
Aggregate threshold          = 0.0
Max coarse size              = 30
Smoothed aggregation damping = 1.5 
Spectral norm calculation    = Anorm
#                              [Anorm, Calc]
# end of sample inputfile
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"
#include "ml_read_utils.h"
#include <math.h>

extern int AZ_using_fortran;
int    parasails_factorized = 0;
int    parasails_sym        = 1;
double parasails_thresh     = 0.01;
int    parasails_nlevels    = 0;
double parasails_filter     = 0.;
double parasails_loadbal    = 0.;
int my_proc_id;


int    *Ke_data_org = NULL, *Kn_data_org = NULL, *global_edge_inds = NULL, *global_node_inds = NULL,
       *global_edge_externs = NULL, *global_node_externs = NULL;
int    *reordered_glob_edges = NULL, *reordered_edge_externs = NULL;
int    *reordered_glob_nodes = NULL, *reordered_node_externs = NULL;
int    *cpntr = NULL, *Ke_bindx = NULL, *Kn_bindx = NULL, Nlocal_edges, iii, *Tmat_bindx;
int    *update, *update_index;
extern void AZ_transform_norowreordering(int proc_config[], int *external[], int bindx[], double val[],
                  int update[], int *update_index[], int *extern_index[],
                  int *data_org[], int N_update, int indx[], int bnptr[],
					 int rnptr[], int *cnptr[], int mat_type);
extern void AZ_input_msr_matrix_nodiag(char datafile[], int update[], double **val, int **bindx, 
				       int N_update, int proc_config[]);
extern void AZ_add_new_row_nodiag(int therow, int *nz_ptr, int *current, double **val,
                    int **bindx, char *input, FILE *dfp, int *msr_len,
				  int *column0);
extern void ML_find_local_indices(int N_update, int bindx[], int update[],
				  int *sorted_ext, int N_external, int map[]);

extern int eye_getrows(void *data, int N_requested_rows, int requested_rows[],
		       int allocated_space, int columns[], double values[], int row_lengths[]);

extern int eye_matvec(void *Amat_in, int ilen, double p[], int olen, double ap[]);
struct reader_context *context;
extern int  ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans);
extern int ML_Operator_Transpose_byrow(ML_Operator *A, ML_Operator *Atrans);

int main(int argc, char *argv[])
{

	int num_PDE_eqns=1, N_levels=3, nsmooth=2;
	int Nglobal_edges, Nglobal_nodes, Nlocal_edges, Nlocal_nodes;
	int level, coarsest_level;

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];

  /* data structure for matrix corresponding to the fine grid */

  double *Ke_val = NULL, *Kn_val, *Tmat_val = NULL, *xxx, *rhs, solve_time, setup_time, start_time;
  AZ_MATRIX *Ke_mat, *Kn_mat;
  ML_Operator *Tmat, *Tmat_trans, *Kn_coarse, *Tcoarse, *Tcoarse_trans, *Pn_coarse, *Rn_coarse, *Pe;
  ML_Operator *Ttrans_byrow;
  AZ_PRECOND *Pmat = NULL;
  ML *ml_edges, *ml_nodes;
  FILE *fp;
  int ch,i, j, Nrigid, *garbage, nblocks, *blocks;
  struct AZ_SCALING *scaling;
  ML_Aggregate *ag;
  double *mode, *rigid, alpha, *newval;
  char filename[80];
  int    one = 1;
  int allocated = 0, *newbindx, offset, current, *block_list = NULL,  k, block;
  int Tmat_cols;
  double *ones, *result;
  int *row_ptr, lower, upper, nz_ptr, Ncols;
  struct ML_CSR_MSRdata *csr_data;
  double *diag;
  int pcount, kk, procs_cols[20], proc_count[20];
  int *bindx = NULL;
  double *val = NULL;
  int counter, row_length, itemp, row_start;
  int *Tcoarse_bindx, *Tcoarse_rowptr, agg1, agg2;
  double *Tcoarse_val, dtemp;
  int Nexterns, *sorted_ext, *map;
  int Nnz;
  char str[80];

#ifdef ML_partition
   FILE *fp2;
   int count, *pcounts, *proc_id;
   int Tmat_size, *proc_assignment, p1, p2, col1, col2;

   nblocks = -1;
   if (argc < 2) {
     printf("Usage: ml_readfile num_processors\n");
     exit(1);
   }
   else sscanf(argv[1],"%d",&nblocks);
   if (nblocks == -1) {
     printf("Usage: ml_readfile num_processors\n");
     exit(1);
   }

#endif

#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif
  my_proc_id = proc_config[AZ_node];

#ifndef ML_partition
   if (proc_config[AZ_node] == 0) {
      ML_Reader_ReadInput("ml_inputfile", &context);
   }
   else context = (struct reader_context *) malloc(sizeof(struct reader_context));
   AZ_broadcast((char *) context,  sizeof(struct reader_context), proc_config,
                AZ_PACK);
   AZ_broadcast((char *) NULL        ,   0          , proc_config, AZ_SEND);

   N_levels = context->N_levels;
   nsmooth   = context->nsmooth;
   num_PDE_eqns = context->N_dofPerNode;
#else
   context = (struct reader_context *) malloc(sizeof(struct reader_context));
   ML_Reader_InitContext(context);
#endif

  /* read in the number of edge matrix equations */
  Nglobal_edges = 0;
  if (proc_config[AZ_node] == 0) {
#    ifdef binary
	fp=fopen("Ke_mat.az","rb");
#    else
	fp=fopen("Ke_mat.az","r");
#    endif
     if (fp==NULL) {
        printf("couldn't open file Ke_mat.az\n");
        exit(1);
     }
#    ifdef binary
        fread(&Nglobal_edges, sizeof(int), 1, fp);
#    else
        fscanf(fp,"%d",&Nglobal_edges);
#    endif
     fclose(fp);
  }
  Nglobal_edges = AZ_gsum_int(Nglobal_edges, proc_config);

  /* read in the number of node matrix equations */
  Nglobal_nodes = 0;
  if (proc_config[AZ_node] == 0) {
#    ifdef binary
	fp=fopen("Kn_mat.az","rb");
#    else
	fp=fopen("Kn_mat.az","r");
#    endif
     if (fp==NULL) {
        printf("couldn't open file Kn_mat.az\n");
        exit(1);
     }
#    ifdef binary
        fread(&Nglobal_nodes, sizeof(int), 1, fp);
#    else
        fscanf(fp,"%d",&Nglobal_nodes);
#    endif
     fclose(fp);
  }
  Nglobal_nodes = AZ_gsum_int(Nglobal_nodes, proc_config);

  /*******************************************************************/
  /* initialize the list of global indices indicating which rows are */
  /* stored on which processor. NOTE: the list of global             */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly.                         */
  /*-----------------------------------------------------------------*/

  if (proc_config[AZ_N_procs] == 1) i = AZ_linear;
  else i = AZ_file;
  AZ_input_update("edge_partition",&Nlocal_edges, &global_edge_inds, 
		  proc_config, Nglobal_edges, num_PDE_eqns,i);

  /********************************************************************/
  /*                      Set up Ke_mat                               */
  /*  1) read it in using Aztec.                                      */
  /*  2) call AZ_transform to convert global column indices to local  */
  /*     indices and to set up Aztec's communication structure.       */
  /*  3) Stuff the arrays into an Aztec matrix.                       */
  /*------------------------------------------------------------------*/
	
  AZ_input_msr_matrix("Ke_mat.az", global_edge_inds, &Ke_val, &Ke_bindx, 
                      Nlocal_edges, proc_config);

  /* For now we add something to the diagonal instead of adding in the */
  /* real mass matrix.                                                 */

  for (i = 0; i < Nlocal_edges; i++) Ke_val[i] += 1.;

  AZ_transform_norowreordering(proc_config, &global_edge_externs, Ke_bindx,
               Ke_val, global_edge_inds, &reordered_glob_edges, 
               &reordered_edge_externs, &Ke_data_org, Nlocal_edges, 0, 0, 0, 
               &cpntr,	       AZ_MSR_MATRIX);
  update = global_edge_inds;
  update_index = reordered_glob_edges;

  Ke_mat = AZ_matrix_create( Nlocal_edges );
  AZ_set_MSR(Ke_mat, Ke_bindx, Ke_val, Ke_data_org, 0, NULL, AZ_LOCAL);
  Ke_mat->matrix_type    = Ke_data_org[AZ_matrix_type];
  Ke_data_org[AZ_N_rows] = Ke_data_org[AZ_N_internal] + 
                           Ke_data_org[AZ_N_border];

  /*******************************************************************/
  /* initialize the list of global indices indicating which rows are */
  /* stored on which processor. NOTE: the list of global             */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly.                         */
  /*-----------------------------------------------------------------*/

  if (proc_config[AZ_N_procs] == 1) i = AZ_linear;
  else i = AZ_file;
  AZ_input_update("node_partition",&Nlocal_nodes, &global_node_inds, 
		  proc_config, Nglobal_nodes, num_PDE_eqns,i);

  /********************************************************************/
  /*                      Set up Kn_mat                               */
  /*  1) read it in using Aztec.                                      */
  /*  2) call AZ_transform to convert global column indices to local  */
  /*     indices and to set up Aztec's communication structure.       */
  /*  3) Stuff the arrays into an Aztec matrix.                       */
  /*------------------------------------------------------------------*/

  AZ_input_msr_matrix("Kn_mat.az", global_node_inds, &Kn_val, &Kn_bindx, 
		      Nlocal_nodes, proc_config);

  AZ_transform_norowreordering(proc_config, &global_node_externs, Kn_bindx, 
               Kn_val, global_node_inds, &reordered_glob_nodes, 
	       &reordered_node_externs, &Kn_data_org, Nlocal_nodes, 0, 0, 0, 
	       &cpntr, AZ_MSR_MATRIX);
  Kn_mat = AZ_matrix_create( Nlocal_nodes );
  AZ_set_MSR(Kn_mat, Kn_bindx, Kn_val, Kn_data_org, 0, NULL, AZ_LOCAL);
  Kn_mat->matrix_type  = Kn_data_org[AZ_matrix_type];
  Kn_data_org[AZ_N_rows]  = Kn_data_org[AZ_N_internal] + Kn_data_org[AZ_N_border];

  /********************************************************************/
  /*                      Set up T_mat                                */
  /*  1) read it in using Aztec.                                      */
  /*  2) convert it to CSR                                            */
  /*      NOTE: we store the diagonal with the offdiagonals (wasting  */
  /*      the diagonal storage). This is so the manual mapping from   */
  /*      global to local indices does not get confusing.             */
  /*  3) Make it into an ML Operator                                  */
  /*     note: since ML_Operator_Create needs a communicator, I call  */
  /*     ML_Create() now to get one. I could have just created a      */
  /*     a communicator and done the ML_Create() later.               */
  /*------------------------------------------------------------------*/
	
  AZ_input_msr_matrix_nodiag("Tmat.az", global_edge_inds, &Tmat_val, 
                      &Tmat_bindx,  Nlocal_edges, proc_config);

  /* compress out any zeros which might occur due to empty rows  */

  lower = Tmat_bindx[0];
  Nnz = Tmat_bindx[Nlocal_edges];
  for (i = 0; i < Nlocal_edges; i++) {
    row_length = 0;
    for (j = lower; j < Tmat_bindx[i+1]; j++) {
      if (Tmat_val[j] != 0.0) row_length++;
    }
    lower = Tmat_bindx[i+1];
    Tmat_bindx[i+1] = Tmat_bindx[i] + row_length;
  }
  lower = Tmat_bindx[0];
  for (i = Tmat_bindx[0]; i < Nnz; i++ ) {
    if (Tmat_val[i] != 0.0) {
      Tmat_val[lower] = Tmat_val[i];
      Tmat_bindx[lower] = Tmat_bindx[i];
      lower++;
    }
  }


  /* Take the global MSR matrix and replace global column indices by */
  /* local column indices                                            */

  Nexterns = Kn_data_org[AZ_N_external];
  sorted_ext = (int *) ML_allocate(sizeof(int)*(Nexterns+1));
  map        = (int *) ML_allocate(sizeof(int)*(Nexterns+1));
  for (i = 0; i < Nexterns; i++) {
    sorted_ext[i] = global_node_externs[i];
    map[i] = reordered_node_externs[i];
  }
  AZ_sort(sorted_ext, Nexterns, map, NULL);

  ML_find_local_indices(Nlocal_nodes,Tmat_bindx, global_node_inds,sorted_ext,
                        Nexterns, map);

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
							  ML_CSR_MSRdata));
  csr_data->columns = Tmat_bindx;
  csr_data->values  = Tmat_val;
  ML_MSR2CSR(csr_data, Nlocal_edges, &Ncols);

  ML_Create(&ml_edges, N_levels);
  ML_Create(&ml_nodes, N_levels);

  Tmat = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Set_ApplyFuncData( Tmat, Nlocal_nodes, Nlocal_edges, 
                                  ML_EMPTY, csr_data, Nlocal_edges, NULL, 0);
  ML_Operator_Set_Getrow(Tmat, ML_EXTERNAL, Nlocal_edges, CSR_getrows);
  ML_Operator_Set_ApplyFunc(Tmat, ML_INTERNAL, CSR_matvec);

  start_time = AZ_second();

  options[AZ_scaling] = AZ_none;
  
  /********************************************************************/
  /* Convert the two discretization matrices (Ke_mat and Kn_mat) to   */
  /* ML matrices and stuff them into the ml grid hierarchies.         */
  /*------------------------------------------------------------------*/

  AZ_ML_Set_Amat(ml_edges, N_levels-1, Nlocal_edges, Nlocal_edges, Ke_mat, 
		 proc_config);
  AZ_ML_Set_Amat(ml_nodes, N_levels-1, Nlocal_nodes, Nlocal_nodes, Kn_mat, 
		 proc_config);

  ML_CommInfoOP_Clone(&(Tmat->getrow->pre_comm), ml_nodes->Amat[N_levels-1].getrow->pre_comm);
  sprintf(str,"P%d",proc_config[AZ_node]);
  /*
  ML_CommInfoOP_Print(Tmat->getrow->pre_comm, str);
  */

  /********************************************************************/
  /*                      Set up Tmat_trans                           */
  /*------------------------------------------------------------------*/

  Tmat_trans = ML_Operator_Create(ml_edges->comm);
  /*
  ML_Operator_Transpose(Tmat, Tmat_trans);
  */
  ML_Operator_Transpose_byrow(Tmat, Tmat_trans);

#ifdef ML_partition

  /* this code is meant to partition the matrices so that things can be */
  /* run in parallel later.                                             */
  /* It is meant to be run on only one processor.                       */
  fp2 = fopen("node_partition","w");

  ML_Operator_AmalgamateAndDropWeak(&(ml_nodes->Amat[N_levels-1]), num_PDE_eqns, 0.0);
  ML_Gen_Blocks_Metis(ml_edges, N_levels-1, &nblocks, &block_list);

  for (i = 0; i < nblocks; i++) {
     count = 0;
     for (j = 0; j < ml_nodes->Amat[N_levels-1].outvec_leng; j++) {
        if (block_list[j] == i) count++;
     }
     fprintf(fp2,"   %d\n",count*num_PDE_eqns);
     for (j = 0; j < ml_nodes->Amat[N_levels-1].outvec_leng; j++) {
        if (block_list[j] == i) {
           for (k = 0; k < num_PDE_eqns; k++)  fprintf(fp2,"%d\n",j*num_PDE_eqns+k);
        }
     }
  }
  fclose(fp2);
  ML_Operator_UnAmalgamateAndDropWeak(&(ml_nodes->Amat[N_levels-1]),num_PDE_eqns,0.0);

  /* Need to partition Tmat given a partitioning of Ke_mat */
  /* Here's how it should work:                          */
  /*    1) if both column points are within processor,   */
  /*       assign row to this processor.                 */
  /*    2) if columns are from different procs,          */
  /*        a) take lower if product of columns is even  */
  /*        b) take upper if product of columns is odd.  */

  csr_data = (struct ML_CSR_MSRdata *) Tmat->data;
  proc_assignment = (int *) AZ_allocate( Tmat->outvec_leng*sizeof(int));

  for (i = 0; i < Tmat->outvec_leng; i++) {
    itemp = (csr_data->rowptr)[i+1] - (csr_data->rowptr)[i];
    row_start = (csr_data->rowptr)[i];
    if (itemp == 2) {
      if ( (csr_data->values)[row_start+1] == 0.0) itemp--;
    }
    if (itemp > 0) {
      if ( (csr_data->values)[row_start] == 0.0) {
         itemp--;
         row_start++;
      }
    }
    if ( itemp > 2) 
      pr_error("Too many nonzeros per row in Tmat   %d\n", itemp);

    if (itemp == 1) {
      col1 = (csr_data->columns)[row_start];
      proc_assignment[i] = block_list[col1];
    }
    else if (itemp == 2) {
      col1 = (csr_data->columns)[row_start];
      col2 = (csr_data->columns)[row_start+1];
      p1   = block_list[col1];
      p2   = block_list[col2];
      if ( (col1*col2)%2 == 0) {
	if (p1 < p2) proc_assignment[i] = p1;
	else proc_assignment[i] = p2;
      }
      else {
	if (p2 < p1) proc_assignment[i] = p1;
	else proc_assignment[i] = p2;
      }
    }
    else proc_assignment[i] = -1;
  }
  pcounts = (int *) ML_allocate( sizeof(int)*nblocks);
  proc_id = (int *) ML_allocate( sizeof(int)*nblocks);
  for (i = 0; i < nblocks; i++) pcounts[i] = 0;
  for (i = 0; i < nblocks; i++) proc_id[i] = i;

  count = 0;
  for (i = 0; i < Tmat->outvec_leng ; i++) {
    if (proc_assignment[i] != -1) pcounts[proc_assignment[i]]++;
    else count++;
  }

  
   ML_az_sort(pcounts, nblocks, proc_id, NULL);

   i = 0; j = 0;
   while ( (i < nblocks) && (pcounts[i] < pcounts[nblocks-1])) {
     if ( proc_assignment[j] == -1) {
       proc_assignment[j] = proc_id[i];
       pcounts[i]++;
       if (pcounts[i] == pcounts[nblocks-1]) i++;
     }
     j++;
   }
   for (i = j; i < Tmat->outvec_leng; i++) {
     if (proc_assignment[i] == -1) proc_assignment[i] = i%nblocks;
   }

  fp2 = fopen("edge_partition","w");

  for (i = 0; i < nblocks; i++) {
     count = 0;
     for (j = 0; j < Tmat->outvec_leng; j++) {
        if (proc_assignment[j] == i) count++;
     }
     fprintf(fp2,"   %d\n",count*num_PDE_eqns);
     for (j = 0; j < Tmat->outvec_leng; j++) {
        if (proc_assignment[j] == i) {
           for (k = 0; k < num_PDE_eqns; k++)  fprintf(fp2,"%d\n",j*num_PDE_eqns+k);
        }
     }
  }
  fclose(fp2);


  exit(1);
#endif

  /********************************************************************/
  /* Set some ML parameters.                                          */
  /*------------------------------------------------------------------*/
	
  ML_Set_ResidualOutputFrequency(ml_edges, context->output);
  ML_Set_Tolerance(ml_edges, context->tol);
  ML_Aggregate_Create( &ag );
  if (ML_strcmp(context->agg_coarsen_scheme,"Mis") == 0) {
     ML_Aggregate_Set_CoarsenScheme_MIS(ag);
  }
  else if (ML_strcmp(context->agg_coarsen_scheme,"Uncoupled") == 0) {
     ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag);
  }
  else if (ML_strcmp(context->agg_coarsen_scheme,"Coupled") == 0) {
     ML_Aggregate_Set_CoarsenScheme_Coupled(ag);
  }
  else {
     printf("ML: Unknown aggregation scheme %s\n",context->agg_coarsen_scheme);
  }
  ML_Aggregate_Set_DampingFactor(ag, context->agg_damping);
  ML_Aggregate_Set_MaxCoarseSize( ag, context->maxcoarsesize);
  ML_Aggregate_Set_Threshold(ag, context->agg_thresh);

  if (ML_strcmp(context->agg_spectral_norm,"Calc") == 0) {
     ML_Aggregate_Set_SpectralNormScheme_Calc(ag);
  }
  else if (ML_strcmp(context->agg_spectral_norm,"Anorm") == 0) {
     ML_Aggregate_Set_SpectralNormScheme_Anorm(ag);
  }
  else {
     printf("ML: Unknown spectral norm scheme %s\n",context->agg_spectral_norm);
  }

  /* read in the rigid body modes */

   Nrigid = 0;
   if (proc_config[AZ_node] == 0) {
      sprintf(filename,"rigid_body_mode%d",Nrigid+1);
      while( (fp = fopen(filename,"r")) != NULL) {
          fclose(fp);
          Nrigid++;
          sprintf(filename,"rigid_body_mode%d",Nrigid+1);
      }
    }
    Nrigid = AZ_gsum_int(Nrigid,proc_config);

    if (Nrigid != 0) {
       rigid = (double *) ML_allocate( sizeof(double)*Nrigid*(Nlocal_edges+1) );
       if (rigid == NULL) {
          printf("Error: Not enough space for rigid body modes\n");
       }
    }

   /* Set rhs */

   rhs=(double *)malloc(Nlocal_edges*sizeof(double));
 
   fp = fopen("rhsfile","r");
   if (fp == NULL) {
      if (proc_config[AZ_node] == 0) printf("taking zero vector for rhs\n");
      for (i = 0; i < Nlocal_edges; i++) rhs[i] = 0.0;
   }
   else {
      fclose(fp);
      if (proc_config[AZ_node] == 0) printf("reading rhs from a file\n");
      AZ_input_msr_matrix("rhsfile", global_edge_inds, &rhs, &garbage, 
			  Nlocal_edges, proc_config);
   }
   AZ_reorder_vec(rhs, Ke_data_org, reordered_glob_edges, NULL);

   for (i = 0; i < Nrigid; i++) {
      sprintf(filename,"rigid_body_mode%d",i+1);
      AZ_input_msr_matrix(filename, global_edge_inds, &mode, &garbage, Nlocal_edges, 
                          proc_config);

   /*
    *  Rescale matrix/rigid body modes and checking 
    *
       AZ_sym_rescale_sl(mode, Ke_mat->Ke_data_org, options, proc_config, scaling);
       Ke_mat->matvec(mode, rigid, Ke_mat, proc_config);
       for (j = 0; j < Nlocal_edges; j++) printf("this is %d %e\n",j,rigid[j]);
    */

    for (j = 0; j < i; j++) {
       alpha = -AZ_gdot(Nlocal_edges, mode, &(rigid[j*Nlocal_edges]), proc_config)/
                  AZ_gdot(Nlocal_edges, &(rigid[j*Nlocal_edges]), &(rigid[j*Nlocal_edges]), 
                               proc_config);
       daxpy_(&Nlocal_edges, &alpha,  &(rigid[j*Nlocal_edges]),  &one, mode, &one);
    }
   
    /* rhs orthogonalization */

    alpha = -AZ_gdot(Nlocal_edges, mode, rhs, proc_config)/
                    AZ_gdot(Nlocal_edges, mode, mode, proc_config);
    daxpy_(&Nlocal_edges, &alpha,  mode,  &one, rhs, &one);

    for (j = 0; j < Nlocal_edges; j++) rigid[i*Nlocal_edges+j] = mode[j];
    free(mode);
    free(garbage);
  }

  for (j = 0; j < Nrigid; j++) {
     alpha = -AZ_gdot(Nlocal_edges, rhs, &(rigid[j*Nlocal_edges]), proc_config)/
              AZ_gdot(Nlocal_edges, &(rigid[j*Nlocal_edges]), &(rigid[j*Nlocal_edges]), 
                      proc_config);
     daxpy_(&Nlocal_edges, &alpha,  &(rigid[j*Nlocal_edges]),  &one, rhs, &one);
  }

  if (Nrigid != 0) {
     ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Nrigid, rigid, Nlocal_edges);
  }

  /********************************************************************/
  /* Set up the operators corresponding to regular unsmoothed         */
  /* aggregation on the nodal matrix.                                 */
  /*------------------------------------------------------------------*/

  coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml_nodes, N_levels-1, 
                                            ML_DECREASING, ag);

#ifdef out
  /*
  printf("ready to built T on the coarse grids\n");
  exit(1);
  */


  /********************************************************************/
  /*                 Build T on the coarse grid.                      */
  /* I'm not completely sure that the space I allocate is sufficient. */
  /* I believe we made a quick calculation and decided that the number*/
  /* of nonzeros in Kn_coarse was an upper bound for the number of    */
  /* nonzeros in T_coarse.                                            */
  /*------------------------------------------------------------------*/

  counter = 0;
  
  Kn_coarse = &(ml_nodes->Amat[N_levels-2]);

  Tcoarse_bindx = (int    *) ML_allocate(Kn_coarse->N_nonzeros*sizeof(int));
  Tcoarse_val   = (double *) ML_allocate(Kn_coarse->N_nonzeros*sizeof(double));
  Tcoarse_rowptr= (int    *) ML_allocate(Kn_coarse->N_nonzeros*sizeof(int));
  Tcoarse_rowptr[0] = 0;
  nz_ptr = 0;
  for (i = 0; i < Kn_coarse->outvec_leng; i++) {
     ML_get_matrix_row(Kn_coarse,1, &i,&allocated,&bindx,&val,&row_length, 0);
     ML_az_sort(bindx, row_length, NULL, NULL);
     for (j = 0; j < row_length; j++) {
       if (bindx[j] > i) {
         Tcoarse_bindx[nz_ptr]  =  bindx[j];
         Tcoarse_val[nz_ptr++]  =  1.;
         Tcoarse_bindx[nz_ptr]  =  i;
         Tcoarse_val[nz_ptr++]  = -1.;
         Tcoarse_rowptr[counter+1] = nz_ptr;
         counter++;
       }
     }
  }

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
							  ML_CSR_MSRdata));
  csr_data->columns = Tcoarse_bindx;
  csr_data->values  = Tcoarse_val;
  csr_data->rowptr  = Tcoarse_rowptr;

  Tcoarse = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Set_ApplyFuncData( Tcoarse, Kn_coarse->outvec_leng, counter, 
                                  ML_EMPTY, csr_data, counter, NULL, 0);
  ML_Operator_Set_Getrow(Tcoarse, ML_EXTERNAL, counter, CSR_getrows);
  ML_Operator_Set_ApplyFunc(Tcoarse, ML_INTERNAL, CSR_matvec);


  /********************************************************************/
  /* Fix P and R so that they are not normalized. This is so that we  */
  /* can use a matrix triple product combined with post-processing to */
  /* generate Pe. The general idea is that the matrix                 */
  /*                        T_h P_n T_H^*                             */
  /* is almost Pe. If we make sure that P_n contains 1's and -1's, the*/
  /* matrix triple product will yield a matrix with +/- 1 and +/- 2's.*/
  /* If we remove all the 1's and divide the 2's by 2. we arrive at Pe*/
  /*------------------------------------------------------------------*/

  Pn_coarse = &(ml_nodes->Pmat[N_levels-2]);
  csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
  for (i = 0; i < Pn_coarse->outvec_leng; i++) {
    if (csr_data->values[i] < 0) csr_data->values[i] = -1.;
    else if (csr_data->values[i] > 0) csr_data->values[i] = 1.;
  }

  Rn_coarse = &(ml_nodes->Rmat[N_levels-1]);
  csr_data = (struct ML_CSR_MSRdata *) Rn_coarse->data;
  for (i = 0; i < csr_data->rowptr[Rn_coarse->outvec_leng]; i++) {
    if (csr_data->values[i] < 0.) csr_data->values[i] = -1.;
    else if (csr_data->values[i] > 0) csr_data->values[i] = 1.;
  }

  /********************************************************************/
  /* Create Tcoarse_trans.                                            */
  /*------------------------------------------------------------------*/

  Tcoarse_trans = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Transpose(Tcoarse, Tcoarse_trans);

  /********************************************************************/
  /* Here is some code that might work someday to generate Pe without */
  /* doing a matrix triple product.                                   */
  /*------------------------------------------------------------------*/
  /*
  csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
  for (i = 0; i < Tmat->outvec_leng; i++) {
    ML_get_matrix_row(Tmat, 1, &i, &allocated, &bindx, &val, &row_length, 0);
    if (row_length == 2) {
     agg1 = csr_data->columns[bindx[0]];
     agg2 = csr_data->columns[bindx[1]];
     printf("%d: agg1 and agg2 %d %d   | %d %d\n",i,agg1,agg2,bindx[0],bindx[1]);
    
     agg1 = aggr_num(bindx[0], PnMat);
     agg2 = aggr_num(bindx[1], PnMat);
     if (agg1 != agg2) {
        printf("Pe(%d,%d) = something;\n",i,thecolumn);

        To get the column suppose we store for each row in Kn_coarse the counter above.
        Thus, stored_counter[i] indicates the lowest coarse grid edge that is 
        assicated with coarse grid node i.
        To get the edge number we compute 
               min_agg = min(agg1,agg2), 
               max_agg = max(agg1,agg2);
        The column that we are looking for is given by 
          stored_counter[min_agg] + number of nonzeros in row min_agg of Kn_coarse
                                    that are greater than min_agg and less than max_agg.


        'something' is either +1 or -1. To get the sign look at
        Tcoarse(thecolumn,) and Tmat(i,). My guess is that it is related to 
        whether agg1 is greater than agg2 (assuming that val[0] =1)
         }
  }
}
  */


  /********************************************************************/
  /* Matrix triple product for Pe where result is stuffed into MG     */
  /* hierarchy of the edge system.                                    */
  /*------------------------------------------------------------------*/

  ML_rap(Tmat, &(ml_nodes->Pmat[N_levels-2]), Tcoarse_trans, 
	 &(ml_edges->Pmat[N_levels-2]),ML_CSR_MATRIX);

  Pe = &(ml_edges->Pmat[N_levels-2]);
  csr_data = (struct ML_CSR_MSRdata *) Pe->data;

  /********************************************************************/
  /* weed out the 1's and divide the 2's by -2.                       */
  /* Note: my signs seem flipped from what we thought the triple      */
  /* matrix product would do. This could be a function of the Pnmat   */
  /* normalization.                                                   */
  /* Finally, convert the matrix to CSR and strip out the zeros.      */
  /* and make sure that it has the level transfer stuff to be in the  */
  /* MG grid hierarchy.                                               */
  /*------------------------------------------------------------------*/

  for (j = 0; j < csr_data->rowptr[Pe->outvec_leng] ; j++) {
    if (csr_data->values[j] == 2) csr_data->values[j] = -1;
    else if (csr_data->values[j] == -2) csr_data->values[j] = 1;
    else if (csr_data->values[j] == -1) csr_data->values[j] = 0;
    else if (csr_data->values[j] ==  1) csr_data->values[j] = 0;
    else if (csr_data->values[j] != 0.0) printf("huh\n");
  }

  /*******************************************************************/
  /* weed out zeros in Pe.                                           */
  /*-----------------------------------------------------------------*/

  lower = csr_data->rowptr[0];
  nz_ptr = 0;
  for (i = 0; i < Pe->outvec_leng; i++) {
    for (j = lower; j < csr_data->rowptr[i+1]; j++) {
      if (csr_data->values[j] != 0.0) nz_ptr++;
    }
    lower = csr_data->rowptr[i+1];
    csr_data->rowptr[i+1] = nz_ptr;
  }
  nz_ptr = 0;
  for (i = 0; i < lower; i++) {
    if (csr_data->values[i] != 0.) {
      csr_data->values[nz_ptr] = csr_data->values[i];
      csr_data->columns[nz_ptr] = csr_data->columns[i];
      nz_ptr++;
    }
  }

  Pe->getrow->external = CSR_getrows;
  Pe->getrow->internal = NULL;
  Pe->getrow->ML_id    = ML_EXTERNAL;
  Pe->matvec->internal = CSR_matvec;
  Pe->matvec->external = NULL;
  Pe->matvec->ML_id = ML_INTERNAL;
  ML_Operator_Set_1Levels(&(ml_edges->Pmat[N_levels-2]),
			  &(ml_edges->SingleLevel[N_levels-2]), 
			  &(ml_edges->SingleLevel[N_levels-1]));

  ML_Gen_Restrictor_TransP(ml_edges, N_levels-1, N_levels-2);
  ML_Gen_AmatrixRAP(ml_edges, N_levels-1, N_levels-2);

  coarsest_level = N_levels - coarsest_level;
  if ( proc_config[AZ_node] == 0 )
	printf("Coarse level = %d \n", coarsest_level);
	
  /* set up smoothers */
  blocks = (int *) ML_allocate(sizeof(int)*Nlocal_edges);
	
  for (level = N_levels-1; level > coarsest_level; level--) {

      num_PDE_eqns = ml_edges->Amat[level].num_PDEs;
      if (proc_config[AZ_node]==0) printf("block size = %d\n",num_PDE_eqns);
		
     /*  Sparse approximate inverse smoother that acutally does both */
     /*  pre and post smoothing.                                     */

     if (ML_strcmp(context->smoother,"Parasails") == 0) {
        ML_Gen_Smoother_ParaSails(ml_edges , level, ML_PRESMOOTHER, nsmooth, 
                                parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
     }

     else if (ML_strcmp(context->smoother,"Hiptmair") == 0) {
       printf("only doing pre smoothing\n");
         ML_Gen_Smoother_Hiptmair(ml_edges , level, ML_PRESMOOTHER, nsmooth,1.,Tmat_trans);
     }
     /* This is the symmetric Gauss-Seidel smoothing that we usually use. */
     /* In parallel, it is not a true Gauss-Seidel in that each processor */
     /* does a Gauss-Seidel on its local submatrix independent of the     */
     /* other processors.                                                 */

     else if (ML_strcmp(context->smoother,"GaussSeidel") == 0) {
       ML_Gen_Smoother_GaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.);
     }
     else if (ML_strcmp(context->smoother,"SymGaussSeidel") == 0) {
       ML_Gen_Smoother_SymGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.);
     }
     else if (ML_strcmp(context->smoother,"BlockGaussSeidel") == 0) {
       ML_Gen_Smoother_BlockGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.,
					 num_PDE_eqns);
     }
     else if (ML_strcmp(context->smoother,"Aggregate") == 0) {
         ML_Gen_Blocks_Aggregates(ag, level, &nblocks, &blocks);
         ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.,
                                        nblocks, blocks);
     }

     /* This is a true Gauss Seidel in parallel. This seems to work for  */
     /* elasticity problems.  However, I don't believe that this is very */
     /* efficient in parallel.                                           */       
     /*
      nblocks = ml_edges->Amat[level].invec_leng;
      for (i =0; i < nblocks; i++) blocks[i] = i;
      ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml_edges , level, ML_PRESMOOTHER,
                                                  nsmooth, 1., nblocks, blocks);
      ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml_edges, level, ML_POSTSMOOTHER,
                                                  nsmooth, 1., nblocks, blocks);
     */

     /* Jacobi Smoothing                                                 */

     else if (ML_strcmp(context->smoother,"Jacobi") == 0) {
        ML_Gen_Smoother_Jacobi(ml_edges , level, ML_PRESMOOTHER, nsmooth,.4);
        ML_Gen_Smoother_Jacobi(ml_edges , level, ML_POSTSMOOTHER, nsmooth,.4);
     }

     /*  This does a block Gauss-Seidel (not true GS in parallel)        */
     /*  where each processor has 'nblocks' blocks.                      */
     /* */

     else if (ML_strcmp(context->smoother,"Metis") == 0) {
         nblocks = 250;
         ML_Gen_Blocks_Metis(ml_edges, level, &nblocks, &blocks);
         ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.,
                                        nblocks, blocks);
     }
     else {
         printf("unknown smoother %s\n",context->smoother);
         exit(1);
     }
   }
   nsmooth   = context->coarse_its;
   /*  Sparse approximate inverse smoother that acutally does both */
   /*  pre and post smoothing.                                     */

   if (ML_strcmp(context->coarse_solve,"Parasails") == 0) {
        ML_Gen_Smoother_ParaSails(ml_edges , coarsest_level, ML_PRESMOOTHER, nsmooth, 
                                parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
   }

   else if (ML_strcmp(context->coarse_solve,"Hiptmair") == 0) {
       printf("only doing pre smoothing\n");
       ML_Gen_Smoother_Hiptmair(ml_edges , coarsest_level, ML_PRESMOOTHER, nsmooth,1.,Tmat_trans);
   }
   else if (ML_strcmp(context->coarse_solve,"GaussSeidel") == 0) {
       ML_Gen_Smoother_GaussSeidel(ml_edges , coarsest_level, ML_BOTH, nsmooth,1.);
   }
   else if (ML_strcmp(context->coarse_solve,"SymGaussSeidel") == 0) {
     printf("this is GS with %d\n",nsmooth);
       ML_Gen_Smoother_SymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, nsmooth,1.);
   }
   else if (ML_strcmp(context->coarse_solve,"BlockGaussSeidel") == 0) {
       ML_Gen_Smoother_BlockGaussSeidel(ml_edges, coarsest_level, ML_BOTH, nsmooth,1.,
					num_PDE_eqns);
   }
   else if (ML_strcmp(context->coarse_solve,"Aggregate") == 0) {
         ML_Gen_Blocks_Aggregates(ag, coarsest_level, &nblocks, &blocks);
         ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, 
                                        nsmooth,1., nblocks, blocks);
   }
   else if (ML_strcmp(context->coarse_solve,"Jacobi") == 0) {
     printf("number of steps is %d\n",nsmooth);
        ML_Gen_Smoother_Jacobi(ml_edges , coarsest_level, ML_BOTH, nsmooth,.5);
   }
   else if (ML_strcmp(context->coarse_solve,"Metis") == 0) {
         nblocks = 250;
         ML_Gen_Blocks_Metis(ml_edges, coarsest_level, &nblocks, &blocks);
         ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, 
                                              nsmooth,1., nblocks, blocks);
   }
   else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0) {
      ML_Gen_CoarseSolverSuperLU( ml_edges, coarsest_level);
   }
   else {
         printf("unknown coarse grid solver %s\n",context->coarse_solve);
         exit(1);
   }
#endif		
   ML_Gen_Solver(ml_edges, ML_MGV, N_levels-1, coarsest_level); 
   AZ_defaults(options, params);
	
   if (ML_strcmp(context->krylov,"Cg") == 0) {
      options[AZ_solver]   = AZ_cg;
   }
   else if (ML_strcmp(context->krylov,"Bicgstab") == 0) {
      options[AZ_solver]   = AZ_bicgstab;
   }
   else if (ML_strcmp(context->krylov,"Tfqmr") == 0) {
      options[AZ_solver]   = AZ_tfqmr;
   }
   else if (ML_strcmp(context->krylov,"Gmres") == 0) {
      options[AZ_solver]   = AZ_gmres;
   }
   else {
      printf("unknown krylov method %s\n",context->krylov);
   }
   options[AZ_scaling]  = AZ_none;
   options[AZ_precond]  = AZ_user_precond;
   options[AZ_conv]     = AZ_noscaled;
   options[AZ_output]   = 1;
   options[AZ_max_iter] = 300;
   options[AZ_poly_ord] = 5;
   options[AZ_kspace]   = 130;
   params[AZ_tol]       = context->tol;
   options[AZ_output]   = context->output;
	
   AZ_set_ML_preconditioner(&Pmat, Ke_mat, ml_edges, options); 
   setup_time = AZ_second() - start_time;
	
   xxx = (double *) malloc( Nlocal_edges*sizeof(double));

   for (iii = 0; iii < Nlocal_edges; iii++) xxx[iii] = 0.0; 
	

   /* Set xxx */

   /*
   fp = fopen("initguessfile","r");
   if (fp != NULL) {
      fclose(fp);
      if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
      AZ_input_msr_matrix("initguessfile", global_edge_inds, &xxx, &garbage, Nlocal_edges, 
                          proc_config);
      options[AZ_conv] = AZ_expected_values;
   }
   else if (proc_config[AZ_node]== 0) printf("taking 0 initial guess \n");

   AZ_reorder_vec(xxx, Ke_data_org, reordered_glob_edges, NULL);
   */
   printf("putting in an node based xxxx\n");
   fp = fopen("initguessfile","r");
   if (fp != NULL) {
      fclose(fp);
      if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
      AZ_input_msr_matrix("initguessfile", global_edge_inds, &xxx, &garbage, Nlocal_edges, 
                          proc_config);
      options[AZ_conv] = AZ_expected_values;
   }
   else if (proc_config[AZ_node]== 0) printf("taking 0 initial guess \n");

   AZ_reorder_vec(xxx, Kn_data_org, reordered_glob_nodes, NULL);



   fp = fopen("AZ_no_multilevel.dat","r");
   scaling = AZ_scaling_create();
   start_time = AZ_second();
   if (fp != NULL) {
      fclose(fp);
      options[AZ_precond] = AZ_none;
      options[AZ_scaling] = AZ_sym_diag;
      options[AZ_ignore_scaling] = AZ_TRUE;

      options[AZ_keep_info] = 1;
      
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, NULL, scaling); 

/*
      options[AZ_pre_calc] = AZ_reuse;
      options[AZ_conv] = AZ_expected_values;
      if (proc_config[AZ_node] == 0) 
              printf("\n-------- Second solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, NULL, scaling); 
      if (proc_config[AZ_node] == 0) 
              printf("\n-------- Third solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, NULL, scaling); 
*/
   }
   else {
      options[AZ_keep_info] = 1;
      options[AZ_conv] = AZ_noscaled;
      options[AZ_output] = 1;
      dtemp = sqrt(ML_gdot(Nlocal_nodes, xxx, xxx, ml_edges->comm));
      printf("norm of x_0 = %e\n",dtemp);

      /*
      ML_Operator_Apply(Tmat, Tmat->invec_leng, xxx, Tmat->outvec_leng,rhs);
      dtemp = sqrt(ML_gdot(Nlocal_edges, rhs, rhs, ml_edges->comm));
      printf("norm of T x_0 = %e\n",dtemp);
      */
      /*
      for (i = 0; i < Nlocal_edges; i++) 
	printf("x(%d)=%10.5e;\n",global_node_inds[i]+1,xxx[reordered_glob_nodes[i]]);
      */
      /*
      ML_Operator_Apply(&(ml_edges->Amat[N_levels-1]),
			ml_edges->Amat[N_levels-1].invec_leng,xxx,
			ml_edges->Amat[N_levels-1].outvec_leng,rhs);
      */

      ML_Operator_Apply(Tmat_trans, Tmat_trans->invec_leng, xxx,
			Tmat_trans->outvec_leng,rhs);
      dtemp = sqrt(ML_gdot(Tmat_trans->outvec_leng, rhs, rhs, ml_edges->comm));
      printf("norm of T^t x_0 = %e\n",dtemp);

      /*
      for (i = 0; i < Tmat_trans->invec_leng; i++) 
	printf("x(%d)=%10.5e; local ind = %d, %d\n",global_edge_inds[i]+1,rhs[reordered_glob_edges[i]],i,proc_config[AZ_node]);
      */
      /*
      for (i = 0; i < Tmat_trans->outvec_leng; i++) 
	printf("Av(%d)=%10.5e; %d\n",global_node_inds[i]+1,xxx[reordered_glob_nodes[i]],i);
      */
      fflush(stdout);


      exit(1);
      /*
      options[AZ_precond] = AZ_none;
     
            Ke_mat->matvec(xxx, rhs, Ke_mat, proc_config);
      for (i = 0; i < Nlocal_edges; i++) printf("%7d     %7d %20.15e %20.15e\n",i+1,i+1,xxx[i],rhs[i]);
      printf("huhhhh %e\n",Ke_mat->val[0]);
     
      printf("the norm is %e\n",sqrt(AZ_gdot(Nlocal_edges, xxx, xxx, proc_config)));
      */


      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, scaling); 
      options[AZ_pre_calc] = AZ_reuse;
      options[AZ_conv] = AZ_expected_values;
/*
      if (proc_config[AZ_node] == 0) 
              printf("\n-------- Second solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, scaling); 
      if (proc_config[AZ_node] == 0) 
              printf("\n-------- Third solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, scaling); 
*/
   }
   solve_time = AZ_second() - start_time;

   if (proc_config[AZ_node] == 0) 
      printf("Solve time = %e, MG Setup time = %e\n", solve_time, setup_time);

   ML_Aggregate_Destroy(&ag);
   ML_Destroy(&ml_edges);
   ML_Destroy(&ml_nodes);
   AZ_free((void *) Ke_mat->data_org);
   AZ_free((void *) Ke_mat->val);
   AZ_free((void *) Ke_mat->bindx);
   AZ_free((void *) global_edge_inds);
   AZ_free((void *) global_edge_externs);
   AZ_free((void *) reordered_edge_externs);
   AZ_free((void *) reordered_glob_edges);
   AZ_scaling_destroy(&scaling);
   if (Ke_mat  != NULL) AZ_matrix_destroy(&Ke_mat);
   if (Pmat  != NULL) AZ_precond_destroy(&Pmat);
   free(xxx);
   free(rhs);


#ifdef ML_MPI
  MPI_Finalize();
#endif
	
  return 0;
	
}

/************************************************************************/
/* Convert the data in csr_data from an MSR matrix to a CSR matrix.     */
/* Also, return the largest column number encountered.                  */
/*----------------------------------------------------------------------*/

int ML_MSR2CSR(struct ML_CSR_MSRdata *csr_data, int Nrows, int *Ncolumns)
{
  int  *row_ptr, *Tmat_bindx, i, j, lower, upper, nz_ptr, Ncols;
  double *Tmat_val, *diag;

  row_ptr = (int *) ML_allocate((Nrows+1)*sizeof(int));  
  csr_data->rowptr  = row_ptr;
  Tmat_bindx = csr_data->columns;
  Tmat_val   = csr_data->values;

  diag    = (double *) ML_allocate(Nrows*sizeof(double));
  for (i = 0; i <= Nrows; i++) row_ptr[i] = Tmat_bindx[i];
  for (i = 0; i < Nrows; i++) diag[i] = Tmat_val[i];
  lower = row_ptr[0];
  row_ptr[0] = 0;
  nz_ptr = 0;
  Ncols = -1;
  for (i = 0; i < Nrows; i++) {

    upper = row_ptr[i+1];
    if ( diag[i] != 0.0) { 
      Tmat_bindx[nz_ptr] = i;
      Tmat_val[nz_ptr++] = diag[i];
      if (i > Ncols) Ncols = i;
    }
    for (j = lower; j < upper; j++) {
      if (Tmat_val[j] != 0.0) {
	Tmat_bindx[nz_ptr] = Tmat_bindx[j];
	Tmat_val[nz_ptr++] = Tmat_val[j];
        if (Tmat_bindx[j] > Ncols) Ncols = Tmat_bindx[j];
      }
    }
    row_ptr[i+1] = nz_ptr;
    lower = upper;
  }
  Ncols++;
  ML_free(diag);
  *Ncolumns = Ncols;
  return 0;
}
/*****************************************************************************/
/*****************************************************************************/

void AZ_transform_norowreordering(int proc_config[], int *external[], int bindx[], double val[],
                  int update[], int *update_index[], int *extern_index[],
                  int *data_org[], int N_update, int indx[], int bnptr[],
                  int rnptr[], int *cnptr[], int mat_type)

/*******************************************************************************

  Convert a global distributed matrix to a parallel local distributed matrix.
  This includes the following steps:
      1) reorder matrix rows so that all the rows corresponding to internal
         points preceed all the rows corresponding to border points.
      2) replace global indicies by local indicies.
      3) make a list of the external unknowns and store them in external[].
      4) make a list of processors which update each external unknown and store
         this list in extern_proc where extern_proc[i] is the processor that
         updates external[i].
      5) make 2 arrays (update_index[], extern_index[]) which define the mapping
         between global and local indicies. In particular, global index
         update[i] corresponds to the locally numbered variable update_index[i]
         and global index external[i] corresponds to the locally numbered
         variable extern_index[i].
      6) Initialize all the quanities in data_org[] to their appropriate values
         (so that communication can properly occur).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc_config:     Processor information corresponding to:
                      proc_config[AZ_node] = name of this processor
                      proc_config[AZ_N_procs] = # of processors used

  external:        On output, list of external blocks

  update:          On input, list of pts to be updated on this node

  update_index,    On output, ordering of update and external locally on this
  extern_index:    processor. For example  'update_index[i]' gives the index
                   location of the block which has the global index 'update[i]'.

  data_org:        On output, indicates how the data is set out on this node.
                   For example, data_org[] contains information on how many
                   unknowns are internal, external and border unknowns as well
                   as which points need to be communicated. See Aztec User's
                   guide for more details.

  N_update         Number of points to be updated on this node.

  val,bindx        On input, global distributed matrix (MSR or VBR) arrays
  indx, bnptr,     holding matrix values. On output, local renumbered matrix
  rnptr, cnptr:    (DMSR or DMVBR).

  mat_type:        Type of matrix (AZ_MSR_MATRIX or AZ_VBR_MATRIX).

*******************************************************************************/

{
  int        i, ii, j;
  static int mat_name = 1;

  int         N_extern;   /* Number of pts needed by this processor for
                             matrix-vector multiply but not updated by this
                             processor.  */
  int         N_internal, /* Number of pts which can be updated without
                             communication */
              N_border;   /* Number of pts to be updated requiring communication
                           */
  int        *extern_proc;
  int        *tcnptr = NULL;

  AZ__MPI_comm_space_ok();
#ifdef AZ_MPI
  if ( proc_config[AZ_Comm_Set] != AZ_Done_by_User) {
      printf("Error: Communicator not set. Use AZ_set_comm()\n");
      printf("       (e.g. AZ_set_comm(proc_config,MPI_COMM_WORLD)).\n");
      exit(1);
  }
#endif

  /*
   * Compute the external points and change the global indices to
   * local indices. That is,
   *   On input:                        On output:
   *      bindx[k] = update[j]      ==>   bindx[k] = j
   *      bindx[k] = external[j]    ==>   bindx[k] = j + N_update
   */

  AZ_find_local_indices(N_update, bindx, update, external, &N_extern, mat_type,
                        bnptr);

  /* compute the cnptr array for VBR matrices */

  if (mat_type == AZ_VBR_MATRIX) {
    if (!AZ_using_fortran) {
      *cnptr = (int *) AZ_allocate((N_update + N_extern + 1)*sizeof(int));
      if (*cnptr == NULL) {
         printf("Out of memory in AZ_transform\n");
         exit(1);
      }
    }

    tcnptr = *cnptr;
    for (i = 0 ; i < N_update + N_extern + 1; i++) tcnptr[i] = 0;

    for (i = 0; i < N_update; i++) tcnptr[i] = rnptr[i+1] - rnptr[i];

    for (i = 0; i < N_update; i++) {
      for (j = bnptr[i]; j < bnptr[i+1]; j++) {
        ii = bindx[j];

        if ((ii >= N_update) && ( tcnptr[ii] == 0)) {
          tcnptr[ii] = (indx[j+1]-indx[j]) / (rnptr[i+1]-rnptr[i]);
        }
      }
    }

    AZ_convert_values_to_ptrs(tcnptr, N_update + N_extern, 0);
  }

  /*
   * Read or compute (and sort) the processor numbers of the processors which
   * update the external points.
   */

  i                = AZ_using_fortran;
  AZ_using_fortran = AZ_FALSE;

  AZ_find_procs_for_externs(N_update, update, *external, N_extern, proc_config,
                            &extern_proc);
  AZ_using_fortran = i;

  /*
   * Determine a new ordering for the points:
   *    a) lowest numbers for internal points,
   *    b) next lowest numbers for border points
   *    c) highest nubers for the external points
   *       NOTE: external points updated by the same processor are consecutively
   *             ordered.
   */

  if (!AZ_using_fortran) {
    *update_index = (int *) AZ_allocate((N_update + 1)*sizeof(int));
    *extern_index = (int *) AZ_allocate((N_extern + 1)*sizeof(int));
  }

  if (*extern_index == NULL)  {
    (void) fprintf(stderr,
                   "Error: Not enough space in main() for extern_index[]\n");
    exit(1);
  }

  AZ_order_ele(*update_index, *extern_index, &N_internal, &N_border, N_update,
               bnptr, bindx, extern_proc, N_extern, AZ_EXTERNS, mat_type);

  /*
   * Permute the matrix using the new ordering.  IMPORTANT: This routine assumes
   * that update_index[] contains 2 sequencies that are ordered but
   * intertwined. See AZ_reorder_matrix().
   */

  AZ_reorder_matrix(N_update, bindx, val, *update_index, *extern_index,
                    indx, rnptr, bnptr, N_extern, tcnptr, AZ_EXTERNS,mat_type);

  /*
   * Initialize 'data_org' so that local information can be exchanged to update
   * the external points.
   */

  AZ_set_message_info(N_extern, *extern_index, N_update, *external, extern_proc,
                      update, *update_index, proc_config, tcnptr, data_org,
                      mat_type);

  (*data_org)[AZ_name]       = mat_name;
  (*data_org)[AZ_N_int_blk]  = N_internal;
  (*data_org)[AZ_N_bord_blk] = N_border;
  (*data_org)[AZ_N_ext_blk]  = N_extern;

  if (mat_type == AZ_VBR_MATRIX) {
    (*data_org)[AZ_N_internal] = rnptr[N_internal];
    (*data_org)[AZ_N_external] = tcnptr[N_update + N_extern] - rnptr[N_update];
    (*data_org)[AZ_N_border]   = rnptr[N_update] - rnptr[N_internal];
  }

  else {
    (*data_org)[AZ_N_internal] = N_internal;
    (*data_org)[AZ_N_external] = N_extern;
    (*data_org)[AZ_N_border]   = N_border;
  }

  mat_name++;
  AZ_free(extern_proc);

} /* AZ_transform */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void AZ_input_msr_matrix_nodiag(char datafile[], int update[], double **val, int **bindx, 
												 int N_update, int proc_config[])

/*******************************************************************************

Exactly the same as AZ_read_msr_matrix except it reads that data in from a 
file specified by the input argument datafile instead from a file called
.data

*******************************************************************************/

{

  /* local variables */

  int    nz_ptr;
  char  *str;
  int    i,j,k, jjj;
  int    current;
  int    st, pnode;
  int    temp, *lil;
  double dtemp;
  int   *requests;
  int    total;
  FILE  *dfp = NULL;
  int    proc, nprocs;
  int    type, type2;
  unsigned int buf_len = 1000;
  int    msr_len;
  int    count = 0;
  int    kkk, m_one = -1, need_request = 0;
  int    column0 = 0;
  MPI_AZRequest request;  /* Message handle */
  int    totalN;

  char   *tchar;
  extern int AZ_sys_msg_type;

  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;
  type2           = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  totalN = AZ_gsum_int(N_update, proc_config);
  str    = (char *) AZ_allocate((buf_len+1)*sizeof(char));
  if (str == NULL) {
    printf("ERROR: NOT enough dynamic memory in AZ_input_msr_matrix_nodiag\n");
    exit(-1);
  }
  msr_len = 8*N_update+2;
  if (!AZ_using_fortran) {
    *bindx = (int *)    AZ_allocate(msr_len*sizeof(int));
    *val   = (double *) AZ_allocate(msr_len*sizeof(double));
  }

  if (*val == NULL) {
    (void) fprintf(stderr,
                   "ERROR: NOT enough dynamic memory in AZ_input_msr_matrix_nodiag\n");
    exit(-1);
  }
  if (!AZ_using_fortran) {
     for (i = 0 ; i < msr_len; i++) (*bindx)[i] = 0;
     for (i = 0 ; i < msr_len; i++) (*val)[i] = 0;
  }

  nz_ptr      = N_update + 1;
  (*bindx)[0] = nz_ptr;
  current     = 0;

  if (proc != 0) {

    /*
     * Send requests to processor 0.  When the response is received add the
     * corresponding row to the matrix and make another request.  When all the
     * requests are done, send -1 as a request to signal processor 0 that we are
     * finished.
     */

    for (i = 0; i < N_update; i++ ) {
      mdwrap_write((void *) &(update[i]), sizeof(int), 0, type, &st);
      pnode = 0;
      mdwrap_iread(str, buf_len, &pnode, &type2, &request);
      j = mdwrap_wait(str, buf_len, &pnode, &type2, &st, &request);
      while (j == sizeof(int)) {
        lil = (int *) str;
        buf_len = (unsigned int) lil[0];
        str = (char *) AZ_realloc(str,buf_len*sizeof(char));
        if (str == 0) {
          (void) fprintf(stderr,
                         "ERROR: Not enough dynamic memory in AZ_input_msr()\n");
          exit(-1);
        }
        mdwrap_iread(str, buf_len, &pnode, &type2, &request);
        j = mdwrap_wait(str, buf_len, &pnode, &type2, &st, &request);
      }

      AZ_add_new_row_nodiag(update[i], &nz_ptr, &current, val, bindx, str, dfp,
                     &msr_len,&column0);
    }

    temp = -1;
    mdwrap_write((void *) &temp, sizeof(int), 0, type, &st);
  }

  else {
#ifdef binary
    dfp = fopen(datafile, "rb");
#else
    dfp = fopen(datafile, "r");
#endif
    if (dfp == NULL) {
      (void) fprintf(stderr, "ERROR: Matrix data file %s not found\n", datafile);
      exit(-1);
    }
    (void) fprintf(stdout, " reading matrix (current version is very slow) .");
    (void) fflush(stdout);

    /* read in number of blks */
    /*
      fscanf(dfp, "%d", &total);
      for (i = 0; i <= total; i++ ) fscanf(dfp, "%d", &temp);
      */

    /* read past cnptr info (not used) */

#ifdef binary
    kkk = fread(&total, sizeof(int), 1, dfp);
#else
    kkk = fscanf(dfp, "%d", &total);  /* read in number of elements */
#endif

    if (kkk <= 0) {
       (void) fprintf(stderr,"data file %s is empty\n", datafile);
       exit(1);
    }

    if (total != totalN) {
       (void) fprintf(stderr,"\nError: data file %s contains %d rows ",datafile, total);
       (void) fprintf(stderr,"while the user's input\n     requires %d rows.\n",
                      totalN);
       exit(1);
    }

    /* get the first requests from all the processors */

    requests    = (int *) AZ_allocate(nprocs*sizeof(int));
    requests[0] = -1;
    for (i = 1; i < nprocs; i++) {
      pnode = -1;
      mdwrap_iread((void *) &temp, sizeof(int), &pnode, &type, &request);
      mdwrap_wait((void *) &temp, sizeof(int), &pnode, &type, &st, &request);
      requests[pnode] = temp;
    }

    /*
     * Go through all the rows, for those rows that we own add them to our local
     * matrix.  Otherwise, read the row into a string, determine which processor
     * has requested the row, send the string to this processor, and receive
     * another request from this processor.
     */

    for (i = 0; i < total; i++) {
      count++;
      if (count%1000 == 0) {
        (void) fprintf(stdout, ".");
        (void) fflush(stdout);
      }
      if ((current < N_update) && (i == update[current])) {
        AZ_add_new_row_nodiag(i, &nz_ptr, &current, val, bindx, 0, dfp, &msr_len,
		       &column0);
      }
      else {
#ifdef binary
        kkk = fread(&temp, sizeof(int), 1, dfp);
#else
        kkk = fscanf(dfp, "%d", &temp);
#endif
        if (kkk <= 0) {
           (void) fprintf(stderr,"\nError: AZ_input_msr(), end-of-file reached");
           (void) fprintf(stderr," while reading row %d.\n",i);
           exit(1);
        }
        if (temp == 0) column0 = 1;

        j = 0;

        while (temp != -1) {
#ifdef binary
          kkk = fread(&dtemp, sizeof(double), 1, dfp);
#else
          kkk = fscanf(dfp, "%lf", &dtemp);
#endif
          if (kkk <= 0) {
           (void) fprintf(stderr,"\nError: AZ_input_msr(), end-of-file reached");
           (void) fprintf(stderr," while reading row %d.\n",i);
           exit(1);
          }

          if (j + 30 > (int) buf_len) {
            buf_len = 2*buf_len + 30;
            str = (char *) AZ_realloc(str,buf_len*sizeof(char));

            if (str == 0) {
              (void) fprintf(stderr,"ERROR: Not Enough dynamic memory in "
                             "AZ_input_msr()\n");
              exit(-1);
            }
            if (need_request != 0)  {
               mdwrap_iread((void *) &(requests[need_request]), 
		        sizeof(int), &need_request,&type,&request);
               mdwrap_wait((void *) &(requests[need_request]), 
		        sizeof(int), &need_request,&type,&st,&request);
               need_request = 0;
            }
            for (jjj = 1; jjj < nprocs; jjj++) {
              if (requests[jjj] != -1) 
                 mdwrap_write((void *) &buf_len, sizeof(int), jjj, type2, &st);
	    }
          }

          /* put 'temp' and 'dtemp' into 'str' so that they can be sent */
          /* to another processor.                                      */

          tchar = (char *) &temp;
          for (kkk = 0 ; kkk < (int)sizeof(int); kkk++) str[j+kkk] = tchar[kkk];
          j += sizeof(int);
          tchar = (char *) &dtemp;
          for (kkk = 0 ; kkk < (int) sizeof(double); kkk++ ) 
             str[j+kkk] = tchar[kkk];
          j += sizeof(double);
#ifdef binary
          kkk = fread(&temp, sizeof(int), 1, dfp);
#else
          kkk = fscanf(dfp, "%d", &temp);
#endif
          if (kkk <= 0) {
           (void) fprintf(stderr,"\nError: AZ_input_msr(), end-of-file reached");
           (void) fprintf(stderr," while reading row %d.\n",i);
           exit(1);
          }
          if (temp == 0) column0 = 1;
        }
        tchar = (char *) &m_one;
        for (kkk = 0 ; kkk < (int)sizeof(int) ; kkk++ ) str[j+kkk] = tchar[kkk];
        j += sizeof(int);

        k = 0;
        if (need_request != 0)  {
           mdwrap_iread((void *) &(requests[need_request]), sizeof(int), 
		    &need_request,&type,&request);
           mdwrap_wait((void *) &(requests[need_request]), sizeof(int), 
		    &need_request,&type,&st, &request);
           need_request = 0;
        }

        while ((k < nprocs) && (requests[k] != i)) k++;
        if (k == nprocs) {
           (void) fprintf(stderr,"\nError: AZ_input_msr(), input file contains");
           (void) fprintf(stderr," a row (%d)\n       that is not ",i);
           (void) fprintf(stderr,"assigned to any processor?\n");
           exit(1);
        }
        mdwrap_write((void *) str, (unsigned int) j, k, type2, &st);
        need_request = k;  /* read is deferred until we will need it */
      }
    }
    if (need_request != 0)  {
       mdwrap_iread((void *) &(requests[need_request]), sizeof(int), 
		&need_request,&type,&request);
       mdwrap_wait((void *) &(requests[need_request]), sizeof(int), 
		&need_request,&type,&st,&request);
       need_request = 0;
    }

    /* at this point all the requests should be -1 */

    for (k = 0 ; k < nprocs ; k++ ) {
       if (requests[k] != -1) {
          (void) fprintf(stderr,"\nError: AZ_input_msr(), processor %d ",k);
          (void) fprintf(stderr,"requested  but never received\n       ");
          (void) fprintf(stderr,"matrix row %d. Check data file.\n",
                         requests[k]);
          exit(1);
       }
    }

    if (column0 == 0) {
       (void) fprintf(stderr,"\nWARNING: AZ_input_msr(), column 0 contains ");
       (void) fprintf(stderr,"no off-diagonal elements.\n         Are you ");
       (void) fprintf(stderr,"sure that you numbered the matrix rows/columns");
       (void) fprintf(stderr," from\n         0 to n-1 (and not 1 to n).\n");
    }


    AZ_free(requests);
    fclose(dfp);
    (void) fprintf(stdout, "\n");
  }

  AZ_free(str);

} /* AZ_input_msr_matrix_nodiag */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_add_new_row_nodiag(int therow, int *nz_ptr, int *current, double **val,
                    int **bindx, char *input, FILE *dfp, int *msr_len,
		    int *column0)

/*******************************************************************************

  Add a new row to the matrix.  If input = 0, the new matrix is read from file
  pointer dfp.  Otherwise, it is read from the string 'input'.  The form of the
  input is as follows:

         col_num1 entry1 col_num2 entry2
         col_num3 entry3 -1

  On output, val[] and  bindx[] are updated to incorporate the new row.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  therow:          The global index of the row being added.

  nz_ptr:          The next available space in val[] and bindx[] for storing
                   nonzero offdiagonals.

  current:         The next available space in a[] to store the matrix diagonal.

  val, bindx:      MSR matrix arrays that will be updated to incorporate new
                   row. See User's Guide.

  input:           Contains the information describing the row to be added (if
                   input == 0, the row information is read from standard input).

*******************************************************************************/

{

  /* local variables */

  int    old_nz;
  double sum = 0.0;
  int    temp;
  double dtemp;
  int    k = 0, kk;
  char   *tchar;

  /**************************** execution begins ******************************/

  old_nz = *nz_ptr;

  if (input == 0) { 
#ifdef binary
    kk  = fread(&temp,sizeof(int),1,dfp);
#else
    kk  = fscanf(dfp, "%d", &temp);
#endif

    if (kk <= 0) {
         (void) fprintf(stderr,"\nError: format error in '.data' file ");
         (void) fprintf(stderr,"on row '%d'\n",*current);
         (void) fprintf(stderr,"      This can be caused if exponents are\n");
         (void) fprintf(stderr,"      given using 'D' instead of 'E'. \n");
       exit(1);
    }
    if (temp == 0) *column0 = 1;
  }
  else {
    tchar = (char *) &temp;
    for (kk = 0 ; kk < (int) sizeof(int) ; kk++ ) tchar[kk] = input[kk];
    k    += sizeof(int);
  }

  while (temp != -1) {
    if (input == 0) {
#ifdef binary
       kk = fread(&dtemp, sizeof(double), 1, dfp);
#else
       kk = fscanf(dfp, "%lf", &dtemp);
#endif
       if (kk <= 0) {
         (void) fprintf(stderr,"\nError: format error in '.data' file ");
         (void) fprintf(stderr,"on row '%d'\n",*current);
         (void) fprintf(stderr,"       This can be caused if exponents are\n");
         (void) fprintf(stderr,"       given using 'D' instead of 'E'. \n");
         exit(1);
       }
    }
    else {
      tchar = (char *) &dtemp;
      for (kk = 0 ; kk < (int) sizeof(double) ; kk++ ) tchar[kk] = input[k+kk];
      k += sizeof(double);
    }

      if (*nz_ptr >= *msr_len) {
        *msr_len = (int) ( 1.5 * (double) *msr_len);
        if (!AZ_using_fortran) {
          *bindx = (int *) AZ_realloc(*bindx,*msr_len*sizeof(int));
          *val   = (double *) AZ_realloc(*val,*msr_len*sizeof(double));
        }
        if (*val == 0) {
          (void) fprintf(stderr,
                         "ERROR: Not enough dynamic memory in AZ_read_msr()\n");
          exit(-1);
        }
      }
      (*bindx)[*nz_ptr] =  temp;
      (*val)[*nz_ptr]   = dtemp;
      (*nz_ptr)++;

    if (input == 0) {
#ifdef binary
       kk  = fread(&temp,sizeof(int),1,dfp);
#else
       kk = fscanf(dfp, "%d", &temp);
#endif
       if (kk <= 0) {
         (void) fprintf(stderr,"\nError: format error in '.data' file ");
         (void) fprintf(stderr,"on row '%d'\n",*current);
         (void) fprintf(stderr,"       This can be caused if exponents are\n");
         (void) fprintf(stderr,"       given using 'D' instead of 'E'. \n");
         exit(1);
       }
       if (temp == 0) *column0 = 1;
    }
    else {
      tchar = (char *) &temp;
      for (kk = 0 ; kk < (int) sizeof(int) ; kk++ ) tchar[kk] = input[kk+k];
      k    += sizeof(int);
    }
  }

  (*val)[*current]     = sum;
  (*bindx)[*current+1] = (*bindx)[*current] + (*nz_ptr - old_nz);
  (*current)++;

} /* AZ_add_new_row_nodiag */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ML_find_local_indices(int N_update, int bindx[], int update[],
                           int *sorted_ext, int N_external, int map[])

/*******************************************************************************

  Given the global column indices for the matrix and a list of elements updated
  on this processor, compute the external elements and store them in the list
  'external' and change the global column indices to local column indices. In
  particular,

  On input, the column index bindx[k] is converted to j on output where

          update[j] = bindx[k]
   or
          external[j - N_update] = bindx[k]

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  bindx:           MSR or VBR column indices. On input, they refer to global
                   column indices. On output, they refer to local column indices
                   as described above. See User's Guide for more information.

  update:          List (global indices) of elements updated on this node.

  external:        List (global indices) of external elements on this node.

  N_external:      Number of external elements on this processor.

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

  bnptr:           'bpntr[N_update]' indicates the location of the
                   last VBR nonzero.

*******************************************************************************/

{

  /* local variables */

  int  j, k;
  int *bins,shift;
  int  start,end;

  /**************************** execution begins ******************************/

  /* set up some bins so that we will be able to use AZ_quick_find() */

  bins = (int *) ML_allocate((N_update / 4 + 10)*sizeof(int));
  if  (bins == NULL) {
    (void) fprintf(stderr, "ERROR: Not enough temp space\n");
    exit(-1);
  }

  AZ_init_quick_find(update, N_update, &shift, bins);

  /*
   * Compute the location of the first and last column index that is stored in
   * the bindx[].
   */

  start = bindx[0]; end = bindx[bindx[0]-1]; 
  
  /*
   * Estimate the amount of space we will need by counting the number of
   * references to columns not found among 'update[]'.  At the same time replace
   * column indices found in update[] by the appropriate index into update[].
   * Add N_update to columns not found in 'update[]' (this effectively marks
   * them as external elements).
   *
   * Note the space estimate will be an over-estimate as we do not take into
   * account that there will be duplicates in the external list.
   */

  for (j = start; j < end; j++) {
    k = AZ_quick_find(bindx[j], update, N_update,shift,bins);

    if (k != -1) bindx[j] = k;
    else {
       k = AZ_find_index(bindx[j], sorted_ext,N_external);
       if (k != -1) bindx[j] = map[k];
       else {
        (void) fprintf(stderr, "%d: ERROR: column number not found %d\n",
		       my_proc_id,
                       bindx[j]);
        exit(-1);
      }
    }
  }

  ML_free((char *) bins);

} /* ML_find_local_indices */

/************************************************************************/
/* Take a matrix that is effectively partitioned by columns and         */
/* transform it into one that is partitioned by rows. The original      */
/* matrix was most likely created by transposing a matrix partitioned   */
/* by row.                                                              */
/*----------------------------------------------------------------------*/


int ML_Operator_ColPartition2RowPartition(ML_Operator *A, ML_Operator *Atrans) {


  ML_Operator *eye1, *eye2;

  eye1 = ML_Operator_Create(A->comm);
  eye2 = ML_Operator_Create(A->comm);

  ML_Operator_Set_ApplyFuncData(eye1, A->invec_leng, A->invec_leng,
			ML_EXTERNAL,NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye1, ML_EXTERNAL, A->invec_leng, eye_getrows);

  ML_Operator_Set_ApplyFuncData(eye2, A->invec_leng, A->invec_leng,
			ML_EXTERNAL,NULL, A->invec_leng, eye_matvec, 0);
  ML_Operator_Set_Getrow(eye2, ML_EXTERNAL, A->invec_leng, eye_getrows);
  ML_2matmult(A, eye1, Atrans);

  /*
  if (A->getrow->use_loc_glob_map == ML_YES) 
     pr_error("ML_Operator_ColPartition2RowPartition: Matrix already has local column indices mapped to global indices\n");
  if (A->getrow->pre_comm != NULL) 
     pr_error("ML_Operator_ColPartition2RowPartiion: Matrix has a pre-communication structure?\n");

  ML_create_unique_col_id(A->invec_leng, &(A->getrow->loc_glob_map),
                           NULL, &max_per_proc, A->comm);

  if (A->getrow->post_comm != NULL)
      ML_exchange_rows( A, &Acomm, A->getrow->post_comm);
  else Acomm = A;

  ML_back_to_csrlocal(Acomm, Atrans, max_per_proc);

  ML_free(A->getrow->loc_glob_map); A->getrow->loc_glob_map = NULL;

  A->getrow->use_loc_glob_map = ML_NO;

  if (A->getrow->post_comm != NULL) {
      tptr = Acomm;
      while ( (tptr!= NULL) && (tptr->sub_matrix != A))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(Acomm);
      ML_Operator_Destroy(Acomm);
   }
  */

  return 1;
}

/************************************************************************/
/* Getrow function for the identity matrix.                             */
/*----------------------------------------------------------------------*/

int eye_getrows(void *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   double *temp;
   int    i;

   temp = (double *) data;

   if (allocated_space < N_requested_rows) return(0);

   for (i = 0; i < N_requested_rows; i++) {
      row_lengths[i] = 1;
      columns[i]     = requested_rows[i];
      values[i]      = 1.;
   }
   return(1);
}

/************************************************************************/
/* Matvec function for the identity matrix.                             */
/*----------------------------------------------------------------------*/

int eye_matvec(void *Amat_in, int ilen, double p[], int olen, double ap[])
{
  int i;

  for (i = 0; i < olen; i++) ap[i] = p[i];

  return(1);
}

/************************************************************************/
/* Take the transpose of an ML_Operator and realign resulting matrix    */
/* so that it is partitioned by rows.                                   */
/*----------------------------------------------------------------------*/

int ML_Operator_Transpose_byrow(ML_Operator *A, ML_Operator *Atrans) {
  ML_Operator *temp;

  temp = ML_Operator_Create(A->comm);
  ML_Operator_Transpose(A, temp);
  ML_Operator_ColPartition2RowPartition(temp, Atrans);
  return 1;
}
