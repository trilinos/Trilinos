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

/*******************************************************************************
Here is a sample input file:
#
#  Test input file ML 
#       
# Notes: 
#   1) Capitalization should not matter
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
*******************************************************************************/

/*******************************************************************************
Output files:
    PPn_1 (n = 1,...,numproc)   Matrix-vector product.  Matrix is the
                                prolongator P, and vector is the vector whose
                                i_th entry is the aggregate that i belongs to. 

    aggn_0 (n=1,...,numproc)    (Node,aggregate) pairs.
*******************************************************************************/
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


int *Ke_data_org = NULL, *Kn_data_org = NULL, *global_edge_inds = NULL,
  *global_node_inds = NULL, *global_edge_externs = NULL,
  *global_node_externs = NULL;
int    *reordered_glob_edges = NULL, *reordered_edge_externs = NULL;
int    *reordered_glob_nodes = NULL, *reordered_node_externs = NULL;
int    *cpntr = NULL, *Ke_bindx = NULL, *Kn_bindx = NULL, Nlocal_edges, iii, *Tmat_bindx;
int    *update, *update_index;
int *external, *extern_index;
struct reader_context *context;
#ifdef ALEGRA
extern void AZ_input_msr_matrix(char datafile[], int update[], double **val, 
				int **bindx, int N_update, int proc_config[]);
extern void AZ_input_update(char datafile[], int *N_update, int *update[], 
			    int proc_config[], int N, int chunk, 
			    int input_option);
#endif

#ifdef debugSmoother
double *xxx;
#endif
/******************************************************************************/

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

  double *Ke_val = NULL, *Kn_val, *Tmat_val = NULL, *rhs, solve_time,
    setup_time, start_time, *yyy, *vvv, *zzz, *fido, *cvec, *nodalvec;
  AZ_MATRIX *Ke_mat, *Kn_mat;
  ML_Operator *Tmat, *Tmat_trans, *Kn_coarse, *Tcoarse, *Tcoarse_trans,
    *Pn_coarse, *Rn_coarse, *Pe, *Amat;
  /*, *Tmat_array[25],
   *Tmat_trans_array[25]; */
  ML_Operator **Tmat_array, **Tmat_trans_array;
  ML_Operator *junkmat;
  ML_Operator *Ttrans_byrow;
  AZ_PRECOND *Pmat = NULL;
  ML *ml_edges, *ml_nodes;
  FILE *fp;
  int ch,i, j, Nrigid, *garbage, nblocks, *blocks;
  struct AZ_SCALING *scaling;
  ML_Aggregate *ag;
  double *mode, *rigid, alpha, *newval;
  char filename[80],vecname[128];
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
  double *Tcoarse_val, dtemp, dtemp2, *node2proc;
  int Nexterns, *sorted_ext, *map;
  int Nnz;
  char str[80];
  int Nghost;
  int mg_cycle_type;
  double omega;
#ifndef debugSmoother
  double *xxx;
#endif

#ifdef ML_partition
  FILE *fp2;
  int count, *pcounts, *proc_id;
  int Tmat_size, *proc_assignment, p1, p2, col1, col2;
  int evenodd[2]; int Iwin[4];
  int tiebreaker;

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
#endif /* ifdef ML_partition */

#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif
  my_proc_id = proc_config[AZ_node];

#ifndef ML_partition
  if (proc_config[AZ_node] == 0)
  {
    printf("Reading settings from 'ml_inputfile'\n");
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
#endif /* ifndef ML_partition */

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

  if (proc_config[AZ_node] == 0)
  {
     printf("Reading edge matrix.\n"); fflush(stdout);
  }

  AZ_input_msr_matrix("Ke_mat.az", global_edge_inds, &Ke_val, &Ke_bindx, 
                      Nlocal_edges, proc_config);

  if (proc_config[AZ_node] == 0)
  {
     printf("Done reading edge matrix.\n"); fflush(stdout);
  }

  /* For now we add something to the diagonal instead of adding in the */
  /* real mass matrix.                                                 */

  /*
    for (i = 0; i < Nlocal_edges; i++) Ke_val[i] += 1.;
    */

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

  /* Set rhs */

  /*rhs=(double *) ML_allocate(Nlocal_edges*sizeof(double));*/
  rhs = (double *) ML_allocate((Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                               *sizeof(double)); 
 
  fp = fopen("rhsfile","r");
  if (fp == NULL)
  {
    printf("%d: rhsfile file pointer is NULL\n",proc_config[AZ_node]); fflush(stdout);
    if (proc_config[AZ_node] == 0) printf("taking zero vector for rhs\n");
    fflush(stdout);
    for (i = 0; i < Nlocal_edges; i++) rhs[i] = 0.0;
  }
  else
  {
    fclose(fp);
    if (proc_config[AZ_node] == 0)
    {
       printf("%d: reading rhs from a file\n",proc_config[AZ_node]);
       fflush(stdout);
    }
    AZ_input_msr_matrix("rhsfile", global_edge_inds, &rhs, &garbage, 
			            Nlocal_edges, proc_config);
    if (proc_config[AZ_node] == 0)
    {
       printf("%d: Done reading rhs from a file\n",proc_config[AZ_node]);
       fflush(stdout);
    }
  }

#ifdef ZEROOUTDIRICHLET
  if (proc_config[AZ_node] == 0) printf("Zeroing out Dirichlet columns\n");
  AZ_zeroDirichletcolumns(Ke_mat, rhs, proc_config);
#else
  if (proc_config[AZ_node] == 0) printf("Not zeroing out Dirichlet columns\n");
#endif

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

  if (proc_config[AZ_node] == 0)
  {
     printf("Reading node matrix.\n"); fflush(stdout);
  }
  AZ_input_msr_matrix("Kn_mat.az", global_node_inds, &Kn_val, &Kn_bindx, 
		      Nlocal_nodes, proc_config);
  if (proc_config[AZ_node] == 0)
  {
     printf("Done reading node matrix.\n"); fflush(stdout);
  }

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
	
  if (proc_config[AZ_node] == 0)
  {
     printf("Reading T matrix\n"); fflush(stdout);
  }
  AZ_input_msr_matrix_nodiag("Tmat.az", global_edge_inds, &Tmat_val, 
			     &Tmat_bindx,  Nlocal_edges, proc_config);
  if (proc_config[AZ_node] == 0)
  {
     printf("Done reading T matrix\n"); fflush(stdout);
  }

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
  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
							  ML_CSR_MSRdata));
  csr_data->columns = Tmat_bindx;
  csr_data->values  = Tmat_val;
  ML_MSR2CSR(csr_data, Nlocal_edges, &Ncols);
  ML_Create(&ml_edges, N_levels);
  ML_Create(&ml_nodes, N_levels);
  Nexterns = Kn_data_org[AZ_N_external];

  AZ_Tmat_transform2ml(Nexterns, global_node_externs, reordered_node_externs,
		       Tmat_bindx, Tmat_val, csr_data->rowptr, Nlocal_nodes,
		       global_node_inds, ml_edges->comm, Nlocal_edges, &Tmat);
  ML_free(csr_data);

  Tmat->data_destroy = ML_CSR_MSRdata_Destroy;
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

/**** check symmetry of Ke ****/
    Amat = &(ml_edges->Amat[N_levels-1]);
    yyy = (double *) ML_allocate( Amat->invec_leng * sizeof(double) );
    vvv = (double *) ML_allocate( Amat->invec_leng * sizeof(double) );
    zzz = (double *) ML_allocate( Amat->invec_leng * sizeof(double) );

    AZ_random_vector(yyy, Ke_data_org, proc_config);
    AZ_random_vector(vvv, Ke_data_org, proc_config);

    ML_Operator_Apply(Amat, Amat->invec_leng, yyy,Amat->outvec_leng,zzz);
    dtemp = sqrt(abs(ML_gdot(Amat->outvec_leng, vvv, zzz, ml_edges->comm)));
    ML_Operator_Apply(Amat, Amat->invec_leng, vvv,Amat->outvec_leng,zzz);
    dtemp2 =  sqrt(abs(ML_gdot(Amat->outvec_leng, yyy, zzz, ml_edges->comm)));

    if ( (abs(dtemp-dtemp2) > 1e-15) && proc_config[AZ_node]== 0)
    {
       fprintf(stderr,"\a\n*****************\n"
                      "WARNING: Edge matrix may not be symmetric.\n");
       fprintf(stderr,"\n              ||vvv^{t} * Ke_mat * yyy|| = %20.15e\n",
               dtemp);
       fprintf(stderr,"\n              ||yyy^{t} * Ke_mat * vvv|| = %20.15e\n",
               dtemp2);
       fprintf(stderr,"               (vvv and yyy are random)\n");
       fprintf(stderr,"*****************\n\n");
       fflush(stdout);
    }      
    else if (proc_config[AZ_node]== 0)
       printf("\nEdge matrix passed symmetry check.\n\n");
    ML_free(yyy); ML_free(zzz); ML_free(vvv);
/**** end of Ke symmetry check ****/

  ML_CommInfoOP_Clone(&(Tmat->getrow->pre_comm),
                      ml_nodes->Amat[N_levels-1].getrow->pre_comm);

  /*
    sprintf(str,"P%d",proc_config[AZ_node]);
    ML_CommInfoOP_Print(Tmat->getrow->pre_comm, str);
    */

  /********************************************************************/
  /*                      Set up Tmat_trans                           */
  /*------------------------------------------------------------------*/

  Tmat_trans = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Transpose_byrow(Tmat, Tmat_trans);

#ifdef ML_partition

  evenodd[0] = 0; evenodd[1] = 0;
  for (i=0; i<nblocks; i++) Iwin[i] = 0;

  /* this code is meant to partition the matrices so that things can be */
  /* run in parallel later.                                             */
  /* It is meant to be run on only one processor.                       */

  ML_Operator_AmalgamateAndDropWeak(&(ml_nodes->Amat[N_levels-1]),
                                    num_PDE_eqns, 0.0);
  ML_Gen_Blocks_Metis(ml_nodes, N_levels-1, &nblocks, &block_list);

  fp2 = fopen("node_partition","w");
  for (i = 0; i < nblocks; i++)
  {
    count = 0;
    for (j = 0; j < ml_nodes->Amat[N_levels-1].outvec_leng; j++)
    {
      if (block_list[j] == i) count++;
    }
    fprintf(fp2,"   %d\n",count*num_PDE_eqns);
    for (j = 0; j < ml_nodes->Amat[N_levels-1].outvec_leng; j++)
    {
      if (block_list[j] == i)
      {
	     for (k = 0; k < num_PDE_eqns; k++)
            fprintf(fp2,"%d\n",j*num_PDE_eqns+k);
      }
    }
  }
  fclose(fp2);
  ML_Operator_UnAmalgamateAndDropWeak(&(ml_nodes->Amat[N_levels-1]),
                                      num_PDE_eqns,0.0);

  /* Need to partition Tmat given a partitioning of Ke_mat */
  /* Here's how it should work:                          */
  /*    1) if both column points are within processor,   */
  /*       assign row to this processor.                 */
  /*    2) if columns are from different procs,          */
  /*        a) take lower if product of columns is even  */
  /*        b) take upper if product of columns is odd.  */

  csr_data = (struct ML_CSR_MSRdata *) Tmat->data;
  proc_assignment = (int *) AZ_allocate( Tmat->outvec_leng*sizeof(int));

  j = 0;
  for (i = 0; i < Tmat->outvec_leng; i++)
  {
    /* Calculate the actual number of nonzeros (<=2) in the row. */
    itemp = (csr_data->rowptr)[i+1] - (csr_data->rowptr)[i];
    row_start = (csr_data->rowptr)[i];
    if (itemp == 2)
    {
      /* Entry could be zero if node is a Dirichlet point. */
      if ( (csr_data->values)[row_start+1] == 0.0) itemp--;
    }
    if (itemp > 0)
    {
      /* Entry could be zero if node is a Dirichlet point. */
       if ( (csr_data->values)[row_start] == 0.0)
       {
          itemp--;
          row_start++;
       }
    }
    if ( itemp > 2) 
      pr_error("Too many nonzeros per row in Tmat   %d\n", itemp);

    if (itemp == 1)
    {
      col1 = (csr_data->columns)[row_start];
      proc_assignment[i] = block_list[col1];
    }
    else if (itemp == 2)
    {
      col1 = (csr_data->columns)[row_start];
      col2 = (csr_data->columns)[row_start+1];
      p1   = block_list[col1];
      p2   = block_list[col2];
      /* 1+(int) (  2.0*rand()/(RAND_MAX+1.0)); */
      /* Randomly pick an integer between 1 and 10. */
      tiebreaker =  1+(int) (  10.0*rand()/(RAND_MAX+1.0));
      if ( tiebreaker > 5  )
      /*if ( (col1*col2)%2 == 0 ) */
      {
         evenodd[0]++;
         if (p1 < p2)
         {
            proc_assignment[i] = p1;
            Iwin[p1]++;
         }
         else
         {
            proc_assignment[i] = p2;
            Iwin[p2]++;
         }
      }
      else
      {
         evenodd[1]++;
         if (p2 < p1)
         {
            proc_assignment[i] = p1;
            Iwin[p1]++;
         }
         else
         {
            proc_assignment[i] = p2;
            Iwin[p2]++;
         }
      }
    }
    else
    {
       proc_assignment[i] = -1;
       j = 1;
    }
  }
  if (j==0) printf("\aWarning: no Dirichlet edges found in Tmat.\n");
  printf("even tie breakers = %d\n",evenodd[0]);
  printf("odd tie breakers = %d\n",evenodd[1]);
  for (i=0; i<nblocks; i++)
     printf("proc. %d won %d times\n",i,Iwin[i]);
  pcounts = (int *) ML_allocate( sizeof(int)*nblocks);
  proc_id = (int *) ML_allocate( sizeof(int)*nblocks);
  for (i = 0; i < nblocks; i++) pcounts[i] = 0;
  for (i = 0; i < nblocks; i++) proc_id[i] = i;

  count = 0;
  for (i = 0; i < Tmat->outvec_leng ; i++)
  {
    if (proc_assignment[i] != -1) pcounts[proc_assignment[i]]++;
    else count++;
  }

  printf("%d:Tmat->invec_leng = %d\n",Tmat->comm->ML_mypid,Tmat->invec_leng);
  printf("%d:Tmat->outvec_leng = %d\n",Tmat->comm->ML_mypid,Tmat->outvec_leng);
  fflush(stdout);
  
  ML_az_sort(pcounts, nblocks, proc_id, NULL);

  i = 0; j = 0;
  while ( (i < nblocks) && (pcounts[i] < pcounts[nblocks-1]) &&
          (j < Tmat->outvec_leng) )
  {
    if ( proc_assignment[j] == -1)
    {
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
#endif /*ifdef ML_partition*/

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

#ifdef HierarchyCheck
  xxx = (double *) ML_allocate((Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                               *sizeof(double)); 


  for (iii = 0; iii < Nlocal_edges; iii++) xxx[iii] = 0.0; 

  /* Set xxx */

  if (proc_config[AZ_node]== 0)
     printf("putting in an edge based initial guess\n");
  fp = fopen("initguessfile","r");
  if (fp != NULL)
  {
    fclose(fp);
    free(xxx);
    if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
    AZ_input_msr_matrix("initguessfile", global_edge_inds, &xxx, &garbage,
			Nlocal_edges, proc_config);
    options[AZ_conv] = AZ_expected_values;
    printf("done reading initial guess\n");
  }
#endif


  coarsest_level = ML_Gen_MGHierarchy_UsingReitzinger(ml_edges, ml_nodes,
						      N_levels-1, ML_DECREASING, ag, Tmat, Tmat_trans,
						      &Tmat_array, &Tmat_trans_array);


  coarsest_level = N_levels - coarsest_level;
  if ( proc_config[AZ_node] == 0 )
    printf("Coarse level = %d \n", coarsest_level);
	

  /* set up smoothers for all levels but the coarsest */
  for (level = N_levels-1; level > coarsest_level; level--)
  {

      num_PDE_eqns = ml_edges->Amat[level].num_PDEs;
		
      /*  Sparse approximate inverse smoother that acutally does both */
      /*  pre and post smoothing.                                     */

      if (ML_strcmp(context->smoother,"Parasails") == 0)
	{
	  ML_Gen_Smoother_ParaSails(ml_edges , level, ML_PRESMOOTHER, nsmooth, 
				    parasails_sym, parasails_thresh, 
				    parasails_nlevels, parasails_filter,
				    parasails_loadbal, parasails_factorized);
	}

      else if (ML_strcmp(context->smoother,"Hiptmair") == 0)
	{
       omega = 1.0;
       if ( proc_config[AZ_node] == 0 )
          printf("Damping parameter = %lf\n",omega);
	  ML_Gen_Smoother_Hiptmair(ml_edges, level, ML_BOTH, nsmooth,
				   omega,Tmat_array, Tmat_trans_array);
	}
      /* This is the symmetric Gauss-Seidel smoothing that we usually use. */
      /* In parallel, it is not a true Gauss-Seidel in that each processor */
      /* does a Gauss-Seidel on its local submatrix independent of the     */
      /* other processors.                                                 */

      else if (ML_strcmp(context->smoother,"GaussSeidel") == 0)
	{
	  ML_Gen_Smoother_GaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.);
	}
      else if (ML_strcmp(context->smoother,"SymGaussSeidel") == 0)
	{
	  ML_Gen_Smoother_SymGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.);
	}
      else if (ML_strcmp(context->smoother,"BlockGaussSeidel") == 0)
	{
	  ML_Gen_Smoother_BlockGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.,
					   num_PDE_eqns);
	}
      else if (ML_strcmp(context->smoother,"Aggregate") == 0)
	{
	  ML_Gen_Blocks_Aggregates(ag, level, &nblocks, &blocks);
	  ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , level, ML_BOTH,
					       nsmooth,1., nblocks, blocks);
	}

      /* This is a true Gauss Seidel in parallel. This seems to work for  */
      /* elasticity problems.  However, I don't believe that this is very */
      /* efficient in parallel.                                           */       
      /*
         blocks = (int *) ML_allocate(sizeof(int)*ml_edges->Amat[level].invec_leng);
	nblocks = ml_edges->Amat[level].invec_leng;
	for (i =0; i < nblocks; i++) blocks[i] = i;
	ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml_edges , level,
	ML_PRESMOOTHER, nsmooth, 1., nblocks, blocks);
	ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml_edges, level,
	ML_POSTSMOOTHER, nsmooth, 1., nblocks, blocks);
	ML_free(blocks);
	*/

      /* Jacobi Smoothing                                                 */

      else if (ML_strcmp(context->smoother,"Jacobi") == 0)
	{
	  ML_Gen_Smoother_Jacobi(ml_edges , level, ML_PRESMOOTHER, nsmooth,.67);
	  ML_Gen_Smoother_Jacobi(ml_edges , level, ML_POSTSMOOTHER, nsmooth,.67);
	}

      /*  This does a block Gauss-Seidel (not true GS in parallel)        */
      /*  where each processor has 'nblocks' blocks.                      */
      /* */

      else if (ML_strcmp(context->smoother,"Metis") == 0)
	{
	  nblocks = 250;
	  nblocks = ml_edges->Amat[level].invec_leng/25;
	  nblocks++;
	  ML_Gen_Blocks_Metis(ml_edges, level, &nblocks, &blocks);
	  ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , level, ML_BOTH,
					       nsmooth,1., nblocks, blocks);
	  ML_free(blocks);
	}
      else
	{
	  printf("unknown smoother %s\n",context->smoother);
	  exit(1);
	}
  }
  nsmooth   = context->coarse_its;
  /*  Sparse approximate inverse smoother that actually does both */
  /*  pre and post smoothing.                                     */

  if (ML_strcmp(context->coarse_solve,"Parasails") == 0)
  {
      ML_Gen_Smoother_ParaSails(ml_edges , coarsest_level, ML_PRESMOOTHER,
				nsmooth, parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
  }

  else if (ML_strcmp(context->coarse_solve,"Hiptmair") == 0)
    {
       omega = 1.0;
       if ( proc_config[AZ_node] == 0 )
          printf("Damping parameter = %lf\n",omega);
      ML_Gen_Smoother_Hiptmair(ml_edges , coarsest_level, ML_BOTH,
			       nsmooth,omega,Tmat_array, Tmat_trans_array);
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
    nblocks = ml_edges->Amat[coarsest_level].invec_leng/25;
    nblocks++;
    ML_Gen_Blocks_Metis(ml_edges, coarsest_level, &nblocks, &blocks);
    ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, 
					 nsmooth,1., nblocks, blocks);
    if (blocks != NULL) ML_free(blocks);
  }
  else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0) {
    ML_Gen_CoarseSolverSuperLU( ml_edges, coarsest_level);
  }
  else {
    printf("unknown coarse grid solver %s\n",context->coarse_solve);
    exit(1);
  }
		
  mg_cycle_type = ML_MGV;
  ML_Gen_Solver(ml_edges, mg_cycle_type, N_levels-1, coarsest_level); 
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
  else
  {
    printf("unknown krylov method %s\n",context->krylov);
  }
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_user_precond;
  /*options[AZ_conv]     = AZ_noscaled;*/
  options[AZ_conv]     = AZ_r0;
  options[AZ_output]   = 1;
  options[AZ_max_iter] = 300;
  options[AZ_poly_ord] = 5;
  options[AZ_kspace]   = 130;
  params[AZ_tol]       = context->tol;
  options[AZ_output]   = context->output;
	
  AZ_set_ML_preconditioner(&Pmat, Ke_mat, ml_edges, options); 
  setup_time = AZ_second() - start_time;
	
  xxx = (double *) ML_allocate((Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                               *sizeof(double)); 

  for (iii = 0; iii < Nlocal_edges; iii++) xxx[iii] = 0.0; 

  /* Set xxx */

  if (proc_config[AZ_node]== 0)
     printf("putting in an edge based xxx\n");
  fp = fopen("initguessfile","r");
  if (fp != NULL) {
    fclose(fp);
    free(xxx);
    if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
    AZ_input_msr_matrix("initguessfile", global_edge_inds, &xxx, &garbage,
			Nlocal_edges, proc_config);
    options[AZ_conv] = AZ_expected_values;
    printf("done reading initial guess\n");
  }
  else if (proc_config[AZ_node]== 0) printf("taking 0 initial guess \n");

  AZ_reorder_vec(xxx, Ke_data_org, reordered_glob_edges, NULL);

  dtemp = sqrt(ML_gdot(Nlocal_edges, xxx, xxx, ml_edges->comm));
  if (proc_config[AZ_node]== 0)
  {
    printf("length of initial guess = %d\n",Nlocal_edges);
    printf("||xxx|| = %e\n",dtemp);
  }
  dtemp = sqrt(ML_gdot(Nlocal_edges, rhs, rhs, ml_edges->comm));
  if (proc_config[AZ_node]== 0)
  {
  printf("||rhs|| = %e\n",dtemp);
  fflush(stdout);
  }

  /*
    printf("putting in an node based xxxx\n");
    fp = fopen("initguessfile","r");
    if (fp != NULL) {
    fclose(fp);
    if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
    AZ_input_msr_matrix("initguessfile", global_node_inds, &xxx, &garbage, Nlocal_nodes, 
    proc_config);
    options[AZ_conv] = AZ_expected_values;
    }
    else if (proc_config[AZ_node]== 0) printf("taking 0 initial guess \n");

    dtemp = sqrt(ML_gdot(Nlocal_nodes, xxx, xxx, ml_edges->comm));
    printf("length of initial guess = %d\n",Nlocal_nodes);
    printf("||xxx|| = %e\n",dtemp);

    AZ_reorder_vec(xxx, Kn_data_org, reordered_glob_nodes, NULL);
    */

  fp = fopen("AZ_no_multilevel.dat","r");
  scaling = AZ_scaling_create();
  start_time = AZ_second();
  if (fp != NULL)
  {
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
  else
  {
    options[AZ_keep_info] = 1;
    /*options[AZ_conv] = AZ_noscaled;*/
    options[AZ_conv] = AZ_r0;
    options[AZ_output] = 1;

    /*
      options[AZ_precond] = AZ_none;
     
      Ke_mat->matvec(xxx, rhs, Ke_mat, proc_config);
      for (i = 0; i < Nlocal_edges; i++)
         printf("%7d     %7d %20.15e %20.15e\n",i+1,i+1,xxx[i],rhs[i]);
      printf("huhhhh %e\n",Ke_mat->val[0]);
      */

    /*options[AZ_conv] = AZ_noscaled;*/
    options[AZ_conv] = AZ_r0;
    /*options[AZ_conv] = AZ_expected_values;*/

#ifdef CHECKOPERATORS
/**** check various operators and vectors ****/
    if (N_levels > 1)
    {
       Amat = &(ml_edges->Rmat[N_levels-1]);
       yyy = (double *) malloc( Amat->outvec_leng * sizeof(double) );
       ML_Operator_Apply(Amat, Amat->invec_leng, rhs,Amat->outvec_leng,yyy);
       dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
       printf("||R_e * rhs|| = %10.7e\n",dtemp);

       /*
       if (Amat->comm->ML_mypid == 0)
          ML_Operator_Print(Amat,"0:Rmat");
       */

       ML_Operator_Apply(Amat, Amat->invec_leng, xxx,Amat->outvec_leng,yyy);
       dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
       printf("||R_e * xxx|| = %10.7e\n",dtemp);
       ML_free(yyy);
    }


    Amat = &(ml_edges->Amat[N_levels-1]);
    yyy = (double *) malloc( Amat->outvec_leng * sizeof(double) );
    ML_Operator_Apply(Amat, Amat->invec_leng, rhs,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Ke_mat * rhs|| = %10.7e\n",dtemp);

    ML_Operator_Apply(Amat, Amat->invec_leng, xxx,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Ke_mat * xxx|| = %10.7e\n",dtemp);
    ML_free(yyy);

    Amat = Tmat_trans;
    yyy = (double *) malloc( Amat->outvec_leng * sizeof(double) );
    ML_Operator_Apply(Amat, Amat->invec_leng, rhs,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Tmat_trans * rhs|| = %10.7e\n",dtemp);

    ML_Operator_Apply(Amat, Amat->invec_leng, xxx,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Tmat_trans * xxx|| = %10.7e\n",dtemp);

    if (N_levels > 1)
    {
       Amat = &(ml_nodes->Rmat[N_levels-1]);
       vvv = (double *) malloc( Amat->outvec_leng * sizeof(double) );
       ML_Operator_Apply(Amat, Amat->invec_leng, yyy,Amat->outvec_leng,vvv);
       dtemp = sqrt(ML_gdot(Amat->outvec_leng, vvv, vvv, ml_edges->comm));
       printf("||R_n * Tmat_trans * xxx|| = %10.7e\n",dtemp);
       ML_free(vvv);
    }

    ML_free(yyy);
    fflush(stdout);
/**** end of check ****/

/**** check prolongation operator ****
    Amat = &(ml_edges->Pmat[N_levels-1]);
    vvv = (double *) malloc( Amat->outvec_leng * sizeof(double) );
    ML_Operator_Apply(Amat, Amat->invec_leng, yyy,Amat->outvec_leng,vvv);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, vvv, vvv, ml_edges->comm));
    printf("||P * R * xxx|| = %10.7e\n",dtemp);
    free(yyy);
    free(vvv);
**** end of check ****/

#endif /*ifdef CHECKOPERATORS*/

    if (proc_config[AZ_node]== 0)
    {
       if (mg_cycle_type == ML_MGV)
          printf("Cycle type = MGV\n");
       else if (mg_cycle_type == ML_MGW)
          printf("Cycle type = MGW\n");
       else
          printf("Cycle type = other\n");
    }
    fflush(stdout);
    AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, scaling); 

    options[AZ_pre_calc] = AZ_reuse;
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
  if (Kn_mat != NULL) AZ_matrix_destroy(&Kn_mat);
  free(garbage);
  free(Kn_val);
  free(Kn_bindx);
  ML_Operator_Destroy(Tmat);
  ML_Operator_Destroy(Tmat_trans);
  ML_MGHierarchy_ReitzingerDestroy(N_levels-2, coarsest_level, &Tmat_array, &Tmat_trans_array);


#ifdef ML_MPI
  MPI_Finalize();
#endif
		
  return 0;
		
}

#ifdef ALEGRA

/********************************** end of main *******************************/
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
	 * $Name$
	 *====================================================================*/
	#ifndef lint
	static char rcsid[] = "$Id$";
	#endif


	/*******************************************************************************
	 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
	 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
	 * Export of this program may require a license from the United States         *
	 * Government.                                                                 *
	 ******************************************************************************/

	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <malloc.h>
	#include <float.h>


	extern int AZ_sys_msg_type;
	extern int AZ_using_fortran;




	/*****************************************************************************/
	/*****************************************************************************/
	void AZ_input_update(char datafile[], int *N_update, int *update[], int proc_config[],
			    int N, int chunk, int input_option)

	/*******************************************************************************

	Exactly the same as AZ_read_update except it reads the update information from
	a file speficied by the input argument datafile instead of .update

	*******************************************************************************/

	{

	  /* local variables */

	  int   t1, t2, i;
	  int   ii, j;
	  int   allocated, length;
	  int   cflag;
	  int   partner;
	  int   proc_x, proc_y, proc_z;
	  int   pts_x, pts_y, pts_z;
	  int   total_pts_x, total_pts_y;
	  int   px, py, pz, k;
	  int   start_x, start_y, start_z;
	  int   end_x, end_y, end_z;
	  int   pt_number;
	  int   count, check;
	  int   proc, nprocs;
	  int   type, type2;
	  FILE *fp = NULL;
	  MPI_AZRequest request;  /* Message handle */


	  /**************************** execution begins ******************************/

	  AZ_MPI_comm_space_ok();
	  proc   = proc_config[AZ_node];
	  nprocs = proc_config[AZ_N_procs];

	  /*
	   * Figure out which chunks should be assigned to this processor using a box
	   * decomposition.  That is, it is assumed that all the chunks are ordered
	   * naturally corresponding to an m x m x m box where m = N^(1/3).  Boxes of
	   * chunks are assigned to processors.
	   *
	   * NOTE: it is assumed that nprocs = 2^power and that the number of chunks in
	   * each direction is divisible by the number of processors in each direction.
	   */

	  if (input_option == AZ_box) {

	    /* determine the number of processors in each direction */

	    if (proc == 0) {
	       (void) fprintf(stdout,"Input the dimensions of the processor cube\n\n");
	       (void) fprintf(stdout,"Enter the number of processors along x axis>");
	       (void) fflush(stdout);
	       scanf("%d",&proc_x);
	       (void) fprintf(stdout,"Enter the number of processors along y axis>");
	       (void) fflush(stdout);
	       scanf("%d",&proc_y);
	       (void) fprintf(stdout,"Enter the number of processors along z axis>");
	       (void) fflush(stdout);
	       scanf("%d",&proc_z);

	       (void) fprintf(stdout,"Input the grid dimensions\n\n");
	       (void) fprintf(stdout,"Enter the number of grid points along x axis>");
	       (void) fflush(stdout);
	       scanf("%d",&pts_x);
	       (void) fprintf(stdout,"Enter the number of grid points along y axis>");
	       (void) fflush(stdout);
	       scanf("%d",&pts_y);
	       (void) fprintf(stdout,"Enter the number of grid points along z axis>");
	       (void) fflush(stdout);
	       scanf("%d",&pts_z);
	    }
	    AZ_broadcast((char *) &proc_x, sizeof(int), proc_config, AZ_PACK);
	    AZ_broadcast((char *) &proc_y, sizeof(int), proc_config, AZ_PACK);
	    AZ_broadcast((char *) &proc_z, sizeof(int), proc_config, AZ_PACK);
	    AZ_broadcast((char *) &pts_x , sizeof(int), proc_config, AZ_PACK);
	    AZ_broadcast((char *) &pts_y , sizeof(int), proc_config, AZ_PACK);
	    AZ_broadcast((char *) &pts_z , sizeof(int), proc_config, AZ_PACK);
	    AZ_broadcast((char *) NULL   , 0          , proc_config, AZ_SEND);
	 
	    total_pts_x = pts_x;
	    total_pts_y = pts_y;


	    if ( proc_x*proc_y*proc_z != nprocs) {
		if (proc == 0) {
		  (void) fprintf(stdout,"Error: %d x %d x %d != %d ",
				 proc_x, proc_y, proc_z, nprocs);
		  (void) fprintf(stdout," (total number of processors)\n");
		}
		exit(1);
	    }

	    if ( pts_x * pts_y * pts_z != N ) {
		if (proc == 0) {
		  (void) fprintf(stdout,"Error: %d x %d x %d != %d ",
				 pts_x, pts_y, pts_z, N);
		  (void) fprintf(stdout," (total number of grid points)\n");
		}
		exit(1);
	    }
	    if ( pts_x%proc_x != 0 ) {
		if (proc == 0) {
		  (void) fprintf(stdout,"Error: grid points along x axis are not an ");
		  (void) fprintf(stdout,"even multiple of processors\n");
		  (void) fprintf(stdout,"       along x axis.");
		}
		exit(1);
	    }
	    if ( pts_y%proc_y != 0 ) {
		if (proc == 0) {
		  (void) fprintf(stdout,"Error: grid points along y axis is not an ");
		  (void) fprintf(stdout,"even multiple of processors\n");
		  (void) fprintf(stdout,"       along y axis.");
		}
		exit(1);
	    }
	    if ( pts_z%proc_z != 0 ) {
		if (proc == 0) {
		  (void) fprintf(stdout,"Error: grid points along z axis is not an ");
		  (void) fprintf(stdout,"even multiple of processors\n");
		  (void) fprintf(stdout,"       along z axis.");
		}
		exit(1);
	    }
	    pts_x = pts_x/proc_x;
	    pts_y = pts_y/proc_y;
	    pts_z = pts_z/proc_z;

	    /* compute the number of elements per processor in each direction */

	    *N_update = pts_x * pts_y * pts_z * chunk;
	    if (!AZ_using_fortran) 
	       *update     = (int *) AZ_allocate((*N_update)*sizeof(int));

	    /* compute the lower left corner and the upper right corner */

	    px = proc % proc_x;
	    pz = (proc-px) / proc_x;
	    py = pz % proc_y;
	    pz = (pz-py) / proc_y;

	    start_x = px * pts_x;
	    end_x   = start_x + pts_x;
	    start_y = py * pts_y;
	    end_y   = start_y + pts_y;
	    start_z = pz * pts_z;
	    end_z   = start_z + pts_z;

	    /* set update[] */

	    count = 0;
	    for (k = start_z; k < end_z; k++ ) {
	      for (j = start_y; j < end_y; j++ ) {
		for (i = start_x; i < end_x; i++ ) {
		  for (ii = 0; ii < chunk; ii++ ) {
		    pt_number = (i + j * total_pts_x + k * total_pts_x * total_pts_y) * 
				    chunk + ii;
		    (*update)[count++] = pt_number;
		  }
		}
	      }
	    }
	  }

	  else if (input_option == AZ_linear) {

	    /*
	     * Figure out which chunks should be assigned to this processor for linear
	     * partitioning.  This means that processor 0 is assigned the chunks
	     * approximately corresponding to 0, ... , N/nprocs and processor 1 is
	     * approximately assigned the chunks 1+N/nprocs to 2*N/nprocs.
	     */

	    t1 = N/nprocs;
	    t2 = N - t1 * nprocs;

	    if ( proc >= t2) t2 += (proc * t1);
	    else {
	      t1++;
	      t2    = proc*t1;
	    }
	    *N_update = t1*chunk;
	    t2   *= chunk;

	    if (!AZ_using_fortran) 
	       *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));

	    if (*update == NULL) {
	      (void) fprintf (stderr, "Not enough space to allocate 'update'\n");
	      exit(-1);
	    }

	    for (i = 0; i < *N_update; i++) (*update)[i] = i + t2;
	  }

	  else if (input_option == AZ_file) {

	    /* read the update elements from the data file */

	    type            = AZ_sys_msg_type;
	    AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;
	    type2           = AZ_sys_msg_type;
	    AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;

	    /*
	     * Processor 0 reads the data file and distributes the lists to the other
	     * processors.
	     */

	    t1 = 0;            /* total number of points distributed */
	    if (proc == 0) {
	      (void) printf("reading from file %s\n", datafile); fflush(stdout);

	      if ( (fp = fopen(datafile, "r")) == NULL) {
		(void) fprintf(stderr, "ERROR: file %s not found\n", datafile);
		exit(-1);
	      }

	      if (!AZ_using_fortran) *update = 0;
	      allocated   = 0;
	      for (i = nprocs - 1; i >= 0; i-- ) {

		/* read in list length and send to processor i  */

		fscanf(fp, "%d", &length);
		t1 += length;
		if (i != 0)
		  mdwrap_write((void *) &length, sizeof(int), i, type, &cflag);

		/*
		 * If this is the last list, we allocate the right amount of space and
		 * keep the list instead of sending it off
		 */

		if (i == 0) {
		  *N_update = length;
		  allocated       = 0;
		}

		/* allocate enough space for list */

		if (length > allocated ) {
		  if ((*update != NULL) && (!AZ_using_fortran)) AZ_free(*update);
		  allocated = length + 1;

		  if (!AZ_using_fortran)
		    *update = (int *) AZ_allocate(allocated*sizeof(int));
		  if (*update == NULL) {
		    (void) fprintf(stderr,
				   "Not enough space to allocate 'update'\n");
		    exit(-1);
		  }
		}

		/* read a list and send it off to proc i (if not last list) */

		for (j = 0; j < length; j++ ) fscanf(fp, "%d", *update + j);
		if (i != 0)
		  mdwrap_write((void *) *update, length*sizeof(int), i, type2, &cflag);
	      }
	      fclose(fp);

	      if (t1 != N*chunk) {
		(void) fprintf(stderr,"AZ_read_update() found %d points in file\n", t1);
		(void) fprintf(stderr,"%s instead of the requested %d\n", datafile,
			       N*chunk);
		exit(-1);
	      }
	    }

	    else {

	      /* read the update list from processor 0 */

	      partner = 0;
	      mdwrap_iread((void *) N_update, sizeof(int), &partner, &type, &request);
	      mdwrap_wait((void *) N_update, sizeof(int), &partner, &type, &cflag, &request);

	      if (!AZ_using_fortran)
		*update = (int *) AZ_allocate((*N_update+1)*sizeof(int));
	      if (*update == NULL)  {
		(void) fprintf(stderr, "Not enough space to allocate 'update'\n");
		exit(-1);
	      }

	      partner = 0;
	      mdwrap_iread((void *) *update, *N_update * sizeof(int), &partner, &type2,
		    &request);
	      mdwrap_wait((void *) *update, *N_update * sizeof(int), &partner, &type2,
		    &cflag, &request);
	    }

	    AZ_sort(*update, *N_update, NULL, NULL);

	    /* check that element '0' is contained on 1 processor. That is,  */
	    /* make sure the user has numbered from 0 to n-1 instead of from */
	    /* 1 to n                                                        */
	    check = 0;
	    if ( (*N_update > 0) && ((*update)[0] == 0) ) check = 1;
	    check = AZ_gsum_int(check, proc_config);
	    if (check != 1) {
	       if (proc == 0) {
		  (void) fprintf(stderr,"Error: In AZ_read_update(), the %s", datafile);
		  (void) fprintf(stderr,"file does not contain\n       one ");
		  (void) fprintf(stderr,"occurance of row 0. Make sure that rows are");
		  (void) fprintf(stderr," numbered\n       from 0 to n-1.\n");
	       }
	       exit(1);
	    }


	  }
	  else {
	    (void) fprintf(stderr,"Unknown input option (%d) in AZ_read_update()\n",
			   input_option);
	    exit(1);
	  }


	} /* AZ_input_update */



	/******************************************************************************/
	/******************************************************************************/
	void AZ_input_msr_matrix(char datafile[], int update[], double **val, int **bindx, 
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

	  /**************************** execution begins ******************************/

	  AZ_MPI_comm_space_ok();
	  type            = AZ_sys_msg_type;
	  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;
	  type2           = AZ_sys_msg_type;
	  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

	  proc   = proc_config[AZ_node];
	  nprocs = proc_config[AZ_N_procs];

	  totalN = AZ_gsum_int(N_update, proc_config);
	  str    = (char *) AZ_allocate((buf_len+1)*sizeof(char));
	  if (str == NULL) {
	    printf("ERROR: NOT enough dynamic memory in AZ_input_msr_matrix\n");
	    exit(-1);
	  }
	  msr_len = 8*N_update+2;
	  if (!AZ_using_fortran) {
	    *bindx = (int *)    AZ_allocate(msr_len*sizeof(int));
	    *val   = (double *) AZ_allocate(msr_len*sizeof(double));
	  }

	  if (*val == NULL) {
	    (void) fprintf(stderr,
			   "ERROR: NOT enough dynamic memory in AZ_input_msr_matrix\n");
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

	      AZ_add_new_row(update[i], &nz_ptr, &current, val, bindx, str, dfp,
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
		AZ_add_new_row(i, &nz_ptr, &current, val, bindx, 0, dfp, &msr_len,
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

	} /* AZ_input_msr_matrix */


#endif
