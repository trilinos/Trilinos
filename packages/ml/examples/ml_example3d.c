/*****************************************************************************/
/* Copyright 1998, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Sample driver for AZTEC/ML package. The software is tested by setting up  */
/* a system of equations and a right hand side and then solving the system   */
/* of equations using AZTECs iterative solvers and ML preconditoner.         */
/*                                                                           */
/* Author:       Ray Tuminaro, Div 1422, Sandia National Labs                */
/* Modified by : Charles Tong (8950)                                         */
/* date:         5/10/99                                                      */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_include.h"

#if defined(HAVE_ML_AZTEC) || defined(HAVE_ML_AZTECOO) || defined(HAVE_ML_AZTEC2_1)

#include "az_aztec.h"
#include "ml_read_utils.h"
extern int AZ_using_fortran;

/*****************************************************************************/
/* embarassingly, we have global variables                                   */
/*****************************************************************************/

int N_grid_pts;
int refine_factor;
int choice;
int num_PDE_eqns=1;
int nsmooth = 1;

int    parasails_factorized = 0;
int    parasails_sym        = 0;
double parasails_thresh     = 0.;
int    parasails_nlevels    = 0;
double parasails_filter     = 0.;
double parasails_loadbal    = 0.;

/*****************************************************************************/
/* Set up and solve a test problem defined in the subroutine                 */
/* init_matrix_vector_structures().                                          */
/*****************************************************************************/
int coarse_iterations = 0, use_cg = 0, num_levels = 2; 
  int    *update = NULL, *external = NULL;
  int    *update_index = NULL, *extern_index = NULL;

/* -------------  external function declarations -------------------------*/

extern void init_options(int options[], double params[]);

extern void init_guess_and_rhs(int update_index[],int update[], double *x[],
                   double *b[], int data_org[], double val[], int indx[], 
                   int bindx[], int rpntr[], int cpntr[], int bpntr[], 
                   int proc_config[]);

extern int construct_ml_grids(int, int *, AZ_MATRIX **, int **, int **, 
                   int **, int **, int, int **, int **, int **, int **, 
                   ML **, int *, double *, AZ_PRECOND **, int);

extern void create_msr_matrix(int*, double **, int **, int);
extern int ml_recv(void*, unsigned int, int *, int *, USR_COMM, USR_REQ *);
extern int ml_wait(void*, unsigned int, int *, int *, USR_COMM, USR_REQ *);
extern int ml_send(void*, unsigned int, int, int, USR_COMM );
extern void Generate_mesh(int, ML_GridAGX **, int *, int **, int *);
extern int ML_mapper0f( void *grid, double *din, double *dout);
extern int ML_mapper0b( void *ml, double *din, double *dout);
extern int ML_mapper1f( void *grid, double *din, double *dout);
extern int ML_mapper1b( void *ml, double *din, double *dout);
extern void add_row_3D(int row,int location,double val[],int bindx[], int *n);

/* -------------  end external function declarations ----------------------*/

#ifdef ML_BENCHMARK
  struct reader_context *context;
#endif

int main(int argc, char *argv[])
{
  int    i, input_option, precon_flag, N_elements_coarse;
  double *b, *x;
#ifdef ML_BENCHMARK
  int    j;
#endif

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];

  /* data structure for matrix corresponding to the fine grid */

  int    *rpntr = NULL,*cpntr = NULL, *indx = NULL, *bpntr = NULL;
  int    proc_factor;

  /* data structures for multilevel grid, discrete operators */
  /* and preconditioners                                     */

  AZ_PRECOND  *Pmat = NULL;
  AZ_MATRIX   *Amat = NULL;
  ML          *ml;
#ifdef ML_BENCHMARK
  char input[MAX_INPUT_STR_LN];
  FILE *ifp;
#endif

  /* ----------------------- execution begins --------------------------------*/

#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  /* Read and broadcast: problem choice, problem size, etc.  */

#ifdef ML_BENCHMARK
   if (proc_config[AZ_node] == 0) {
      ML_Reader_ReadInput("ml_inputfile", &context);
      ifp = fopen("ml_inputfile", "r");
      if (!ML_Reader_LookFor(ifp, "number of x nodes per proc", input, '='))
	proc_factor = 20;   /* Default */
      else {
	ML_Reader_ReadString(ifp, input, '\n');
	if (sscanf(input, "%d", &(proc_factor)) != 1) {
	  fprintf(stderr, "ERROR: can\'t interp int while looking for \"%s\"\n",
		  "number of x nodes per proc");
	  exit(-1);
	}
      }
   }
   else context = (struct reader_context *) malloc(sizeof(struct reader_context)
);
   AZ_broadcast((char *) &proc_factor,   sizeof(int), proc_config, AZ_PACK);
   AZ_broadcast((char *) context,  sizeof(struct reader_context), proc_config,
                AZ_PACK);
   AZ_broadcast((char *) NULL        ,   0          , proc_config, AZ_SEND);
   precon_flag = 1;
   num_levels = context->N_levels;
   num_PDE_eqns = context->N_dofPerNode;
   nsmooth   = context->nsmooth;
   use_cg = 1;
   coarse_iterations = context->coarse_its;
   refine_factor = 1;
   choice = 1;


#else

  if (proc_config[AZ_node] == 0) 
  {
     precon_flag = 1;
     printf("n x n grid in each processor, n = ");
     scanf("%d", &proc_factor);
     printf("proc_factor = %d\n", proc_factor);
     refine_factor = 1;

     choice = 1;
     printf("Coarse solver (0 for SuperLU, > 0 for GS) : \n");
     scanf("%d", &coarse_iterations);
     printf("Use CG (0 for no , 1 for yes) : \n");
     scanf("%d", &use_cg);
     printf("Number of levels : \n");
     scanf("%d", &num_levels);
     printf("number of smoothing step : \n");
     scanf("%d", &nsmooth);
     printf("number of smoothing step = %d\n", nsmooth);

/*
     printf("ParaSails factorized, (1=true or 0=false) : \n");
     scanf("%d",  &parasails_factorized);
     printf("ParaSails symmetric smoother, (1=true or 0=false) : \n");
     scanf("%d",  &parasails_sym);
     printf("ParaSails threshold, (0.01 or 0.1) : \n");
     scanf("%lf", &parasails_thresh);
     printf("ParaSails number of levels, (0 or higher) : \n");
     scanf("%d",  &parasails_nlevels);
     printf("ParaSails filter value, (0. or 0.05 or 0.10) : \n");
     scanf("%lf", &parasails_filter);
     printf("ParaSails loadbal param, (0. or 0.9) : \n");
     scanf("%lf", &parasails_loadbal);
*/
  }
  AZ_broadcast((char *) &coarse_iterations,  sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &use_cg,  sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &num_levels,  sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &num_PDE_eqns,  sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &precon_flag,   sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &proc_factor,   sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &refine_factor, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &choice       , sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &nsmooth      , sizeof(int), proc_config, AZ_PACK);
/*
  AZ_broadcast((char *) &parasails_factorized, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &parasails_sym,     sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &parasails_thresh,  sizeof(double), proc_config, AZ_PACK);
  AZ_broadcast((char *) &parasails_nlevels, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &parasails_filter,  sizeof(double), proc_config, AZ_PACK);
  AZ_broadcast((char *) &parasails_loadbal, sizeof(double), proc_config, AZ_PACK);
*/
  AZ_broadcast((char *) NULL        ,   0          , proc_config, AZ_SEND);
#endif
  i = proc_config[AZ_N_procs];
  i = (int) pow( (double) i, 0.50001 );
  i = i * proc_factor;
  N_elements_coarse = i * i * i;
  if (proc_config[AZ_node] == 0) {
     i = i * refine_factor;
     printf("Fine grid = %d x %d x %d \n", i+1 , i+1 , i+1 ); 
  }
  input_option = 1;

  /* set AZTEC options */

  init_options(options,params);

  /* create coarse grid discretization */

  construct_ml_grids(N_elements_coarse, proc_config,
                     &Amat, &update, &update_index, 
                     &external, &extern_index, input_option, &indx,
                     &bpntr, &rpntr, &cpntr, &ml, options, params, 
                     &Pmat, precon_flag);

  /* set up ML preconditioning */

  AZ_set_ML_preconditioner(&Pmat, Amat, ml, options);

  /* solve the system of equations on the fine grid */

  for ( i = 0; i < 1; i++ ) 
  {
     init_guess_and_rhs(update_index, update, &x, &b, Amat->data_org, 
                        Amat->val, indx, Amat->bindx, rpntr, cpntr, 
                        bpntr, proc_config);
     if (use_cg == 0) {
        ML_Iterate(ml, x, b);
     }
     else
     AZ_iterate(x, b, options, params, status, proc_config, Amat, Pmat, 
#ifdef AZ_ver2_1_0_5
                    (struct AZ_SCALING *)
#else
                    (struct scale*)
#endif
                                    NULL);

  }
  
  /*
  for (i=0; i<Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border]; i++) {
     diff = (x[update_index[i]] - (double) update[i]);
     if (x[i] > 1.0E-8) diff = diff / (double) update[i];
     if (diff < 0.0) diff = -diff;
     if ( diff > 1.0E-5)
        printf("%d : data %d = %e %e \n", proc_config[AZ_node], i, 
               x[update_index[i]], (double) update[i]);
  }
  */
#ifdef ML_BENCHMARK
   if (proc_config[AZ_node] == 0) 
     printf("Printing out a few entries of the solution ...\n");

   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 7) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],x[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 23) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],x[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 47) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],x[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 101) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],x[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 171) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],x[update_index[j]]); fflush(stdout);}
#endif


  ML_Destroy(&ml);
  free( (void *) Amat->data_org);
  free( (void *) Amat->val);
  free( (void *) Amat->bindx);
  if (Amat  != NULL) AZ_matrix_destroy(&Amat);
  if (Pmat != NULL) AZ_precond_destroy(&Pmat); 

  /* Free allocated memory */

  free((void *) update);    free((void *) update_index);
  free((void *) external);  free((void *) extern_index);
  free((void *) x);         free((void *) b);           

#ifdef ML_MPI
  MPI_Finalize();
#endif
  return 0;
}

/*****************************************************************************/
/* more global variables                                                     */
/*****************************************************************************/

int local_index_f_leng, *local_index_f;
int local_index_c_leng, *local_index_c;
int *ML_update;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void init_options(int options[], double params[])
{

  /* Choose among AZTEC options (see User's Guide). */

  AZ_defaults(options, params);

  options[AZ_solver]   = AZ_cg;
#ifdef ML_BENCHMARK
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
#endif
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_user_precond;
  options[AZ_conv]     = AZ_r0;
  options[AZ_output]   = 1;
  options[AZ_pre_calc] = AZ_calc;
  options[AZ_max_iter] = 1550;
  options[AZ_poly_ord] = 5;
  options[AZ_overlap]  = AZ_none;
  options[AZ_kspace]   = 130;
  options[AZ_orthog]   = AZ_modified;
  options[AZ_aux_vec]  = AZ_resid;

  params[AZ_tol]       = 1.00e-8;
  params[AZ_drop]      = 0.0;
  params[AZ_ilut_fill] = 1;
  params[AZ_omega]     = 1.;

} /* init_options */

/******************************************************************************/
/******************************************************************************/
/* Set the initial guess and the right hand side where the right hand side
 * is obtained by doing a matrix-vector multiplication.
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date :  3/15/95
 *
 * Parameters
 *
 *    update_index   ==      On input, ordering of update and external
 *                           locally on this processor. For example
 *                           'update_index[i]' gives the index location
 *                           of the block which has the global index
 *                           'update[i]'.
 *    update         ==      On input, list of pts to be updated on this node
 *    data_org       ==      On input, indicates how data is set on this node.
 *                           For example, data_org[] contains information on
 *                           how many unknowns are internal, external and
 *                           border unknowns as well as which points need
 *                           to be communicated. See User's Guide for more
 *                           details.
 *    val, indx,     ==      On input, holds matrix nonzeros. See User's Guide
 *    bindx, rpntr,          for more details.
 *    cpntr, bpntr
 *    x              ==      On output, 'x' is allocated and set to all zeros.
 *    b              ==      On output, 'b' is allocated and is set to the
 *                           result of a matrix-vector product.
 ******************************************************************************/
void init_guess_and_rhs(int update_index[], int update[], double *x[],double
                        *b[],int data_org[], double val[], int indx[], int
                        bindx[], int rpntr[], int cpntr[], int bpntr[], int
                        proc_config[])
{

  int    i,j;
  int    temp,num;
  AZ_MATRIX Amat;

  temp = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  num  = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* allocate vectors */

  i       = num + data_org[AZ_N_external];
  *x      = (double *) calloc(i, sizeof(double));
  *b      = (double *) calloc(i, sizeof(double));
  if ((*b == NULL) && (i != 0)) {
    (void) fprintf(stderr, "Not enough space in init_guess_and_rhs() for ax\n");
    exit(1);
  }

  /* initialize 'x' to a function which will be used in matrix-vector product */

  for (i = 0; i < temp; i++) 
     for (j = 0; j < num_PDE_eqns; j++) 
        (*x)[i*num_PDE_eqns+j] = (double) (update[i]);
  /* for (i = 0; i < temp; i++) (*x)[i] = (double) 1; */

  /* Reorder 'x' so that it conforms to the transformed matrix */
 
  AZ_reorder_vec(*x,data_org,update_index,rpntr);

  /* take out the constant vector. Used for the */
  /* finite element problem because it is singular */

  Amat.rpntr      = rpntr;   Amat.cpntr    = cpntr;
  Amat.bpntr      = bpntr;   Amat.bindx    = bindx;
  Amat.indx       = indx;    Amat.val      = val;
  Amat.data_org   = data_org;
  Amat.matrix_type = data_org[AZ_matrix_type];
  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
      Amat.matvec = AZ_MSR_matvec_mult;
  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      Amat.matvec = AZ_VBR_matvec_mult;

  Amat.matvec(*x, *b, &Amat, proc_config);

  for (i = 0; i < num; i++) (*x)[i] = 0.0;
  ML_update = update;

} /* init_guess_and_rhs */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int construct_ml_grids(int N_elements, int *proc_config, AZ_MATRIX **Amat_f, 
              int *update_f[], int *update_index_f[], int *external_f[], 
              int *extern_index_f[], int input_option, int *indx[], 
              int *bpntr[], int *rpntr[], int *cpntr[], ML **ml_ptr, 
              int *options, double *params, AZ_PRECOND **Pmat, int flag)
{
  int         MSRorVBR = AZ_MSR_MATRIX, N_update_f;
  int         N2, *bindx, *data_org, leng;
  double      *val;

  ML           *ml;
  ML_GridAGX   *f_mesh;
  int N_levels = 5, level, coarsest_level;
  ML_Aggregate *ag;
  /*
  int          nblocks, *blocks;
  */
  /*  hacker stuff
  ML          *subml;
  struct aztec_context *context;
  */
#ifdef ML_BENCHMARK
  int nblocks, *blocks;
#endif

  /**************************************************************************/

  N_levels = num_levels;

  /* generate the fine grid and its partitionings */

  N2 = N_elements * refine_factor * refine_factor * refine_factor;
  Generate_mesh(N2, &f_mesh, &N_update_f, update_f, proc_config);

  /* use the partitioning for the fine grid to initialize matrix */

  N_grid_pts = AZ_gsum_int(N_update_f, proc_config);
  create_msr_matrix(*update_f,&val,&bindx,N_update_f);
  AZ_transform(proc_config, external_f, bindx, val, *update_f, 
               update_index_f, extern_index_f, &data_org, N_update_f, 
               *indx, *bpntr, *rpntr, cpntr, MSRorVBR);

  *Amat_f = AZ_matrix_create(data_org[AZ_N_internal] + data_org[AZ_N_border]);
  AZ_set_MSR(*Amat_f, bindx, val, data_org, 0, NULL, AZ_LOCAL);      
  (*Amat_f)->matrix_type  = data_org[AZ_matrix_type];
  data_org[AZ_N_rows]  = data_org[AZ_N_internal] + data_org[AZ_N_border];
  ML_GridAGX_Destroy(&f_mesh);

  if (flag == 1) {

    /* initialization and set up scheme */

    ML_Create(&ml, N_levels);
    ML_Set_PrintLevel(3);
    (*ml_ptr) = ml;

    /* set up processor information */

    ML_Set_Comm_MyRank(ml, proc_config[AZ_node]);
    ML_Set_Comm_Nprocs(ml, proc_config[AZ_N_procs]);
#ifdef ML_MPI
    ML_Set_Comm_Communicator(ml, MPI_COMM_WORLD);
#else
    ML_Set_Comm_Communicator(ml, 91);
#endif
    ML_Set_Comm_Send(ml, ml_send);
    ML_Set_Comm_Recv(ml, ml_recv);
    ML_Set_Comm_Wait(ml, ml_wait);

    /* set up discretization matrix and matrix vector function */

    leng = (*Amat_f)->data_org[AZ_N_internal] + 
           (*Amat_f)->data_org[AZ_N_border];
    AZ_ML_Set_Amat(ml, N_levels-1, leng, leng, *Amat_f, proc_config);

    ML_Aggregate_Create( &ag );
    ML_Aggregate_Set_MaxCoarseSize(ag, 10);
    ML_Aggregate_Set_Threshold(ag, 0.0);
    coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, N_levels-1, ML_DECREASING,ag);
    coarsest_level = N_levels - coarsest_level;
    if ( proc_config[AZ_node] == 0 )
       printf("Coarse level = %d \n", coarsest_level);

    /* set up smoothers */

    for (level = N_levels-1; level > coarsest_level; level--) {
       /* Block Gauss-Seidel using all variables within an aggregate */
       /* as a block.                                                */
       /*

       nblocks = ML_Aggregate_Get_AggrCount( ag, level );
       ML_Aggregate_Get_AggrMap( ag, level, &blocks);
       ML_Gen_Smoother_VBlockSymGaussSeidel( ml , level, ML_PRESMOOTHER, 
                                             nsmooth, 1.0, 
                                             nblocks, blocks);
       ML_Gen_Smoother_VBlockSymGaussSeidel( ml , level, ML_POSTSMOOTHER, 
 					     nsmooth, 1.0, 
                                             nblocks, blocks);
       */


       /* This is the symmetric Gauss-Seidel smoothing. In parallel,    */
       /* it is not a true Gauss-Seidel in that each processor          */
       /* does a Gauss-Seidel on its local submatrix independent of the */
       /* other processors.                                             */

#ifndef ML_BENCHMARK
       ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_POSTSMOOTHER, nsmooth,1.);
       ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_PRESMOOTHER, nsmooth,1.);

#else
	 if (ML_strcmp(context->smoother,"Parasails") == 0) {
	   ML_Gen_Smoother_ParaSails(ml , level, ML_PRESMOOTHER, nsmooth,
				     parasails_sym, parasails_thresh,
				     parasails_nlevels, parasails_filter,
				     parasails_loadbal, parasails_factorized);
	 }

	 /* This is the symmetric Gauss-Seidel smoothing that we usually use. */
	 /* In parallel, it is not a true Gauss-Seidel in that each processor */
	 /* does a Gauss-Seidel on its local submatrix independent of the     */
	 /* other processors.                                                 */

	 else if (ML_strcmp(context->smoother,"MLS") == 0) {
	    ML_Gen_Smoother_MLS(ml, level, ML_BOTH, 10.,nsmooth);
	 }
	 else if (ML_strcmp(context->smoother,"GaussSeidel") == 0) {
	   ML_Gen_Smoother_GaussSeidel(ml , level, ML_BOTH, nsmooth,1.);
	 }
	 else if (ML_strcmp(context->smoother,"SymGaussSeidel") == 0) {
	   ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_BOTH, nsmooth,1.);
	 }
	 else if (ML_strcmp(context->smoother,"BlockGaussSeidel") == 0) {
	   ML_Gen_Smoother_BlockGaussSeidel(ml , level, ML_BOTH, nsmooth,1.,
					    num_PDE_eqns);
	 }
	 else if (ML_strcmp(context->smoother,"Aggregate") == 0) {
	   ML_Gen_Blocks_Aggregates(ag, level, &nblocks, &blocks);
	   ML_Gen_Smoother_VBlockSymGaussSeidel(ml , level, ML_BOTH, nsmooth,1.,
						nblocks, blocks);
	 }

	 /* This is a true Gauss Seidel in parallel. This seems to work for  */
	 /* elasticity problems.  However, I don't believe that this is very */
	 /* efficient in parallel.                                           */
	 /*
	   nblocks = ml->Amat[level].invec_leng;
	   for (i =0; i < nblocks; i++) blocks[i] = i;
	   ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml , level, ML_PRESMOOTHER,
	   nsmooth, 1., nblocks, blocks);
	   ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml, level, ML_POSTSMOOTHER,
	   nsmooth, 1., nblocks, blocks);
	 */

	 /* Jacobi Smoothing                                                 */

	 else if (ML_strcmp(context->smoother,"Jacobi") == 0) {
	   ML_Gen_Smoother_Jacobi(ml , level, ML_PRESMOOTHER, nsmooth,.4);
	   ML_Gen_Smoother_Jacobi(ml , level, ML_POSTSMOOTHER, nsmooth,.4);
	 }

	 /*  This does a block Gauss-Seidel (not true GS in parallel)        */
	 /*  where each processor has 'nblocks' blocks.                      */

	 else if (ML_strcmp(context->smoother,"Metis") == 0) {
	   nblocks = 250;
	   ML_Gen_Blocks_Metis(ml, level, &nblocks, &blocks);
	   ML_Gen_Smoother_VBlockSymGaussSeidel(ml , level, ML_BOTH, nsmooth,1.,
						nblocks, blocks);
	 }
	 else {
	   printf("unknown smoother %s\n",context->smoother);
	   exit(1);
	 }
#endif
       /*  Sparse approximate inverse smoother that acutally does both */
       /*  pre and post smoothing.                                     */

       /*
       ML_Gen_Smoother_ParaSails(ml , level, ML_PRESMOOTHER, nsmooth,
                                parasails_sym, parasails_thresh,
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
       */

       /*
       ML_Gen_Smoother_Jacobi(ml , level, ML_PRESMOOTHER, 1, 0.67 );
       ML_Gen_Smoother_Jacobi(ml , level, ML_POSTSMOOTHER, 2, 0.67 );
       */
    }

#ifndef ML_BENCHMARK
    if (coarse_iterations == 0) ML_Gen_CoarseSolverSuperLU(ml,coarsest_level);
    else {
       ML_Gen_Smoother_SymGaussSeidel(ml, coarsest_level, ML_PRESMOOTHER, 
                                    coarse_iterations,1.);
    /* Some advanced options for speeding up coarse grid solves  */
    /*                                                           */
    /* This one is used to block the matrix (throw away nonzeros */
    /* between blocks) before calling superlu.                   */

    /* ml->Amat[coarsest_level].invec_leng = -7;
       ML_Gen_CoarseSolverSuperLU(ml,coarsest_level);
       ml->Amat[coarsest_level].invec_leng=ml->Amat[coarsest_level].outvec_leng;
     */

    /* This one is used to use Aztec on a global matrix. That is */
    /* replicate global matrix on each processor and then apply  */
    /* Aztec to solve it.                                        */

    /* oldvalue = options[AZ_max_iter];
       options[AZ_max_iter] = -42;
       ML_Gen_SmootherAztec(ml, coarsest_level, options, params, 
            proc_config, status, 2, ML_PRESMOOTHER,NULL);
       options[AZ_max_iter] = oldvalue;
     */

    /* This one is used to use Aztec with superlu as a preconditoner */
    /* when we have dropped values from coarse grid matrix. To do    */
    /* this, we                                                      */
    /*      1) replicate matrix within a new ml structure (subml).   */
    /*      2) invoke superlu on matrix within new structure.        */
    /*      3) Set Aztec smoother for original ml structure.         */
    /*      4) Hard wire ML function and data pointer into Aztec     */
    /*         data structure.                                       */
    /*
       ML_Create(&subml,1);
       ML_Init_Amatrix(subml, 0, ml->Amat[coarsest_level].invec_leng,
                ml->Amat[coarsest_level].invec_leng,
                ml->Amat[coarsest_level].data);
       ML_Operator_Set_Getrow( &(subml->Amat[0]), ML_EXTERNAL, 
                        subml->Amat[0].outvec_leng, 
			ml->Amat[coarsest_level].getrow->external);
       ML_CommInfoOP_Clone( &(subml->Amat[0].getrow->pre_comm), 
                            ml->Amat[coarsest_level].getrow->pre_comm);

       if (ml->Amat[coarsest_level].matvec->ML_id == ML_EXTERNAL)
            ML_Operator_Set_ApplyFunc(&(subml->Amat[0]),ML_EXTERNAL,
                                 ml->Amat[coarsest_level].matvec->external);
       else ML_Operator_Set_ApplyFunc(&(subml->Amat[0]),ML_INTERNAL,
                                 ml->Amat[coarsest_level].matvec->internal);

       subml->Amat[0].invec_leng = -7;
       ML_Gen_CoarseSolverSuperLU(subml,0);
       subml->Amat[0].invec_leng = subml->Amat[0].outvec_leng;
       ML_Gen_Solver( subml, ML_MGV, 0, 0);

       options[AZ_precond]    = AZ_user_precond;
       ML_Gen_SmootherAztec(ml, coarsest_level, options, params, proc_config,
                   status, 2, ML_PRESMOOTHER,ML_precondition);

       context = (struct aztec_context *) 
                         ml->pre_smoother[coarsest_level].smoother->data;
       context->Prec->ml_ptr = subml;
       context->Prec->precond_data = subml;
       context->Prec->prec_function = ML_precondition;
     */
    }
#else
      nsmooth   = context->coarse_its;
      /*  Sparse approximate inverse smoother that acutally does both */
      /*  pre and post smoothing.                                     */

      if (ML_strcmp(context->coarse_solve,"Parasails") == 0) {
        ML_Gen_Smoother_ParaSails(ml , coarsest_level, ML_PRESMOOTHER, nsmooth,
				  parasails_sym, parasails_thresh,
				  parasails_nlevels, parasails_filter,
				  parasails_loadbal, parasails_factorized);
      }

      else if (ML_strcmp(context->coarse_solve,"MLS") == 0) {
	ML_Gen_Smoother_MLS(ml, coarsest_level, ML_BOTH, 20.,nsmooth);
      }
      else if (ML_strcmp(context->coarse_solve,"GaussSeidel") == 0) {
	ML_Gen_Smoother_GaussSeidel(ml , coarsest_level, ML_BOTH, nsmooth,1.);
      }
      else if (ML_strcmp(context->coarse_solve,"SymGaussSeidel") == 0) {
	ML_Gen_Smoother_SymGaussSeidel(ml , coarsest_level, ML_BOTH, nsmooth,1.);
      }
      else if (ML_strcmp(context->coarse_solve,"BlockGaussSeidel") == 0) {
	ML_Gen_Smoother_BlockGaussSeidel(ml, coarsest_level, ML_BOTH, nsmooth,1.,
					 num_PDE_eqns);
      }
      else if (ML_strcmp(context->coarse_solve,"Aggregate") == 0) {
	ML_Gen_Blocks_Aggregates(ag, coarsest_level, &nblocks, &blocks);
	ML_Gen_Smoother_VBlockSymGaussSeidel(ml , coarsest_level, ML_BOTH,
					     nsmooth,1., nblocks, blocks);
      }
      else if (ML_strcmp(context->coarse_solve,"Jacobi") == 0) {
        ML_Gen_Smoother_Jacobi(ml , coarsest_level, ML_BOTH, nsmooth,.5);
      }
      else if (ML_strcmp(context->coarse_solve,"Metis") == 0) {
	nblocks = 250;
	ML_Gen_Blocks_Metis(ml, coarsest_level, &nblocks, &blocks);
	ML_Gen_Smoother_VBlockSymGaussSeidel(ml , coarsest_level, ML_BOTH,
					     nsmooth,1., nblocks, blocks);
      }
      else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0) {
	ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
      }
      else {
	printf("unknown coarse grid solver %s\n",context->coarse_solve);
	exit(1);
      }
#endif


    ML_Gen_Solver( ml, ML_MGV, N_levels-1, coarsest_level);
    ML_Aggregate_Destroy(&ag);

  }
  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Generate_mesh(int N_elements, ML_GridAGX **meshp, int *N_update, 
                   int *update[], int proc_config[])
{
   int    i, j, k, nprocs, nprocs_1d, mypid, mypid_y, mypid_z;
   int    icnt, *intlist, offset, nelmnt_1d, nelmnt_part_yz, estart_y; 
   int    estart_z, nelmnt_local, nstart_y, nstart_z, index, eend_y, eend_z;
   int    nelmnt_2d, nnode_1d, nnode_local_y, nnode_local_z, nnode_local;
   int    nnode_2d, nnode_expand, nnode_x, nnode_y, nnode_z, ***cube;
   int    y_begin, z_begin;
   double *coord, h;
   ML_GridAGX *mesh;

   /* --------------------------------------------------------------- */
   /* box in the processors                                           */
   /* --------------------------------------------------------------- */

   mypid    = proc_config[AZ_node];
   nprocs   = proc_config[AZ_N_procs];
   
   nprocs_1d = (int) pow( (double) nprocs, 0.50001 );
   if ( nprocs_1d * nprocs_1d != nprocs ) {
      printf("Error: nprocs should be a square (%d).\n",nprocs_1d);
      exit(1);
   }
   mypid_y = mypid % nprocs_1d;
   mypid_z = mypid / nprocs_1d;

   /* --------------------------------------------------------------- */
   /* process the grid                                                */
   /* --------------------------------------------------------------- */

   nelmnt_1d = (int) pow( (double) N_elements, 0.333334 );

   /* partition the elements along the yz-axis */

   nelmnt_part_yz = nelmnt_1d / nprocs_1d;
   if (nelmnt_part_yz * nprocs_1d != nelmnt_1d) {
      printf("Error: nelmnt_part not good. %d %d\n",nelmnt_part_yz,nelmnt_1d);
      exit(-1);
   }
   nelmnt_local = nelmnt_part_yz * nelmnt_part_yz * nelmnt_1d;
  
   /* initialize the mesh */

   ML_GridAGX_Create( meshp );
   mesh = (*meshp);
   h = 1.0 / (double) nelmnt_1d;
   ML_GridAGX_Set_Dimension(mesh, 3);
   ML_GridAGX_Set_NElmnts(mesh, nelmnt_local, 8);

   /* generate the global element numbers of the local elements */

   ML_memory_alloc((void**) &intlist, nelmnt_local * sizeof(int), "AP1");
   estart_y  = nelmnt_part_yz * mypid_y;
   estart_z  = nelmnt_part_yz * mypid_z;
   eend_y    = estart_y + nelmnt_part_yz - 1;
   eend_z    = estart_z + nelmnt_part_yz - 1;
   icnt      = 0;
   nelmnt_2d = nelmnt_1d * nelmnt_1d;
   for (k = estart_z; k <= eend_z; k++) 
      for (j = estart_y; j <= eend_y; j++) 
         for (i = 0; i < nelmnt_1d; i++)
            intlist[icnt++] = i + j * nelmnt_1d + k * nelmnt_2d;
   ML_GridAGX_Load_ElmntGlobalNum(mesh, nelmnt_local, intlist);
   ML_memory_free( (void **) &intlist);

   /* generate the local computational grid */

   nnode_x       = nelmnt_1d + 1;
   nnode_y       = nelmnt_part_yz + 1;
   nnode_z       = nelmnt_part_yz + 1;
   ML_memory_alloc((void**) &cube, nnode_x * sizeof(int**), "AP2");
   for ( i = 0; i < nnode_x; i++ )
      ML_memory_alloc((void**) &(cube[i]), nnode_y * sizeof(int*), "AP3");
   for ( i = 0; i < nnode_x; i++ )
      for ( j = 0; j < nnode_y; j++ )
         ML_memory_alloc((void**) &(cube[i][j]),nnode_z*sizeof(int),"AP4");
   icnt = 0;
   y_begin = 1;
   z_begin = 1;
   if ( estart_y == 0 ) y_begin--;
   if ( estart_z == 0 ) z_begin--;
   for ( k = z_begin; k < nnode_z; k++ )
      for ( j = y_begin; j < nnode_y; j++ )
         for ( i = 0; i < nnode_x; i++ ) cube[i][j][k] = icnt++;
   if ( y_begin == 1 || z_begin == 1 ) {
      for ( i = 0; i < nnode_x; i++ ) cube[i][0][0] = icnt++;
   }
   if ( y_begin == 1 ) {
      for ( k = 1; k < nnode_z; k++ )
         for ( i = 0; i < nnode_x; i++ )
            cube[i][0][k] = icnt++;
   }
   if ( z_begin == 1 ) {
      for ( j = 1; j < nnode_y; j++ )
         for ( i = 0; i < nnode_x; i++ )
            cube[i][j][0] = icnt++;
   }
   
   /* generate the global node IDs */

   nnode_local_y = nnode_y - 1;
   nnode_local_z = nnode_z - 1;
   if (mypid_y == 0 ) nnode_local_y++;
   if (mypid_z == 0 ) nnode_local_z++;
   nnode_1d = nelmnt_1d + 1;
   nnode_local  = nnode_local_y * nnode_local_z * nnode_1d;
   (*N_update)  = nnode_local;
   nnode_expand = nnode_y * nnode_z * nnode_x;
   (*update)   = (int*) malloc(nnode_expand * sizeof(int));

   nstart_y    = mypid_y * nnode_local_y;
   nstart_z    = mypid_z * nnode_local_z;
   nnode_2d    = nnode_1d * nnode_1d;
   for (k = 0; k < nnode_z; k++) 
      for (j = 0; j < nnode_y; j++) 
         for (i = 0; i < nnode_x; i++) { 
            index = cube[i][j][k];
            offset  = i + (j+nstart_y)*nnode_1d + (k+nstart_z)*nnode_2d;
            (*update)[index] = offset;
         }

   ML_GridAGX_Set_NVert(mesh, *N_update);
   ML_GridAGX_Load_VertGlobalNum(mesh, nnode_expand, *update);

   /* generate coordinates of local nodes */

   ML_memory_alloc((void**) &coord, 3*nnode_expand*sizeof(double),"AP6");
   for (k = 0; k < nnode_z; k++) 
      for (j = 0; j < nnode_y; j++) 
         for (i = 0; i < nnode_x; i++) { 
            index = cube[i][j][k];
            offset  = i + (j+nstart_y)*nnode_1d + (k+nstart_z)*nnode_2d;
            coord[index*3]   = (double) (h * i);
            coord[index*3+1] = (double) (h * (j + nstart_y));
            coord[index*3+2] = (double) (h * (k + nstart_z));
         } 
   ML_GridAGX_Load_AllVertCoordinates(mesh, nnode_expand, coord);
   ML_memory_free( (void **) &coord);

   /* finally construct the element-to-node lists */

   ML_memory_alloc((void**) &intlist, 8 * sizeof(int), "AP7");

   k = (eend_z - estart_z + 1) * (eend_y - estart_y + 1) * nelmnt_1d;
   for (k = estart_z; k <= eend_z; k++) 
      for (j = estart_y; j <= eend_y; j++) 
         for (i = 0; i < nelmnt_1d; i++) {
            intlist[0] = cube[i][j-estart_y][k-estart_z];
            intlist[1] = cube[i+1][j-estart_y][k-estart_z];
            intlist[2] = cube[i][j-estart_y+1][k-estart_z];
            intlist[3] = cube[i+1][j-estart_y+1][k-estart_z];
            intlist[4] = cube[i][j-estart_y][k-estart_z+1];
            intlist[5] = cube[i+1][j-estart_y][k-estart_z+1];
            intlist[6] = cube[i][j-estart_y+1][k-estart_z+1];
            intlist[7] = cube[i+1][j-estart_y+1][k-estart_z+1];
            ML_IntList_Load_Sublist(mesh->ele_nodes, 8, intlist);
         }  
   ML_memory_free( (void **) &intlist);
   for ( i = 0; i < nnode_x; i++ )
      for ( j = 0; j < nnode_y; j++ )
         ML_memory_free( (void**) &(cube[i][j]) );
   for ( i = 0; i < nnode_x; i++ )
      ML_memory_free( (void**) &(cube[i]));
   ML_memory_free( (void**) &cube);
}


int ML_mapper0f( void *grid, double *din, double *dout)
{
   int        i, k, leng, index, step, istep, jstep;
   ML_GridAGX *mlgrid;
   
   step = num_PDE_eqns;
   mlgrid = (ML_GridAGX *) grid;
   leng = ML_GridAGX_Get_NVert(mlgrid);
   for (i = 0; i < leng*step; i++) dout[i] = 0.0;
   for (i = 0; i < local_index_c_leng; i++) {
      index = local_index_c[i];
      istep = i * step;
      jstep = index * step;
      for ( k = 0; k < step; k++ ) dout[istep+k] = din[jstep+k];
   }
   return 0;
} 

int ML_mapper0b( void *ml, double *din, double *dout)
{
   int     i, k, index, step, istep, jstep;
   
   k = (int) ml;
   step = num_PDE_eqns;
   for (i = 0; i < local_index_c_leng; i++) {
      index = local_index_c[i];
      istep = i * step;
      jstep = index * step;
      for ( k = 0; k < step; k++ ) dout[jstep+k] = din[istep+k];
   }
   return 0;
}

int ML_mapper1f( void *grid, double *din, double *dout)
{
   int        i, k, leng, index, step, istep, jstep;
   ML_GridAGX *mlgrid;
   
   step = num_PDE_eqns;
   mlgrid = (ML_GridAGX *) grid;
   leng = ML_GridAGX_Get_NVert(mlgrid);
   for (i = 0; i < leng*step; i++) dout[i] = 0.0;
   for (i = 0; i < local_index_f_leng; i++) {
      index = local_index_f[i];
      istep = i * step;
      jstep = index * step;
      for ( k = 0; k < step; k++ ) dout[istep+k] = din[jstep+k];
   }
   return 0;
}

int ML_mapper1b( void *ml, double *din, double *dout)
{
   int     i, k, index, step, istep, jstep;
   
   k = (int) ml;
   step = num_PDE_eqns;
   for (i = 0; i < local_index_f_leng; i++) {
      index = local_index_f[i];
      istep = i * step;
      jstep = index * step;
      for ( k = 0; k < step; k++ ) dout[jstep+k] = din[istep+k];
   }
   return 0;
}

/**************************************************************************/
/* communication subroutines                                              */
/*------------------------------------------------------------------------*/
#ifndef ML_MPI
#define MPI_Request int
#endif
int ml_recv(void* buf, unsigned int count, int *src, int *mid, 
            USR_COMM comm, USR_REQ *request )
{
   return(md_wrap_iread( buf, count, src, mid, (MPI_Request *) request ));
}

int ml_wait(void* buf, unsigned int count, int *src, int *mid, 
            USR_COMM comm, USR_REQ *request )
{
   int status;
   return(md_wrap_wait( buf, count, src, mid, &status, (MPI_Request *) request ));
}

int ml_send(void* buf, unsigned int count, int dest, int mid, 
            USR_COMM comm )
{
   int status;
   return(md_wrap_write( buf, count, dest, mid, &status));
}

/*****************************************************************************
 *                             MSR                                           *
 *****************************************************************************/

void create_msr_matrix(int update[], double **val, int **bindx, int N_update)

{
  int i,total_nz, n = -1;
  int avg_nonzeros_per_row = 27;

  total_nz = N_update*avg_nonzeros_per_row + 1;
  *bindx   = (int *) AZ_allocate(total_nz*sizeof(int));
  *val     = (double *) AZ_allocate(total_nz*sizeof(double));
  if ((*val == NULL) && (total_nz != 0) ) {
    (void) fprintf(stderr, "Error: Not enough space to create matrix\n");
    (void) fprintf(stderr,
                   "      Try reducing the variable 'avg_nonzeros_per_row'\n");
    exit(1);
  }
  for (i = 0 ; i < total_nz ; i++ ) (*bindx)[i] = 0;

  (*bindx)[0] = N_update+1;
  for (i = 0; i < N_update; i++) {
    add_row_3D(update[i], i, *val, *bindx,&n);
    if ( (*bindx)[i+1] > total_nz) {
      (void) fprintf(stderr, "Error:total_nz not large enough to accomodate");
      (void) fprintf(stderr, " nonzeros\n       Try increasing the variable");
      (void) fprintf(stderr, " 'avg_nonzeros_per_row'\n       Finished \n");
      (void) fprintf(stderr, "first %d rows out of %d\n", i, N_update);
      exit(1);
    }
  }

} /* create_msr_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void add_row_3D(int row, int location, double val[], int bindx[], int *n)

/*
 * Add one row to an MSR matrix corresponding to a 7pt discrete approximation
 * to the 3D Poisson operator on an n x n x n square.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/15/95
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    val,bindx    == (see Aztec User's guide). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{

  int        m;
  int        k, NP;
int nx,ny,nz;

  /* determine grid dimensions */

  m = *n;
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) pow( ((double) m )+0.01, 0.33334);
    if (m == 1) {

      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }

    if (m*m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^3 where k is an integer\n");
      exit(1);
    }
  }
  *n = m;
  k = bindx[location];
  NP = num_PDE_eqns;

  /*
   * Check neighboring points in each direction and add nonzero entry if
   * neighbor exists.
   */

  bindx[k] = row + NP;     if ((row/NP)%m      != m-1) val[k++] = -1.0;
  bindx[k] = row - NP;     if ((row/NP)%m      != 0  ) val[k++] = -1.0;
  bindx[k] = row + m*NP;   if ((row/(NP*m))%m  != m-1) val[k++] = -1.0;
  bindx[k] = row - m*NP;   if ((row/(NP*m))%m  != 0  ) val[k++] = -1.0;
  bindx[k] = row + m*m*NP; if (bindx[k]    < m*m*m*NP) val[k++] = -1.0;
  bindx[k] = row - m*m*NP; if (bindx[k]          >= 0) val[k++] = -1.0;

nx = (row/NP)%m;
ny = (row/(NP*m))%m;
nz = (row/(NP*m*m));

/*
  if (nx != m-1) {
     bindx[k] = row + NP                ; val[k++] = -1.0;
     if (ny != m-1) {
        bindx[k] = row + NP + m*NP         ; val[k++] = -.5;
        if (nz != m-1) {
           bindx[k] = row + NP + m*NP + m*m*NP; val[k++] = -.25;
        }
        if (nz != 0) {
           bindx[k] = row + NP + m*NP - m*m*NP; val[k++] = -.25;
        }
     }
     if (nz != m-1) {
        bindx[k] = row + NP        + m*m*NP; val[k++] = -0.5;
     }
     if (nz != 0) {
        bindx[k] = row + NP        - m*m*NP; val[k++] = -0.5;
     }
     if (ny != 0) {
        bindx[k] = row + NP - m*NP         ; val[k++] = -0.5;
        if (nz != m-1) {
           bindx[k] = row + NP - m*NP + m*m*NP; val[k++] = -0.25;
        }
        if (nz != 0) {
           bindx[k] = row + NP - m*NP - m*m*NP; val[k++] = -0.25;
        }
     }
  }
  if (nx != 0) {
     bindx[k] = row - NP                ; val[k++] = -1.0;
     if (ny != m-1) {
        bindx[k] = row - NP + m*NP         ; val[k++] = -0.5;
        if (nz != m-1) {
           bindx[k] = row - NP + m*NP + m*m*NP; val[k++] = -0.25;
        }
        if (nz != 0) {
           bindx[k] = row - NP + m*NP - m*m*NP; val[k++] = -0.25;
        }
     }
     if (nz != m-1) {
        bindx[k] = row - NP        + m*m*NP; val[k++] = -0.5;
     }
     if (nz != 0) {
        bindx[k] = row - NP        - m*m*NP; val[k++] = -0.5;
     }
     if (ny != 0) {
        bindx[k] = row - NP - m*NP         ; val[k++] = -0.5;
        if (nz != m-1) {
           bindx[k] = row - NP - m*NP + m*m*NP; val[k++] = -0.25;
        }
        if (nz != 0) {
           bindx[k] = row - NP - m*NP - m*m*NP; val[k++] = -0.25;
        }
     }
  }
     if (ny != m-1) {
        bindx[k] = row      + m*NP         ; val[k++] = -1.0;
        if (nz != m-1) {
           bindx[k] = row      + m*NP + m*m*NP; val[k++] = -0.5;
        }
        if (nz != 0) {
           bindx[k] = row      + m*NP - m*m*NP; val[k++] = -0.5;
        }
     }
     if (nz != m-1) {
        bindx[k] = row             + m*m*NP; val[k++] = -1.0;
     }
     if (nz != 0) {
        bindx[k] = row             - m*m*NP; val[k++] = -1.0;
     }
     if (ny != 0) {
        bindx[k] = row      - m*NP         ; val[k++] = -1.0;
        if (nz != m-1) {
           bindx[k] = row      - m*NP + m*m*NP; val[k++] = -0.5;
        }
        if (nz != 0) {
           bindx[k] = row      - m*NP - m*m*NP; val[k++] = -0.5;
        }
     }
*/


  bindx[location+1] = k;
  val[location]     = 6.0;
/*
  val[location]     = 14.0;
*/

} /* add_row_3D */

#else

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])
{

  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("This test requires Aztec");

#ifdef ML_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif

