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
/* date:         5/10/99                                                     */
/*****************************************************************************/

/*
#define MATLAB
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_include.h"

#if defined(HAVE_ML_AZTEC) || defined(HAVE_ML_AZTECOO) || defined(HAVE_ML_AZTEC2_1)

#include "az_aztec.h"
#include "ml_read_utils.h"

extern int AZ_using_fortran;
int    parasails_factorized = 0;
int    parasails_sym        = 1;
double parasails_thresh     = 0.01;
int    parasails_nlevels    = 0;
double parasails_filter     = 0.;
double parasails_loadbal    = 0.;
double *scaling_vect = NULL;
  /* -------------  external function declarations -------------------------*/

extern void add_row_5pt(int, int , double val[], int bindx[], int *n);

extern void init_options(int options[], double params[]);

extern void init_guess_and_rhs(int update_index[],int update[], double *x[],
			       double *b[], int data_org[], double val[], int indx[], 
			       int bindx[], int rpntr[], int cpntr[], int bpntr[], 
			       int proc_config[]);

extern int construct_ml_grids(int, int *, AZ_MATRIX **, int **, int **, 
			      int **, int **, int **, int **, int **, int **, 
			      ML **, int, ML_Aggregate **);
extern void create_msr_matrix(int*, double **, int **, int);
extern void create_rot_msr_matrix(int*, double **, int **, int);
extern void create_msr_matrix_conv_diff(int*, double **, int **, int,
					double, double);
extern int  ml_recv(void*, unsigned int, int *, int *, USR_COMM, USR_REQ *);
extern int  ml_wait(void*, unsigned int, int *, int *, USR_COMM, USR_REQ *);
extern int  ml_send(void*, unsigned int, int, int, USR_COMM );
extern void Generate_mesh(int, ML_GridAGX **, int *, int **, int *);
extern void ml_print_global_element(int element,int *update,int *data_org,
				    int *update_index, int *vector,int *proc_config, 
				    int *external, int *extern_index);
extern int  ml_find_global_row(int i, int proc, int global_nrows, 
			       int **whole_glob_map);

extern void add_row_rot(int, int , double val[], int bindx[], int *n);

extern void add_row_7pt(int, int , double val[], int bindx[], int *n, 
			double alpha, double beta);
extern void AZ_print_out2(int update_index[], int extern_index[], int update[],
                      int external[], double val[], int indx[],  int
                      bindx[], int rpntr[], int cpntr[], int bpntr[], int
                      proc_config[], int choice, int matrix, int N_update,
                      int N_external, int off_set );
extern double ml_variable_coef(double x, double y, int m);

/****************************************************************************/
/* To Run with a rotated null space do the following:                       */
/*     1) change num_PDE_eqns to 2                                          */
/*     2) uncomment 'null_vect' code just before ML_Aggregate_Set_NullSpace */
/****************************************************************************/

int N_grid_pts;
int refine_factor;
int choice;
int num_PDE_eqns=1;
int nsmooth=1;
/* for rotated nullspace test: */
double c, s, theta=0.39710678;
double agg_thresh;

#ifdef PARTEST
int **global_mapping;
int global_nrows, global_ncoarse;
#endif


/*****************************************************************************/
/* Set up and solve a test problem defined in the subroutine                 */
/* init_matrix_vector_structures().                                          */
/*****************************************************************************/
int coarse_iterations = 0, use_cg = 0, num_levels = 2; 
  int    *update = NULL, *external = NULL;
  int    *update_index = NULL, *extern_index = NULL;
#ifdef ML_BENCHMARK
  struct reader_context *context;
#endif

int main(int argc, char *argv[])
{
  int    i, j, input_option, precon_flag, N_elements_coarse;
  double *b, *x;

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];

  /* data structure for matrix corresponding to the fine grid */

  int    *rpntr = NULL,*cpntr = NULL, *indx = NULL, *bpntr = NULL;
  int    proc_factor;

  /* NOTE: alpha and beta are only used in convection diffusion problems */
  /* (alpha is the strength of diffusion and beta is the strength of     */
  /* convection)   */

  /* data structures for multilevel grid, discrete operators */
  /* and preconditioners                                     */

  AZ_PRECOND  *Pmat = NULL;
  AZ_MATRIX   *Amat = NULL;
  ML          *ml   = NULL;
  ML_Aggregate *ml_ag = NULL;
#ifdef ML_BENCHMARK
  char input[MAX_INPUT_STR_LN];
  FILE *ifp;
#endif
#ifdef MATLAB
	FILE *xfile, *rhsfile;
#endif

  /* ----------------------- execution begins --------------------------------*/

#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif
  if (argc == -7) ML_avoid_unused_param((void *) argv);

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
	  pr_error("ERROR: can\'t interp int while looking for \"%s\"\n",
		  "number of x nodes per proc");
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
   agg_thresh = context->agg_thresh;


#else

  if (proc_config[AZ_node] == 0) 
  {
     printf("Preconditioning (0 for none, 1 for ML, 2 for ilu) : \n");
     scanf("%d", &precon_flag);
/*
     precon_flag = 1;
*/
     printf("n x n grid in each processor, n = ");
     scanf("%d", &proc_factor);
     refine_factor = 1;
     choice = 1;
     printf("Coarse solver (0 for SuperLU, >= 1 for GS) : \n");
     scanf("%d", &coarse_iterations);
     printf("Use CG (0 for no , 1 for yes) : \n");
     scanf("%d", &use_cg);
     printf("Number of levels : \n");
     scanf("%d", &num_levels);
     printf("number of smoothing steps : \n");
     scanf("%d", &nsmooth);
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
  AZ_broadcast((char *) NULL        ,   0          , proc_config, AZ_SEND);
#endif
  i = proc_config[AZ_N_procs];
  i = (int) pow( (double) i, 0.50001 );
  i = i * proc_factor;
  N_elements_coarse = i * i;
  if (proc_config[AZ_node] == 0) {
     i = i * refine_factor;
     printf("Fine grid = %d x %d \n", i+1 , i+1 ); 
  }
  input_option = 1;

#ifdef PARTEST
	global_nrows=(i+1)*(i+1)*num_PDE_eqns;
	global_ncoarse=global_nrows;
#endif

  /* set AZTEC options */

  init_options(options,params);

  /* create coarse grid discretization */
  construct_ml_grids(N_elements_coarse, proc_config,
                     &Amat, &update, &update_index, 
                    &external, &extern_index, &indx,
                     &bpntr, &rpntr, &cpntr, &ml, precon_flag, &ml_ag);
  /* set up ML preconditioning */

  if ( precon_flag == 1 ) {
     AZ_set_ML_preconditioner(&Pmat, Amat, ml, options);
  }
  else if ( precon_flag == 2 ) { 
     options[AZ_precond]         = AZ_dom_decomp;
     options[AZ_subdomain_solve] = AZ_ilut; 
  } else options[AZ_precond] = AZ_none;

  /*printf("\n\n----------------\nBEFORE SOLVE\n----------------\n\n");
  AZ_check_input(Amat->data_org,options,params,proc_config);
  printf("\n\n----------------\nBEFORE SOLVE\n----------------\n\n");*/


  /* solve the system of equations on the fine grid */

  for ( i = 0; i < 1; i++ ) 
  {
     init_guess_and_rhs(update_index, update, &x, &b, Amat->data_org, 
                        Amat->val, indx, Amat->bindx, rpntr, cpntr, 
                        bpntr, proc_config);

#ifdef MATLAB
	xfile=fopen("x0.txt","w");
	for (i=0 ; i<(proc_factor+1)*(proc_factor+1) ; i++) {
		fprintf(xfile,"%lf\n",x[i]);
	}
	fclose(xfile);
	
#endif
#ifdef SCALE_ME
	ML_MSR_scalerhs(b, scaling_vect, Amat->data_org[AZ_N_internal] +
                    Amat->data_org[AZ_N_border]);
	ML_MSR_scalesol(x, scaling_vect, Amat->data_org[AZ_N_internal] +
			Amat->data_org[AZ_N_border]);
#endif

      if (use_cg == 0) {
         ML_Iterate(ml, x, b); 
      }
      else {
         AZ_iterate(x, b, options, params, status, proc_config, Amat, Pmat, 
#ifdef AZ_ver2_1_0_5
                    (struct AZ_SCALING *) 
#else
                    (struct scale*) 
#endif
                                    NULL);
      }

#ifdef MATLAB
	rhsfile=fopen("rhs.txt","w");
	xfile=fopen("x.txt","w");
	if ((rhsfile==NULL) || (xfile==NULL))
		printf("couldn't open output files\n");
	printf("proc_factor=%d\n",proc_factor);
	for (i=0 ; i<(proc_factor+1)*(proc_factor+1) ; i++) {
		fprintf(rhsfile,"%lf\n",b[i]);
		fprintf(xfile,"%lf\n",x[i]);
	}
	fclose(rhsfile);
	fclose(xfile);
#endif


  }

  /*printf("\n\n----------------\nAFTER SOLVE\n----------------\n\n");
  AZ_check_input(Amat->data_org,options,params,proc_config);
  printf("\n\n----------------\nAFTER SOLVE\n----------------\n\n");*/
  
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



ML_Aggregate_Destroy(&ml_ag);
   ML_Destroy(&ml);
   AZ_free((void *) Amat->data_org);
   AZ_free((void *) Amat->val);
   AZ_free((void *) Amat->bindx);
   if (Amat  != NULL) AZ_matrix_destroy(&Amat);
   if (Pmat  != NULL) AZ_precond_destroy(&Pmat);


  /* Free allocated memory */

  free((void *) update);    AZ_free((void *) update_index);
  AZ_free((void *) external);  AZ_free((void *) extern_index);
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
/*
  options[AZ_overlap]  = AZ_none;
*/
  options[AZ_overlap]  = 1;
  options[AZ_kspace]   = 130;
  options[AZ_orthog]   = AZ_modified;
  options[AZ_aux_vec]  = AZ_resid;

  params[AZ_tol]       = 1.00e-8;
  params[AZ_drop]      = 0.0;
  params[AZ_ilut_fill] = 10;
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

  int    i, j;
  int    temp,num;
  AZ_MATRIX Amat;

  temp = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  num  = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* allocate vectors */

  i       = num + data_org[AZ_N_external];
  *x      = (double *) calloc(i, sizeof(double));
  *b      = (double *) calloc(i, sizeof(double));
  if ((*b == NULL) && (i != 0)) {
    pr_error("Not enough space in init_guess_and_rhs() for ax\n");
  }

  /*initialize 'x' to a function which will be used in matrix-vector product */

  for (i = 0; i < temp; i++) 
     (*x)[i] = (double) (update[i])/(double)(num_PDE_eqns);
  /* for (i = 0; i < temp; i++) (*x)[i] = (double) 1; */

  for (i = 0; i < temp/num_PDE_eqns; i++)
     for (j = 0; j < num_PDE_eqns; j++)
        (*x)[i*num_PDE_eqns+j] = (*x)[i*num_PDE_eqns];

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
              int *extern_index_f[], int *indx[], 
              int *bpntr[], int *rpntr[], int *cpntr[], ML **ml_ptr, 
              int flag, ML_Aggregate **ml_ag)
{
  int         MSRorVBR = AZ_MSR_MATRIX, N_update_f;
  int         i, j, N2, *bindx, *data_org, leng;
  double      *val;
  double *null_vect = NULL;

  ML           *ml;
  ML_GridAGX   *f_mesh;
  int          N_levels = 2, level;
  int          coarsest_level;
  int          old_guy, ndim;
#ifdef ML_BENCHMARK
  int nblocks, *blocks;
#endif

   /**************************************************************************/


   N_levels = num_levels;
   /* generate the fine grid and its partitionings */

   N2 = N_elements * refine_factor * refine_factor;
   Generate_mesh(N2, &f_mesh, &N_update_f, update_f, proc_config);

   /* use the partitioning for the fine grid to initialize matrix */

   N_grid_pts = AZ_gsum_int(N_update_f, proc_config);

   /* Change N_update_f and update_f to reflect num_pdes */

   N_update_f *= num_PDE_eqns;
   *update_f = (int *) realloc(*update_f, sizeof(int)*N_update_f);
   for (i = N_update_f/num_PDE_eqns-1; i >= 0; i--) {
      old_guy = (*update_f)[i];
      for (j = 0; j < num_PDE_eqns; j++) {
         (*update_f)[i*num_PDE_eqns+j] = old_guy*num_PDE_eqns + j;
      }
   }

/*
   create_msr_matrix_conv_diff(*update_f, &val, &bindx, N_update_f, 1.0, 1000.0); 
*/
   create_msr_matrix(*update_f, &val, &bindx, N_update_f);

   AZ_transform(proc_config, external_f, bindx, val, *update_f, update_index_f, 
                extern_index_f, &data_org, N_update_f, *indx, *bpntr, *rpntr, 
                cpntr, MSRorVBR);

#ifdef PARTEST

   global_mapping=(int **) malloc(global_nrows*sizeof(int *));
	 
   for (i=0; i < global_nrows; i++) {
      global_mapping[i]=(int *)malloc(2*sizeof(int));
      ml_print_global_element(i,*update_f, data_org, *update_index_f, 
                      global_mapping[i], proc_config, *external_f, *extern_index_f);
   }

#endif
   *Amat_f = AZ_matrix_create(data_org[AZ_N_internal] + data_org[AZ_N_border]);
   AZ_set_MSR(*Amat_f, bindx, val, data_org, 0, NULL, AZ_LOCAL);

   (*Amat_f)->matrix_type  = data_org[AZ_matrix_type];
   data_org[AZ_N_rows]  = data_org[AZ_N_internal] + data_org[AZ_N_border];
   ML_GridAGX_Destroy(&f_mesh); 
#ifdef SCALE_ME
   ML_MSR_sym_diagonal_scaling(*Amat_f, proc_config, &scaling_vect);
#endif



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

ML_Aggregate_Create(ml_ag);
ML_Aggregate_Set_Threshold( *ml_ag, agg_thresh );
ndim = num_PDE_eqns;
/*

ndim = 2;
null_vect=(double *) calloc(leng*ndim, sizeof(double));
c=cos(theta);
s=sin(theta);

for (i=0; i<leng/ndim; i++) {
null_vect[ i*ndim]=c;
null_vect[ i*ndim+1 ]=-s;
null_vect[ i*ndim+ leng ]=s;
null_vect[ i*ndim+ leng + 1 ]=c;
null_vect[ i*ndim]=1.;
null_vect[ i*ndim+1 ]=1.;
null_vect[ i*ndim+ leng ]=1.;
null_vect[ i*ndim+ leng + 1 ]=-1.;
}
*/

      ML_Aggregate_Set_NullSpace(*ml_ag, num_PDE_eqns, ndim, null_vect, leng); 
      if (null_vect != NULL) free(null_vect);
#ifdef SCALE_ME
      ML_Aggregate_Scale_NullSpace(*ml_ag, scaling_vect,leng);
#endif

      coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, N_levels-1, 
                                                          ML_DECREASING, *ml_ag);
      coarsest_level = N_levels - coarsest_level;
      if ( proc_config[AZ_node] == 0 )
         printf("Coarse level = %d \n", coarsest_level);

      /* set up smoothers */

      for (level = N_levels-1; level > coarsest_level; level--) {
         /*  Sparse approximate inverse smoother that acutally does both */
         /*  pre and post smoothing.                                     */

         /*

         ML_Gen_Smoother_ParaSails(ml , level, ML_PRESMOOTHER, nsmooth,
                                parasails_sym, parasails_thresh,
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
         */

         /* This is the symmetric Gauss-Seidel smoothing. In parallel,    */
         /* it is not a true Gauss-Seidel in that each processor          */
         /* does a Gauss-Seidel on its local submatrix independent of the */
         /* other processors.                                             */
#define MB_MODIF
#define MB_MODIF2
#ifndef ML_BENCHMARK
         ML_Gen_Smoother_SymGaussSeidel(ml, level, ML_BOTH,  nsmooth,1.);
	 /* ml->pre_smoother[level].smoother->internal = 
	                                  ML_Smoother_MSR_GSforwardnodamping;
	    ml->post_smoother[level].smoother->internal = 
	                                  ML_Smoother_MSR_GSbackwardnodamping;

	    temp1[0] = ml->Amat[level].outvec_leng;
	    if (level != 0)
	       temp1[1] = ml->Amat[level-1].outvec_leng;
	    else temp1[1] = 0;
	    ML_gsum_vec_int(temp1, temp2, 2, ml->comm);
	    eig_ratio = 20.;
	    if (temp1[1] != 0)
	       eig_ratio = ((double) temp1[0])/ ((double) temp1[1]);
	    if (eig_ratio < 4.) eig_ratio = 4.;
	    ML_Gen_Smoother_MLS(ml, level, ML_BOTH, eig_ratio,nsmooth);
	 */
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
	   ML_Gen_Blocks_Aggregates(*ml_ag, level, &nblocks, &blocks);
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
	   pr_error("unknown smoother %s\n",context->smoother);
	 }
#endif
      }
#ifndef ML_BENCHMARK
      if (coarse_iterations == 0) ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
      else ML_Gen_Smoother_SymGaussSeidel(ml, coarsest_level, ML_PRESMOOTHER, 
					  coarse_iterations,1.);
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
	ML_Gen_Blocks_Aggregates(*ml_ag, coarsest_level, &nblocks, &blocks);
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
#ifdef HAVE_ML_SUPERLU
      else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0) {
	ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
      }
#elif defined(HAVE_ML_AMESOS)
      else if (ML_strcmp(context->coarse_solve,"KLU") == 0) {
	ML_Gen_Smoother_Amesos( ml, coarsest_level, ML_AMESOS_KLU,-1);
      }
#endif
      else {
	pr_error("unknown coarse grid solver %s\n",context->coarse_solve);
      }
#endif


      ML_Gen_Solver(ml, ML_MGV, N_levels-1, coarsest_level); 

   }
   return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Generate_mesh(int N_elements, ML_GridAGX **meshp, int *N_update, 
                   int *update[], int proc_config[])
{
   int    i, j, k, nprocs, nprocs_1d, mypid, mypid_x, mypid_y;
   int    icnt, *intlist, offset, nelmnt_1d, nelmnt_part_xy, estart_x; 
   int    estart_y, nelmnt_local, nstart_x, nstart_y, index, eend_x, eend_y;
   int    nnode_1d, nnode_local_x, nnode_local_y, nnode_local;
   int    nnode_expand, nnode_x, nnode_y, **square;
   int    x_begin, y_begin;
   double *coord, h;
   ML_GridAGX *mesh;

   /* --------------------------------------------------------------- */
   /* box in the processors                                           */
   /* --------------------------------------------------------------- */

   mypid    = proc_config[AZ_node];
   nprocs   = proc_config[AZ_N_procs];
   
   nprocs_1d = (int) pow( (double) nprocs, 0.50001 );
   if ( nprocs_1d * nprocs_1d != nprocs ) {
      pr_error("Error: nprocs should be a square (%d).\n",nprocs_1d);
   }
   mypid_x = mypid % nprocs_1d;
   mypid_y = mypid / nprocs_1d;

   /* --------------------------------------------------------------- */
   /* process the grid                                                */
   /* --------------------------------------------------------------- */

   nelmnt_1d = (int) pow( (double) N_elements, 0.500001 );

   /* partition the elements along the xy-axis */

   nelmnt_part_xy = nelmnt_1d / nprocs_1d;
   if (nelmnt_part_xy * nprocs_1d != nelmnt_1d) {
      pr_error("Error: nelmnt_part not good. %d %d\n",nelmnt_part_xy,nelmnt_1d);
   }
   nelmnt_local = nelmnt_part_xy * nelmnt_part_xy;
  
   /* initialize the mesh */

   ML_GridAGX_Create( meshp );
   mesh = (*meshp);
   h = 1.0 / (double) nelmnt_1d;
   ML_GridAGX_Set_Dimension(mesh, 2);
   ML_GridAGX_Set_NElmnts(mesh, nelmnt_local, 4);

   /* generate the global element numbers of the local elements */

   ML_memory_alloc((void**) &intlist, nelmnt_local * sizeof(int), "AP1");
   estart_x  = nelmnt_part_xy * mypid_x;
   estart_y  = nelmnt_part_xy * mypid_y;
   eend_x    = estart_x + nelmnt_part_xy - 1;
   eend_y    = estart_y + nelmnt_part_xy - 1;
   icnt      = 0;

   for (k = estart_y; k <= eend_y; k++) 
      for (j = estart_x; j <= eend_x; j++) 
         intlist[icnt++] = j + k * nelmnt_1d;
   ML_GridAGX_Load_ElmntGlobalNum(mesh, nelmnt_local, intlist);
   ML_memory_free( (void **) &intlist);

   /* generate the local computational grid */

   nnode_x       = nelmnt_part_xy + 1;
   nnode_y       = nelmnt_part_xy + 1;
   ML_memory_alloc((void**) &square, nnode_x * sizeof(int*), "AP2");
   for ( i = 0; i < nnode_x; i++ )
      ML_memory_alloc((void**) &(square[i]), nnode_y * sizeof(int), "AP3");
   icnt = 0;
   x_begin = 1;
   y_begin = 1;
   if ( estart_x == 0 ) x_begin--;
   if ( estart_y == 0 ) y_begin--;
   for ( k = y_begin; k < nnode_y; k++ )
      for ( j = x_begin; j < nnode_x; j++ )
         square[j][k] = icnt++;
   if ( x_begin == 1 || y_begin == 1 ) square[0][0] = icnt++;
   if ( x_begin == 1 ) {
      for ( k = 1; k < nnode_y; k++ ) square[0][k] = icnt++;
   }
   if ( y_begin == 1 ) {
      for ( j = 1; j < nnode_x; j++ ) square[j][0] = icnt++;
   }
   
   /* generate the global node IDs */

   nnode_local_x = nnode_x - 1;
   nnode_local_y = nnode_y - 1;
   if (mypid_x == 0 ) nnode_local_x++;
   if (mypid_y == 0 ) nnode_local_y++;
   nnode_1d = nelmnt_1d + 1;
   nnode_local  = nnode_local_x * nnode_local_y;
   (*N_update)  = nnode_local;
   nnode_expand = nnode_y * nnode_x;
   (*update)   = (int*) malloc(nnode_expand * sizeof(int));

   nstart_x    = mypid_x * nnode_local_x;
   nstart_y    = mypid_y * nnode_local_y;
   for (k = 0; k < nnode_y; k++) 
      for (j = 0; j < nnode_x; j++) {
            index = square[j][k];
            offset  = j + nstart_x + (k+nstart_y)*nnode_1d;
            (*update)[index] = offset;
      }

   ML_GridAGX_Set_NVert(mesh, *N_update);
   ML_GridAGX_Load_VertGlobalNum(mesh, nnode_expand, *update);

   /* generate coordinates of local nodes */

   ML_memory_alloc((void**) &coord, 2*nnode_expand*sizeof(double),"AP6");
   for (k = 0; k < nnode_y; k++) 
      for (j = 0; j < nnode_x; j++) {
         index = square[j][k];
         offset  = j + nstart_x + (k+nstart_y)*nnode_1d;
         coord[index*2+0] = (double) (h * (j + nstart_x));
         coord[index*2+1] = (double) (h * (k + nstart_y));
      } 
   ML_GridAGX_Load_AllVertCoordinates(mesh, nnode_expand, coord);

   /* finally construct the element-to-node lists */

   ML_memory_alloc((void**) &intlist, 4 * sizeof(int), "AP7");

   k = (eend_y - estart_y + 1) * (eend_x - estart_x + 1);
   for (k = estart_y; k <= eend_y; k++) 
      for (j = estart_x; j <= eend_x; j++) {
            intlist[0] = square[j-estart_x][k-estart_y];
            intlist[1] = square[j-estart_x+1][k-estart_y];
            intlist[2] = square[j-estart_x][k-estart_y+1];
            intlist[3] = square[j-estart_x+1][k-estart_y+1];
            ML_IntList_Load_Sublist(mesh->ele_nodes, 4, intlist);
      }  
   ML_memory_free( (void **) &intlist);
   for ( i = 0; i < nnode_x; i++ )
      ML_memory_free( (void**) &(square[i]) );
   ML_memory_free( (void**) &square);
   ML_memory_free( (void **) &coord);
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
  if (count == -57) ML_avoid_unused_param((void *) comm);
  return(md_wrap_iread( buf, count, src, mid, (MPI_Request *) request));
}

int ml_wait(void* buf, unsigned int count, int *src, int *mid, 
            USR_COMM comm, USR_REQ *request )
{
   int status;
  if (count == -57) ML_avoid_unused_param((void *) comm);
   return(md_wrap_wait( buf, count, src, mid, &status, (MPI_Request *) request));
}

int ml_send(void* buf, unsigned int count, int dest, int mid, 
            USR_COMM comm )
{
   int status;
   if (count == -57) ML_avoid_unused_param((void *) comm);
   return(md_wrap_write( buf, count, dest, mid, &status));
}

/*****************************************************************************
 *                             MSR                                           *
 *****************************************************************************/

void create_msr_matrix(int update[], double **val, int **bindx, int N_update)

{
  int i, total_nz, n = -1;
  int avg_nonzeros_per_row = 9;

  total_nz = /*num_PDE_eqns**/N_update*avg_nonzeros_per_row + 1;
  *bindx   = (int *) AZ_allocate(total_nz*sizeof(int));
  *val     = (double *) AZ_allocate(total_nz*sizeof(double));
  if ((*val == NULL) && (total_nz != 0) ) {
    (void) fprintf(stderr, "Error: Not enough space to create matrix\n");
    (void) fprintf(stderr,
                   "      Try reducing the variable 'avg_nonzeros_per_row'\n");
    exit(1);
  }
  for (i = 0 ; i < total_nz ; i++ ) (*bindx)[i] = 0;

  (*bindx)[0] = /*num_PDE_eqns**/N_update+1;
  for (i = 0; i < N_update/**num_PDE_eqns*/; i++) {
		add_row_5pt(update[i], i, *val, *bindx,&n);
		if ( (*bindx)[i+1] > total_nz) {
			(void) fprintf(stderr, "Error:total_nz not large enough to accomodate");
			(void) fprintf(stderr, " nonzeros\n       Try increasing the variable");
			(void) fprintf(stderr, " 'avg_nonzeros_per_row'\n       Finished \n");
			(void) fprintf(stderr, "first %d rows out of %d\n", i, N_update);
			exit(1);
		}
	}

} /* create_msr_matrix */

void create_rot_msr_matrix(int update[], double **val, int **bindx, int N_update)

{
  int i, total_nz, n = -1;
  int avg_nonzeros_per_row = 10;

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
		add_row_rot(update[i], i, *val, *bindx, &n);
		if ( (*bindx)[i+1] > total_nz) {
			(void) fprintf(stderr, "Error:total_nz not large enough to accomodate");
			(void) fprintf(stderr, " nonzeros\n       Try increasing the variable");
			(void) fprintf(stderr, " 'avg_nonzeros_per_row'\n       Finished \n");
			(void) fprintf(stderr, "first %d rows out of %d\n", i, N_update);
			exit(1);
		}
	}

} /* create_rot_msr_matrix */


void create_msr_matrix_conv_diff(int update[], double **val, int **bindx, 
                                 int N_update, double alpha, double beta)
{
  int i,total_nz, n = -1;
  int avg_nonzeros_per_row = 7;

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
    add_row_7pt(update[i], i, *val, *bindx,&n, alpha, beta);
    if ( (*bindx)[i+1] > total_nz) {
      (void) fprintf(stderr, "Error:total_nz not large enough to accomodate");
      (void) fprintf(stderr, " nonzeros\n       Try increasing the variable");
      (void) fprintf(stderr, " 'avg_nonzeros_per_row'\n       Finished \n");
      (void) fprintf(stderr, "first %d rows out of %d\n", i, N_update);
      exit(1);
    }
  }

} /* create_msr_matrix_conv_diff */

/******************************************************************************/
/******************************************************************************/
double *biggie = NULL;
/* @@@@ */
double ml_variable_coef(double x, double y, int m)
{
  double mpi = 3.14159;
  static int seed = 493217272;
  int i;

  return(1.);
  if (biggie == NULL) {
    biggie = (double *) malloc(sizeof(double)*2*m*2*m + 1000);
    for (i = 0; i < 4*m*m + 1; i++) biggie[i] = pow(10.,AZ_srandom1(&seed));
  }
  return( biggie[   (int) (2.*x*m + 2.*y*m + .00001) ] );
  printf("%e\n",2.*x*m);
  mpi *= (double) m;
  /*
  printf("a few %e %e   %e %e\n",
	 mpi*x,mpi*y,sin(mpi*x),sin(mpi*y));
  */
  return( sin(mpi*x) * sin(mpi*y) + 1.);
  return( pow(sin(mpi*x),2.) + pow(sin(mpi*y),2.) );
  if ((x < .315) && (x > .2) ) return(1.0e-5);
  if ((x < .515) && (x > .4) ) return(1.0e+5);
  if ((x < .315) && (x > .2) ) return(1.0e-5);
}
void add_row_5pt(int row, int location, double val[], 
								 int bindx[], int *n)

/*
 * Add one row to an MSR matrix corresponding to a 5pt discrete approximation
 * to the 2D Poisson operator on an n x n square.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/15/95
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    val,bindx    == (see user's guide). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{

/*
  static int n = -1;
*/
  int m;
  int        k;
  int        NP;
  double x,y,h;

  /* determine grid dimensions */

  m = *n;
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) sqrt( ((double) m)+0.01);
    if ( m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^2 where k is an integer\n");
      exit(1);
    }
    if (m == 1) {

      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }
  }
  *n = m;

  k  = bindx[location];
  NP = num_PDE_eqns;

  x = (double )( (row/NP)%m);
  y = (double )( (row/(NP*m)));
  h = 1./((double) m);

  x = h/2. + x/((double) m);
  y = h/2. + y/((double) m);
  /*
   * Check neighboring points in each direction and add nonzero entry if
   * neighbor exists.
   */

#ifndef MB_MODIF
  bindx[k] = row + NP;   if ((row/NP)%m !=     m-1) val[k++] = -1.00;
  bindx[k] = row - NP;   if ((row/NP)%m !=       0) val[k++] = -1.00;
  bindx[k] = row + m*NP; if ((row/(NP*m))%m != m-1) val[k++] = -1.00;
  bindx[k] = row - m*NP; if ((row/(NP*m))%m !=   0) val[k++] = -1.00;
  val[location]     = 4.0;
#else

  /* hardwired stencil stuff ...

  bindx[k] = row + NP;   if ((row/NP)%m !=     m-1) val[k++] = sten_f;
  bindx[k] = row - NP;   if ((row/NP)%m !=       0) val[k++] = sten_d;
  bindx[k] = row + m*NP-NP; if (((row/(NP*m))%m != m-1) &&((row/NP)%m != 0))  val[k++] = sten_a;
  bindx[k] = row + m*NP; if ((row/(NP*m))%m != m-1) val[k++] = sten_b;
  bindx[k] = row + m*NP+NP; if (((row/(NP*m))%m != m-1) &&((row/NP)%m != m-1)) val[k++] = sten_c;
  bindx[k] = row - m*NP-NP; if (((row/(NP*m))%m !=   0) &&((row/NP)%m != 0)) val[k++] = sten_g;
  bindx[k] = row - m*NP; if ((row/(NP*m))%m !=   0) val[k++] = sten_h;
  bindx[k] = row - m*NP+NP; if (((row/(NP*m))%m !=   0) && ((row/NP)%m != m-1)) val[k++] = sten_i;
  val[location]     = sten_e;
  */

  bindx[k] = row + NP;   if ((row/NP)%m !=     m-1) val[k++] = -ml_variable_coef(x+h/2.,y,m);
  bindx[k] = row - NP;   if ((row/NP)%m !=       0) val[k++] = -ml_variable_coef(x-h/2.,y,m);
  bindx[k] = row + m*NP; if ((row/(NP*m))%m != m-1) val[k++] = -ml_variable_coef(x,y+h/2.,m);
  bindx[k] = row - m*NP; if ((row/(NP*m))%m !=   0) val[k++] = -ml_variable_coef(x,y-h/2.,m);
  val[location]     = ml_variable_coef(x+h/2.,y,m) + ml_variable_coef(x-h/2.,y,m)+
                      ml_variable_coef(x,y+h/2.,m) 
                       + ml_variable_coef(x,y-h/2.,m);
#endif

  bindx[location+1] = k;

} /* add_row_5pt */


void add_row_rot(int row, int location, double val[], 
								 int bindx[], int *n)

/*
 * same as add_row_5pt except a rotation matrix is added and NP is assumed to be 2
 */

{

/*
  static int n = -1;
*/
  int m;
  int        k;
  int        NP, oddeven;
	double c, s;
double new_theta = 0.170796;

	c=cos(new_theta);
	s=sin(new_theta);

  /* determine grid dimensions */

  m = *n;
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) sqrt( ((double) m)+0.01);
    if ( m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^2 where k is an integer\n");
      exit(1);
    }
    if (m == 1) {

      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }
  }
  *n = m;

  k  = bindx[location];
  NP = num_PDE_eqns;

	if (NP != 2)
		pr_error("error in setting A matrix - the rotated version is only implemented for 2 PDEs\n");
	
	oddeven=row%2;

  /*
   * Check neighboring points in each direction and add nonzero entry if
   * neighbor exists.
   */

  bindx[k] = row + NP;   if ((row/NP)%m !=     m-1) val[k++] = -c;
  bindx[k] = row - NP;   if ((row/NP)%m !=       0) val[k++] = -c;
  bindx[k] = row + m*NP; if ((row/(NP*m))%m != m-1) val[k++] = -c;
  bindx[k] = row - m*NP; if ((row/(NP*m))%m !=   0) val[k++] = -c;

	if (oddeven==1) {
		bindx[k] = row + 1;   if ((row/NP)%m !=     m-1) val[k++] =      s;
		bindx[k] = row - 1;                              val[k++] =   -4*s;
		bindx[k] = row - 3;   if ((row/NP)%m !=       0) val[k++] =      s;
		bindx[k] = row + m*NP - 1; if ((row/(NP*m))%m != m-1) val[k++] = s;
		bindx[k] = row - m*NP - 1; if ((row/(NP*m))%m !=   0) val[k++] = s;
	}
	else {
		bindx[k] = row + 1;                              val[k++] =    4*s;
		bindx[k] = row - 1;   if ((row/NP)%m !=       0) val[k++] =     -s;
		bindx[k] = row + 3; if ((row/(NP))%m != m-1)   val[k++] =     -s;
		bindx[k] = row + m*NP+1; if ((row/(NP*m))%m != m-1) val[k++] =  -s;
		bindx[k] = row - m*NP+1; if ((row/(NP*m))%m !=   0) val[k++] =  -s;
	}


  bindx[location+1] = k;
  val[location]     = 4.0*c;

} /* add_row_rot */


void add_row_7pt(int row, int location, double val[], int bindx[], int *n, 
								 double alpha, double beta)

/*
 * Add one row to an MSR matrix corresponding to a 7pt discrete approximation
 * to the 2D Convection Diffusion operator with diffusion coeffcient alpha and 
 * convection coefficient beta on an n x n square.
 *
 * Author: Dawn Chamberlain
 * Date:   7/22/99
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    val,bindx    == (see user's guide). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{

/*
  static int n = -1;
*/
  int m;
  int        k;
  int        NP;
	double h;

  /* determine grid dimensions */
	/* it appears to me (dmc) that m is supposed to be the number of grid points per
		 row, but it doesn't appear to be commented anywhere */

  m = *n;
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) sqrt( ((double) m)+0.01);
    if ( m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^2 where k is an integer\n");
      exit(1);
    }
    if (m == 1) {
			
      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }
  }
  *n = m;
	
	h=1.0/(m-1);  /* grid spacing */

  k  = bindx[location];
  NP = num_PDE_eqns;
	
  /*
   * Check neighboring points in each direction and add nonzero entry if
   * neighbor exists.
   */
		
  bindx[k] = row + NP;   
	if ((row/NP)%m !=     m-1) 
		val[k++] = -alpha - beta*h/6.0;
	
  bindx[k] = row - NP;   
	if ((row/NP)%m !=       0) 
		val[k++] = -alpha + beta*h/6.0;
	
  bindx[k] = row + m*NP; 
	if ((row/(NP*m))%m != m-1) 
		val[k++] = -alpha - beta*h/6.0;
	
  bindx[k] = row - m*NP; 
	if ((row/(NP*m))%m !=   0) 
		val[k++] = -alpha + beta*h/6.0;

  bindx[k] = row - (m+1)*NP; 
	if ((row/(NP*m))%m >   1) 
		val[k++] = beta*h/3.0;

  bindx[k] = row + (m+1)*NP; 
	if ((row/(NP*m))%m <   m-2) 
		val[k++] = -beta*h/3.0;

  bindx[location+1] = k;
  val[location]     = 4.0;

} /* add_row_7pt */


void AZ_print_out2(int update_index[], int extern_index[], int update[], 
	int external[], double val[], int indx[], int bindx[], int rpntr[], 
	int cpntr[], int bpntr[], int proc_config[], int choice, int matrix, 
	int N_update, int N_external, int of)
{
/*******************************************************************************

  Print out the matrix in 1 of 3 formats described below.
  starting point of the matrix.

  Author:          Ray Tuminaro, SNL, 9222
  =======

  Return code:     none.
  ============ 

  Parameter list:
  ===============

  update_index,    On input, ordering of update and external locally on this
  extern_index:    processor. For example  'update_index[i]' gives the index
                   location of the block which has the global index 'update[i]'.
                   (AZ_global_mat only).
 
  update:          On input, blks updated on this node (AZ_global_mat only).

  external:        On input, list of external blocks (AZ_global_mat only).

  val,bindx        On input, matrix (MSR or VBR) arrays holding matrix values.
  indx, bnptr,     When using either AZ_input_form or AZ_explicit, these can
  rnptr, cnptr:    be either pre or post-AZ_transform() values depending on what
                   the user wants to see. When using AZ_global_form, these must
                   be post-AZ_transform() values.

  proc_config:     On input, processor information corresponding to:
                      proc_config[AZ_node] = name of this processor
                      proc_config[AZ_N_procs] = # of processors used
 
  choice:          On input, 'choice' determines the output to be printed
		   as described above.
 
  matrix:          On input, type of matrix (AZ_MSR_MATRIX or AZ_VBR_MATRIX).

  N_update:        On input, number of points/blks to be updated on this node.

  N_external:      On input, number of external points/blks needed by this node.

  of:              On input, an offset used with AZ_global_matrix and 
		   AZ_explicit. In particular, a(of,of) is printed for 
                   the matrix element stored as a(0,0).

*******************************************************************************/

   int type, neighbor, cflag;
   int ii,i = 1,j,tj;
   int iblk, jblk, m1, n1, ival, new_iblk, new_jblk;
   MPI_Request request;  /* Message handle */
#ifdef MATLAB
	 FILE *outfp;
#endif

 
	 /*   type            = AZ_sys_msg_type;
   AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % 
                     AZ_NUM_MSGS +AZ_MSG_TYPE;*/
   /* Synchronize things so that processor 0 prints first and then */
   /* sends a message to processor 1 so that he prints, etc.      */

#ifdef MATLAB
	 outfp=fopen("setAmat.m","w");
	 fprintf(outfp,"A=sparse(%d,%d);\n",N_grid_pts, N_grid_pts);
#endif

   neighbor = proc_config[AZ_node] - 1;
   if ( proc_config[AZ_node] != 0) {
      mdwrap_iread((void *) &i, 0, &neighbor, &type, &request);
      mdwrap_wait((void *) &i, 0, &neighbor, &type, &cflag, &request);
   }
   printf("proc %d:\n",proc_config[AZ_node]);

   if (choice == AZ_input_form ) {
     if ( update != (int *) NULL) {
        printf("  N_update: %5d\n",N_update); printf("  update:");
        AZ_list_print(update, N_update, (double *) NULL , 0);
     }

     if (matrix == AZ_MSR_MATRIX) {
        printf("  bindx: ");
        AZ_list_print(bindx, bindx[N_update], (double *) NULL , 0);

        printf("  val:   ");
        AZ_list_print((int *) NULL , N_update, val , bindx[N_update]);
     }
     else if (matrix == AZ_VBR_MATRIX) {
        printf("  rpntr: ");
        AZ_list_print(rpntr, N_update+1, (double *) NULL , 0);
        if ( cpntr != (int *) NULL ) {
           printf("  cpntr: ");
           AZ_list_print(cpntr, N_update+1+ N_external, (double *) NULL , 0);
        }
        printf("  bpntr: ");
        AZ_list_print(bpntr, N_update+1, (double *) NULL , 0);
        printf("  bindx: ");
        AZ_list_print(bindx, bpntr[N_update], (double *) NULL , 0);
        printf("  indx:  ");
        AZ_list_print(indx, bpntr[N_update]+1, (double *) NULL , 0);
        printf("  val:   ");
        AZ_list_print((int *) NULL, indx[bpntr[N_update]], val, 0);
     }
   }
   else if (choice == AZ_global_mat ) {
     if ( matrix == AZ_MSR_MATRIX) {
        for (i = 0 ; i < N_update; i++ ) {
          ii = update_index[i];
          printf("   a(%d,%d) = %20.13e;\n",update[i]+of,update[i]+of,val[ii]);
          for (j = bindx[ii] ; j < bindx[ii+1] ; j++ ) {
            tj = AZ_find_simple(bindx[j], update_index, N_update, extern_index,
                              N_external,update,external);
            if (tj != -1) 
               printf("   a(%d,%d) = %20.13e;\n",update[i]+of,tj+of,val[j]);
            else (void) fprintf(stderr,"col %d (= bindx[%d]) is undefined\n",
                                tj, j);
          }
        }
     }
     else if (matrix == AZ_VBR_MATRIX) {
        for (iblk= 0; iblk < N_update; iblk++) {
           new_iblk = update_index[iblk];

           m1 = rpntr[new_iblk+1] - rpntr[new_iblk];
 
           /* loop over blocks in the current block-row */
 
           for (ii = bpntr[new_iblk]; ii < bpntr[new_iblk+1]; ii++) {
              new_jblk = AZ_find_simple(bindx[ii], update_index, N_update, 
			 extern_index, N_external,update,external);
              if (new_jblk == -1) {
                 pr_error("local column %d not found\n",new_jblk);
              }
              jblk = bindx[ii];
              ival =  indx[ii];

              n1 = cpntr[jblk+1] - cpntr[jblk];

              for (i = 0; i < m1; i++) {
                 for (j = 0; j < n1; j++)
                    (void) printf("   a(%d(%d),%d(%d)) = %20.13e;\n",
                                  update[iblk]+ of,i+of, new_jblk+of, j+of, 
                                  val[ival+j*m1+i]);
              }
           }
        }
     }
   }
   else if (choice == AZ_explicit) {
     if ( matrix == AZ_MSR_MATRIX) {
        for (i = 0 ; i < N_update; i++ ) {
          if (update == NULL) tj = i+of;
          else tj = update[i] + of;
#ifdef MATLAB
          fprintf(outfp,"A(%d,%d)=%20.13e;\n",tj+1,tj+1,val[i]);
          for (j = bindx[i] ; j < bindx[i+1] ; j++ ) {
               fprintf(outfp,"A(%d,%d) = %20.13e;\n",tj+1,bindx[j]+of+1,val[j]);
          }
#endif
        }
     }
     else if (matrix == AZ_VBR_MATRIX) {
        for (iblk = 0; iblk < N_update; iblk++) {
           if (update == NULL) tj = iblk+of;
           else tj = update[iblk] + of;

           m1 = rpntr[iblk+1] - rpntr[iblk];
 
           /* loop over blocks in the current block-row */
 
           for (ii = bpntr[iblk]; ii < bpntr[iblk+1]; ii++) {
              jblk = bindx[ii];
              ival =  indx[ii];
              n1 = (indx[ii+1]-ival)/m1;

              for (i = 0; i < m1; i++) {
                 for (j = 0; j < n1; j++)
                    (void) printf("   a(%d(%d),%d(%d)) = %20.13e;\n", tj, 
				  i+of, jblk+of, j+of, val[ival+j*m1+i]);
              }
 
           }
        }
     }
   }
   else (void) fprintf(stderr,"AZ_matrix_out: output choice unknown\n");

   neighbor = proc_config[AZ_node] + 1;
   if ( proc_config[AZ_node] != proc_config[AZ_N_procs] - 1) 
      mdwrap_write((char *) &i, 0, neighbor, type, &cflag);

   i = AZ_gsum_int(i,proc_config);
#ifdef MATLAB
	 fclose(outfp);
#endif

}



void ml_print_global_element(int element,int update[],int data_org[],
           int update_index[], int vector[],int proc_config[],
	   int external[], int extern_index[])
{
/*
 * Print out the vector element corresponding to the global numbering
 * 'element'. Note: if the VBR format is used, this routine will print
 * out all the vector elements corresponding to this block.
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date:   6/15/96
 *
 * Parameters
 *
 *    element        ==      On input, global number of vector element that
 *                           will be printed.
 *    update         ==      On input, list of pts updated on this node
 *    data_org       ==      On input, indicates how the data is set out on
 *                           this node. For example, data_org[] contains
 *                           information on how many unknowns are internal,
 *                           external, and border unknowns as well as which
 *                           points need to be communicated. See User's Guide
 *                           for more details.
 *    update_index   ==      On input, ordering of update locally on this
 *                           processor. For example, 'update_index[i]' gives 
 *                           the index location of the block which has the 
 *                           global index 'update[i]'.
 *    rpntr          ==      On input, rpntr[i+1]-rpntr[i] gives the block
 *                           size of the ith local block.
 *    vector         ==      On input, vector to be printed (just one element).
 *    proc_config    ==      On input, processor information:
 *                              proc_config[AZ_node] = name of this processor
 *                              proc_config[AZ_N_procs] = # of processors used
 */
   int i;
   int N_local;

   /* synchronize things */

   i = AZ_gsum_int(1,proc_config);
   N_local = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
 


   i = AZ_find_index(element,update, N_local);

   if (i !=-1) {
      i = update_index[i];
			vector[0]=element;
			vector[1]=i;
   }
	 else {
            i = AZ_find_index(element,external, data_org[AZ_N_external]);
            if (i !=-1) {
               i = extern_index[i];
               vector[0]=element;
               vector[1]=i;
            }
	    else {
		 vector[0]=-1;
		 vector[1]=-1;
	    }
}

   /* synchronize things */
   i = AZ_gsum_int(i,proc_config);

}

int ml_find_global_row(int i, int proc, int global_nrows, int **whole_glob_map)
{
	/* this function takes a local row number and returns the global row number */
	int row, found_it;

	row = 0;
	found_it = 0;
	while ((row < global_nrows*num_PDE_eqns) && (!found_it)) {
		if ((whole_glob_map[row][1] == proc) && (whole_glob_map[row][2] == i)) {
			found_it = 1;
		}
		else
			row ++;
	}
	if (row==global_nrows*num_PDE_eqns)
		printf("error: didn't find local row (proc=%d, local row=%d) in global mapping\n", proc, i);
	
	return(row);
}

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
