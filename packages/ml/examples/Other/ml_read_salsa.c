/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/*****************************************************************************/
/* Sample driver for AZTEC/ML package. The software is tested by reading in  */
/* a matrix stored in a file called .data, using a zero initial guess        */
/* and a random right hand side, and then solving the system of equations    */
/* using AZTECs gmres solver and ML preconditioner                           */
/*                                                                           */
/* NOTE: the file .data must exist on all processors (though it need only    */
/* contain the number of rows in the matrix on all but the first processor)  */
/*                                                                           */
/* Author:       Dawn Chamberlain, Div 9222, Sandia National Labs            */
/* date:         10/21/99                                                    */
/*****************************************************************************/

#include "ml_include.h"
#if defined(HAVE_ML_AZTEC2_1) || defined(HAVE_ML_AZTECOO)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "az_aztec.h"
#include <math.h>

#ifndef GCC_VERSION
#if (defined(__GNUC__) && defined(__GNUC_MINOR__)) && defined(__GNUC_PATCHLEVEL__)
#define GCC_VERSION  (__GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__)
#endif
#endif

extern int AZ_using_fortran;
int    parasails_factorized = 0;
int    parasails_sym        = 0;
double parasails_thresh     = 0.15;
int    parasails_nlevels    = 0;
double parasails_filter     = 0.;
double parasails_loadbal    = 0.;



int main(int argc, char *argv[])
{
#if defined(GCC_VERSION) && GCC_VERSION >= 460
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif
	int num_PDE_eqns=5, N_levels=3;
    /* int nsmooth=1; */

	int    leng, level, N_grid_pts, coarsest_level;

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];

  /* data structure for matrix corresponding to the fine grid */

  int    *data_org = NULL, *update = NULL, *external = NULL;
  int    *update_index = NULL, *extern_index = NULL;
  int    *cpntr = NULL;
  int    *bindx = NULL, N_update, iii;
  double *val = NULL;
	double *xxx, *rhs;

	AZ_MATRIX *Amat;
	AZ_PRECOND *Pmat = NULL;
	ML *ml;
	FILE *fp;
  int ch,i;
   struct AZ_SCALING *scaling;
double solve_time, setup_time, start_time;
ML_Aggregate *ag;
int *ivec;
#ifdef VBR_VERSION
ML_Operator *B, *C, *D;
int *vbr_cnptr, *vbr_rnptr, *vbr_indx, *vbr_bindx, *vbr_bnptr, total_blk_rows;
int total_blk_cols, blk_space, nz_space;
double *vbr_val;
struct ML_CSR_MSRdata *csr_data;
#endif

#if defined(GCC_VERSION) && GCC_VERSION >= 460
#pragma GCC diagnostic pop
#endif


#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

#ifdef binary
	fp=fopen(".data","rb");
#else
	fp=fopen(".data","r");
#endif
	if (fp==NULL)
		{
			printf("couldn't open file .data\n");
			exit(1);
		}
#ifdef binary
        fread(&leng, sizeof(int), 1, fp);
#else
	fscanf(fp,"%d",&leng);
#endif

	fclose(fp);

	N_grid_pts=leng/num_PDE_eqns;



  /* initialize the list of global indices. NOTE: the list of global */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly. */

  AZ_read_update(&N_update, &update, proc_config, N_grid_pts, num_PDE_eqns,
                 AZ_linear);

  AZ_read_msr_matrix(update, &val, &bindx, N_update, proc_config);

  /* This code is to fix things up so that we are sure we have */
  /* all block (including the ghost nodes the same size.       */

  AZ_block_MSR(&bindx, &val, N_update, num_PDE_eqns, update);


  AZ_transform(proc_config, &external, bindx, val,  update, &update_index,
	       &extern_index, &data_org, N_update, 0, 0, 0, &cpntr,
               AZ_MSR_MATRIX);

  Amat = AZ_matrix_create( leng );

#ifndef VBR_VERSION

  AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);

  Amat->matrix_type  = data_org[AZ_matrix_type];

  data_org[AZ_N_rows]  = data_org[AZ_N_internal] + data_org[AZ_N_border];

#else

total_blk_rows = N_update/num_PDE_eqns;
total_blk_cols = total_blk_rows;
blk_space      = total_blk_rows*20;
nz_space       = blk_space*num_PDE_eqns*num_PDE_eqns;

vbr_cnptr = (int    *) ML_allocate(sizeof(int   )*(total_blk_cols+1));
vbr_rnptr = (int    *) ML_allocate(sizeof(int   )*(total_blk_cols+1));
vbr_bnptr = (int    *) ML_allocate(sizeof(int   )*(total_blk_cols+2));
vbr_indx  = (int    *) ML_allocate(sizeof(int   )*(blk_space+1));
vbr_bindx = (int    *) ML_allocate(sizeof(int   )*(blk_space+1));
vbr_val   = (double *) ML_allocate(sizeof(double)*(nz_space+1));

for (i = 0; i <= total_blk_cols; i++) vbr_cnptr[i] = num_PDE_eqns;

  AZ_msr2vbr(vbr_val, vbr_indx, vbr_rnptr,  vbr_cnptr, vbr_bnptr,
                vbr_bindx, bindx, val,
                total_blk_rows, total_blk_cols, blk_space,
                nz_space, -1);

  data_org[AZ_N_rows]  = data_org[AZ_N_internal] + data_org[AZ_N_border];
  data_org[AZ_N_int_blk]  = data_org[AZ_N_internal]/num_PDE_eqns;
  data_org[AZ_N_bord_blk] = data_org[AZ_N_bord_blk]/num_PDE_eqns;
  data_org[AZ_N_ext_blk]  = data_org[AZ_N_ext_blk]/num_PDE_eqns;
  data_org[AZ_matrix_type] = AZ_VBR_MATRIX;


  AZ_set_VBR(Amat, vbr_rnptr, vbr_cnptr, vbr_bnptr, vbr_indx, vbr_bindx,
             vbr_val, data_org, 0, NULL, AZ_LOCAL);

  Amat->matrix_type  = data_org[AZ_matrix_type];
#endif

  start_time = AZ_second();

  ML_Create(&ml, N_levels);
  ML_Set_PrintLevel(3);


  /* set up discretization matrix and matrix vector function */

  AZ_ML_Set_Amat(ml, N_levels-1, N_update, N_update, Amat, proc_config);

  ML_Aggregate_Create( &ag );
  ML_Aggregate_Set_Threshold(ag,0.0);
  ML_Set_SpectralNormScheme_PowerMethod(ml);
/*
   To run SA:
     a) set damping factor to 1 and use power method
        ML_Aggregate_Set_DampingFactor(ag, 4./3.);
   To run NSA:
     a) set damping factor to 0
        ML_Aggregate_Set_DampingFactor(ag, 0.);
   To run NSR
     a) set damping factor to 1 and use power method
        ML_Aggregate_Set_DampingFactor(ag, 1.);
        ag->Restriction_smoothagg_transpose = ML_FALSE;
        ag->keep_agg_information=1;
        ag->keep_P_tentative=1;
     b) hack code so it calls the energy minimizing restriction
          line 2973 of ml_agg_genP.c
     c) turn on the NSR flag in ml_agg_energy_min.cpp
   To run Emin
     a) set min_eneryg = 2 and keep_agg_info = 1;
      ag->minimizing_energy=2;
      ag->keep_agg_information=1;
      ag->cheap_minimizing_energy = 0;
      ag->block_scaled_SA = 1;
*/
  ag->minimizing_energy=2;
  ag->keep_agg_information=1;
  ag->block_scaled_SA = 1;

  ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, num_PDE_eqns, NULL, N_update);
  ML_Aggregate_Set_MaxCoarseSize( ag, 20);

/*
ML_Aggregate_Set_RandomOrdering( ag );
ML_Aggregate_Set_DampingFactor(ag, .1);
ag->drop_tol_for_smoothing = 1.0e-3;
ML_Aggregate_Set_Threshold(ag, 1.0e-3);
ML_Aggregate_Set_MaxCoarseSize( ag, 300);
*/


	coarsest_level = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml, N_levels-1, ML_DECREASING, ag);
	coarsest_level = N_levels - coarsest_level;
	if ( proc_config[AZ_node] == 0 )
		printf("Coarse level = %d \n", coarsest_level);

	/* set up smoothers */

        AZ_defaults(options, params);

	for (level = N_levels-1; level > coarsest_level; level--) {
          /* This is the Aztec domain decomp/ilu smoother that we */
          /* usually use for this problem.                        */

/*
          options[AZ_precond] = AZ_dom_decomp;
          options[AZ_subdomain_solve] = AZ_ilut;
          params[AZ_ilut_fill] = 1.0;
          options[AZ_reorder] = 1;
          ML_Gen_SmootherAztec(ml, level, options, params,
                        proc_config, status, AZ_ONLY_PRECONDITIONER,
                        ML_PRESMOOTHER,NULL);
*/

          /*  Sparse approximate inverse smoother that acutally does both */
          /*  pre and post smoothing.                                     */
          /*

          ML_Gen_Smoother_ParaSails(ml , level, ML_PRESMOOTHER, nsmooth,
                                parasails_sym, parasails_thresh,
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);

          parasails_thresh /= 4.;
          */


          /* This is the symmetric Gauss-Seidel smoothing. In parallel,    */
          /* it is not a true Gauss-Seidel in that each processor          */
          /* does a Gauss-Seidel on its local submatrix independent of the */
          /* other processors.                                             */
          /*
	  ML_Gen_Smoother_SymGaussSeidel(ml,level,ML_PRESMOOTHER, nsmooth,1.);
	  ML_Gen_Smoother_SymGaussSeidel(ml,level,ML_POSTSMOOTHER,nsmooth,1.);
          */

          /* Block Gauss-Seidel with block size equal to #DOF per node.    */
          /* Not a true Gauss-Seidel in that each processor does a         */
          /* Gauss-Seidel on its local submatrix independent of the other  */
          /* processors.                                                   */
          /*

	  ML_Gen_Smoother_BlockGaussSeidel(ml,level,ML_PRESMOOTHER,
                                           nsmooth,0.67, num_PDE_eqns);
	  ML_Gen_Smoother_BlockGaussSeidel(ml,level,ML_POSTSMOOTHER,
                                           nsmooth, 0.67, num_PDE_eqns);
          */


  	  ML_Gen_Smoother_SymBlockGaussSeidel(ml,level,ML_POSTSMOOTHER,
                                                1, 1.0, num_PDE_eqns);
	}

        ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
	ML_Gen_Solver(ml, ML_MGW, N_levels-1, coarsest_level);
	AZ_defaults(options, params);

        options[AZ_solver]   = AZ_gmres;
        options[AZ_scaling]  = AZ_none;
        options[AZ_precond]  = AZ_user_precond;
/*
        options[AZ_conv]     = AZ_r0;
*/
        options[AZ_output]   = 1;
        options[AZ_max_iter] = 1500;
        options[AZ_poly_ord] = 5;
        options[AZ_kspace]   = 130;
        params[AZ_tol]       = 1.0e-8;
/*
options[AZ_precond] = AZ_dom_decomp;
options[AZ_subdomain_solve] = AZ_ilut;
params[AZ_ilut_fill] = 2.0;
*/

	AZ_set_ML_preconditioner(&Pmat, Amat, ml, options);
setup_time = AZ_second() - start_time;

	xxx = (double *) malloc( leng*sizeof(double));
	rhs=(double *)malloc(leng*sizeof(double));

	for (iii = 0; iii < leng; iii++) xxx[iii] = 0.0;

        /* Set rhs */

        fp = fopen("AZ_capture_rhs.mat","r");
        if (fp == NULL) {
           if (proc_config[AZ_node] == 0) printf("taking random vector for rhs\n");
           AZ_random_vector(rhs, data_org, proc_config);
           AZ_reorder_vec(rhs, data_org, update_index, NULL);
        }
        else {
           fclose(fp);
	   ivec =(int *)malloc((leng+1)*sizeof(int));
           AZ_input_msr_matrix("AZ_capture_rhs.mat", update, &rhs, &ivec,
                                N_update, proc_config);
           free(ivec);
           AZ_reorder_vec(rhs, data_org, update_index, NULL);
        }

        /* Set x */

        fp = fopen("AZ_capture_init_guess.mat","r");
        if (fp != NULL) {
           fclose(fp);
	   ivec =(int *)malloc((leng+1)*sizeof(int));
           AZ_input_msr_matrix("AZ_capture_init_guess.mat",update, &xxx, &ivec,
                                N_update, proc_config);
           free(ivec);
           AZ_reorder_vec(xxx, data_org, update_index, NULL);
        }

        /* if Dirichlet BC ... put the answer in */

        for (i = 0; i < data_org[AZ_N_internal]+data_org[AZ_N_border]; i++) {
           if ( (val[i] > .99999999) && (val[i] < 1.0000001))
              xxx[i] = rhs[i];
        }

        fp = fopen("AZ_no_multilevel.dat","r");
        scaling = AZ_scaling_create();
start_time = AZ_second();
        if (fp != NULL) {
           fclose(fp);
           options[AZ_precond] = AZ_none;
           options[AZ_scaling] = AZ_sym_diag;
           options[AZ_ignore_scaling] = AZ_TRUE;

           options[AZ_keep_info] = 1;
           AZ_iterate(xxx, rhs, options, params, status, proc_config, Amat, NULL, scaling);

/*
           options[AZ_pre_calc] = AZ_reuse;
           options[AZ_conv] = AZ_expected_values;
           if (proc_config[AZ_node] == 0)
              printf("\n-------- Second solve with improved convergence test -----\n");
           AZ_iterate(xxx, rhs, options, params, status, proc_config, Amat, NULL, scaling);
           if (proc_config[AZ_node] == 0)
              printf("\n-------- Third solve with improved convergence test -----\n");
           AZ_iterate(xxx, rhs, options, params, status, proc_config, Amat, NULL, scaling);
*/
        }
        else {
           options[AZ_keep_info] = 1;
           AZ_iterate(xxx, rhs, options, params, status, proc_config, Amat, Pmat, scaling);
           options[AZ_pre_calc] = AZ_reuse;
           options[AZ_conv] = AZ_expected_values;
/*
           if (proc_config[AZ_node] == 0)
              printf("\n-------- Second solve with improved convergence test -----\n");
           AZ_iterate(xxx, rhs, options, params, status, proc_config, Amat, Pmat, scaling);
           if (proc_config[AZ_node] == 0)
              printf("\n-------- Third solve with improved convergence test -----\n");
           AZ_iterate(xxx, rhs, options, params, status, proc_config, Amat, Pmat, scaling);
*/
        }
   solve_time = AZ_second() - start_time;

   if (proc_config[AZ_node] == 0)
      printf("Solve time = %e, MG Setup time = %e\n", solve_time, setup_time);

   ML_Aggregate_Destroy(&ag);
   ML_Destroy(&ml);
   AZ_free((void *) Amat->data_org);
   AZ_free((void *) Amat->val);
   AZ_free((void *) Amat->bindx);
   AZ_free((void *) update);
   AZ_free((void *) external);
   AZ_free((void *) extern_index);
   AZ_free((void *) update_index);
   AZ_scaling_destroy(&scaling);
   if (Amat  != NULL) AZ_matrix_destroy(&Amat);
   if (Pmat  != NULL) AZ_precond_destroy(&Pmat);
   free(xxx);
   free(rhs);


#ifdef ML_MPI
  MPI_Finalize();
#endif

  return 0;

}

#else

int main(int argc, char *argv[])
{
#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
  puts("Please configure ML with --enable-aztecoo to run this example");
#ifdef ML_MPI
  MPI_Finalize();
#endif
  /* returns ok not to break the test harness */
  return(EXIT_SUCCESS);
}

#endif /* HAVE_ML_AZTECOO || HAVE_ML_AZTEC2_1 */

