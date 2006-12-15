/*****************************************************************************/
/* Copyright 1999, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

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
#if (defined(HAVE_ML_AZTEC2_1) || defined(HAVE_ML_AZTECOO))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"
#include <math.h>


extern int AZ_using_fortran;


int main(int argc, char *argv[])
{
	int num_PDE_eqns=1, N_levels=10, nsmooth=1;

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

#ifdef HAVE_MPI
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
			printf("Couldn't open file .data\nYou may need to gunzip it and copy it from the example source directory.\n");
#ifdef HAVE_MPI
                        MPI_Finalize();
#endif
			exit(EXIT_SUCCESS);
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

	AZ_transform( proc_config, &external, bindx,
								val,  update, &update_index,
								&extern_index, &data_org, N_update,
								0, 0, 0, &cpntr, AZ_MSR_MATRIX);
	
	Amat = AZ_matrix_create( leng );
	AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);

	Amat->matrix_type  = data_org[AZ_matrix_type];
	
	data_org[AZ_N_rows]  = data_org[AZ_N_internal] + data_org[AZ_N_border];
			
start_time = AZ_second();
	
	ML_Create(&ml, N_levels);
			
			
	/* set up discretization matrix and matrix vector function */
	
	AZ_ML_Set_Amat(ml, N_levels-1, N_update, N_update, Amat, proc_config);
	
        ML_Aggregate_Create( &ag );
        ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag);
        ML_Aggregate_Set_Threshold(ag,0.0);


	coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, N_levels-1, ML_DECREASING, ag);
	coarsest_level = N_levels - coarsest_level;
	if ( proc_config[AZ_node] == 0 )
		printf("Coarse level = %d \n", coarsest_level);
	
	/* set up smoothers */
	
	for (level = N_levels-1; level > coarsest_level; level--) {
		
          /* This is the symmetric Gauss-Seidel smoothing. In parallel,    */
          /* it is not a true Gauss-Seidel in that each processor          */
          /* does a Gauss-Seidel on its local submatrix independent of the */
          /* other processors.                                             */

	  ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_PRESMOOTHER, nsmooth,1.);
	  ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_POSTSMOOTHER, nsmooth,1.);

          /*
	  ML_Gen_Smoother_Jacobi(ml , level, ML_PRESMOOTHER, nsmooth, .67);
	  ML_Gen_Smoother_Jacobi(ml , level, ML_POSTSMOOTHER, nsmooth, .67 );
          */
		
	}
	
	ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
		
	
	ML_Gen_Solver(ml, ML_MGV, N_levels-1, coarsest_level); 
	AZ_defaults(options, params);
	
        options[AZ_solver]   = AZ_cg;
        options[AZ_scaling]  = AZ_none;
        options[AZ_precond]  = AZ_user_precond;
        options[AZ_conv]     = AZ_r0;
        options[AZ_output]   = 1;
        options[AZ_max_iter] = 1500;
        options[AZ_poly_ord] = 5;
        options[AZ_kspace]   = 130;
        params[AZ_tol]       = 1.0e-8;
	
	AZ_set_ML_preconditioner(&Pmat, Amat, ml, options); 
setup_time = AZ_second() - start_time;
	
	xxx = (double *) malloc( leng*sizeof(double));
	rhs=(double *)malloc(leng*sizeof(double));

	for (iii = 0; iii < leng; iii++) xxx[iii] = 0.0; 
	
        /* Set rhs */
 
        fp = fopen("AZ_capture_rhs.dat","r");
        if (fp == NULL) {
           if (proc_config[AZ_node] == 0) printf("taking random vector for rhs\n");
           AZ_random_vector(rhs, data_org, proc_config);
           AZ_reorder_vec(rhs, data_org, update_index, NULL);
        }
        else {
           ch = getc(fp);
           if (ch == 'S') {
              while ( (ch = getc(fp)) != '\n') ;
           }
           else ungetc(ch,fp);
           for (i = 0; i < data_org[AZ_N_internal]+data_org[AZ_N_border]; i++) 
              fscanf(fp,"%lf",&(rhs[i]));
           fclose(fp);
        }

        /* Set x */

        fp = fopen("AZ_capture_init_guess.dat","r");
        if (fp != NULL) {
           ch = getc(fp);
           if (ch == 'S') {
              while ( (ch = getc(fp)) != '\n') ;
           }
           else ungetc(ch,fp);
           for (i = 0; i < data_org[AZ_N_internal]+data_org[AZ_N_border]; i++)
              fscanf(fp,"%lf",&(xxx[i]));
           fclose(fp);
           options[AZ_conv] = AZ_expected_values;
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


#ifdef HAVE_MPI
  MPI_Finalize();
#endif
	
  return 0;
	
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{

  /* still need to deal with MPI, some architecture don't like
     an exit(0) without MPI_Finalize() */
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("This test requires Aztec.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif
