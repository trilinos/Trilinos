#ifndef MB_MODIF
#define MB_MODIF
#endif
/*
#define ML_partition 
*/
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"
#include <math.h>


extern int AZ_using_fortran;
int    parasails_factorized = 0;
int    parasails_sym        = 1;
double parasails_thresh     = 0.01;
int    parasails_nlevels    = 0;
double parasails_filter     = 0.;
double parasails_loadbal    = 0.;
int    which_filename = 0;


  int    *data_org = NULL, *update = NULL, *external = NULL;
  int    *update_index = NULL, *extern_index = NULL;
  int    *cpntr = NULL, *bindx = NULL, N_update, iii;

double *scaling_vect = NULL;
#define SCXLE_ME

int main(int argc, char *argv[])
{
	int num_PDE_eqns=6, N_levels=4, nsmooth=2;

	int    leng, level, N_grid_pts, coarsest_level;

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];

  /* data structure for matrix corresponding to the fine grid */

  double *val = NULL, *xxx, *rhs, solve_time, setup_time, start_time;
  AZ_MATRIX *Amat;
  AZ_PRECOND *Pmat = NULL;
  ML *ml;
  FILE *fp;
  int ch,i, j, Nrigid, *garbage = NULL, nblocks, *blocks;
  struct AZ_SCALING *scaling;
  ML_Aggregate *ag;
double *mode, *rigid;
char filename[80];
double alpha;
int    one = 1;
int allocated = 0, *newbindx, offset, current, *block_list = NULL,  k, block;
double *newval;
int old_prec, old_sol;
double old_tol;
double *Amode, beta, biggest;
ML_Operator *Amatrix;
int *rowi_col = NULL, rowi_N, count2, ccc;
double *rowi_val = NULL;
int max_nz_row, big_ind = -1, ii;
double max_diag, min_diag, max_sum, sum;
 int nBlocks, Nper_block, *blockIndices, Ndof;
#ifdef ML_partition
   FILE *fp2;
   int count;

   if (argc != 2) {
     printf("Usage: ml_read_elas num_processors\n");
     exit(1);
   }
   else sscanf(argv[1],"%d",&nblocks);
#endif

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  /* read in the number of matrix equations */
  leng = 0;
  if (proc_config[AZ_node] == 0) {
#    ifdef binary
	fp=fopen(".data","rb");
#    else
	fp=fopen(".data","r");
#    endif
     if (fp==NULL) {
        printf("couldn't open file .data\n");
        exit(1);
     }
#    ifdef binary
        fread(&leng, sizeof(int), 1, fp);
#    else
        fscanf(fp,"%d",&leng);
#    endif
     fclose(fp);
  }
  leng = AZ_gsum_int(leng, proc_config);

  N_grid_pts=leng/num_PDE_eqns;

  /* initialize the list of global indices. NOTE: the list of global */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly. */

  if (proc_config[AZ_N_procs] == 1) i = AZ_linear;
  else i = AZ_file;
  AZ_read_update(&N_update, &update, proc_config, N_grid_pts, num_PDE_eqns,i);

  AZ_read_msr_matrix(update, &val, &bindx, N_update, proc_config);


  /* This code is to fix things up so that we are sure we have */ 
  /* all block (including the ghost nodes the same size.       */

  AZ_block_MSR(&bindx, &val, N_update, num_PDE_eqns, update);

  AZ_transform_norowreordering(proc_config, &external, bindx, val,  update, &update_index,
	       &extern_index, &data_org, N_update, 0, 0, 0, &cpntr,
	       AZ_MSR_MATRIX);
	
  Amat = AZ_matrix_create( leng );
  AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);

  Amat->matrix_type  = data_org[AZ_matrix_type];
	
  data_org[AZ_N_rows]  = data_org[AZ_N_internal] + data_org[AZ_N_border];

#ifdef SCALE_ME
  ML_MSR_sym_diagonal_scaling(Amat, proc_config, &scaling_vect);  
#endif
			
  start_time = AZ_second();

  options[AZ_scaling] = AZ_none;
  ML_Create(&ml, N_levels);
  ML_Set_PrintLevel(10);
			
			
  /* set up discretization matrix and matrix vector function */
	
  AZ_ML_Set_Amat(ml, N_levels-1, N_update, N_update, Amat, proc_config);

#ifdef ML_partition

  /* this code is meant to partition the matrices so that things can be */
  /* run in parallel later.                                             */
  /* It is meant to be run on only one processor.                       */
#ifdef	MB_MODIF
  fp2 = fopen(".update","w");
#else
  fp2 = fopen("partition_file","w");
#endif

  ML_Operator_AmalgamateAndDropWeak(&(ml->Amat[N_levels-1]), num_PDE_eqns, 0.0);
  ML_Gen_Blocks_Metis(ml, N_levels-1, &nblocks, &block_list);

  for (i = 0; i < nblocks; i++) {
     count = 0;
     for (j = 0; j < ml->Amat[N_levels-1].outvec_leng; j++) {
        if (block_list[j] == i) count++;
     }
     fprintf(fp2,"   %d\n",count*num_PDE_eqns);
     for (j = 0; j < ml->Amat[N_levels-1].outvec_leng; j++) {
        if (block_list[j] == i) {
           for (k = 0; k < num_PDE_eqns; k++)  fprintf(fp2,"%d\n",j*num_PDE_eqns+k);
        }
     }
  }
  fclose(fp2);
  ML_Operator_UnAmalgamateAndDropWeak(&(ml->Amat[N_levels-1]),num_PDE_eqns,0.0);
#ifdef	MB_MODIF
  printf(" partition file dumped in .update\n");
#endif
  exit(1);
#endif
	
  ML_Aggregate_Create( &ag );
/*
  ML_Aggregate_Set_CoarsenScheme_MIS(ag);
*/
#ifdef MB_MODIF
  ML_Aggregate_Set_DampingFactor(ag,1.50);
#else
  ML_Aggregate_Set_DampingFactor(ag,1.5);
#endif
  ML_Aggregate_Set_CoarsenScheme_METIS(ag);
  ML_Aggregate_Set_NodesPerAggr( ml, ag, -1, 35);
  /*
  ML_Aggregate_Set_Phase3AggregateCreationAggressiveness(ag, 10.001);
  */


  ML_Aggregate_Set_Threshold(ag, 0.0);
  ML_Aggregate_Set_MaxCoarseSize( ag, 300);


  /* read in the rigid body modes */

   Nrigid = 0;

  /* to ensure compatibility with RBM dumping software */
   if (proc_config[AZ_node] == 0) {

      sprintf(filename,"rigid_body_mode%02d",Nrigid+1);
      while( (fp = fopen(filename,"r")) != NULL) {
	which_filename = 1;
          fclose(fp);
          Nrigid++;
          sprintf(filename,"rigid_body_mode%02d",Nrigid+1);
      }
      sprintf(filename,"rigid_body_mode%d",Nrigid+1);
      while( (fp = fopen(filename,"r")) != NULL) {
          fclose(fp);
          Nrigid++;
          sprintf(filename,"rigid_body_mode%d",Nrigid+1);
      }
    }

    Nrigid = AZ_gsum_int(Nrigid,proc_config);

    if (Nrigid != 0) {
       rigid = (double *) ML_allocate( sizeof(double)*Nrigid*(N_update+1) );
       if (rigid == NULL) {
          printf("Error: Not enough space for rigid body modes\n");
       }
    }

    rhs   = (double *) malloc(leng*sizeof(double));
    xxx   = (double *) malloc(leng*sizeof(double));

    for (iii = 0; iii < leng; iii++) xxx[iii] = 0.0; 
	


    for (i = 0; i < Nrigid; i++) {
       if (which_filename == 1) sprintf(filename,"rigid_body_mode%02d",i+1);
       else sprintf(filename,"rigid_body_mode%d",i+1);
       AZ_input_msr_matrix(filename,update,&mode,&garbage,N_update,proc_config);
       AZ_reorder_vec(mode, data_org, update_index, NULL);
       /* here is something to stick a rigid body mode as the initial */
       /* The idea is to solve A x = 0 without smoothing with a two   */
       /* level method. If everything is done properly, we should     */
       /* converge in 2 iterations.                                   */
       /* Note: we must also zero out components of the rigid body    */
       /* mode that correspond to Dirichlet bcs.                      */

       if (i == -4) {
          for (iii = 0; iii < leng; iii++) xxx[iii] = mode[iii];

          ccc = 0;
          Amatrix = &(ml->Amat[N_levels-1]);
          for (iii = 0; iii < Amatrix->outvec_leng; iii++) {
             ML_get_matrix_row(Amatrix,1,&iii,&allocated,&rowi_col,&rowi_val,
                               &rowi_N, 0);
             count2 = 0;
             for (j = 0; j < rowi_N; j++) if (rowi_val[j] != 0.) count2++;
             if (count2 <= 1) { xxx[iii] = 0.; ccc++; }
          }
          free(rowi_col); free(rowi_val);
          allocated = 0; rowi_col = NULL; rowi_val = NULL;
       }

       /*
        *  Rescale matrix/rigid body modes and checking 
        *
        AZ_sym_rescale_sl(mode, Amat->data_org, options, proc_config, scaling);
        Amat->matvec(mode, rigid, Amat, proc_config);
        for (j = 0; j < N_update; j++) printf("this is %d %e\n",j,rigid[j]);
        */

        /* Here is some code to check that the rigid body modes are  */
        /* really rigid body modes. The idea is to multiply by A and */
        /* then to zero out things that we "think" are boundaries.   */
        /* In this hardwired example, things near boundaries         */
        /* correspond to matrix rows that do not have 81 nonzeros.   */
        /*

        Amode = (double *) malloc(leng*sizeof(double));
        Amat->matvec(mode, Amode, Amat, proc_config);
        j = 0;
        biggest = 0.0;
        for (ii = 0; ii < N_update; ii++) {
           if ( Amat->bindx[ii+1] - Amat->bindx[ii] != 80) {
              Amode[ii] = 0.; j++;
           }
           else { 
              if ( fabs(Amode[ii]) > biggest) {
                 biggest=fabs(Amode[ii]); big_ind = ii;
              }
           }
        }
        printf("%d entries zeroed out of %d elements\n",j,N_update);
        alpha = AZ_gdot(N_update, Amode, Amode, proc_config);
        beta  = AZ_gdot(N_update,  mode,  mode, proc_config);
        printf("||A r||^2 =%e, ||r||^2 = %e, ratio = %e\n",
               alpha,beta,alpha/beta);
        printf("the biggest is %e at row %d\n",biggest,big_ind);
        free(Amode);

        */

        /* orthogonalize mode with respect to previous modes. */

        for (j = 0; j < i; j++) {
           alpha = -AZ_gdot(N_update, mode, &(rigid[j*N_update]), proc_config)/
                    AZ_gdot(N_update, &(rigid[j*N_update]), 
                               &(rigid[j*N_update]), proc_config);
	   /*           daxpy_(&N_update,&alpha,&(rigid[j*N_update]),  &one, mode, &one); */
        }
#ifndef	MB_MODIF
       printf(" after mb %e %e %e\n",mode[0],mode[1],mode[2]);
#endif

        for (j = 0; j < N_update; j++) rigid[i*N_update+j] = mode[j];
        free(mode);
        free(garbage); garbage = NULL;

    }

    if (Nrigid != 0) {
             ML_Aggregate_Set_BlockDiagScaling(ag);
       ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Nrigid, rigid, N_update);
       free(rigid);
    }
#ifdef SCALE_ME
    ML_Aggregate_Scale_NullSpace(ag, scaling_vect, N_update);
#endif

    coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, N_levels-1, 
				ML_DECREASING, ag);
   AZ_defaults(options, params);
   coarsest_level = N_levels - coarsest_level;
   if ( proc_config[AZ_node] == 0 )
	printf("Coarse level = %d \n", coarsest_level);
	
   /* set up smoothers */
	
   for (level = N_levels-1; level > coarsest_level; level--) {
		
/*
      ML_Gen_Smoother_BlockGaussSeidel(ml, level,ML_BOTH, 1, 1., num_PDE_eqns);
*/

    /*  Sparse approximate inverse smoother that acutally does both */
    /*  pre and post smoothing.                                     */
    /*
      ML_Gen_Smoother_ParaSails(ml , level, ML_PRESMOOTHER, nsmooth, 
                                parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
     */

     /* This is the symmetric Gauss-Seidel smoothing that we usually use. */
     /* In parallel, it is not a true Gauss-Seidel in that each processor */
     /* does a Gauss-Seidel on its local submatrix independent of the     */
     /* other processors.                                                 */

     /* ML_Gen_Smoother_MLS(ml, level, ML_BOTH, 30., nsmooth); */
     Ndof = ml->Amat[level].invec_leng;

     ML_Gen_Blocks_Aggregates(ag, level, &nBlocks, &blockIndices);

     ML_Gen_Smoother_BlockDiagScaledCheby(ml, level, ML_BOTH, 30.,nsmooth,
					  nBlocks, blockIndices);

     /*
      ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_BOTH, nsmooth,1.);
     */
     

      /* This is a true Gauss Seidel in parallel. This seems to work for  */
      /* elasticity problems.  However, I don't believe that this is very */
      /* efficient in parallel.                                           */       
     /*
      nblocks = ml->Amat[level].invec_leng/num_PDE_eqns;
      blocks = (int *) ML_allocate(sizeof(int)*N_update);
      for (i =0; i < ml->Amat[level].invec_leng; i++) 
         blocks[i] = i/num_PDE_eqns;

      ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml , level, ML_PRESMOOTHER,
                                                  nsmooth, 1., nblocks, blocks);
      ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml, level, ML_POSTSMOOTHER,
                                                  nsmooth, 1., nblocks, blocks);
      free(blocks);
*/

      /* Block Jacobi Smoothing */
      /*
      nblocks = ml->Amat[level].invec_leng/num_PDE_eqns;
      blocks = (int *) ML_allocate(sizeof(int)*N_update);
      for (i =0; i < ml->Amat[level].invec_leng; i++) 
         blocks[i] = i/num_PDE_eqns;

      ML_Gen_Smoother_VBlockJacobi(ml , level, ML_BOTH, nsmooth, 
                                   ML_ONE_STEP_CG, nblocks, blocks);
      free(blocks);
      */
   
      /* Jacobi Smoothing                                                 */
     /*
     
      ML_Gen_Smoother_Jacobi(ml , level, ML_PRESMOOTHER, nsmooth, ML_ONE_STEP_CG);
      ML_Gen_Smoother_Jacobi(ml , level, ML_POSTSMOOTHER, nsmooth,ML_ONE_STEP_CG);
     */
     


      /*  This does a block Gauss-Seidel (not true GS in parallel)        */
      /*  where each processor has 'nblocks' blocks.                      */
      /*
      nblocks = 250;
      ML_Gen_Blocks_Metis(ml, level, &nblocks, &blocks);
      ML_Gen_Smoother_VBlockJacobi(ml , level, ML_BOTH, nsmooth,ML_ONE_STEP_CG,
                                        nblocks, blocks);
      free(blocks);
      */
      num_PDE_eqns = 6;
   }
   /* Choose coarse grid solver: mls, superlu, symGS, or Aztec */

   /*
   ML_Gen_Smoother_MLS(ml, coarsest_level, ML_BOTH, 30., nsmooth); 	   
   ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
   */
   /*
   ML_Gen_Smoother_SymGaussSeidel(ml , coarsest_level, ML_BOTH, nsmooth,1.);
   */

   old_prec = options[AZ_precond];
   old_sol  = options[AZ_solver];
   old_tol  = params[AZ_tol];
   params[AZ_tol] = 1.0e-9;
   params[AZ_tol] = 1.0e-5;
   options[AZ_precond] = AZ_Jacobi;
   options[AZ_solver]  = AZ_cg;
   options[AZ_poly_ord] = 1;
   options[AZ_conv] = AZ_r0;
   options[AZ_orth_kvecs] = AZ_TRUE;

   j = AZ_gsum_int(ml->Amat[coarsest_level].outvec_leng, proc_config); 

   options[AZ_keep_kvecs] = j - 6;
   options[AZ_max_iter] =  options[AZ_keep_kvecs];

   ML_Gen_SmootherAztec(ml, coarsest_level, options, params,
            proc_config, status, options[AZ_keep_kvecs], ML_PRESMOOTHER, NULL);

   options[AZ_conv] = AZ_noscaled;
   options[AZ_keep_kvecs] = 0;
   options[AZ_orth_kvecs] = 0;
   options[AZ_precond] = old_prec;
   options[AZ_solver] = old_sol;
   params[AZ_tol] = old_tol;

   /*   */


#ifdef RST_MODIF
   ML_Gen_Solver(ml, ML_MGV, N_levels-1, coarsest_level); 
#else
#ifdef	MB_MODIF
   ML_Gen_Solver(ml, ML_SAAMG,   N_levels-1, coarsest_level); 
#else
   ML_Gen_Solver(ml, ML_MGFULLV, N_levels-1, coarsest_level); 
#endif
#endif
	
   options[AZ_solver]   = AZ_GMRESR;
         options[AZ_solver]   = AZ_cg;
   options[AZ_scaling]  = AZ_none;
   options[AZ_precond]  = AZ_user_precond;
   options[AZ_conv]     = AZ_r0;
   options[AZ_conv] = AZ_noscaled;
   options[AZ_output]   = 1;
   options[AZ_max_iter] = 500;
   options[AZ_poly_ord] = 5;
   options[AZ_kspace]   = 40;
   params[AZ_tol]       = 4.8e-6;

   AZ_set_ML_preconditioner(&Pmat, Amat, ml, options); 
   setup_time = AZ_second() - start_time;
	
   /* Set rhs */

   fp = fopen("AZ_capture_rhs.dat","r");
   if (fp == NULL) {
      AZ_random_vector(rhs, data_org, proc_config);
      if (proc_config[AZ_node] == 0) printf("taking random vector for rhs\n");
      for (i = 0; i < -N_update; i++) {
        rhs[i] = (double) update[i]; rhs[i] = 7.;
      }
   }
   else {
      if (proc_config[AZ_node]== 0) printf("reading rhs guess from file\n");
      AZ_input_msr_matrix("AZ_capture_rhs.dat", update, &rhs, &garbage,
			  N_update, proc_config);
      free(garbage);
   }
   AZ_reorder_vec(rhs, data_org, update_index, NULL);

   printf("changing rhs by multiplying with A\n");
  Amat->matvec(rhs, xxx, Amat, proc_config);
  for (i = 0; i < N_update; i++) rhs[i] = xxx[i];

   fp = fopen("AZ_capture_init_guess.dat","r");
   if (fp != NULL) {
      fclose(fp);
      if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
      AZ_input_msr_matrix("AZ_capture_init_guess.dat", update, &xxx, &garbage,
      			  N_update, proc_config);
      free(garbage);


      xxx = (double *) realloc(xxx, sizeof(double)*(
					 Amat->data_org[AZ_N_internal]+
					 Amat->data_org[AZ_N_border] +
					 Amat->data_org[AZ_N_external]));
   }
   AZ_reorder_vec(xxx, data_org, update_index, NULL);

   /* if Dirichlet BC ... put the answer in */

/*
   for (i = 0; i < data_org[AZ_N_internal]+data_org[AZ_N_border]; i++) {
      if ( (val[i] > .99999999) && (val[i] < 1.0000001))
         xxx[i] = rhs[i];      
   }
*/

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
      options[AZ_conv] = AZ_noscaled;
      options[AZ_conv] = AZ_r0;
      params[AZ_tol] = 1.0e-7;
      /* ML_Iterate(ml, xxx, rhs); */
alpha = sqrt(AZ_gdot(N_update, xxx, xxx, proc_config));
printf("init guess = %e\n",alpha);
alpha = sqrt(AZ_gdot(N_update, rhs, rhs, proc_config));
printf("rhs = %e\n",alpha);
#ifdef SCALE_ME
	ML_MSR_scalerhs(rhs, scaling_vect, data_org[AZ_N_internal] +
                    data_org[AZ_N_border]);
	ML_MSR_scalesol(xxx, scaling_vect, data_org[AZ_N_internal] +
			data_org[AZ_N_border]);
#endif

max_diag = 0.;
min_diag = 1.e30;
max_sum  = 0.;
for (i = 0; i < N_update; i++) {
   if (Amat->val[i] < 0.) printf("woops negative diagonal A(%d,%d) = %e\n",
				 i,Amat->val[i]);
   if (Amat->val[i] > max_diag) max_diag = Amat->val[i];
   if (Amat->val[i] < min_diag) min_diag = Amat->val[i];
   sum = fabs(Amat->val[i]);
   for (j = Amat->bindx[i]; j < Amat->bindx[i+1]; j++) {
      sum += fabs(Amat->val[j]);
   }
   if (sum > max_sum) max_sum = sum;
}
printf("Largest diagonal = %e, min diag = %e large abs row sum = %e\n",
max_diag, min_diag, max_sum);

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
   if (proc_config[AZ_node] == 0)
     printf("Printing out a few entries of the solution ...\n");

   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 7) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],xxx[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 23) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],xxx[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 47) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],xxx[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 101) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],xxx[update_index[j]]); fflush(stdout);}
   j = AZ_gsum_int(7, proc_config); /* sync processors */
   for (j=0;j<Amat->data_org[AZ_N_internal]+ Amat->data_org[AZ_N_border];j++)
     if (update[j] == 171) {printf("solution(gid = %d) = %10.4e\n",
			      update[j],xxx[update_index[j]]); fflush(stdout);}


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

