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
#                              [Parasails, GaussSeidel, SymGaussSeidel, Poly, 
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis]
Smoother steps per level     = 7
Coarse grid solver           = SuperLU
#                              [Parasails, GaussSeidel, SymGaussSeidel, Poly,
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


  int    *data_org = NULL, *update = NULL, *external = NULL;
  int    *update_index = NULL, *extern_index = NULL;
  int    *cpntr = NULL, *bindx = NULL, N_update, iii;

struct reader_context *context;

int main(int argc, char *argv[])
{
	int num_PDE_eqns=1, N_levels=3, nsmooth=2;

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
  int i, j, Nrigid, *garbage, nblocks, *blocks;
  struct AZ_SCALING *scaling;
  ML_Aggregate *ag;
  double *mode, *rigid, alpha;
  char filename[80];
  int    one = 1;
#ifdef ML_partition
  int *block_list = NULL, k;
#endif

#ifdef ML_partition
   FILE *fp2;
   int count;

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
  ML_Set_PrintLevel(context->output_level);

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
  AZ_input_update(context->partition_file,&N_update, &update, proc_config, 
                   N_grid_pts, num_PDE_eqns,i);
	
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
			
  start_time = AZ_second();

  options[AZ_scaling] = AZ_none;
  ML_Create(&ml, N_levels);
			
			
  /* set up discretization matrix and matrix vector function */
	
  AZ_ML_Set_Amat(ml, N_levels-1, N_update, N_update, Amat, proc_config);

#ifdef ML_partition

  /* this code is meant to partition the matrices so that things can be */
  /* run in parallel later.                                             */
  /* It is meant to be run on only one processor.                       */
  fp2 = fopen("partition_file","w");

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
  exit(1);
#endif

	
  ML_Set_ResidualOutputFrequency(ml, context->output);
  ML_Set_Tolerance(ml, context->tol);
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
       rigid = (double *) ML_allocate( sizeof(double)*Nrigid*(N_update+1) );
       if (rigid == NULL) {
          printf("Error: Not enough space for rigid body modes\n");
       }
    }

   /* Set rhs */

   rhs=(double *)malloc(leng*sizeof(double));
 
   fp = fopen("rhsfile","r");
   if (fp == NULL) {
      if (proc_config[AZ_node] == 0) printf("taking linear vector for rhs\n");
      for (i = 0; i < N_update; i++) rhs[i] = (double) update[i];
   }
   else {
      fclose(fp);
      if (proc_config[AZ_node] == 0) printf("reading rhs from a file\n");
      AZ_input_msr_matrix("rhsfile", update, &rhs, &garbage, N_update, 
                          proc_config);
   }
   AZ_reorder_vec(rhs, data_org, update_index, NULL);

   for (i = 0; i < Nrigid; i++) {
      sprintf(filename,"rigid_body_mode%d",i+1);
      AZ_input_msr_matrix(filename, update, &mode, &garbage, N_update, 
                          proc_config);
      AZ_reorder_vec(mode, data_org, update_index, NULL);

   /*
    *  Rescale matrix/rigid body modes and checking 
    *
       AZ_sym_rescale_sl(mode, Amat->data_org, options, proc_config, scaling);
       Amat->matvec(mode, rigid, Amat, proc_config);
       for (j = 0; j < N_update; j++) printf("this is %d %e\n",j,rigid[j]);
    */

    for (j = 0; j < i; j++) {
       alpha = -AZ_gdot(N_update, mode, &(rigid[j*N_update]), proc_config)/
                  AZ_gdot(N_update, &(rigid[j*N_update]), &(rigid[j*N_update]), 
                               proc_config);
       daxpy_(&N_update, &alpha,  &(rigid[j*N_update]),  &one, mode, &one);
    }
   
    /* rhs orthogonalization */

    alpha = -AZ_gdot(N_update, mode, rhs, proc_config)/
                    AZ_gdot(N_update, mode, mode, proc_config);
    daxpy_(&N_update, &alpha,  mode,  &one, rhs, &one);

    for (j = 0; j < N_update; j++) rigid[i*N_update+j] = mode[j];
    free(mode);
    free(garbage);
  }

  for (j = 0; j < Nrigid; j++) {
     alpha = -AZ_gdot(N_update, rhs, &(rigid[j*N_update]), proc_config)/
              AZ_gdot(N_update, &(rigid[j*N_update]), &(rigid[j*N_update]), 
                      proc_config);
     daxpy_(&N_update, &alpha,  &(rigid[j*N_update]),  &one, rhs, &one);
  }

  if (Nrigid != 0) {
     ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Nrigid, rigid, N_update);
  }


  coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, N_levels-1, 
                                            ML_DECREASING, ag);
  coarsest_level = N_levels - coarsest_level;
  if ( proc_config[AZ_node] == 0 )
	printf("Coarse level = %d \n", coarsest_level);
	
  /* set up smoothers */
  blocks = (int *) ML_allocate(sizeof(int)*N_update);
	
  for (level = N_levels-1; level > coarsest_level; level--) {

      num_PDE_eqns = ml->Amat[level].num_PDEs;
      if (proc_config[AZ_node]==0) printf("block size = %d\n",num_PDE_eqns);
		
     /*  Sparse approximate inverse smoother that acutally does both */
     /*  pre and post smoothing.                                     */

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

     else if (ML_strcmp(context->smoother,"GaussSeidel") == 0) {
       ML_Gen_Smoother_GaussSeidel(ml , level, ML_BOTH, nsmooth,1.);
     }
     else if (ML_strcmp(context->smoother,"SymGaussSeidel") == 0) {
       ML_Gen_Smoother_SymGaussSeidel(ml , level, ML_BOTH, nsmooth,1.);
     }
     else if (ML_strcmp(context->smoother,"Poly") == 0) {
       ML_Gen_Smoother_MLS(ml, level, ML_BOTH, 30., nsmooth);
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
     /* */

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
   }
	
   nsmooth   = context->coarse_its;
   /*  Sparse approximate inverse smoother that acutally does both */
   /*  pre and post smoothing.                                     */

   if (ML_strcmp(context->coarse_solve,"Parasails") == 0) {
        ML_Gen_Smoother_ParaSails(ml , coarsest_level, ML_PRESMOOTHER, nsmooth, 
                                parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
   }

   else if (ML_strcmp(context->coarse_solve,"GaussSeidel") == 0) {
       ML_Gen_Smoother_GaussSeidel(ml , coarsest_level, ML_BOTH, nsmooth,1.);
   }
   else if (ML_strcmp(context->coarse_solve,"Poly") == 0) {
     ML_Gen_Smoother_MLS(ml, coarsest_level, ML_BOTH, 30., nsmooth);
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
		
   ML_Gen_Solver(ml, ML_MGV, N_levels-1, coarsest_level); 
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
   options[AZ_conv]     = AZ_r0;
   options[AZ_output]   = 1;
   options[AZ_max_iter] = 500;
   options[AZ_poly_ord] = 5;
   options[AZ_kspace]   = 130;
   params[AZ_tol]       = context->tol;
   options[AZ_output]   = context->output;
	
   AZ_set_ML_preconditioner(&Pmat, Amat, ml, options); 
   setup_time = AZ_second() - start_time;
	
   xxx = (double *) malloc( leng*sizeof(double));

   for (iii = 0; iii < leng; iii++) xxx[iii] = 0.0; 
	

   /* Set x */

   fp = fopen("initguessfile","r");
   if (fp != NULL) {
      fclose(fp);
      if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
      AZ_input_msr_matrix("initguessfile", update, &xxx, &garbage, N_update, 
                          proc_config);
/*
      ch = getc(fp);
      if (ch == 'S') {
         while ( (ch = getc(fp)) != '\n') ;
      }
      else ungetc(ch,fp);
      for (i = 0; i < data_org[AZ_N_internal]+data_org[AZ_N_border]; i++)
         fscanf(fp,"%lf",&(xxx[i]));
*/
      options[AZ_conv] = AZ_expected_values;
   }
   else if (proc_config[AZ_node]== 0) printf("taking 0 initial guess \n");

   AZ_reorder_vec(xxx, data_org, update_index, NULL);

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


#ifdef ML_MPI
  MPI_Finalize();
#endif
	
  return 0;
	
}
