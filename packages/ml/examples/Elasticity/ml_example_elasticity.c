/*****************************************************************************/
/* Copyright 1999, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************
   Sample driver for AZTEC/ML package. 
   this driver comes together with the examples in the directories

   ExampleMatrices/shell_vb_intersection/ (<-under construction)
   ExampleMatrices/shell_vb_lineal/       
   
   Each directory holds a file "inputfile" specifying ML-options and
   a couple of data*.txt files specifying the variable block system of equations

   usage is :
   ml_example_elasticity.exe ExampleMatrices/shell_vb_lineal
   or:
   ml_example_elasticity.exe ExampleMatrices/shell_vb_intersection

   These are small examples derived elasticity problems using finite element
   shell elements with 6 degrees of freedom per node
   The *_vb_* (Variable Block) directories contain matrices with condensed
   Dirichlet-boundary conditions,
   so the blocks size (default 6x6) is NOT constant.
   They are supposed to be used applying the type of aggregation
   "VBMetis" which can handle variable AND constant blocked matrices
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ml_include.h"
#include "ml_agg_VBMETIS.h"
#include "ml_read_utils.h"
#include "ml_lapack.h"
#include <math.h>
#if defined(HAVE_ML_AZTEC2_1) || defined(HAVE_ML_AZTECOO)
#include "az_aztec.h"
#if 1
#undef ML_partition 
#endif
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

  int leng, level, N_grid_pts, coarsest_level;
  int leng1,leng2;
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
  int i, j, Nrigid, *garbage, nblocks=0, *blocks = NULL, *block_pde=NULL;
  struct AZ_SCALING *scaling;
  ML_Aggregate *ag;
  double *mode, *rigid=NULL, alpha; 
  char filename[80];
  int    one = 1;
  int    proc,nprocs;
  char pathfilename[100];

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
  /* get number of processors and the name of this processor */
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
  proc   = 0;
  nprocs = 1;
#endif

   if (proc_config[AZ_node] == 0) {
      sprintf(pathfilename,"%s/inputfile",argv[1]);
      ML_Reader_ReadInput(pathfilename, &context);
   }
   else context = (struct reader_context *) ML_allocate(sizeof(struct reader_context));
   AZ_broadcast((char *) context,  sizeof(struct reader_context), proc_config,
                AZ_PACK);
   AZ_broadcast((char *) NULL        ,   0          , proc_config, AZ_SEND);

   N_levels = context->N_levels;
   printf("N_levels %d\n",N_levels);
   nsmooth   = context->nsmooth;
   num_PDE_eqns = context->N_dofPerNode;
   printf("num_PDE_eqns %d\n",num_PDE_eqns);

   ML_Set_PrintLevel(context->output_level);

  /* read in the number of matrix equations */
  leng = 0;
  if (proc_config[AZ_node] == 0) {
        sprintf(pathfilename,"%s/data_matrix.txt",argv[1]);
        fp=fopen(pathfilename,"r");
     if (fp==NULL) {
        printf("**ERR** couldn't open file data_matrix.txt\n");
        exit(1);
     }
        fscanf(fp,"%d",&leng);
     fclose(fp);
  }
  leng = AZ_gsum_int(leng, proc_config);

  N_grid_pts=leng/num_PDE_eqns;


  /* initialize the list of global indices. NOTE: the list of global */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly. */
#if 0
  if (proc_config[AZ_N_procs] == 1) i = AZ_linear;
  else i = AZ_file;
#endif
  i = AZ_linear;

  /* cannot use AZ_input_update for variable blocks (forgot why, but debugged through it)*/
  /* make a linear distribution of the matrix       */
  /* if the linear distribution does not align with the blocks, */
  /* this is corrected in ML_AZ_Reader_ReadVariableBlocks */
  leng1 = leng/nprocs;
  leng2 = leng-leng1*nprocs;
  if (proc >= leng2)
  {
     leng2 += (proc*leng1);
  }
  else
  {
     leng1++;
     leng2 = proc*leng1;
  }
  N_update = leng1;
  update = (int*)AZ_allocate((N_update+1)*sizeof(int));
  if (update==NULL)
  {
      (void) fprintf (stderr, "Not enough space to allocate 'update'\n");
      fflush(stderr); exit(EXIT_FAILURE);
  }
  for (i=0; i<N_update; i++) update[i] = i+leng2;
  
#if 0 /* debug */
  printf("proc %d N_update %d\n",proc_config[AZ_node],N_update);
  fflush(stdout);                   
#endif
  sprintf(pathfilename,"%s/data_vblocks.txt",argv[1]);
  ML_AZ_Reader_ReadVariableBlocks(pathfilename,&nblocks,&blocks,&block_pde,
                                  &N_update,&update,proc_config);
#if 0 /* debug */
  printf("proc %d N_update %d\n",proc_config[AZ_node],N_update);
  fflush(stdout);                   
#endif

  sprintf(pathfilename,"%s/data_matrix.txt",argv[1]);
  AZ_input_msr_matrix(pathfilename,update, &val, &bindx, N_update, proc_config);

  /* This code is to fix things up so that we are sure we have   */ 
  /* all blocks (including the ghost nodes) the same size.       */
  /* not sure, whether this is a good idea with variable blocks  */
  /* the examples inpufiles (see top of this file) don't need it */
  /* anyway                                                      */
  /*
  AZ_block_MSR(&bindx, &val, N_update, num_PDE_eqns, update);
  */
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
  AZ_ML_Set_Amat(ml, 0, N_update, N_update, Amat, proc_config);

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
  else if (ML_strcmp(context->agg_coarsen_scheme,"Metis") == 0) {
     ML_Aggregate_Set_CoarsenScheme_METIS(ag);
     for (i=0; i<N_levels; i++)
        ML_Aggregate_Set_NodesPerAggr(ml,ag,i,9);
  }
  else if (ML_strcmp(context->agg_coarsen_scheme,"VBMetis") == 0) {
     /* when no blocks read, use standard metis assuming constant block sizes */
     if (!blocks) 
        ML_Aggregate_Set_CoarsenScheme_METIS(ag);
     else {
        ML_Aggregate_Set_CoarsenScheme_VBMETIS(ag);
        ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS(ag,0,N_levels,nblocks,
                                                       blocks,block_pde,N_update);
     }
     for (i=0; i<N_levels; i++)
        ML_Aggregate_Set_NodesPerAggr(ml,ag,i,9);
  }
  else {
     printf("**ERR** ML: Unknown aggregation scheme %s\n",context->agg_coarsen_scheme);
     exit(-1);
  }
  ML_Aggregate_Set_DampingFactor(ag, context->agg_damping);
  ML_Aggregate_Set_MaxCoarseSize( ag, context->maxcoarsesize);
  ML_Aggregate_Set_Threshold(ag, context->agg_thresh);

  if (ML_strcmp(context->agg_spectral_norm,"Calc") == 0) {
     ML_Set_SpectralNormScheme_Calc(ml);
  }
  else if (ML_strcmp(context->agg_spectral_norm,"Anorm") == 0) {
     ML_Set_SpectralNormScheme_Anorm(ml);
  }
  else {
     printf("**WRN** ML: Unknown spectral norm scheme %s\n",context->agg_spectral_norm);
  }

  /* read in the rigid body modes */

   Nrigid = 0;
   if (proc_config[AZ_node] == 0) {
      sprintf(filename,"data_nullsp%d.txt",Nrigid);
      sprintf(pathfilename,"%s/%s",argv[1],filename);
      while( (fp = fopen(pathfilename,"r")) != NULL) {
          fclose(fp);
          Nrigid++;
          sprintf(filename,"data_nullsp%d.txt",Nrigid);
          sprintf(pathfilename,"%s/%s",argv[1],filename);
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
   sprintf(pathfilename,"%s/data_rhs.txt",argv[1]);
   fp = fopen(pathfilename,"r");
   if (fp == NULL) {
      rhs=(double *)ML_allocate(leng*sizeof(double));
      if (proc_config[AZ_node] == 0) printf("taking linear vector for rhs\n");
      for (i = 0; i < N_update; i++) rhs[i] = (double) update[i];
   }
   else {
      fclose(fp);
      if (proc_config[AZ_node] == 0) printf("reading rhs from a file\n");
      AZ_input_msr_matrix(pathfilename, update, &rhs, &garbage, N_update, 
                          proc_config);
   }
   AZ_reorder_vec(rhs, data_org, update_index, NULL);

   for (i = 0; i < Nrigid; i++) {
      sprintf(filename,"data_nullsp%d.txt",i);
      sprintf(pathfilename,"%s/%s",argv[1],filename);
      AZ_input_msr_matrix(pathfilename, update, &mode, &garbage, N_update, 
                          proc_config);
      AZ_reorder_vec(mode, data_org, update_index, NULL);

#if 0 /* test the given rigid body mode, output-vector should be ~0 */
       Amat->matvec(mode, rigid, Amat, proc_config);
       for (j = 0; j < N_update; j++) printf("this is %d %e\n",j,rigid[j]);
#endif

    for (j = 0; j < i; j++) {
       alpha = -AZ_gdot(N_update, mode, &(rigid[j*N_update]), proc_config)/
                  AZ_gdot(N_update, &(rigid[j*N_update]), &(rigid[j*N_update]), 
                               proc_config);
       DAXPY_F77(&N_update, &alpha,  &(rigid[j*N_update]),  &one, mode, &one);
    }
   
    /* rhs orthogonalization */

    alpha = -AZ_gdot(N_update, mode, rhs, proc_config)/
                    AZ_gdot(N_update, mode, mode, proc_config);
    DAXPY_F77(&N_update, &alpha,  mode,  &one, rhs, &one);

    for (j = 0; j < N_update; j++) rigid[i*N_update+j] = mode[j];
    free(mode);
    free(garbage);
  }

  for (j = 0; j < Nrigid; j++) {
     alpha = -AZ_gdot(N_update, rhs, &(rigid[j*N_update]), proc_config)/
              AZ_gdot(N_update, &(rigid[j*N_update]), &(rigid[j*N_update]), 
                      proc_config);
     DAXPY_F77(&N_update, &alpha,  &(rigid[j*N_update]),  &one, rhs, &one);
  }

#if 0 /* for testing the default nullsp */
  ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, 6, NULL, N_update);
#else
  if (Nrigid != 0) {
     ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Nrigid, rigid, N_update);
  }
#endif
  if (rigid) ML_free(rigid);

  ag->keep_agg_information = 1;
  coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, 0, 
                                            ML_INCREASING, ag);
  coarsest_level--;                                            

  if ( proc_config[AZ_node] == 0 )
	printf("Coarse level = %d \n", coarsest_level);
	
#if 0
  /* set up smoothers */
  if (!blocks)
     blocks = (int *) ML_allocate(sizeof(int)*N_update);
#endif

  for (level = 0; level < coarsest_level; level++) {

      num_PDE_eqns = ml->Amat[level].num_PDEs;
		
     /*  Sparse approximate inverse smoother that acutally does both */
     /*  pre and post smoothing.                                     */

     if (ML_strcmp(context->smoother,"Parasails") == 0) {
        ML_Gen_Smoother_ParaSails(ml , level, ML_PRESMOOTHER, nsmooth, 
                                parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                (int) parasails_loadbal, parasails_factorized);
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
       ML_Gen_Smoother_Cheby(ml, level, ML_BOTH, 30., nsmooth);
     }
     else if (ML_strcmp(context->smoother,"BlockGaussSeidel") == 0) {
       ML_Gen_Smoother_BlockGaussSeidel(ml , level, ML_BOTH, nsmooth,1.,
					 num_PDE_eqns);
     }
     else if (ML_strcmp(context->smoother,"VBSymGaussSeidel") == 0) {
         if (blocks)    ML_free(blocks);
         if (block_pde) ML_free(block_pde);
         blocks    = NULL;
         block_pde = NULL;
         nblocks   = 0;
         ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag,level,N_levels,&nblocks,
                                                        &blocks,&block_pde);
         if (blocks==NULL) ML_Gen_Blocks_Aggregates(ag, level, &nblocks, &blocks);
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
         if (blocks)    ML_free(blocks);
         if (block_pde) ML_free(block_pde);
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
	
   /* set coarse level solver */
   nsmooth   = context->coarse_its;
   /*  Sparse approximate inverse smoother that acutally does both */
   /*  pre and post smoothing.                                     */

   if (ML_strcmp(context->coarse_solve,"Parasails") == 0) {
        ML_Gen_Smoother_ParaSails(ml , coarsest_level, ML_PRESMOOTHER, nsmooth, 
                                parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                (int) parasails_loadbal, parasails_factorized);
   }

   else if (ML_strcmp(context->coarse_solve,"GaussSeidel") == 0) {
       ML_Gen_Smoother_GaussSeidel(ml , coarsest_level, ML_BOTH, nsmooth,1.);
   }
   else if (ML_strcmp(context->coarse_solve,"Poly") == 0) {
     ML_Gen_Smoother_Cheby(ml, coarsest_level, ML_BOTH, 30., nsmooth);
   }
   else if (ML_strcmp(context->coarse_solve,"SymGaussSeidel") == 0) {
       ML_Gen_Smoother_SymGaussSeidel(ml , coarsest_level, ML_BOTH, nsmooth,1.);
   }
   else if (ML_strcmp(context->coarse_solve,"BlockGaussSeidel") == 0) {
       ML_Gen_Smoother_BlockGaussSeidel(ml, coarsest_level, ML_BOTH, nsmooth,1.,
					num_PDE_eqns);
   }
   else if (ML_strcmp(context->coarse_solve,"Aggregate") == 0) {
         if (blocks)    ML_free(blocks);
         if (block_pde) ML_free(block_pde);
         ML_Gen_Blocks_Aggregates(ag, coarsest_level, &nblocks, &blocks);
         ML_Gen_Smoother_VBlockSymGaussSeidel(ml , coarsest_level, ML_BOTH, 
                                        nsmooth,1., nblocks, blocks);
   }
   else if (ML_strcmp(context->coarse_solve,"Jacobi") == 0) {
        ML_Gen_Smoother_Jacobi(ml , coarsest_level, ML_BOTH, nsmooth,.5);
   }
   else if (ML_strcmp(context->coarse_solve,"Metis") == 0) {
         if (blocks)    ML_free(blocks);
         if (block_pde) ML_free(block_pde);
         nblocks = 250;
         ML_Gen_Blocks_Metis(ml, coarsest_level, &nblocks, &blocks);
         ML_Gen_Smoother_VBlockSymGaussSeidel(ml , coarsest_level, ML_BOTH, 
                                              nsmooth,1., nblocks, blocks);
   }
   else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0) {
      ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
   }
   else if (ML_strcmp(context->coarse_solve,"Amesos") == 0) {
      ML_Gen_Smoother_Amesos(ml,coarsest_level,ML_AMESOS_KLU,-1, 0.0);
   }
   else {
         printf("unknown coarse grid solver %s\n",context->coarse_solve);
         exit(1);
   }
		
   ML_Gen_Solver(ml, ML_MGV, 0, coarsest_level); 

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
   if (blocks)            ML_free(blocks);
   if (block_pde)         ML_free(block_pde);
   options[AZ_scaling]  = AZ_none;
   options[AZ_precond]  = AZ_user_precond;
   options[AZ_conv]     = AZ_r0;
   options[AZ_output]   = 1;
   options[AZ_max_iter] = context->max_outer_its;
   options[AZ_poly_ord] = 5;
   options[AZ_kspace]   = 130;
   params[AZ_tol]       = context->tol;
   options[AZ_output]   = context->output;
   ML_free(context);
	
   AZ_set_ML_preconditioner(&Pmat, Amat, ml, options); 
   setup_time = AZ_second() - start_time;
	
   xxx = (double *) malloc( leng*sizeof(double));

   for (iii = 0; iii < leng; iii++) xxx[iii] = 0.0; 
	

   /* Set x */
   /*
   there is no initguess supplied with these examples for the moment....
   */
   fp = fopen("initguessfile","r");
   if (fp != NULL) {
      fclose(fp);
      if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
      AZ_input_msr_matrix("data_initguess.txt", update, &xxx, &garbage, N_update, 
                          proc_config);

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


#else /* do not have aztec */
int main(int argc, char *argv[]) 
{
#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
  fprintf(stdout,"Please configure ML with --enable-aztecoo to run this example\n");
  fflush(stdout);
#ifdef ML_MPI
  MPI_Finalize();
#endif
  /* returns ok not to break the test harness */
  return(EXIT_SUCCESS);
}
#endif /* HAVE_ML_AZTECOO || HAVE_ML_AZTEC2_1 */
