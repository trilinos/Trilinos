/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/*****************************************************************************/
/* Sample driver for segregated solver package. The software is tested by    */
/* reading in matrices stored in 2 files (specified in the body of the code),*/ 
/* and setting up a block system from the first matrix.  It uses ML to       */
/* precondition one diagonal block and preconditions the other diagonal      */
/* block by replacing that block with the other matrix, and then doing one   */
/* step of Jacobi on that block.                                             */
/*                                                                           */
/* This program also tests new AZTEC functions for creating block matrices   */
/* and submatrices from other AZTEC matrices.                                */
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


extern int AZ_using_fortran;


void main(int argc, char *argv[])
{
	int num_PDE_eqns=1, N_levels=30, nsmooth=1;

	int    leng, lengMp, level, N_grid_pts, coarsest_level;

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];
	int    *data_org2;
	int    *update2;
	
  /* data structure for matrix corresponding to the fine grid */

  int    *data_org = NULL, *update = NULL, *external = NULL, *external2=NULL;
  int    *update_index = NULL, *extern_index = NULL;
  int    *update_index2 = NULL, *extern_index2 = NULL;
  int    *cpntr = NULL, *cpntr2 = NULL, i, ctr;
  int    *bindx = NULL, *bindx2 = NULL, N_update;
  double *val = NULL, *val2 = NULL;
	double *xxx, *rhs;

	/* data structures for splitting A into a 2x2 block matrix */
	int Nsubrows[2];
	int Nsubcols[2];
	int **subrows, **subcols;
	int Nblock_rows=2, Nblock_cols=2, **submat_locs;

	/* datafiles containing the main matrix (A) and the matrix which 
		 will replace the (2,2) block of A in the preconditioner (Mp) */
	char tmp, *Amatfile = "../examples/AmatrixAZ.dat\0";
	char *Mpmatfile = "../examples/MpmatrixAZ.dat\0";
	char *rhsFile = "../examples/bvector.dat\0";

	/* structures for storing block and submatrices that we will create
		 with the new AZTEC functions */
	AZ_MATRIX **subAmat_list, *blockAmat;

	AZ_MATRIX *Amat, *Mpmat;
	AZ_PRECOND *Pmat = NULL, *blockPmat;
	ML *ml;
	FILE *fp;

	/* seg is the object that will store the segregated matrix and 
		 preconditioner */
	struct SEG_Struct *seg;


#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, NULL);
#endif


	fp=fopen(Amatfile,"r");
	if (fp==NULL)
		{
			printf("couldn't open A file\n");
			exit(1);
		}
	fscanf(fp,"%d",&leng);
	fclose(fp);

	N_grid_pts=leng/num_PDE_eqns;

	fp=fopen(Mpmatfile,"r");
	if (fp==NULL)
		{
			printf("couldn't open Mp file\n");
			exit(1);
		}
	fscanf(fp,"%d",&lengMp);
	fclose(fp);

	N_grid_pts=leng/num_PDE_eqns;




	
	/* read in A and convert it to an AZTEC MSR matrix */
  AZ_read_update(&N_update, &update, proc_config, N_grid_pts, num_PDE_eqns,
                 AZ_linear);
	
		AZ_input_msr_matrix(Amatfile, update, &val, &bindx, N_update, proc_config);

	AZ_transform( proc_config, &external, bindx,
								val,  update, &update_index,
								&extern_index, &data_org, N_update,
								0, 0, 0, &cpntr, AZ_MSR_MATRIX);

	Amat = AZ_matrix_create( leng );
	AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);


	/* read in Mp and convert it to an AZTEC MSR matrix */
  AZ_read_update(&N_update, &update2, proc_config, N_grid_pts, 1,
                 AZ_linear);
	
	AZ_input_msr_matrix(Mpmatfile, update2, &val2, &bindx2, lengMp, proc_config);

	AZ_transform( proc_config, &external2, bindx2,
								val2,  update2, &update_index2,
								&extern_index2, &data_org2, lengMp,
								0, 0, 0, &cpntr2, AZ_MSR_MATRIX);
	
 	Mpmat = AZ_matrix_create( lengMp );
	AZ_set_MSR(Mpmat, bindx2, val2, data_org2, 0, NULL, AZ_LOCAL);
	

	AZ_defaults(options,params);
 

	/* now set up the arrays that will define the blocks of A.  In this case, 
		 we're taking a matrix in which the rows are grouped by node in the
		 mesh, and we want to separate out the 2 velocity equations from the
		 pressure equation, so the "first" block will be the first 2 of each
		 group of 3 rows in the matrix and the "second" block will be the 
		 remaining 1/3 of the equations */

	/* Note: the later function calls require that the row and column
		 lists be sorted */

	subrows=(int **) malloc(2 * sizeof(int *));
	Nsubrows[0]=2*leng/3;
	Nsubrows[1]=leng/3;
	subrows[0]=(int *) malloc(Nsubrows[0]*sizeof(int));
	subrows[1]=(int *) malloc(Nsubrows[1]*sizeof(int));
	
	subcols=(int **) malloc(2 * sizeof(int *));
	Nsubcols[0]=2*leng/3;
	Nsubcols[1]=leng/3;
	subcols[0]=(int *) malloc(Nsubcols[0]*sizeof(int));
	subcols[1]=(int *) malloc(Nsubcols[1]*sizeof(int));

	
	
	for (ctr=0; ctr < Nsubrows[0]/2; ctr++) {
		subrows[0][2*ctr]=3*ctr;
		subrows[0][2*ctr+1]=3*ctr+1;
	}
	for (ctr=0; ctr < Nsubcols[0]; ctr++)
		subcols[0][ctr]=subrows[0][ctr];
	
	for (ctr=0; ctr < Nsubrows[1]; ctr++) 
		subrows[1][ctr]=3*ctr+2;
	for (ctr=0; ctr < Nsubcols[1]; ctr++)
		subcols[1][ctr]=subrows[1][ctr];
	

	/* first play with submatrices and block matrices.  
		 AZ_blockmatrix_create requires the submatrices to
		 be created before passing them in in order to allow
		 you to create a block matrix from unrelated matrices */
	subAmat_list= (AZ_MATRIX **) malloc(4*sizeof(AZ_MATRIX *));
	
	subAmat_list[0] = AZ_submatrix_create(Amat, Nsubrows[0], 
																				subrows[0], Nsubcols[0], 
																				subcols[0], proc_config);
	subAmat_list[1] = AZ_submatrix_create(Amat, Nsubrows[0], subrows[0], 
																				Nsubcols[1], subcols[1], 
																				proc_config);
	subAmat_list[2] = AZ_submatrix_create(Amat, Nsubrows[1], 
																				subrows[1], Nsubcols[0], 
																				subcols[0], proc_config);
	subAmat_list[3] = AZ_submatrix_create(Amat, Nsubrows[1], 
																				subrows[1], Nsubcols[1], 
																				subcols[1], proc_config);
	
	/* submat_locs give the block row and block col where each submatrix
		 goes */
	submat_locs=(int **) malloc(4*sizeof(int *));
	for (ctr=0; ctr<4; ctr++)
		submat_locs[ctr]=(int *) malloc(2*sizeof(int));
	
	submat_locs[0][0]=0;
	submat_locs[0][1]=0;
	
	submat_locs[1][0]=0;
	submat_locs[1][1]=1;
	
	submat_locs[2][0]=1;
	submat_locs[2][1]=0;
	
	submat_locs[3][0]=1;
	submat_locs[3][1]=1;
	
	blockAmat = AZ_blockmatrix_create(subAmat_list, 4, submat_locs, Nblock_rows, 
																		Nblock_cols, Nsubrows, subrows, Nsubcols, 
																		subcols, proc_config);

	blockPmat = AZ_precond_create(blockAmat, AZ_precondition, NULL); 
	
	/* now create a segregated object from A.  It will create its own
		 versions of the submatrices of A */
	SEG_Create(&seg, Amat, 2, subrows, Nsubrows, SEG_LOWER_TRIANGULAR, proc_config);
	
	
	/* in order to use an ML preconditioner on the (0,0) block of A, we
		 need to create an ML object just as we would to use ML on a 
		 regular AZTEC matrix */
	ML_Create(&ml, N_levels);
	
	
	/* since the submatrices of A are internal to the seg object, have
		 to call this function instead of the usual ML AZ_ML_Set_Amat */
	SEG_ML_Set_Amat(seg, ml, N_levels-1, 0, proc_config);
	
	coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, N_levels-1, 
																											 ML_DECREASING, NULL);
	coarsest_level = N_levels - coarsest_level;
	printf("proc_config[AZ_node]=%d\n",proc_config[AZ_node]);
	if ( proc_config[AZ_node] == 0 )
		printf("Coarse level = %d \n", coarsest_level);
		
	/* set up smoothers */
	for (level = N_levels-1; level > coarsest_level; level--) {
				
		ML_Gen_SmootherGaussSeidel(ml , level, ML_PRESMOOTHER, nsmooth);
		ML_Gen_SmootherGaussSeidel(ml , level, ML_POSTSMOOTHER, nsmooth);
		
	}
		
	ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
	printf("coarsestlevel=%d\n", coarsest_level);
	
	ML_Gen_Solver(ml, ML_MGV, N_levels-1, coarsest_level); 
	

	/* now that we have ml set up fully, we can set it to be the
		 preconditioner for the (0,0) block */

	SEG_Set_ML_Precond(seg, ml, 0, params, options, proc_config);
	

	/* Replace the (1,1) block with Mp since that's the matrix we want to
		 precondition with there */
	SEG_Replace_Submat(seg, 1, 1, Mpmat, proc_config);
	
	/* this simply shows that you can replace off-diagonal block as well.
		 In this case, subAmat_list[2] is identical to the (1,0) block of A
		 so the replacement makes no difference */
	SEG_Replace_Submat(seg, 1, 0, subAmat_list[2], proc_config);
	
	/* I want to precondition the new (1,1) block with Jacobi.  I don't
		 really care what any of the options are except for AZ_precond, so
		 I won't bother to make a new copy here.  However, the following 
		 functions will make a full copy of the options I send in and
		 associate them with this preconditioner */
	options[AZ_precond] = AZ_Jacobi;	
	SEG_Set_AZ_Precond(seg, 1, params, options, proc_config);
	

	/* now create and set the overall preconditioner */	
	Pmat = AZ_precond_create(Amat, AZ_precondition, NULL);
	AZ_Set_SEG_Preconditioner(Pmat, seg, options);
	
		
	/* set various options for the AZTEC solve */
	options[AZ_kspace]=75;
	options[AZ_solver]=AZ_gmres;
	options[AZ_max_iter] = 400;
	params[AZ_tol]=1e-6;
	
	/* read in a right hand side for this solve */
	rhs=(double *)malloc(leng*sizeof(double));
	
	fp=fopen(rhsFile,"r");
	if ((fp==NULL) || (rhs==NULL))
		printf("damn\n");
	fscanf(fp,"%c",&tmp);
	i=0;
	while (tmp != ']'){
		fscanf(fp,"%lf%c",rhs+i,&tmp);
		i++;
	}
	
	/* set the initial guess to zero */
	xxx = (double *) malloc( leng*sizeof(double));
	for (i=0; i < leng; i++)
		xxx[i] = 0.0;
	
	
	/* call AZ_iterate exactly as usual */
	AZ_iterate(xxx, rhs, options, params,
						 status, proc_config, Amat, Pmat, 0);
	
	
	/* now go through and free all the memory we have allocated */

	/* I believe it is important to destroy ml before destroying seg - otherwise
		 it will try to free the same memory twice */
	ML_Destroy(&ml);
	
	SEG_Destroy(&seg);
	AZ_blockmatrix_destroy(&blockAmat);
	
	for (i=0; i < 4; i++)
		AZ_submatrix_destroy(&(subAmat_list[i]));
		
	AZ_precond_destroy(&Pmat);
	AZ_precond_destroy(&blockPmat);
	AZ_matrix_destroy(&Amat);
	
	AZ_matrix_destroy(&Mpmat);
	
	for (i=0; i < 4; i++)
		free(submat_locs[i]);
	
	free(submat_locs);
	free(subAmat_list);
	free(rhs);
	free(xxx);
	free(update);
	free(update2);
	free(update_index);
	free(update_index2);
	free(bindx);
	free(val);
	free(bindx2);
	free(val2);
	free(data_org);
	free(data_org2);
	free(external);
	free(external2);
	free(extern_index);
	free(extern_index2);
	free(subrows[0]);
	free(subrows[1]);
	free(subcols[0]);
	free(subcols[1]);
	free(subrows);
	free(subcols);
	
	
#ifdef ML_MPI
  MPI_Finalize();
#endif
	
	
}


