/*****************************************************************************/
/* Copyright 2002, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Sample driver for Poisson equation AMG smoothed aggregation solver in the */
/* ML package using Aztec. This test is on a square using a 2-dimensional    */
/* uniform grid with Dirichlet boundary conditions.                          */
/* Note: Dirichlet nodes are not stored in the matrix.                       */
/* The node numbering in this example is lexicographical. That is,           */
/*                                                                           */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*     ----------(12)----------(13)----------(14)----------(15)----------    */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*     ----------( 8)----------( 9)----------(10)----------(11)----------    */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*     ----------( 4)----------( 5)----------( 6)----------( 7)----------    */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*     ----------( 0)----------( 1)----------( 2)----------( 3)----------    */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                 |             |             |             |               */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"
#include "ml_viz_stats.h"
#include "ml_agg_info.h"

/*****************************************************************************/
/* All functions/structures starting with the word 'user' denote things that */
/* are not in ML and are specific to this particular example.                */
/*                                                                           */
/* User defined structure holding information on how the PDE is partitioned  */
/* over the processor system.                                                */
/*****************************************************************************/
struct user_partition_data {               
  int *my_global_ids;      /* my_global_ids[i]: id of ith local unknown.     */
  int *needed_external_ids;/* global ids of ghost unknowns.                  */
  int Nlocal;              /* Number of local unknowns.                      */
  int Nglobal;             /* Number of global unknowns.                     */
  int *my_local_ids;       /* my_local_ids[i]: local id of ith global unknown*/
  int mypid;               /* processor id                                   */
  int nprocs;              /* total number of processors.                    */
  int Nghost;              /* number of ghost variables on processor.        */
};

/*****************************************************************************/
/* Function definitions.                                                     */
/*****************************************************************************/
extern void        user_partition(struct user_partition_data *Partition);
extern AZ_MATRIX   *user_Kn_build(struct user_partition_data *,
				  double **, double **);

#ifdef ML_MPI
#define COMMUNICATOR   MPI_COMM_WORLD
#else
#define COMMUNICATOR   AZ_NOT_MPI
#endif


int main(int argc, char *argv[])
{
  int    Nnodes=128*128;          /* Total number of nodes in the problem.*/
                                    /* 'Nnodes' must be a perfect square.   */
  int    MaxMgLevels=6;             /* Maximum number of Multigrid Levels   */
  double tolerance = 1.0e-8;        /* At convergence:                      */
                                    /*   ||r_k||_2 < tolerance ||r_0||_2    */

  ML           *ml;          
  ML_Aggregate *ag;
  int          Nlevels;
  int          level, coarsest_level;
  double       *rhs, *xxx;
  struct       user_partition_data Partition = {NULL, NULL,0,0,NULL,0,0,0};

  /* See Aztec User's Guide for information on these variables */

  AZ_MATRIX    *Kn_mat;
  AZ_PRECOND   *Pmat = NULL;
  int          proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double       params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];
  double       *x = NULL, *y = NULL;
  

  /* get processor information (id & # of procs) and set ML's printlevel. */

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
  AZ_set_proc_config(proc_config, COMMUNICATOR);
  
  ML_Set_PrintLevel(10);   /* set ML's output level: 0 gives least output */

  /* Set the # of global nodes and partition over the processors.         */

  Partition.Nglobal = Nnodes;
  user_partition(&Partition);

  /* Create an empty multigrid hierarchy and set the 'MaxMGLevels-1'th      */
  /* level discretization within this hierarchy to the ML matrix            */
  /* representing Kn (Poisson discretization).                              */

  ML_Create(&ml, MaxMgLevels);

  /* Build Kn as Aztec matrices. Use built-in function AZ_ML_Set_Amat()    */
  /* to convert to an ML matrix and put in hierarchy.                      */

  Kn_mat = user_Kn_build( &Partition,&x,&y);
  AZ_ML_Set_Amat(ml, MaxMgLevels-1, Partition.Nlocal, 
		 Partition.Nlocal, Kn_mat, proc_config);
 

  /********************************************************************/
  /* Set some ML parameters.                                          */
  /* Some details:                                                    */
  /* - ML_Aggregate_VizAndStats_Setup tells ML to keep trace of some  */
  /*   information about the aggregates, so that it will be possible  */
  /*   to visualize them (after the construction of the ML hierarchy) */
  /*   with ML_Aggregate_VizAndStats_Compute (this required the nodal */
  /*   coordinates). Memory used for visualization is free using      */
  /*   ML_Aggregate_VizAndStats_Destoy.                               */
  /* - ML_Aggregate_Set_GlobalNumber defined the global number of     */
  /*   aggregates. Other options are available as well. For instance, */
  /*   one can define the local number of agggregate, as              */
  /*      ML_Aggregate_Set_LocalNumber( ml, ag, -1, 128 );            */
  /*   (-1 means for all levels), or fix the number of nodes per each */
  /*   aggregte, as                                                   */
  /*      ML_Aggregate_Set_NodesPerAggr( ml, ag, MaxMgLevels-2, 100 );*/
  /*   With METIS aggregation, one can use reordering for the next    */
  /* level:                                                           */
  /*      ML_Aggregate_Set_ReorderingFlag( ml, ag,-1, ML_YES );       */
  /*------------------------------------------------------------------*/

  ML_Aggregate_Create( &ag );  

  ML_Aggregate_VizAndStats_Setup( ag, MaxMgLevels );
  ML_Aggregate_Set_CoarsenScheme_METIS(ag);
  ML_Aggregate_Set_MaxCoarseSize(ag, 30);
  ML_Aggregate_Set_Threshold(ag, 0.0);
  
  /* example of setting global number of aggregates */
  ML_Aggregate_Set_GlobalNumber( ml, ag, MaxMgLevels-1, 128 );
  ML_Aggregate_Set_GlobalNumber( ml, ag, MaxMgLevels-2, 16);
  ML_Aggregate_Set_GlobalNumber( ml, ag, MaxMgLevels-3, 4 );

  /********************************************************************/
  /* Build hierarchy using smoothed aggregation.                      */
  /*------------------------------------------------------------------*/

  Nlevels = ML_Gen_MGHierarchy_UsingAggregation(ml, MaxMgLevels-1,
						ML_DECREASING, ag);
  
  coarsest_level = MaxMgLevels - Nlevels;

  /* ********************************************************************** */
  /* Set up smoother for all levels (here, MLS polynomial smoother)         */
  /* ********************************************************************** */


  for (level = MaxMgLevels-1; level > coarsest_level; level--)
    ML_Gen_Smoother_MLS(ml, level, ML_BOTH, 30., 3);

#ifdef HAVE_ML_AMESOS
  ML_Gen_Smoother_Amesos(ml,coarsest_level, ML_AMESOS_KLU, -1);
#else
  ML_Gen_CoarseSolverSuperLU( ml, coarsest_level);
#endif

  /* Must be called before invoking the preconditioner */
  ML_Gen_Solver(ml, ML_MGV, MaxMgLevels-1, coarsest_level); 

  /* Set initial guess and right hand side. */

  xxx = (double *) ML_allocate((Partition.Nlocal+
				Partition.Nghost)*sizeof(double)); 
  rhs = (double *) ML_allocate(Partition.Nlocal*sizeof(double)); 
  ML_random_vec(xxx, Partition.Nlocal, ml->comm);
  ML_random_vec(rhs, Partition.Nlocal, ml->comm);

  /* Invoke solver */

  AZ_defaults(options, params);
  options[AZ_solver]   = AZ_cg;
  params[AZ_tol]       = tolerance;
  options[AZ_max_iter] = 1550;
  options[AZ_output] = 32;
  
  AZ_set_ML_preconditioner(&Pmat, Kn_mat, ml, options); 
  AZ_iterate(xxx, rhs, options, params, status, proc_config, Kn_mat, Pmat, NULL);
  
  ML_Aggregate_VizAndStats_Compute( ml, ag, MaxMgLevels,
				    x, y, NULL, 2, NULL);

  ML_Aggregate_VizAndStats_Clean( ag, MaxMgLevels );
  
  /* clean up. */
  
  if (Partition.my_local_ids != NULL) free(Partition.my_local_ids);
  if (Partition.my_global_ids != NULL) free(Partition.my_global_ids);
  if (Partition.needed_external_ids != NULL) 
    free(Partition.needed_external_ids);

  ML_Aggregate_Destroy(&ag);
  ML_Destroy(&ml);
  if (Pmat  != NULL) AZ_precond_destroy(&Pmat);
  if (Kn_mat != NULL) {
    AZ_free(Kn_mat->bindx);
    AZ_free(Kn_mat->val);
    AZ_free(Kn_mat->data_org);
    AZ_matrix_destroy(&Kn_mat);
  }
  free(xxx);
  free(rhs);

  if( x != NULL ) ML_free( x );
  if( y != NULL ) ML_free( y );
  
#ifdef ML_MPI
  MPI_Finalize();
#endif
  return 0;
}

/* Assign unknowns to processors */
void user_partition(struct user_partition_data *Partition)
{
  int    proc_config[AZ_PROC_SIZE];

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_input_update(NULL,&(Partition->Nlocal), &(Partition->my_global_ids),
		  proc_config, Partition->Nglobal, 1, AZ_linear);
  Partition->Nghost = 0;   /* will be computed later */
}

/********************************************************************/
/* Set up Kn_mat                                                    */
/*  1) First set it up as an Aztec DMSR matrix                      */
/*     This stores the diagonal of row j in Ke_val[j] and stores a  */
/*     pointer to row j's off-diagonal nonzeros in Ke_bindx[j] with */
/*     column/nonzero entries in Ke_bindx[Ke_bindx[j]:Ke_bindx[j+1]-1] */
/*     and Ke_val[Ke_bindx[j]:Ke_bindx[j+1]-1].                     */
/*  2) call AZ_transform to convert global column indices to local  */
/*     indices and to set up Aztec's communication structure.       */
/*  3) Stuff the arrays into an Aztec matrix.                       */
/*------------------------------------------------------------------*/

AZ_MATRIX *user_Kn_build(struct user_partition_data *Partition,
			 double **x, double **y )

{
  int *Kn_bindx;
  double *Kn_val;
  int    proc_config[AZ_PROC_SIZE];
  AZ_MATRIX *Kn_mat;
  int    *reordered_glob = NULL, *cpntr = NULL, *Kn_data_org = NULL;
  int i, ii, jj, nx, gid, Nlocal, nz_ptr;
  int *reordered_externs = NULL;  /* Aztec thing */
  size_t size;

  Nlocal = Partition->Nlocal;
  Kn_bindx = (int    *) malloc((27*Nlocal+5)*sizeof(int));
  Kn_val   = (double *) malloc((27*Nlocal+5)*sizeof(double));
  Kn_bindx[0] = Nlocal+1;

  nx = (int) sqrt( ((double) Partition->Nglobal) + .00001);

  for (i = 0; i < Nlocal; i++) {
    gid = (Partition->my_global_ids)[i];

    nz_ptr = Kn_bindx[i];
    ii = gid%nx;
    jj = (gid - ii)/nx;


    if (ii != nx-1) { Kn_bindx[nz_ptr] = gid+ 1; Kn_val[nz_ptr++] = -1.;}
    if (jj != nx-1) { Kn_bindx[nz_ptr] = gid+nx; Kn_val[nz_ptr++] = -1.;}
    if (jj !=    0) { Kn_bindx[nz_ptr] = gid-nx; Kn_val[nz_ptr++] = -1.;}
    if (ii !=    0) { Kn_bindx[nz_ptr] = gid- 1; Kn_val[nz_ptr++] = -1.;}

    Kn_val[i] = 4.;
    Kn_bindx[i+1] = nz_ptr;
  }

  /* Transform the global Aztec matrix into a local Aztec matrix. That is,   */
  /* replace global column indices by local indices and set up communication */
  /* data structure 'Ke_data_org' that will be used for matvec's.            */

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_transform(proc_config,&(Partition->needed_external_ids),
			       Kn_bindx, Kn_val, Partition->my_global_ids,
			       &reordered_glob, &reordered_externs, 
			       &Kn_data_org, Nlocal, 0, 0, 0, 
			       &cpntr, AZ_MSR_MATRIX);
  Partition->Nghost = Kn_data_org[AZ_N_external];

  /* Convert old style Aztec matrix to newer style Aztec matrix */

  Kn_mat = AZ_matrix_create( Nlocal );
  AZ_set_MSR(Kn_mat, Kn_bindx, Kn_val, Kn_data_org, 0, NULL, AZ_LOCAL);

  size = sizeof(double)*(Nlocal+Partition->Nghost);

  *x = (double *) ML_allocate( size );
  *y = (double *) ML_allocate( size );
  
  AZ_ML_Build_NodalCoordinates( Partition->Nglobal,
			  Partition->Nlocal, Kn_data_org[AZ_N_external],
			  Partition->my_global_ids,
			  Partition->needed_external_ids,
			  reordered_glob, reordered_externs,
			  *x, *y, NULL, 2 );
  AZ_free(reordered_glob);
  AZ_free(reordered_externs);
  
  return(Kn_mat);
} 
