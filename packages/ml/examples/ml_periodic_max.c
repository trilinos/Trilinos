/*****************************************************************************/
/* Copyright 2002, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Sample driver for Maxwell equation AMG solver in the ML package. The      */
/* software is tested by setting up a 2-dimensional uniform grid example on  */
/* a square. Two different ways of creating matrices are shown:              */
/*     -DAZTEC : Matrices are first set up as Aztec matrices and then        */
/*               converted to ML matrices.                                   */
/*     default : Matrix-free matvec() and getrow() functions are created to  */
/*               represent Ke, Kn, Tmat. These are then used to define ML    */
/*               matrices.                                                   */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define AxTEC
#ifdef AZTEC
#undef AZTEC
#endif
#ifdef AZTEC
#include "az_aztec.h"
#endif
#include "ml_include.h"

#define EDGE 0
#define NODE 1
/*****************************************************************************/
/* User defined structure holding information on how the PDE is partitioned  */
/* over the processor system.                                                */
/*****************************************************************************/
struct user_partition {               
  int *my_global_ids;      /* my_global_ids[i]: id of ith local unknown.     */
  int *needed_external_ids;/* global ids of ghost unknowns.                  */
  int Nlocal;              /* Number of local unknowns.                      */
  int Nglobal;             /* Number of global unknowns.                     */
  int *my_local_ids;       /* my_local_ids[i]: local id of ith global unknown*/
  int type;                /* EDGE or NODE                                   */
  int mypid;               /* processor id                                   */
  int nprocs;              /* total number of processors.                    */
  int Nghost;
};

/******************************************************************************/
/* User defined structure for performing a matrix-vector product and for      */
/* getting a row of the T matrix (null space).                                */
/******************************************************************************/
struct Tmat_data {
  struct user_partition *edge;
  struct user_partition *node;
  ML_Operator *Kn;
};


/******************************************************************************/
/* Function definitions.                                                      */
/******************************************************************************/
extern void partition_edges(struct user_partition *Partition);
extern void partition_nodes(struct user_partition *Partition);
extern int update_ghost_edges(double vector[], void *data);
extern int update_ghost_nodes(double vector[], void *data);
#ifdef AZTEC
extern AZ_MATRIX   *user_Ke_build(struct user_partition *);
extern AZ_MATRIX   *user_Kn_build(struct user_partition *);
extern ML_Operator *user_T_build (struct user_partition *, struct user_partition *, 
			      ML_Operator *);
#else
extern int Ke_matvec(void *, int, double *, int , double *);
extern int Ke_getrow(void *, int, int *, int, int *, double *, int *);
extern int Kn_getrow(void *, int, int *, int, int *, double *, int *);
extern int Tmat_getrow(void *, int, int *, int, int *, double *, int *);
extern int Tmat_matvec(void *, int, double *, int, double *);
#endif


#ifdef ML_MPI
#define COMMUNICATOR   MPI_COMM_WORLD
#else
#define COMMUNICATOR   AZ_NOT_MPI
#endif


int main(int argc, char *argv[])
{
  int    Nnodes=4*4;              /* Total number of nodes in the problem.*/
                                    /* 'Nnodes' must be a perfect square.   */
  int    MaxMgLevels=6;             /* Maximum number of Multigrid Levels   */
  int    Nits_per_presmooth=1;      /* # of pre & post smoothings per level */
  double tolerance = 1.0e-8;        /* At convergence:                      */
                                    /*   ||r_k||_2 < tolerance ||r_0||_2    */
  int smoothPe_flag = ML_YES;       /* ML_YES: smooth tentative prolongator */
                                    /* ML_NO: don't smooth prolongator      */

  /***************************************************************************/
  /* Select Hiptmair relaxation subsmoothers for the nodal and edge problems */
  /* Choices include                                                         */
  /*   1) ML_Gen_Smoother_SymGaussSeidel: this corresponds to a processor    */
  /*      local version of symmetric Gauss-Seidel/SOR. The number of sweeps  */
  /*      can be set via either 'edge_its' or 'nodal_its'. The damping can   */
  /*      be set via 'edge_omega' or 'nodal_omega'. When set to ML_DDEFAULT, */
  /*      the damping is set to '1' on one processor. On multiple processors */
  /*      a lower damping value is set. This is needed to converge processor */
  /*      local SOR.                                                         */
  /*   2) ML_Gen_Smoother_MLS: this corresponds to polynomial relaxation.    */
  /*      The degree of the polynomial is set via 'edge_its' or 'nodal_its'. */
  /*      If the degree is '-1', Marian Brezina's MLS polynomial is chosen.  */
  /*      Otherwise, a Chebyshev polynomial is used over high frequencies    */
  /*      [ lambda_max/alpha , lambda_max]. Lambda_max is computed. 'alpha'  */
  /*      is hardwired in this example to correspond to twice the ratio of   */
  /*      unknowns in the fine and coarse meshes.                            */
  /*                                                                         */
  /* Using 'hiptmair_type' (see comments below) it is also possible to choose*/
  /* when edge and nodal problems are relaxed within the Hiptmair smoother.  */
  /***************************************************************************/

  void  *edge_smoother=(void *)     /* Edge relaxation:                     */
               ML_Gen_Smoother_MLS; /*     ML_Gen_Smoother_MLS              */
                                    /*     ML_Gen_Smoother_SymGaussSeidel   */
  void *nodal_smoother=(void *)     /* Nodal relaxation                     */
               ML_Gen_Smoother_MLS; /*     ML_Gen_Smoother_MLS              */
                                    /*     ML_Gen_Smoother_SymGaussSeidel   */

  int  edge_its = 3;                /* Iterations or polynomial degree for  */
  int  nodal_its = 3;               /* edge/nodal subsmoothers.             */
  double nodal_omega = ML_DDEFAULT, /* SOR damping parameter for noda/edge  */
         edge_omega  = ML_DDEFAULT; /* subsmoothers (see comments above).   */
  int   hiptmair_type=HALF_HIPTMAIR;/* FULL_HIPTMAIR: each invokation       */
                                    /*     smoothes on edges, then nodes,   */
                                    /*     and then once again on edges.    */
                                    /* HALF_HIPTMAIR: each pre-invokation   */
                                    /*     smoothes on edges, then nodes.   */
                                    /*     Each post-invokation smoothes    */
                                    /*     on nodes then edges. .           */


  ML_Operator  *Tmat, *Tmat_trans, **Tmat_array, **Tmat_trans_array;
  ML           *ml_edges, *ml_nodes;
  ML_Aggregate *ag;
  int          Nfine_edge, Ncoarse_edge, Nfine_node, Ncoarse_node, Nlevels;
  int          level, coarsest_level, itmp;
  double       edge_coarsening_rate, node_coarsening_rate, *rhs, *xxx;
  void         **edge_args, **nodal_args;
  struct       user_partition Edge_Partition = {NULL, NULL,0,0}, 
                                Node_Partition = {NULL, NULL,0,0};
  struct Tmat_data Tmat_data;
int i, Ntotal;
 ML_Comm *comm;

  /* See Aztec User's Guide for information on these variables */

#ifdef AZTEC
  AZ_MATRIX    *Ke_mat, *Kn_mat;
  AZ_PRECOND   *Pmat = NULL;
  int          proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double       params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];
#endif


  /* get processor information (proc id & # of procs) and set ML's printlevel. */

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
#ifdef AZTEC
  AZ_set_proc_config(proc_config, COMMUNICATOR);
#endif
  ML_Set_PrintLevel(10);   /* set ML's output level: 0 gives least output */

  /* Set the # of global nodes/edges and partition both the edges and the */
  /* nodes over the processors. NOTE: I believe we assume that if an edge */
  /* is assigned to a processor at least one of its nodes must be also    */
  /* assigned to that processor.                                          */

  Node_Partition.Nglobal = Nnodes;
  Edge_Partition.Nglobal = Node_Partition.Nglobal*2;
  Node_Partition.type = NODE;
  Edge_Partition.type = EDGE;
#define perxodic
#ifdef periodic
Node_Partition.Nglobal += 2; 
#endif
  partition_edges(&Edge_Partition);
  partition_nodes(&Node_Partition);
xxx = (double *) ML_allocate((Edge_Partition.Nlocal+100)*sizeof(double)); 
rhs = (double *) ML_allocate((Edge_Partition.Nlocal+100)*sizeof(double)); 
 for (i = 0; i < Edge_Partition.Nlocal + 100; i++) xxx[i] = -1.;
 for (i = 0; i < Edge_Partition.Nlocal; i++) xxx[i] = (double) 
        Edge_Partition.my_global_ids[i];
 printf("before update\n"); fflush(stdout);
update_ghost_edges(xxx, (void *) &Edge_Partition);
 printf("after update\n"); fflush(stdout);


  //  exit(1);
  /* Create an empty multigrid hierarchy and set the 'MaxMGLevels-1'th   */
  /* level discretization within this hierarchy to the ML matrix         */
  /* representing Ke (Maxwell edge discretization).                      */

  ML_Create(&ml_edges, MaxMgLevels);
#ifdef AZTEC
  /* Build Ke as an Aztec matrix. Use built-in function AZ_ML_Set_Amat() */
  /* to convert to an ML matrix and put in hierarchy.                    */

  Ke_mat = user_Ke_build(&Edge_Partition);
  AZ_ML_Set_Amat(ml_edges, MaxMgLevels-1, Edge_Partition.Nlocal,
      		 Edge_Partition.Nlocal, Ke_mat, proc_config);
#else
  /* Build Ke directly as an ML matrix.                                  */

  ML_Init_Amatrix      (ml_edges, MaxMgLevels-1, Edge_Partition.Nlocal,
			Edge_Partition.Nlocal, &Edge_Partition);
Ke_matvec(&Edge_Partition, Edge_Partition.Nlocal, xxx, Edge_Partition.Nlocal, rhs);


  printf("before getrow\n"); fflush(stdout);
  Ntotal = Edge_Partition.Nlocal;
  if (Edge_Partition.nprocs == 2) Ntotal += Edge_Partition.Nghost;
  ML_Set_Amatrix_Getrow(ml_edges, MaxMgLevels-1,  Ke_getrow, update_ghost_edges, Ntotal);
  printf("after getrow\n"); fflush(stdout);
  ML_Set_Amatrix_Matvec(ml_edges, MaxMgLevels-1,  Ke_matvec);
  ML_Operator_Check_Getrow(
		    &(ml_edges->Amat[MaxMgLevels-1]), 7,"Amat");
printf("end of program\n"); fflush(stdout);
#endif

  //  ML_Operator_Print( &(ml_edges->Amat[MaxMgLevels-1]),"Ke");


  /* Build an Aztec matrix representing an auxiliary nodal PDE problem.  */
  /* This should be a variable coefficient Poisson problem (with unknowns*/
  /* at the nodes). The coefficients should be chosen to reflect the     */
  /* conductivity of the original edge problems.                         */
  /* Create an empty multigrid hierarchy. Convert the Aztec matrix to an */
  /* ML matrix and put it in the 'MaxMGLevels-1' level of the hierarchy. */
  /* Note it is possible to multiply T'*T for get this matrix though this*/
  /* will not incorporate material properties.                           */

  ML_Create(&ml_nodes, MaxMgLevels);

#ifdef AZTEC
  Kn_mat = user_Kn_build( &Node_Partition);
  AZ_ML_Set_Amat(ml_nodes, MaxMgLevels-1, Node_Partition.Nlocal, 
		 Node_Partition.Nlocal, Kn_mat, proc_config);
#else
  ML_Init_Amatrix      (ml_nodes, MaxMgLevels-1 , Node_Partition.Nlocal,
			Node_Partition.Nlocal, &Node_Partition);
  Ntotal = Node_Partition.Nlocal;
  if (Node_Partition.nprocs == 2) Ntotal += Node_Partition.Nghost;
  ML_Set_Amatrix_Getrow(ml_nodes, MaxMgLevels-1,  Kn_getrow, update_ghost_nodes, Ntotal);
#endif
  //  ML_Operator_Print( &(ml_nodes->Amat[MaxMgLevels-1]),"Kn");

  /* Build an ML matrix representing the null space of the PDE problem. */
  /* This should be a discrete gradient (nodes to edges).               */

#ifdef AZTEC
    Tmat = user_T_build (&Edge_Partition, &Node_Partition, 
  		   &(ml_nodes->Amat[MaxMgLevels-1]));
#else
    Tmat = ML_Operator_Create(ml_nodes->comm);
    Tmat_data.edge = &Edge_Partition;
    Tmat_data.node = &Node_Partition;
    Tmat_data.Kn   = &(ml_nodes->Amat[MaxMgLevels-1]);

    ML_Operator_Set_ApplyFuncData( Tmat,	Node_Partition.Nlocal,
				   Edge_Partition.Nlocal, ML_EMPTY, (void *) &Tmat_data, 
				   Edge_Partition.Nlocal, NULL, 0);
    //    ML_Operator_Set_ApplyFuncData( Tmat,	Node_Partition.Nlocal,
    //				   Edge_Partition.Nlocal, ML_EXTERNAL, (void *) &Tmat_data, 
    //				   Edge_Partition.Nlocal, Tmat_getrow, 0);
    ML_Operator_Set_Getrow( Tmat, ML_EXTERNAL, Edge_Partition.Nlocal,Tmat_getrow);
    ML_Operator_Set_ApplyFunc(Tmat, ML_EXTERNAL, Tmat_matvec);
  ML_Comm_Create( &comm);
  printf("before comminfo generate %d\n",
	 Node_Partition.Nghost); fflush(stdout);


 ML_CommInfoOP_Generate( &(Tmat->getrow->pre_comm), update_ghost_nodes, 
&Node_Partition,
                              comm, Tmat->invec_leng, 
Node_Partition.Nghost);

#endif
    //    ML_Operator_Print( Tmat, "Tmat");

 for (i = 0; i < Node_Partition.Nlocal; i++) xxx[i] = (double) 
        Node_Partition.my_global_ids[i];
 ML_Operator_Apply(Tmat, Tmat->invec_leng, xxx, Tmat->outvec_leng, rhs);
for (i = 0; i < Edge_Partition.Nlocal ; i++) 
  printf("%d: %d %e \n", Edge_Partition.mypid, Edge_Partition.my_global_ids[i], rhs[i]);

 fflush(stdout);


  /********************************************************************/
  /* Set some ML parameters.                                          */
  /*------------------------------------------------------------------*/
	
  ML_Set_ResidualOutputFrequency(ml_edges, 1);
  ML_Set_Tolerance(ml_edges, 1.0e-8);
  ML_Aggregate_Create( &ag );
  ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag);
  ML_Aggregate_Set_DampingFactor(ag, 0.0); /* must use 0 for maxwell */
  ML_Aggregate_Set_MaxCoarseSize(ag, 30);
  ML_Aggregate_Set_Threshold(ag, 0.0);


  /********************************************************************/
  /*                      Set up Tmat_trans                           */
  /*------------------------------------------------------------------*/

  Tmat_trans = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Transpose_byrow(Tmat, Tmat_trans);


  Nlevels=ML_Gen_MGHierarchy_UsingReitzinger(ml_edges, &ml_nodes,MaxMgLevels-1,
					     ML_DECREASING,ag,Tmat,Tmat_trans, 
					     &Tmat_array,&Tmat_trans_array, 
					     smoothPe_flag, 1.5);

  /* Set the Hiptmair subsmoothers */

  if (nodal_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    ML_Smoother_Arglist_Set(nodal_args, 1, &nodal_omega);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    ML_Smoother_Arglist_Set(edge_args, 1, &edge_omega);
  }
  if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    Nfine_node = Tmat_array[MaxMgLevels-1]->invec_leng;
    ML_gsum_vec_int(&Nfine_node, &itmp, 1, ml_edges->comm);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    Nfine_edge = Tmat_array[MaxMgLevels-1]->outvec_leng;
    ML_gsum_vec_int(&Nfine_edge, &itmp, 1, ml_edges->comm);
  }

  /****************************************************
  * Set up smoothers for all levels but the coarsest. *
  ****************************************************/
  coarsest_level = MaxMgLevels - Nlevels;

  for (level = MaxMgLevels-1; level > coarsest_level; level--)
    {
      if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
	Ncoarse_edge = Tmat_array[level-1]->outvec_leng;
	ML_gsum_vec_int(&Ncoarse_edge, &itmp, 1, ml_edges->comm);
	edge_coarsening_rate =  2.*((double) Nfine_edge)/ ((double) Ncoarse_edge);
	ML_Smoother_Arglist_Set(edge_args, 1, &edge_coarsening_rate);
	Nfine_edge = Ncoarse_edge;
      }
      if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
	Ncoarse_node = Tmat_array[level-1]->invec_leng;
	ML_gsum_vec_int(&Ncoarse_node, &itmp, 1, ml_edges->comm);
	node_coarsening_rate =  2.*((double) Nfine_node)/ ((double) Ncoarse_node);
	ML_Smoother_Arglist_Set(nodal_args, 1, &node_coarsening_rate);
	Nfine_node = Ncoarse_node;
      }
      ML_Gen_Smoother_Hiptmair(ml_edges, level, ML_BOTH, Nits_per_presmooth,
			       Tmat_array, Tmat_trans_array, NULL, edge_smoother,
			       edge_args, nodal_smoother,nodal_args, hiptmair_type);
    }

  /*******************************************
  * Set up coarsest level smoother
  *******************************************/

  if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
    edge_coarsening_rate = (double) Nfine_edge;
    ML_Smoother_Arglist_Set(edge_args, 1, &edge_coarsening_rate);
  }
  if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
    node_coarsening_rate = (double) Nfine_node;
    ML_Smoother_Arglist_Set(nodal_args,1,&node_coarsening_rate);
  }
  //  ML_Gen_Smoother_Hiptmair(ml_edges, coarsest_level, ML_BOTH, Nits_per_presmooth,
  //  			   Tmat_array, Tmat_trans_array, NULL, edge_smoother,
  //  			   edge_args, nodal_smoother,nodal_args, hiptmair_type);
       ML_Gen_CoarseSolverSuperLU( ml_edges, coarsest_level);
  

  /* Must be called before invoking the preconditioner */
  ML_Gen_Solver(ml_edges, ML_MGV, MaxMgLevels-1, coarsest_level); 



  /* Set the initial guess and the right hand side. Invoke solver */	

  xxx = (double *) ML_allocate(Edge_Partition.Nlocal*sizeof(double)); 
  ML_random_vec(xxx, Edge_Partition.Nlocal, ml_edges->comm);
  rhs = (double *) ML_allocate(Edge_Partition.Nlocal*sizeof(double)); 
  ML_random_vec(rhs, Edge_Partition.Nlocal, ml_edges->comm);

#ifdef AZTEC
  /* Choose the Aztec solver and criteria. Also tell Aztec that */
  /* ML will be supplying the preconditioner.                   */

  AZ_defaults(options, params);
  options[AZ_solver]   = AZ_fixed_pt;
  options[AZ_solver]   = AZ_gmres;
  options[AZ_kspace]   = 80;
  params[AZ_tol]       = tolerance;
  AZ_set_ML_preconditioner(&Pmat, Ke_mat, ml_edges, options); 
  options[AZ_conv] = AZ_noscaled;
  AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, NULL);
  AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, NULL);
#else
  ML_Iterate(ml_edges, xxx, rhs);
#endif


  /* clean up. */

  ML_Smoother_Arglist_Delete(&nodal_args);
  ML_Smoother_Arglist_Delete(&edge_args);
  ML_Aggregate_Destroy(&ag);
  ML_Destroy(&ml_edges);
  ML_Destroy(&ml_nodes);
#ifdef AZTEC
  AZ_free((void *) Ke_mat->data_org);
  AZ_free((void *) Ke_mat->val);
  AZ_free((void *) Ke_mat->bindx);
  if (Ke_mat  != NULL) AZ_matrix_destroy(&Ke_mat);
  if (Pmat  != NULL) AZ_precond_destroy(&Pmat);
  if (Kn_mat != NULL) AZ_matrix_destroy(&Kn_mat);
#endif
  free(xxx);
  free(rhs);
  ML_Operator_Destroy(Tmat);
  ML_Operator_Destroy(Tmat_trans);
  ML_MGHierarchy_ReitzingerDestroy(MaxMgLevels-2, coarsest_level, &Tmat_array, &Tmat_trans_array);

#ifdef ML_MPI
  MPI_Finalize();
#endif
		
  return 0;
		
}

/* Example specific stuff to setup Ke, Kn, T                */
/* This example corresponds to periodic boundary conditions */
/* The numbering I think is like this                       */
/*                                                          */
/*               |----6----|----7----|----8----|            */
/*               |         |         |         |            */
/*               12        13        14        12           */
/*               |         |         |         |            */
/*               |----3----|----4----|----5----|            */
/*               |         |         |         |            */
/*               9         10        11        9            */
/*               |         |         |         |            */
/*               |----0----|----1----|----2----|            */
/*               |         |         |         |            */
/*               15        16        17        15           */
/*               |         |         |         |            */
/*               |----6----|----7----|----8----|            */
/*                                                          */
#define HORIZONTAL 1
#define VERTICAL   2
extern int inv2dindex(int index, int *i, int *j, int n, int *hor_or_vert);
extern int inv2dnodeindex(int index, int *i, int *j, int n);
extern int northeast2d(int i, int j, int n);
extern int northwest2d(int i, int j, int n);
extern int southeast2d(int i, int j, int n);
extern int southwest2d(int i, int j, int n);
extern int north2d(int i, int j, int n);
extern int south2d(int i, int j, int n);
extern int east2d(int i, int j, int n);
extern int west2d(int i, int j, int n);
extern void compress_matrix(double val[], int bindx[], int N_update);

/********************************************************************/
/* Given the (i,j)th cell return the node indices defining the cell */
/*                                                                  */
/*                northwest . ------- . northeast                   */
/*                          |         |                             */
/*                          |         |                             */
/*                          |         |                             */
/*                southeast . ------- . southeast                   */
/********************************************************************/
int northwest2d(int i, int j, int n) { return(i + (j)*n); }
int southwest2d(int i, int j, int n) { if (j == 0) j = n;    /* periodic */
                                       return(i + (j-1)*n);}
int northeast2d(int i, int j, int n) { if (i == n-1) i = -1; /* periodic */
                                       return(i+1 + (j)*n);}
int southeast2d(int i, int j, int n) { if (j == 0) j = n;    /* periodic */
                                       if (i == n-1) i = -1; /* periodic */
                                       return(i+1 + (j-1)*n);}

/********************************************************************/
/* Given the (i,j)th cell return the edge indices defining the cell */
/*                                                                  */
/*                          .--north--.                             */
/*                          |         |                             */
/*                        west       east                           */
/*                          |         |                             */
/*                          .--south--.                             */
/*                                                                  */
/********************************************************************/

int north2d(int i, int j, int n) { return(i + j*n);}
int south2d(int i, int j, int n) { if (j == 0) j = n;       /* periodic */
                                   return(i + (j-1)*n);}
int west2d(int i, int j, int n)  { if (j == 0) j = n;       /* periodic */
                                   return(i+(j-1)*n+n*n);}
int east2d(int i, int j, int n)  { if (j == 0) j = n;       /* periodic */
                                   if (i == n-1) i = -1;    /* periodic */
                                   return(i+1+(j-1)*n+n*n);}

/* Given the index of either a 'south' or 'west' edge, return the */
/* cell coordinates (i,j) as well as the orientation (horizontal  */
/* or vertical) of the edge.                                      */
int inv2dindex(int index, int *i, int *j, int n, int *hor_or_vert)
{
  *hor_or_vert = HORIZONTAL;
  if (index >= n*n) {
    *hor_or_vert = VERTICAL;
    index -= n*n;
  }
  *i = (index%n);
  *j = (index - (*i))/n;
  (*j) =  ((*j)+1)%n;
  return 1;
}
/* Given the index of either a 'southwest' node return the */
/* cell coordinates (i,j).                                 */
int inv2dnodeindex(int index, int *i, int *j, int n)
{
  *i = (index%n);
  *j = (index - (*i))/n;
  (*j) =  ((*j)+1)%n;

  return 1;
}


#ifdef AZTEC
/********************************************************************/
/* Set up Ke_mat                                                    */
/*  1) First set it up as an Aztec DMSR matrix                      */
/*     This stores the diagonal of row j in Ke_val[j] and stores a  */
/*     pointer to row j's off-diagonal nonzeros in Ke_bindx[j] with */
/*     column/nonzero entries in Ke_bindx[Ke_bindx[j]:Ke_bindx[j+1]-1] */
/*     and Ke_val[Ke_bindx[j]:Ke_bindx[j+1]-1].                     */
/*  2) call AZ_transform to convert global column indices to local  */
/*     indices and to set up Aztec's communication structure.       */
/*  3) Stuff the arrays into an Aztec matrix.                       */
/*------------------------------------------------------------------*/

AZ_MATRIX *user_Ke_build(struct user_partition *Edge_Partition)
{
  double dcenter, doffdiag, sigma = .0001;
  int ii,jj, horv, i, nx, global_id, nz_ptr, Nlocal_edges;

  /* Aztec matrix and temp variables */

  int       *Ke_bindx, *Ke_data_org = NULL;
  double    *Ke_val;
  AZ_MATRIX *Ke_mat;
  int       proc_config[AZ_PROC_SIZE], *cpntr = NULL;
  int       *reordered_glob_edges = NULL, *reordered_edge_externs = NULL;

  Nlocal_edges = Edge_Partition->Nlocal;
  nx = (int) sqrt( ((double) Edge_Partition->Nglobal/2) + .00001);

  Ke_bindx = (int    *) malloc((7*Nlocal_edges+1)*sizeof(int));
  Ke_val   = (double *) malloc((7*Nlocal_edges+1)*sizeof(double));
  Ke_bindx[0] = Nlocal_edges+1;

  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx));
  doffdiag = -1 + sigma/((double) ( 6 * nx * nx));

  /* Create a DMSR matrix with global column indices */

  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    Ke_val[i] = dcenter;
    Ke_bindx[i+1] = Ke_bindx[i] + 6;
    inv2dindex(global_id, &ii, &jj, nx, &horv);
    nz_ptr = Ke_bindx[global_id];
    if (horv == HORIZONTAL) {
      Ke_bindx[nz_ptr] = north2d(ii,jj,nx);     Ke_val[nz_ptr++] = doffdiag;
      Ke_bindx[nz_ptr] = west2d(ii,jj,nx);      Ke_val[nz_ptr++] = 1.;
      Ke_bindx[nz_ptr] = east2d(ii,jj,nx);      Ke_val[nz_ptr++] = -1.;
      if (jj == 0) jj = nx-1;
      else jj--;
      Ke_bindx[nz_ptr] = west2d(ii,jj,nx);      Ke_val[nz_ptr++] = -1.;
      Ke_bindx[nz_ptr] = south2d(ii,jj,nx);     Ke_val[nz_ptr++] = doffdiag;
      Ke_bindx[nz_ptr] = east2d(ii,jj,nx);      Ke_val[nz_ptr++] = 1.;
    }
    else {
      Ke_bindx[nz_ptr] = north2d(ii,jj,nx);     Ke_val[nz_ptr++] = -1.;
      Ke_bindx[nz_ptr] = east2d(ii,jj,nx);      Ke_val[nz_ptr++] = doffdiag;
      Ke_bindx[nz_ptr] = south2d(ii,jj,nx);     Ke_val[nz_ptr++] = 1.;
      if (ii == 0) ii = nx-1;
      else ii--;
      Ke_bindx[nz_ptr] = west2d(ii,jj,nx);     Ke_val[nz_ptr++] = doffdiag;
      Ke_bindx[nz_ptr] = south2d(ii,jj,nx);     Ke_val[nz_ptr++] = -1.;
      Ke_bindx[nz_ptr] = north2d(ii,jj,nx);     Ke_val[nz_ptr++] = 1.;
    }
  }

  /* Transform the global Aztec matrix into a local Aztec matrix. That is,   */
  /* replace global column indices by local indices and set up communication */
  /* data structure 'Ke_data_org' that will be used for matvec's.            */

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_transform_norowreordering(proc_config, &(Edge_Partition->needed_external_ids),
			       Ke_bindx, Ke_val, Edge_Partition->my_global_ids,
			       &reordered_glob_edges, &reordered_edge_externs, 
			       &Ke_data_org, Nlocal_edges, 0, 0, 0, 
			       &cpntr,	       AZ_MSR_MATRIX);

  /* Convert old style Aztec matrix to newer style Aztec matrix */

  Ke_mat = AZ_matrix_create( Nlocal_edges );
  AZ_set_MSR(Ke_mat, Ke_bindx, Ke_val, Ke_data_org, 0, NULL, AZ_LOCAL);

  return(Ke_mat);
}

int *reordered_node_externs = NULL;  /* Aztec thing */

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

AZ_MATRIX *user_Kn_build(struct user_partition *Node_Partition)

{
  int *Kn_bindx;
  double *Kn_val;
  int    proc_config[AZ_PROC_SIZE];
  AZ_MATRIX *Kn_mat;
  int    *reordered_glob_nodes = NULL, *cpntr = NULL, *Kn_data_org = NULL;
  int i, ii, jj, nx, global_id, Nlocal_nodes, nz_ptr;


  Nlocal_nodes = Node_Partition->Nlocal;
 Kn_bindx = (int    *) malloc((27*Nlocal_nodes+5)*sizeof(int));
  Kn_val   = (double *) malloc((27*Nlocal_nodes+5)*sizeof(double));
  Kn_bindx[0] = Nlocal_nodes+1;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);
#ifdef periodic
  nx = (int) sqrt( ((double) Node_Partition->Nglobal - 2) + .00001); //periodic stuff
Nlocal_nodes -= 2; // periodic stuff
#endif

  for (i = 0; i < Nlocal_nodes; i++) {
    global_id = (Node_Partition->my_global_ids)[i];
    Kn_bindx[i+1] = Kn_bindx[i] + 8;
    Kn_val[i] = 8.;
    inv2dnodeindex(global_id, &ii, &jj, nx);
    nz_ptr = Kn_bindx[i];
    Kn_bindx[nz_ptr] = southeast2d(ii,jj,nx);  Kn_val[nz_ptr++] = -1.;
    Kn_bindx[nz_ptr] = northwest2d(ii,jj,nx);  Kn_val[nz_ptr++] = -1.;
    Kn_bindx[nz_ptr] = northeast2d(ii,jj,nx);  Kn_val[nz_ptr++] = -.00000001;
    if (ii == 0) ii = nx-1;
    else ii--;
    Kn_bindx[nz_ptr] = northwest2d(ii,jj,nx);  Kn_val[nz_ptr++] = -.00000001;
    if (jj == 0) jj = nx-1;
    else jj--;
    Kn_bindx[nz_ptr] = southeast2d(ii,jj,nx);  Kn_val[nz_ptr++] = -1.;
    Kn_bindx[nz_ptr] = northwest2d(ii,jj,nx);  Kn_val[nz_ptr++] = -1.;
    Kn_bindx[nz_ptr] = southwest2d(ii,jj,nx);  Kn_val[nz_ptr++] = -.00000001;
    if (ii == nx-1) ii = 0;
    else ii++;
    Kn_bindx[nz_ptr] = southeast2d(ii,jj,nx);  Kn_val[nz_ptr++] = -.00000001;
  }
#ifdef periodic
i = Nlocal_nodes;            //periodic stuff
Kn_bindx[i+1] = Kn_bindx[i]; //periodic stuff
Kn_val[i] = 1.;              //periodic stuff
i++;                         //periodic stuff
Kn_bindx[i+1] = Kn_bindx[i]; //periodic stuff
Kn_val[i] = 1.;              //periodic stuff
Nlocal_nodes += 2;           //periodic stuff
#endif

  /* Transform the global Aztec matrix into a local Aztec matrix. That is,   */
  /* replace global column indices by local indices and set up communication */
  /* data structure 'Ke_data_org' that will be used for matvec's.            */

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_transform_norowreordering(proc_config,&(Node_Partition->needed_external_ids),
			       Kn_bindx, Kn_val, Node_Partition->my_global_ids,
			       &reordered_glob_nodes, &reordered_node_externs, 
			       &Kn_data_org, Nlocal_nodes, 0, 0, 0, 
			       &cpntr, AZ_MSR_MATRIX);

  /* Convert old style Aztec matrix to newer style Aztec matrix */

  Kn_mat = AZ_matrix_create( Nlocal_nodes );
  AZ_set_MSR(Kn_mat, Kn_bindx, Kn_val, Kn_data_org, 0, NULL, AZ_LOCAL);

  return(Kn_mat);
}

/********************************************************************/
/* Set up T_mat                                                     */
/*  1) First set it up as a DMSR matrix. This is a bit dumb because */
/*     a rectangular matrix doesn't have a diagonal. Anyway, DMSR   */
/*     This stores the diagonal of row j in Ke_val[j] and stores a  */
/*     pointer to row j's off-diagonal nonzeros in Ke_bindx[j] with */
/*     column/nonzero entries in Ke_bindx[Ke_bindx[j]:Ke_bindx[j+1]-1] */
/*     and Ke_val[Ke_bindx[j]:Ke_bindx[j+1]-1].                     */
/*  2) Convert the matrix to CSR format.                            */
/*  3) call a modified Aztec routine to convert global columns to   */
/*     local indices and to make an ML matrix out of it.            */
/*  4) Since the above routine does not compute a communication data*/
/*     structure, we clone Knmat's communication structure (i.e. we */
/*     assume that Tmat and Knmat have the same communication.      */
/*------------------------------------------------------------------*/

ML_Operator *user_T_build(struct user_partition *Edge_Partition,
		      struct user_partition *Node_Partition,
		      ML_Operator *Kn_mat)
{

  int nx, i, ii, jj, horv, Ncols, Nexterns;
  int *Tmat_bindx;
  double *Tmat_val;
  ML_Operator *Tmat;
  struct ML_CSR_MSRdata *csr_data;
  ML_Comm *comm;
  struct aztec_context *aztec_context;
  int global_id;
  int Nlocal_nodes, Nlocal_edges; 
  int nz_ptr;

  Nlocal_nodes = Node_Partition->Nlocal;
  Nlocal_edges = Edge_Partition->Nlocal;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);
#ifdef periodic
  nx = (int) sqrt( ((double) Node_Partition->Nglobal - 2) + .00001); //periodic stuff
#endif

  Tmat_bindx = (int    *) malloc((6*Nlocal_edges+5)*sizeof(int));
  Tmat_val   = (double *) malloc((6*Nlocal_edges+5)*sizeof(double));
  // periodic stuff above
  Tmat_bindx[0] = Nlocal_edges + 1;
  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    Tmat_bindx[i+1] = Tmat_bindx[i] + 2;
#ifdef periodic
Tmat_bindx[i+1] += 1; // periodic stuff
#endif
    Tmat_val[i] = 0.0;

    inv2dindex(global_id, &ii, &jj, nx, &horv);
    nz_ptr = Tmat_bindx[i];
    if (horv == HORIZONTAL) {
      Tmat_bindx[nz_ptr] = southwest2d(ii,jj,nx);  Tmat_val[nz_ptr++] = -1.;
      Tmat_bindx[nz_ptr] = southeast2d(ii,jj,nx);  Tmat_val[nz_ptr++] = 1.;
#ifdef periodic
      Tmat_bindx[nz_ptr] = Nlocal_nodes-2;         Tmat_val[nz_ptr++] = 1.;
#endif
    }
    else {
      Tmat_bindx[nz_ptr] = northwest2d(ii,jj,nx);  Tmat_val[nz_ptr++] = -1.;
      Tmat_bindx[nz_ptr] = southwest2d(ii,jj,nx);  Tmat_val[nz_ptr++] =  1.;
#ifdef periodic
      Tmat_bindx[nz_ptr] = Nlocal_nodes - 1;       Tmat_val[nz_ptr++] = 1.;
#endif
    }
  }

  /* Convert the MSR matrix to a CSR matrix. Then use a modified Aztec */
  /* routine to convert the global CSR matrix to a local ML matrix     */
  /* Since this routine does not compute the communication structure,  */
  /* we assume that it is identical to Kn's and just clone it.         */

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  csr_data->columns = Tmat_bindx;
  csr_data->values  = Tmat_val;
  ML_MSR2CSR(csr_data, Nlocal_edges, &Ncols);
  aztec_context = (struct aztec_context *) Kn_mat->data;
  Nexterns = (aztec_context->Amat->data_org)[AZ_N_external];

  Nexterns = 0;

  ML_Comm_Create( &comm);

  AZ_Tmat_transform2ml(Nexterns,   Node_Partition->needed_external_ids,
		       reordered_node_externs,
		       Tmat_bindx, Tmat_val, csr_data->rowptr, Nlocal_nodes,
		       Node_Partition->my_global_ids,
		       comm, Nlocal_edges, &Tmat);
  ML_free(csr_data);
  Tmat->data_destroy = ML_CSR_MSRdata_Destroy;
  ML_CommInfoOP_Clone(&(Tmat->getrow->pre_comm), Kn_mat->getrow->pre_comm);

  return(Tmat);
}
#endif

/* Assign unknowns to processors */
void partition_nodes(struct user_partition *Partition)
{
  int i, nx, Nloc;
  ML_Comm *comm;
  int *my_local_ids;

  Partition->Nghost = 0;
#ifdef AZTEC
  int    proc_config[AZ_PROC_SIZE];

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_input_update(NULL,&(Partition->Nlocal),
		  &(Partition->my_global_ids),
		    proc_config,     Partition->Nglobal, 1, AZ_linear);
#else
 // nonaztec
  nx = (int) sqrt( ((double) Partition->Nglobal) + .00001);

  ML_Comm_Create( &comm);
  Partition->mypid = comm->ML_mypid;
  Partition->nprocs= comm->ML_nprocs;
  Partition->Nghost = 2*nx;

  if (comm->ML_nprocs > 2) {
    if (comm->ML_mypid == 0)
      printf("This example only works for one or two processors\n");
    exit(1);
  }
  if ((comm->ML_nprocs == 2) && (Partition->Nglobal%2 == 1)) {
    if (comm->ML_mypid == 0)
      printf("Two processor case requires an even number of points.\n");
    exit(1);
  }
  if (comm->ML_nprocs == 1) {
    Partition->Nlocal = Partition->Nglobal; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    for (i = 0; i < Partition->Nlocal; i++) 
      (Partition->my_global_ids)[i] = i;
    
    Partition->my_local_ids = (int *) malloc(sizeof(int)*Partition->Nglobal);
    for (i = 0; i < Partition->Nglobal; i++) 
      (Partition->my_local_ids)[i] = i;
  }
  else {
    Partition->Nlocal = Partition->Nglobal/2; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    my_local_ids =  (int *) malloc(sizeof(int)*Partition->Nglobal);
    Partition->my_local_ids = my_local_ids;

    for (i = 0; i < Partition->Nglobal; i++) 
      (Partition->my_local_ids)[i] = -1;

    for (i = 0; i < Partition->Nlocal; i++) {
      (Partition->my_global_ids)[i] = i + comm->ML_mypid*Partition->Nlocal;
      my_local_ids[(Partition->my_global_ids)[i]] = i;
    }

    Nloc = Partition->Nlocal;
    if (comm->ML_mypid == 0) {
      for (i = 0; i < nx; i++) {
	my_local_ids[Nloc + i ] = Nloc + i;
	my_local_ids[2*Nloc - nx + i ] = Nloc + i + nx;
      }
    }
    else {
      for (i = 0; i < nx; i++) {
	my_local_ids[Nloc - nx + i ] = Nloc + i;
	my_local_ids[            i ] = Nloc + i + nx;
      }
    }

  }
#endif
}
/* Assign unknowns to processors */
void partition_edges(struct user_partition *Partition)
{
  int i, nx, Nloc;
  ML_Comm *comm;
  int *my_local_ids;

#ifdef AZTEC
  int    proc_config[AZ_PROC_SIZE];

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_input_update(NULL,&(Partition->Nlocal),
		  &(Partition->my_global_ids),
		    proc_config,     Partition->Nglobal, 1, AZ_linear);
#else
  Partition->Nghost = 0;
 // nonaztec
  nx = (int) sqrt( ((double) Partition->Nglobal/2) + .00001);

  ML_Comm_Create( &comm);
  Partition->mypid = comm->ML_mypid;
  Partition->nprocs= comm->ML_nprocs;

  if (comm->ML_nprocs > 2) {
    if (comm->ML_mypid == 0)
      printf("This example only works for one or two processors\n");
    exit(1);
  }
  if ((comm->ML_nprocs == 2) && (Partition->Nglobal%2 == 1)) {
    if (comm->ML_mypid == 0)
      printf("Two processor case requires an even number of points.\n");
    exit(1);
  }


  if (comm->ML_nprocs == 1) {  /* 1 processor case is simple */
    Partition->Nlocal = Partition->Nglobal; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    for (i = 0; i < Partition->Nlocal; i++) 
      (Partition->my_global_ids)[i] = i;
    
    Partition->my_local_ids = (int *) malloc(sizeof(int)*Partition->Nglobal);
    for (i = 0; i < Partition->Nglobal; i++) 
      (Partition->my_local_ids)[i] = i;
  }
  else {

    /* allocate space */

    Partition->Nlocal = Partition->Nglobal/2; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    my_local_ids =  (int *) malloc(sizeof(int)*Partition->Nglobal);
    Partition->my_local_ids = my_local_ids;

    /* initialize local ids to '-1' (not owned by me) */
    for (i = 0; i < Partition->Nglobal; i++) 
      (Partition->my_local_ids)[i] = -1;

    /* set global ids */
    for (i = 0; i < Partition->Nlocal/2; i++) {
      (Partition->my_global_ids)[i] = i + comm->ML_mypid*Partition->Nlocal/2;
    }
    for (i = 0; i < Partition->Nlocal/2; i++) {
      (Partition->my_global_ids)[i + Partition->Nlocal/2] = 
	Partition->Nlocal+
	i + comm->ML_mypid*Partition->Nlocal/2;
    }

    /* set local ids of nonghost unknowns */

    for (i = 0; i < Partition->Nlocal; i++) 
      my_local_ids[(Partition->my_global_ids)[i]] = i;

    /* set the ghost unknowns */
    Partition->Nghost = 3*nx;

    Nloc = Partition->Nlocal;
    if (comm->ML_mypid == 0) {
      for (i = 0; i < nx; i++) {
	my_local_ids[Nloc/2 + i ] = Nloc + i;
	my_local_ids[Nloc - nx + i ] = Nloc + i + nx;
	my_local_ids[2*Nloc - nx + i ] = Nloc + i + 2*nx;
      }
    }
    else {
      for (i = 0; i < nx; i++) {
	my_local_ids[Nloc/2 - nx + i ] = Nloc + i;
	my_local_ids[              i ] = Nloc + i + nx;
	my_local_ids[3*Nloc/2 - nx + i ] = Nloc + i + 2*nx;
      }
    }

  }

#endif
}


int Ke_matvec(void *data, int Nlocal_edges, double p[], int N_out, double Ap[])
{
  int i, global_id, ii, jj, nx, horv;
  double dcenter, doffdiag, sigma = .0001, *temp;
  int *my_local_ids;
  struct user_partition *Edge_Partition;

  Edge_Partition = (struct user_partition *) data;
  my_local_ids = Edge_Partition->my_local_ids;
  nx = (int) sqrt( ((double) Edge_Partition->Nglobal/2) + .00001);
  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx));
  doffdiag = -1 + sigma/((double) ( 6 * nx * nx));
  temp = (double *) malloc(sizeof(double)*(3*nx + Nlocal_edges));

  for (i = 0; i < Nlocal_edges; i++) temp[i] = p[i];
  update_ghost_edges(temp, (void *) Edge_Partition);

  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    Ap[i] = dcenter*temp[i];
    inv2dindex(global_id, &ii, &jj, nx, &horv);
    if (horv == HORIZONTAL) {
      Ap[i] += doffdiag*temp[my_local_ids[north2d(ii,jj,nx)]];
      Ap[i] +=       1.*temp[my_local_ids[west2d(ii,jj,nx)]];
      Ap[i] +=      -1.*temp[my_local_ids[east2d(ii,jj,nx)]];

      if (jj == 0) jj = nx-1;
      else jj--;

      Ap[i] +=      -1.*temp[my_local_ids[west2d(ii,jj,nx)]];
      Ap[i] += doffdiag*temp[my_local_ids[south2d(ii,jj,nx)]];
      Ap[i] +=       1.*temp[my_local_ids[east2d(ii,jj,nx)]];
    }
    else {
      Ap[i] += -1.*temp[my_local_ids[north2d(ii,jj,nx)]];
      Ap[i] += doffdiag*temp[my_local_ids[east2d(ii,jj,nx)]];
      Ap[i] +=  1.*temp[my_local_ids[south2d(ii,jj,nx)]];

      if (ii == 0) ii = nx-1;
      else ii--;

      Ap[i] += doffdiag*temp[my_local_ids[west2d(ii,jj,nx)]];
      Ap[i] += -1.*temp[my_local_ids[south2d(ii,jj,nx)]];
      Ap[i] +=  1.*temp[my_local_ids[north2d(ii,jj,nx)]];
    }
  }
  free(temp);
  return 1;
}

int Ke_getrow(void *data, int N_requested_rows, int requested_rows[],
	      int allocated_space, int columns[], double values[], 
	      int row_lengths[])
{
  int i, global_id, ii, jj, nx, horv;
  double dcenter, doffdiag, sigma = .001;
  int *my_local_ids;
  struct user_partition *Edge_Partition;

  if (allocated_space < 7) return 0;
  Edge_Partition = (struct user_partition *) data;
  my_local_ids = Edge_Partition->my_local_ids;
  nx = (int) sqrt( ((double) Edge_Partition->Nglobal/2) + .00001);
  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx));
  doffdiag = -1 + sigma/((double) ( 6 * nx * nx));
  global_id = (Edge_Partition->my_global_ids)[*requested_rows];
  i = 0;
  columns[i] = *requested_rows; values[i++] = dcenter;
  inv2dindex(global_id, &ii, &jj, nx, &horv);
  if (horv == HORIZONTAL) {
    columns[i] = my_local_ids[north2d(ii,jj,nx)]; values[i++] = doffdiag;
    columns[i] = my_local_ids[ west2d(ii,jj,nx)]; values[i++] =  1.;
    columns[i] = my_local_ids[ east2d(ii,jj,nx)]; values[i++] = -1.;

    if (jj == 0) jj = nx-1;
    else jj--;

    columns[i] = my_local_ids[ west2d(ii,jj,nx)]; values[i++] = -1.;
    columns[i] = my_local_ids[south2d(ii,jj,nx)]; values[i++] = doffdiag;
    columns[i] = my_local_ids[ east2d(ii,jj,nx)]; values[i++] =  1.;

  }
  else {
    columns[i] = my_local_ids[north2d(ii,jj,nx)]; values[i++] = -1.;
    columns[i] = my_local_ids[ east2d(ii,jj,nx)]; values[i++] = doffdiag;
    columns[i] = my_local_ids[south2d(ii,jj,nx)]; values[i++] =  1.;

    if (ii == 0) ii = nx-1;
    else ii--;

    columns[i] = my_local_ids[ west2d(ii,jj,nx)]; values[i++] = doffdiag;
    columns[i] = my_local_ids[south2d(ii,jj,nx)]; values[i++] = -1.;
    columns[i] = my_local_ids[north2d(ii,jj,nx)]; values[i++] = 1.;

  }
  *row_lengths = i;

  return 1;
}

int Kn_getrow(void *data, int N_requested_rows, int requested_rows[],
	      int allocated_space, int columns[], double values[], 
	      int row_lengths[])

{
  int i, ii, jj, nx, global_id, Nlocal_nodes;

  int *my_local_id;
  struct user_partition *Node_Partition;

  if (allocated_space < 9) return 0;
  Node_Partition = (struct user_partition *) data;
  my_local_id = Node_Partition->my_local_ids;

  Nlocal_nodes = Node_Partition->Nlocal;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);
#ifdef periodic
  nx = (int) sqrt( ((double) Node_Partition->Nglobal - 2) + .00001); //periodic stuff
Nlocal_nodes -= 2; // periodic stuff
#endif
 i = 0;

 global_id = (Node_Partition->my_global_ids)[*requested_rows];
 if (global_id < Node_Partition->Nglobal - 2) {
   inv2dnodeindex(global_id, &ii, &jj, nx);
   columns[i] = *requested_rows; values[i++] = 8;

   columns[i] = my_local_id[ southeast2d(ii,jj,nx)]; values[i++] = -1.;
   columns[i] = my_local_id[ northwest2d(ii,jj,nx)]; values[i++] = -1.;
   columns[i] = my_local_id[ northeast2d(ii,jj,nx)]; values[i++] = -1.;

   if (ii == 0) ii = nx-1;
   else ii--;

   columns[i] = my_local_id[ northwest2d(ii,jj,nx)]; values[i++] = -1.;

   if (jj == 0) jj = nx-1;
   else jj--;

   columns[i] = my_local_id[ southeast2d(ii,jj,nx)]; values[i++] = -1.;
   columns[i] = my_local_id[ northwest2d(ii,jj,nx)]; values[i++] = -1.;
   columns[i] = my_local_id[ southwest2d(ii,jj,nx)]; values[i++] = -.00000001;

   if (ii == nx-1) ii = 0;
   else ii++;

   columns[i] = my_local_id[ southeast2d(ii,jj,nx)]; values[i++] = -.00000001;
  }
#ifdef periodic
 else {
   columns[i] = *requested_rows; values[i++] = 8;
 }
#endif

  *row_lengths = i;

  return(1);
}

int Tmat_getrow(void *data, int N_requested_rows, int requested_rows[],
	      int allocated_space, int columns[], double values[], 
	      int row_lengths[])

{

  int nx, i, ii, jj, horv;
  int global_id;
  int Nlocal_nodes, Nlocal_edges; 
     struct user_partition *Edge_Partition;
     struct user_partition *Node_Partition;
     ML_Operator *Kn_mat;
     struct Tmat_data *Tmat_data;
     int *my_local_id;



  if (allocated_space < 3) return 0;
  Tmat_data = (struct Tmat_data *) data;
  Node_Partition = Tmat_data->node;
  Edge_Partition = Tmat_data->edge;
  Kn_mat              = Tmat_data->Kn;
  my_local_id = Node_Partition->my_local_ids;

  Nlocal_nodes = Node_Partition->Nlocal;
  Nlocal_edges = Edge_Partition->Nlocal;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);
#ifdef periodic
  nx = (int) sqrt( ((double) Node_Partition->Nglobal - 2) + .00001); //periodic stuff
#endif

  global_id = (Edge_Partition->my_global_ids)[*requested_rows];
  inv2dindex(global_id, &ii, &jj, nx, &horv);
  i = 0;

  if (horv == HORIZONTAL) {
    columns[i] = my_local_id[southwest2d(ii,jj,nx)]; values[i++] = -1.;
    columns[i] = my_local_id[southeast2d(ii,jj,nx)]; values[i++] =  1.;
#ifdef periodic
    columns[i] = my_local_id[Node_Partition->Nglobal-2]; values[i++] =  1.;
#endif
  }
  else {
    columns[i] = my_local_id[northwest2d(ii,jj,nx)]; values[i++] = -1.;
    columns[i] = my_local_id[southwest2d(ii,jj,nx)]; values[i++] =  1.;

#ifdef periodic
    columns[i] = my_local_id[Node_Partition->Nglobal-1]; values[i++] =  1.;
#endif
  }
  *row_lengths = i;

  return 1;
}
int Tmat_matvec(void *data, int Nlocal_nodes, double p[], int Nlocal_edges, double Ap[])
{
  int i, global_id, ii, jj, nx, horv;
  struct user_partition *Edge_Partition;
  struct user_partition *Node_Partition;
  ML_Operator *Kn_mat;
  struct Tmat_data *Tmat_data;
  int *my_local_id;
  double *temp;

  Tmat_data = (struct Tmat_data *) data;
  Node_Partition = Tmat_data->node;
  Edge_Partition = Tmat_data->edge;
  Kn_mat              = Tmat_data->Kn;
  my_local_id = Node_Partition->my_local_ids;
  temp = (double *) malloc(sizeof(double)*(Nlocal_nodes + 
					   Node_Partition->Nghost));

  for (i = 0; i < Nlocal_nodes; i++) temp[i] = p[i];
  update_ghost_nodes(temp, Node_Partition);



  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);
#ifdef periodic
  nx = (int) sqrt( ((double) Node_Partition->Nglobal - 2) + .00001); //periodic stuff
#endif

  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    inv2dindex(global_id, &ii, &jj, nx, &horv);
    Ap[i] = 0.;


    if (horv == HORIZONTAL) {
      Ap[i] += -1.*temp[my_local_id[southwest2d(ii,jj,nx)]];
      Ap[i] +=  1.*temp[my_local_id[southeast2d(ii,jj,nx)]];
#ifdef periodic
      Ap[i] +=  1.*temp[Nlocal_nodes-2];
#endif
    }
    else {
      Ap[i] += -1.*temp[my_local_id[northwest2d(ii,jj,nx)]];
      Ap[i] +=  1.*temp[my_local_id[southwest2d(ii,jj,nx)]];
#ifdef periodic
      Ap[i] +=  1.*temp[Nlocal_nodes-1];
#endif
    }
  }
  free(temp);
  return 1;

}

int update_ghost_edges(double vector[], void *data)
{

  struct user_partition *Partition;
  int *my_local_ids, i, nx, Nloc, type = 7103;
  MPI_Request  request;
  MPI_Status status;
  double *send_buf;

  printf("in side comm\n"); fflush(stdout);
  Partition = (struct user_partition *) data;
  if (Partition->nprocs == 1) return 1;

  my_local_ids = Partition->my_local_ids;
  
  nx = (int) sqrt( ((double) Partition->Nglobal/2) + .00001);
  Nloc = Partition->Nlocal;
  send_buf = (double *) malloc( (3*nx + 1)*sizeof(double));


  /* post receive */

  MPI_Irecv(&(vector[Nloc]), 3*nx, MPI_DOUBLE, 1 - Partition->mypid, type, 
	    MPI_COMM_WORLD, &request);

  /* send data */

  if (Partition->mypid == 1) {
      for (i = 0; i < nx; i++) {
	send_buf[i     ] = vector[my_local_ids[Nloc/2 + i    ]];
	send_buf[i + nx] = vector[my_local_ids[Nloc - nx + i ]];
	send_buf[i+2*nx] = vector[my_local_ids[2*Nloc - nx + i ]];
      }
  }
  else {
    for (i = 0; i < nx; i++) {
      send_buf[i     ] = vector[my_local_ids[Nloc/2 - nx + i ]];
      send_buf[i + nx] = vector[my_local_ids[              i ]];
      send_buf[i+2*nx] = vector[my_local_ids[3*Nloc/2 - nx + i ]];
    }
  }

  MPI_Send(send_buf, 3*nx, MPI_DOUBLE, 1 - Partition->mypid, type,
	   MPI_COMM_WORLD);


  /* receive data */

  MPI_Wait(&request, &status);
  printf("out side comm\n"); fflush(stdout);

  return 1;
}

int update_ghost_nodes(double vector[], void *data)
{

  struct user_partition *Partition;
  int *my_local_ids, i, nx, Nloc, type = 7107;
  MPI_Request  request;
  MPI_Status status;
  double *send_buf;

  Partition = (struct user_partition *) data;
  if (Partition->nprocs == 1) return 1;
  printf("nprocs = %d, proc id = %d\n",Partition->nprocs,
Partition->mypid); fflush(stdout);

  my_local_ids = Partition->my_local_ids;
  
  nx = (int) sqrt( ((double) Partition->Nglobal) + .00001);
  Nloc = Partition->Nlocal;
  send_buf = (double *) malloc( (2*nx + 1)*sizeof(double));


  /* post receive */

  MPI_Irecv(&(vector[Nloc]), 2*nx, MPI_DOUBLE, 1 - Partition->mypid, type, 
	    MPI_COMM_WORLD, &request);

  /* send data */

  if (Partition->mypid == 1) {
    for (i = 0; i < nx; i++) {
      send_buf[i   ] = vector[my_local_ids[Nloc + i ]];
      send_buf[i+nx] = vector[my_local_ids[2*Nloc - nx + i ]];
    }
  }
  else {
    for (i = 0; i < nx; i++) {
      send_buf[i   ] = vector[my_local_ids[Nloc - nx + i ]];
      send_buf[i+nx] = vector[my_local_ids[            i ]];
    }
  }

  MPI_Send(send_buf, 2*nx, MPI_DOUBLE, 1 - Partition->mypid, type,
	   MPI_COMM_WORLD);


  /* receive data */

  MPI_Wait(&request, &status);
  printf("out side comm\n"); fflush(stdout);

  return 1;
}
