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
/*                                                                           */
/* This particular example corresponds to Dirichlet boundary conditions on   */
/* the west and south boundaries and Neumman boundary conditions on the east */
/* and north boundaries. Dirichlet edges and Dirichlet nodes are not stored  */
/* in the matrix. The number of unknowns in this example is as follows:      */
/*                                                                           */
/*         ----12----(12)----13----(13)----14----(14)----15----(15)          */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*                    28            29            30            31           */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*         ---- 8----( 8)---- 9----( 9)----10----(10)----11----(11)          */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*                    24            25            26            27           */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*         ---- 4----( 4)---- 5----( 5)---- 6----( 6)---- 7----( 7)          */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*                    20            21            22            23           */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*         ---- 0----( 0)---- 1----( 1)---- 2----( 2)---- 3----( 3)          */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*                    16            17            18            19           */
/*                     |             |             |             |           */
/*                     |             |             |             |           */
/*                                                                           */
/* Quantities in () are nodes and the rest are edges.                        */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef AZTEC
#undef AxTEC
#endif
#ifdef AZTEC
#include "az_aztec.h"
#endif
#include "ml_include.h"

/*****************************************************************************/
/* All functions/structures starting with the word 'user' denote things that */
/* are not in ML and are specific to this particular example.                */
/* User defined structure holding information on how the PDE is partitioned  */
/* over the processor system.                                                */
/*****************************************************************************/
struct user_partition {               
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
/* User defined structure for performing a matrix-vector product and for     */
/* getting a row of the T matrix (null space).                               */
/*****************************************************************************/
struct user_Tmat_data {
  struct user_partition *edge;
  struct user_partition *node;
  ML_Operator *Kn;
};

/*****************************************************************************/
/* Function definitions.                                                     */
/*****************************************************************************/
extern void user_partition_edges(struct user_partition *,
				 struct user_partition *);
extern void user_partition_nodes(struct user_partition *Partition);
#ifdef AZTEC
extern AZ_MATRIX   *user_Ke_build(struct user_partition *);
extern AZ_MATRIX   *user_Kn_build(struct user_partition *);
extern ML_Operator *user_T_build (struct user_partition *, 
                                  struct user_partition *, ML_Operator *, ML_Comm *);
#ifdef ML_MPI
#define COMMUNICATOR   MPI_COMM_WORLD
#else
#define COMMUNICATOR   AZ_NOT_MPI
#endif

#else
extern  int user_update_ghost_edges(double vector[], void *data);
extern  int user_update_ghost_nodes(double vector[], void *data);
extern  int user_Ke_matvec(void *, int, double *, int , double *);
extern  int user_Ke_getrow(void *, int, int *, int, int *, double *, int *);
extern  int user_T_matvec(void *, int, double *, int, double *);
extern  int user_T_getrow(void *, int, int *, int, int *, double *, int *);
extern  int user_Kn_getrow(void *, int, int *, int, int *, double *, int *);
#endif


int main(int argc, char *argv[])
{
  int    Nnodes=32*32;              /* Total number of nodes in the problem.*/
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
  struct       user_partition Edge_Partition = {NULL, NULL,0,0,NULL,0,0,0}, 
                                Node_Partition = {NULL, NULL,0,0,NULL,0,0,0};

  /* See Aztec User's Guide for information on these variables */

#ifdef AZTEC
  AZ_MATRIX    *Ke_mat, *Kn_mat;
  AZ_PRECOND   *Pmat = NULL;
  int          proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double       params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];
#else
  struct user_Tmat_data user_Tmat_data;
  ML_Comm *comm;
  ML_Krylov *kdata;
#endif


  /* get processor information (id & # of procs) and set ML's printlevel. */

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
#ifdef AZTEC
  AZ_set_proc_config(proc_config, COMMUNICATOR);
#endif
  ML_Set_PrintLevel(15);   /* set ML's output level: 0 gives least output */

  /* Set the # of global nodes/edges and partition both the edges and the */
  /* nodes over the processors. NOTE: I believe we assume that if an edge */
  /* is assigned to a processor at least one of its nodes must be also    */
  /* assigned to that processor.                                          */

  Node_Partition.Nglobal = Nnodes;
  Edge_Partition.Nglobal = Node_Partition.Nglobal*2;
          /* in this simple example the number of edges is 2X the */
          /* number of nodes. This will not generally be true.    */

  user_partition_nodes(&Node_Partition);
  user_partition_edges(&Edge_Partition, &Node_Partition);
       /* IMPORTANT: An edge can be assigned to a processor only if at  */
       /*            least one its nodes is assigned to that processor. */

  /* 1) Create an empty multigrid hierarchy and set the 'MaxMGLevels-1'th   */
  /* level discretization within this hierarchy to the ML matrix            */
  /* representing Ke (Maxwell edge discretization).                         */
  /* 2) Create a second empty multigrid hierarchy. Build an auxiliary nodal */
  /* PDE problem. This should be a variable coefficient Poisson problem     */
  /* (with unknowns at the nodes). The coefficients should be chosen to     */
  /* reflect the conductivity of the original edge problems. Set the        */
  /* 'MaxMGLevels-1'th of the second hierarchy to this nodal matrix.        */
  /* 3) Build an ML matrix representing the null space of the PDE problem.  */
  /* This should be a discrete gradient (nodes to edges).                   */

  /* Note it is possible to multiply T'*T for get a nodal discretization    */
  /* though this will not incorportate material properties.                 */

  ML_Create(&ml_edges, MaxMgLevels);
  ML_Create(&ml_nodes, MaxMgLevels);

#ifdef AZTEC
  /* Build Ke/Kn as Aztec matrices. Use built-in function AZ_ML_Set_Amat() */
  /* to convert to an ML matrix and put in hierarchy.                      */

  Ke_mat = user_Ke_build(&Edge_Partition);
  AZ_ML_Set_Amat(ml_edges, MaxMgLevels-1, Edge_Partition.Nlocal,
      		 Edge_Partition.Nlocal, Ke_mat, proc_config);

  Kn_mat = user_Kn_build( &Node_Partition);
  AZ_ML_Set_Amat(ml_nodes, MaxMgLevels-1, Node_Partition.Nlocal, 
		 Node_Partition.Nlocal, Kn_mat, proc_config);

  /* Use Aztec to build ML matrix representing Tmat.                      */

  Tmat = user_T_build (&Edge_Partition, &Node_Partition, 
  		   &(ml_nodes->Amat[MaxMgLevels-1]),ml_nodes->comm);

#else
  /* Build Ke directly as an ML matrix.                                  */

  ML_Init_Amatrix      (ml_edges, MaxMgLevels-1, Edge_Partition.Nlocal,
			Edge_Partition.Nlocal, &Edge_Partition);

  ML_Set_Amatrix_Getrow(ml_edges, MaxMgLevels-1,  user_Ke_getrow, 
			user_update_ghost_edges,  
			Edge_Partition.Nlocal + Edge_Partition.Nghost);

  ML_Set_Amatrix_Matvec(ml_edges, MaxMgLevels-1,  user_Ke_matvec);

  /* Build Kn directly as an ML matrix.                                  */

  ML_Init_Amatrix      (ml_nodes, MaxMgLevels-1 , Node_Partition.Nlocal,
			Node_Partition.Nlocal, &Node_Partition);

  ML_Set_Amatrix_Getrow(ml_nodes, MaxMgLevels-1,  user_Kn_getrow, 
			user_update_ghost_nodes, 
			Node_Partition.Nlocal + Node_Partition.Nghost);

  Tmat = ML_Operator_Create(ml_nodes->comm);
  user_Tmat_data.edge = &Edge_Partition;
  user_Tmat_data.node = &Node_Partition;
  user_Tmat_data.Kn   = &(ml_nodes->Amat[MaxMgLevels-1]);

  ML_Operator_Set_ApplyFuncData(Tmat, Node_Partition.Nlocal,
				Edge_Partition.Nlocal, ML_EXTERNAL, 
				(void *) &user_Tmat_data, Edge_Partition.Nlocal,
				  user_T_matvec, 0);
  ML_Operator_Set_Getrow(Tmat,ML_EXTERNAL,Edge_Partition.Nlocal,user_T_getrow);

  ML_CommInfoOP_Generate( &(Tmat->getrow->pre_comm), user_update_ghost_nodes, 
			    &Node_Partition, ml_edges->comm, Tmat->invec_leng, 
			    Node_Partition.Nghost);

#endif


  /********************************************************************/
  /* Set some ML parameters.                                          */
  /*------------------------------------------------------------------*/
	
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
    ML_gsum_scalar_int(&Nfine_node, &itmp, ml_edges->comm);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    Nfine_edge = Tmat_array[MaxMgLevels-1]->outvec_leng;
    ML_gsum_scalar_int(&Nfine_edge, &itmp, ml_edges->comm);
  }

  /***************************************************************
  * Set up smoothers for all levels. For the MLS polynomial pick 
  * parameters based on the coarsening rate. See paper 'Parallel
  * Multigrid Smoothing: Polynomial Versus Gauss-Seidel' by Adams,
  * Brezina, Hu, Tuminaro
  ****************************************************************/
  coarsest_level = MaxMgLevels - Nlevels;

  for (level = MaxMgLevels-1; level >= coarsest_level; level--) {
    if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
      if (level != coarsest_level) {
	Ncoarse_edge = Tmat_array[level-1]->outvec_leng;
	ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_edges->comm);
	edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                   ((double) Ncoarse_edge);
      }
      else edge_coarsening_rate =  (double) Nfine_edge;

      ML_Smoother_Arglist_Set(edge_args, 1, &edge_coarsening_rate);
      Nfine_edge = Ncoarse_edge;
    }
    if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
      if (level != coarsest_level) {
	Ncoarse_node = Tmat_array[level-1]->invec_leng;
	ML_gsum_scalar_int(&Ncoarse_node, &itmp, ml_edges->comm);
	node_coarsening_rate =  2.*((double) Nfine_node)/ 
                                   ((double) Ncoarse_node);
      }
      else node_coarsening_rate = (double) Nfine_node;

      ML_Smoother_Arglist_Set(nodal_args, 1, &node_coarsening_rate);
      Nfine_node = Ncoarse_node;
    }
    ML_Gen_Smoother_Hiptmair(ml_edges, level, ML_BOTH, Nits_per_presmooth,
			     Tmat_array, Tmat_trans_array, NULL, 
			     edge_smoother, edge_args, nodal_smoother,
			     nodal_args, hiptmair_type);
  }

  /* Must be called before invoking the preconditioner */
  ML_Gen_Solver(ml_edges, ML_MGV, MaxMgLevels-1, coarsest_level); 


  /* Set initial guess and right hand side. */

#ifdef AZTEC
  xxx = (double *) ML_allocate((Edge_Partition.Nlocal+
				Edge_Partition.Nghost)*sizeof(double)); 
#else
  xxx = (double *) ML_allocate(Edge_Partition.Nlocal*sizeof(double)); 
#endif
  rhs = (double *) ML_allocate(Edge_Partition.Nlocal*sizeof(double)); 
  ML_random_vec(xxx, Edge_Partition.Nlocal, ml_edges->comm);
  ML_random_vec(rhs, Edge_Partition.Nlocal, ml_edges->comm);

  /* Invoke solver */
#ifdef AZTEC

  AZ_defaults(options, params);
  options[AZ_solver]   = AZ_cg;
  params[AZ_tol]       = tolerance;
  options[AZ_conv]     = AZ_noscaled;
  AZ_set_ML_preconditioner(&Pmat, Ke_mat, ml_edges, options); 
  AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, NULL);
#else
  kdata = ML_Krylov_Create(ml_edges->comm);
  ML_Krylov_Set_PrintFreq( kdata, 1 );
  ML_Krylov_Set_Method(kdata, ML_CG);
  ML_Krylov_Set_Amatrix(kdata, &(ml_edges->Amat[MaxMgLevels-1]));
  ML_Krylov_Set_PreconFunc(kdata, ML_MGVSolve_Wrapper);
  ML_Krylov_Set_Precon(kdata, ml_edges);
  ML_Krylov_Set_Tolerance(kdata, tolerance);
  ML_Krylov_Solve(kdata, Edge_Partition.Nlocal, rhs, xxx);
  ML_Krylov_Destroy( &kdata );
#endif
  
  /* clean up. */

  if (Edge_Partition.my_local_ids != NULL) free(Edge_Partition.my_local_ids);
  if (Node_Partition.my_local_ids != NULL) free(Node_Partition.my_local_ids);
  if (Node_Partition.my_global_ids != NULL) free(Node_Partition.my_global_ids);
  if (Edge_Partition.my_global_ids != NULL) free(Edge_Partition.my_global_ids);
  if (Node_Partition.needed_external_ids != NULL) 
    free(Node_Partition.needed_external_ids);
  if (Edge_Partition.needed_external_ids != NULL) 
    free(Edge_Partition.needed_external_ids);


  ML_Smoother_Arglist_Delete(&nodal_args);
  ML_Smoother_Arglist_Delete(&edge_args);
  ML_Aggregate_Destroy(&ag);
  ML_Destroy(&ml_edges);
  ML_Destroy(&ml_nodes);
#ifdef AZTEC
  if (Ke_mat  != NULL) {
    AZ_free(Ke_mat->bindx);
    AZ_free(Ke_mat->val);
    AZ_free(Ke_mat->data_org);
    AZ_matrix_destroy(&Ke_mat);
  }
  if (Pmat  != NULL) AZ_precond_destroy(&Pmat);
  if (Kn_mat != NULL) {
    AZ_free(Kn_mat->bindx);
    AZ_free(Kn_mat->val);
    AZ_free(Kn_mat->data_org);
    AZ_matrix_destroy(&Kn_mat);
  }
#endif
  free(xxx);
  free(rhs);
  ML_Operator_Destroy(&Tmat);
  ML_Operator_Destroy(&Tmat_trans);
  ML_MGHierarchy_ReitzingerDestroy(MaxMgLevels-2, &Tmat_array,
				   &Tmat_trans_array);
#ifdef ML_MPI
  MPI_Finalize();
#endif
		
  return 0;
		
}

#define HORIZONTAL 1
#define VERTICAL   2
extern int invindex(int index, int *i, int *j, int n, int *hor_or_vert);
extern int northwest(int i, int j, int n);
extern int southeast(int i, int j, int n);
extern int southwest(int i, int j, int n);
extern int north(int i, int j, int n);
extern int south(int i, int j, int n);
extern int east(int i, int j, int n);
extern int west(int i, int j, int n);

/********************************************************************/
/* Given the (i,j)th cell return the node indices defining the cell */
/*                                                                  */
/*                          NW-------NE                             */
/*                          |         |                             */
/*                          |         |                             */
/*                          |         |                             */
/*                          SW--------SE                            */
/*                                                                  */
/********************************************************************/

int northwest(int i, int j, int n) { return(i + (j)*n); }
int southwest(int i, int j, int n) { if (j == 0) j = n;
                                       return(i + (j-1)*n);}
int southeast(int i, int j, int n)  { if (j == 0) j = n;    
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

int north(int i, int j, int n) { return(i + (j)*n);}
int south(int i, int j, int n) { return(i + (j-1)*n);}
int west(int i, int j, int n)  { return(i+(j)*n+n*n -1);}
int east(int i, int j, int n)  { return(i+1+(j)*n+n*n -1);}

/* Given the index of either a 'south' or 'west' edge, return the */
/* cell coordinates (i,j) as well as the orientation (horizontal  */
/* or vertical) of the edge.                                      */
int invindex(int index, int *i, int *j, int n, int *hor_or_vert)
{
  *hor_or_vert = HORIZONTAL;
  if (index >= n*n) {
    *hor_or_vert = VERTICAL;
    index -= n*n;
  }
  *i = (index%n);
  *j = (index - (*i))/n;
  if ( *hor_or_vert == HORIZONTAL)   (*j) =  ((*j)+1)%n;
  if ( *hor_or_vert == VERTICAL)     (*i) =  ((*i)+1)%n;
  return 1;
}

/* Assign unknowns to processors */
void user_partition_nodes(struct user_partition *Partition)
{
#ifdef AZTEC
  int    proc_config[AZ_PROC_SIZE];
#else
  int i, nx, Nloc;
  ML_Comm *comm;
  int *my_local_ids;
#endif

#ifdef AZTEC
  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_input_update(NULL,&(Partition->Nlocal), &(Partition->my_global_ids),
		  proc_config, Partition->Nglobal, 1, AZ_linear);
  Partition->Nghost = 0;
#else
  ML_Comm_Create( &comm);
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

  Partition->mypid = comm->ML_mypid;
  Partition->nprocs= comm->ML_nprocs;

  /* Fill partition structure */

  if (comm->ML_nprocs == 1) {
    Partition->Nghost = 0;
    Partition->Nlocal = Partition->Nglobal; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    for (i = 0; i < Partition->Nlocal; i++) (Partition->my_global_ids)[i] = i;
    
    Partition->my_local_ids = (int *) malloc(sizeof(int)*Partition->Nglobal);
    for (i = 0; i < Partition->Nglobal; i++)  (Partition->my_local_ids)[i] = i;
  }
  else {
    /* on 2 processors assignment is done via 2 horizontal strips  */
    /* Note: local_ids's for ghost unknowns must follow local_id's */
    /*       for uknowns assigned to this processor. Further,      */
    /*       ghost nodes assigned to the same processor must have  */
    /*       consecutively numbered local_id's.                    */

    nx = (int) sqrt( ((double) Partition->Nglobal) + .00001);
    Partition->Nghost = nx;
    Partition->Nlocal = Partition->Nglobal/2; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    my_local_ids            =  (int *) malloc(sizeof(int)*Partition->Nglobal);
    Partition->my_local_ids = my_local_ids;

    /* set local and global id's of unknowns assigned to this processor */
    for (i = 0; i < Partition->Nlocal; i++) {
      (Partition->my_global_ids)[i] = i + comm->ML_mypid*Partition->Nlocal;
      my_local_ids[(Partition->my_global_ids)[i]] = i;
    }

    /* set local id's of ghost unknowns */
    Nloc = Partition->Nlocal;
    if (comm->ML_mypid == 0) {
      for (i = 0; i < nx; i++) 	my_local_ids[Nloc + i] = Nloc + i;
    }
    else {
      for (i = 0; i < nx; i++) 	my_local_ids[Nloc - nx + i ] = Nloc + i;
    }

  }
  ML_Comm_Destroy(&comm);
#endif
}
/* Assign unknowns to processors */
void user_partition_edges(struct user_partition *Partition,
			  struct user_partition *NodePartition)
{
  int i;
#ifndef AZTEC
  int nx, Nloc, *my_local_ids;
  ML_Comm *comm;
#endif

#ifdef AZTEC
  /* Edge partition is derived from node partition */

  Partition->Nlocal = NodePartition->Nlocal*2;
  Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);

  for (i = 0; i < NodePartition->Nlocal; i++)
    (Partition->my_global_ids)[i] = (NodePartition->my_global_ids)[i];
  for (i = 0; i < NodePartition->Nlocal; i++)
    (Partition->my_global_ids)[i+NodePartition->Nlocal] = 
               (NodePartition->my_global_ids)[i] + NodePartition->Nglobal;

#else
  ML_Comm_Create( &comm);

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

  Partition->mypid = comm->ML_mypid;
  Partition->nprocs= comm->ML_nprocs;

  if (comm->ML_nprocs == 1) { 
    Partition->Nghost = 0;
    Partition->Nlocal = Partition->Nglobal; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    Partition->my_local_ids  = (int *) malloc(sizeof(int)*Partition->Nglobal);

    /* Set local and global id's */
    for (i = 0; i < Partition->Nlocal; i++)  (Partition->my_global_ids)[i] = i;
    for (i = 0; i < Partition->Nglobal; i++) (Partition->my_local_ids)[i] = i;
  }
  else {
    nx = (int) sqrt( ((double) Partition->Nglobal/2) + .00001);
    Partition->Nlocal = Partition->Nglobal/2; 
    Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);
    my_local_ids             =  (int *) malloc(sizeof(int)*Partition->Nglobal);
    Partition->my_local_ids  = my_local_ids;


    /* set global ids */
    for (i = 0; i < Partition->Nlocal/2; i++) {
      (Partition->my_global_ids)[i] = i + comm->ML_mypid*Partition->Nlocal/2;
      (Partition->my_global_ids)[i + Partition->Nlocal/2] = Partition->Nlocal +
	                                i + comm->ML_mypid*Partition->Nlocal/2;
    }

    /* set local ids of nonghost unknowns */

    for (i = 0; i < Partition->Nlocal; i++) 
      my_local_ids[(Partition->my_global_ids)[i]] = i;

    /* set the ghost unknowns */

    Nloc = Partition->Nlocal;
    if (comm->ML_mypid == 0) {
      Partition->Nghost = 2*nx;
      for (i = 0; i < nx; i++) {
	my_local_ids[Nloc/2 + i ] = Nloc + i;
	my_local_ids[3*Nloc/2  + i ] = Nloc + i + nx;
      }
    }
    else {
      Partition->Nghost = nx;
      for (i = 0; i < nx; i++) {
	my_local_ids[Nloc/2 - nx + i ] = Nloc + i;
      }
    }

  }
  ML_Comm_Destroy(&comm);
#endif
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
  double dcenter, doff, sigma = .0001;
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
  doff = -1 + sigma/((double) ( 6 * nx * nx));

  /* Create a DMSR matrix with global column indices */

  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    invindex(global_id, &ii, &jj, nx, &horv);
    nz_ptr = Ke_bindx[i];

    Ke_val[i] = dcenter;

    if (horv == HORIZONTAL) {
      if (jj != 0) {
	Ke_bindx[nz_ptr] = north(ii,jj,nx);     Ke_val[nz_ptr++] = doff;
	Ke_bindx[nz_ptr] = east(ii,jj,nx);      Ke_val[nz_ptr++] = -1.;
	if (ii != 0) {Ke_bindx[nz_ptr]=west(ii,jj,nx); Ke_val[nz_ptr++]= 1.;}
	jj--;
      }
      else {
	Ke_val[i] = 1. +  2.*sigma/((double) ( 3 * nx * nx));
	jj = nx-1;
      }
      Ke_bindx[nz_ptr] = east(ii,jj,nx);      Ke_val[nz_ptr++] = 1.;
      if (ii != 0){ Ke_bindx[nz_ptr]=west(ii,jj,nx);  Ke_val[nz_ptr++]=-1.;}
      if (jj != 0){ Ke_bindx[nz_ptr]=south(ii,jj,nx); Ke_val[nz_ptr++]=doff;}
    }
    else {
      if (ii != 0) {
	Ke_bindx[nz_ptr] = north(ii,jj,nx);     Ke_val[nz_ptr++] = -1.;
	Ke_bindx[nz_ptr] = east(ii,jj,nx);      Ke_val[nz_ptr++] = doff;
	if (jj != 0) {Ke_bindx[nz_ptr]=south(ii,jj,nx); Ke_val[nz_ptr++]=1.;}
	ii--;
      }
      else {
	Ke_val[i]  = 1 +  2.*sigma/((double) ( 3 * nx * nx));
	ii = nx-1;
      }
      Ke_bindx[nz_ptr] = north(ii,jj,nx);     Ke_val[nz_ptr++] = 1.;
      if (ii != 0) {Ke_bindx[nz_ptr]=west(ii,jj,nx);  Ke_val[nz_ptr++]=doff;}
      if (jj != 0) {Ke_bindx[nz_ptr]=south(ii,jj,nx); Ke_val[nz_ptr++]=-1.;}
    }
    Ke_bindx[i+1] = nz_ptr;
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
  AZ_free(reordered_glob_edges);
  AZ_free(reordered_edge_externs);
  Edge_Partition->Nghost = Ke_data_org[AZ_N_external];

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
  int i, ii, jj, nx, gid, Nlocal_nodes, nz_ptr;


  Nlocal_nodes = Node_Partition->Nlocal;
  Kn_bindx = (int    *) malloc((27*Nlocal_nodes+5)*sizeof(int));
  Kn_val   = (double *) malloc((27*Nlocal_nodes+5)*sizeof(double));
  Kn_bindx[0] = Nlocal_nodes+1;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);

  for (i = 0; i < Nlocal_nodes; i++) {
    gid = (Node_Partition->my_global_ids)[i];

    nz_ptr = Kn_bindx[i];
    ii = gid%nx;
    jj = (gid - ii)/nx;


    if (ii != nx-1) { Kn_bindx[nz_ptr] = gid+ 1; Kn_val[nz_ptr++] = -1.;}
    if (jj != nx-1) { Kn_bindx[nz_ptr] = gid+nx; Kn_val[nz_ptr++] = -1.;}
    if (jj !=    0) { Kn_bindx[nz_ptr] = gid-nx; Kn_val[nz_ptr++] = -1.;}
    if (ii !=    0) { Kn_bindx[nz_ptr] = gid- 1; Kn_val[nz_ptr++] = -1.;}

    if ((ii != nx-1) && (jj !=    0)) 
      {Kn_bindx[nz_ptr] = gid-nx+1; Kn_val[nz_ptr++] = -1.;}
    if ((ii != nx-1) && (jj != nx-1)) 
      {Kn_bindx[nz_ptr] = gid+nx+1; Kn_val[nz_ptr++] = -1.;}
    if ((ii !=    0) && (jj != nx-1)) 
      {Kn_bindx[nz_ptr] = gid+nx-1; Kn_val[nz_ptr++] = -1.;}
    if ((ii !=    0) && (jj !=    0)) 
      {Kn_bindx[nz_ptr] = gid-nx-1; Kn_val[nz_ptr++] = -1.;}
    Kn_val[i] = (double) (nz_ptr - Kn_bindx[i]);
    Kn_bindx[i+1] = nz_ptr;
  }

  /* Transform the global Aztec matrix into a local Aztec matrix. That is,   */
  /* replace global column indices by local indices and set up communication */
  /* data structure 'Ke_data_org' that will be used for matvec's.            */

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_transform_norowreordering(proc_config,&(Node_Partition->needed_external_ids),
			       Kn_bindx, Kn_val, Node_Partition->my_global_ids,
			       &reordered_glob_nodes, &reordered_node_externs, 
			       &Kn_data_org, Nlocal_nodes, 0, 0, 0, 
			       &cpntr, AZ_MSR_MATRIX);
  Node_Partition->Nghost = Kn_data_org[AZ_N_external];
  AZ_free(reordered_glob_nodes);

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
		      ML_Operator *Kn_mat, ML_Comm *comm)
{

  int nx, i, ii, jj, horv, Ncols, Nexterns;
  int *Tmat_bindx;
  double *Tmat_val;
  ML_Operator *Tmat;
  struct ML_CSR_MSRdata *csr_data;
  struct aztec_context *aztec_context;
  int global_id;
  int Nlocal_nodes, Nlocal_edges; 
  int nz_ptr;

  Nlocal_nodes = Node_Partition->Nlocal;
  Nlocal_edges = Edge_Partition->Nlocal;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);

  Tmat_bindx = (int    *) malloc((3*Nlocal_edges+1)*sizeof(int));
  Tmat_val   = (double *) malloc((3*Nlocal_edges+1)*sizeof(double));

  Tmat_bindx[0] = Nlocal_edges + 1;
  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    Tmat_val[i] = 0.0;

    invindex(global_id, &ii, &jj, nx, &horv);
    nz_ptr = Tmat_bindx[i];


    ii--;

    if (horv == HORIZONTAL) {
      if(ii != -1) {
	Tmat_bindx[nz_ptr] = southwest(ii,jj,nx); Tmat_val[nz_ptr++] = -1.;
      }
      Tmat_bindx[nz_ptr]   = southeast(ii,jj,nx); Tmat_val[nz_ptr++] =  1.; 
    }
    else {
      if (ii == -1) ii = nx-1;
      Tmat_bindx[nz_ptr] = northwest(ii,jj,nx); Tmat_val[nz_ptr++] = -1.;
      if (jj != 0) {
	Tmat_bindx[nz_ptr] = southwest(ii,jj,nx); Tmat_val[nz_ptr++] =  1.;}
    }
    Tmat_bindx[i+1] = nz_ptr;

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


  /*  ML_Comm_Create( &comm); */

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
#else

int user_Ke_matvec(void *data, int Nlocal_edges, double p[], int N_out, 
		   double Ap[])
{
  int i, global_id, ii, jj, nx, horv;
  double dcenter, doff, sigma = .0001, *temp;
  int *id;
  struct user_partition *Edge_Partition;

  Edge_Partition = (struct user_partition *) data;
  id = Edge_Partition->my_local_ids;
  nx = (int) sqrt( ((double) Edge_Partition->Nglobal/2) + .00001);
  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx));
  doff = -1 + sigma/((double) ( 6 * nx * nx));
  temp = (double *) malloc(sizeof(double)*(Edge_Partition->Nghost+Nlocal_edges));

  /* communicate */
  for (i = 0; i < Nlocal_edges; i++) temp[i] = p[i];
  user_update_ghost_edges(temp, (void *) Edge_Partition);

  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    Ap[i] = dcenter*temp[i];
    invindex(global_id, &ii, &jj, nx, &horv);
    if (horv == HORIZONTAL) {
      if (jj != 0) {
	Ap[i] = dcenter*temp[i];
	Ap[i] += doff*temp[id[north(ii,jj,nx)]];
	Ap[i] +=  -1.*temp[id[ east(ii,jj,nx)]];
	if (ii !=    0) Ap[i] +=  1.*temp[id[ west(ii,jj,nx)]]; 
	jj--;
      }
      else {
	Ap[i] = ( 1. +  2.*sigma/((double) ( 3 * nx * nx)))*temp[i];
	jj = nx-1;
      }
      Ap[i] +=  1.*temp[id[ east(ii,jj,nx)]];
      if (ii !=    0)     Ap[i] +=  -1.*temp[id[west(ii,jj,nx)]]; 
      if (jj !=    0) Ap[i] += doff*temp[id[south(ii,jj,nx)]];
    }
    else {
      if (ii != 0) {
	Ap[i] = dcenter*temp[i];
	Ap[i] += -1.*temp[id[north(ii,jj,nx)]];
	Ap[i] += doff*temp[id[ east(ii,jj,nx)]];
	if (jj !=    0) Ap[i] += 1.*temp[id[south(ii,jj,nx)]];
	ii--;
      }
      else {
	Ap[i] = ( 1. +  2.*sigma/((double) ( 3 * nx * nx)))*temp[i];
	ii = nx-1;
      }
      Ap[i] += 1.*temp[id[north(ii,jj,nx)]]; 
      if (ii !=    0) Ap[i] += doff*temp[id[ west(ii,jj,nx)]];
      if (jj !=    0) Ap[i] += -1.*temp[id[south(ii,jj,nx)]];
    }
  }
  free(temp);
  return 1;
}

int user_Ke_getrow(void *data, int N_requested_rows, int requested_rows[],
	      int allocated_space, int columns[], double values[], 
	      int row_lengths[])
{
  int i, global_id, ii, jj, nx, horv;
  double dcenter, doff, sigma = .0001;
  int *id;
  struct user_partition *Edge_Partition;

  if (allocated_space < 7) return 0;
  Edge_Partition = (struct user_partition *) data;
  id = Edge_Partition->my_local_ids;
  nx = (int) sqrt( ((double) Edge_Partition->Nglobal/2) + .00001);
  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx));
  doff = -1 + sigma/((double) ( 6 * nx * nx));
  global_id = (Edge_Partition->my_global_ids)[*requested_rows];
  i = 0;
  invindex(global_id, &ii, &jj, nx, &horv);

  columns[i] = *requested_rows; 

  if (horv == HORIZONTAL) {
    if (jj != 0) {
      values[i++] = dcenter;
      columns[i] = id[north(ii,jj,nx)]; values[i++] = doff;
      columns[i] = id[ east(ii,jj,nx)]; values[i++] = -1.;
      if (ii !=    0) {columns[i] = id[ west(ii,jj,nx)]; values[i++] =  1.;}
      jj--;
    }
    else {
      values[i++] = 1 +  2.*sigma/((double) ( 3 * nx * nx));
      jj = nx-1;
    }
    columns[i] = id[ east(ii,jj,nx)]; values[i++] =  1.;
    if (ii !=    0) {columns[i] = id[ west(ii,jj,nx)]; values[i++] = -1.;}
    if (jj !=    0) {columns[i] = id[south(ii,jj,nx)]; values[i++] = doff;}
  }
  else {
    if (ii != 0) {
      values[i++] = dcenter;
      columns[i] = id[north(ii,jj,nx)]; values[i++] = -1.;
      columns[i] = id[ east(ii,jj,nx)]; values[i++] = doff;
      if (jj !=    0) {columns[i] = id[south(ii,jj,nx)]; values[i++] =  1.;}
      ii--;
    }
    else {
      values[i++] = 1 +  2.*sigma/((double) ( 3 * nx * nx));
      ii = nx-1;
    }
    columns[i] = id[north(ii,jj,nx)]; values[i++] = 1.;
    if (ii !=    0) {columns[i] = id[ west(ii,jj,nx)]; values[i++] = doff;}
    if (jj !=    0) {columns[i] = id[south(ii,jj,nx)]; values[i++] = -1.;}
  }
  *row_lengths = i;

  return 1;
}

int user_Kn_getrow(void *data, int N_requested_rows, int requested_rows[],
	      int allocated_space, int cols[], double val[], 
	      int row_lengths[])

{
  int i, ii, jj, nx, gid;

  int *ids;
  struct user_partition *Node_Partition;

  if (allocated_space < 9) return 0;
  Node_Partition = (struct user_partition *) data;
  ids = Node_Partition->my_local_ids;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);
  i = 0;

  gid = (Node_Partition->my_global_ids)[*requested_rows];

  ii = gid%nx;
  jj = (gid - ii)/nx;


  if (ii != nx-1) { cols[i] = ids[gid+ 1]; val[i++] = -1.;}
  if (jj != nx-1) { cols[i] = ids[gid+nx]; val[i++] = -1.;}
  if (jj !=    0) { cols[i] = ids[gid-nx]; val[i++] = -1.;}
  if (ii !=    0) { cols[i] = ids[gid- 1]; val[i++] = -1.;}

  if ((ii != nx-1) && (jj !=    0)) {cols[i] = ids[gid-nx+1]; val[i++] = -1.;}
  if ((ii != nx-1) && (jj != nx-1)) {cols[i] = ids[gid+nx+1]; val[i++] = -1.;}
  if ((ii !=    0) && (jj != nx-1)) {cols[i] = ids[gid+nx-1]; val[i++] = -1.;}
  if ((ii !=    0) && (jj !=    0)) {cols[i] = ids[gid-nx-1]; val[i++] = -1.;}
  cols[i] = *requested_rows; val[i++] = (double) i;

  *row_lengths = i;

  return(1);
}

int user_T_getrow(void *data, int N_requested_rows, int requested_rows[],
	      int allocated_space, int columns[], double values[], 
	      int row_lengths[])

{

  int nx, i, ii, jj, horv;
  int global_id;
  struct user_partition *Edge_Partition;
  struct user_partition *Node_Partition;
  struct user_Tmat_data *user_Tmat_data;
  int *my_local_id;



  if (allocated_space < 3) return 0;
  user_Tmat_data = (struct user_Tmat_data *) data;
  Node_Partition = user_Tmat_data->node;
  Edge_Partition = user_Tmat_data->edge;
  my_local_id = Node_Partition->my_local_ids;

  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);

  global_id = (Edge_Partition->my_global_ids)[*requested_rows];
  invindex(global_id, &ii, &jj, nx, &horv);

  i = 0;
  ii--;

  if (horv == HORIZONTAL) {
    if(ii != -1) {
      columns[i] = my_local_id[southwest(ii,jj,nx)]; values[i++] = -1.;}
      columns[i] = my_local_id[southeast(ii,jj,nx)]; values[i++] =  1.; 
  }
  else {
    if (ii == -1) ii = nx-1;
    columns[i] = my_local_id[northwest(ii,jj,nx)]; values[i++] = -1.;
    if (jj != 0) {
      columns[i] = my_local_id[southwest(ii,jj,nx)]; values[i++] =  1.;}
  }
  *row_lengths = i;

  return 1;
}
int user_T_matvec(void *data, int Nlocal_nodes, double p[], int Nlocal_edges, double Ap[])
{
  int i, global_id, ii, jj, nx, horv;
  struct user_partition *Edge_Partition;
  struct user_partition *Node_Partition;
  struct user_Tmat_data *user_Tmat_data;
  int *my_local_id;
  double *temp;

  user_Tmat_data = (struct user_Tmat_data *) data;
  Node_Partition = user_Tmat_data->node;
  Edge_Partition = user_Tmat_data->edge;
  my_local_id = Node_Partition->my_local_ids;
  temp = (double *) malloc(sizeof(double)*(Nlocal_nodes + 
					   Node_Partition->Nghost));

  for (i = 0; i < Nlocal_nodes; i++) temp[i] = p[i];
  user_update_ghost_nodes(temp, Node_Partition);



  nx = (int) sqrt( ((double) Node_Partition->Nglobal) + .00001);


  for (i = 0; i < Nlocal_edges; i++) {
    global_id = (Edge_Partition->my_global_ids)[i];
    invindex(global_id, &ii, &jj, nx, &horv);
    Ap[i] = 0.;

    ii--;

    if (horv == HORIZONTAL) {
      if(ii != -1) Ap[i] += -1.*temp[my_local_id[southwest(ii,jj,nx)]];
      Ap[i] += 1.*temp[my_local_id[southeast(ii,jj,nx)]];
    }
    else {
      if (ii == -1) ii = nx-1;
      Ap[i] += -1.*temp[my_local_id[northwest(ii,jj,nx)]];
      if (jj != 0) Ap[i] += 1.*temp[my_local_id[southwest(ii,jj,nx)]];
    }
  }
  free(temp);
  return 1;

}

int user_update_ghost_edges(double vector[], void *data)
{

  struct user_partition *Partition;
  int *my_local_ids, i, nx, Nloc, type = 7103;
  MPI_Request  request;
  MPI_Status status;
  double *send_buf;

  Partition = (struct user_partition *) data;
  if (Partition->nprocs == 1) return 1;

  my_local_ids = Partition->my_local_ids;
  
  nx = (int) sqrt( ((double) Partition->Nglobal/2) + .00001);
  Nloc = Partition->Nlocal;
  send_buf = (double *) malloc( (3*nx + 1)*sizeof(double));


  /* post receive */

  MPI_Irecv(&(vector[Nloc]), Partition->Nghost, MPI_DOUBLE, 
	    1 - Partition->mypid, type, MPI_COMM_WORLD, &request);

  /* send data */

  if (Partition->mypid == 1) {
      for (i = 0; i < nx; i++) {
	send_buf[i   ] = vector[my_local_ids[Nloc/2   + i ]];
	send_buf[i+nx] = vector[my_local_ids[3*Nloc/2 + i ]];
      }
  }
  else {
    for (i = 0; i < nx; i++) {
      send_buf[i     ] = vector[my_local_ids[Nloc/2 - nx + i ]];
    }
  }

  MPI_Send(send_buf, 3*nx - Partition->Nghost, MPI_DOUBLE, 
	   1 - Partition->mypid, type, MPI_COMM_WORLD);


  /* receive data */

  MPI_Wait(&request, &status);

  free(send_buf);
  return 1;
}

int user_update_ghost_nodes(double vector[], void *data)
{

  struct user_partition *Partition;
  int *my_local_ids, i, nx, Nloc, type = 7107;
  MPI_Request  request;
  MPI_Status status;
  double *send_buf;

  Partition = (struct user_partition *) data;
  if (Partition->nprocs == 1) return 1;

  my_local_ids = Partition->my_local_ids;
  
  nx = (int) sqrt( ((double) Partition->Nglobal) + .00001);
  Nloc = Partition->Nlocal;
  send_buf = (double *) malloc( (nx + 1)*sizeof(double));

  /* post receive */

  MPI_Irecv(&(vector[Nloc]), nx, MPI_DOUBLE, 1 - Partition->mypid, type, 
	    MPI_COMM_WORLD, &request);

  /* send data */

  if (Partition->mypid == 1) {
    for (i = 0; i < nx; i++) send_buf[i] = vector[my_local_ids[Nloc + i ]];
  }
  else {
    for (i = 0; i < nx; i++) send_buf[i] = vector[my_local_ids[Nloc - nx + i]];
  }

  MPI_Send(send_buf, nx, MPI_DOUBLE, 1 - Partition->mypid, type,
	   MPI_COMM_WORLD);


  /* receive data */

  MPI_Wait(&request, &status);
  free(send_buf);

  return 1;
}

#endif
