/*****************************************************************************/
/* Copyright 2002, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************
 * Sample driver for Maxwell equation AMG solver in the ML package. The
 * software is tested by setting up a 2-dimensional uniform grid example on 
 * a square. For details about the problem at hand, please refer to file
 * ml_simple_max.c, of which this file is the C++ counterpart.
 *
 * This file shows how to use the class
 * ML_Epetra::MultiLevelPreconditioner to solve this formulation of the
 * Maxwell equations. The class takes care of building the node and edge
 * hierarchy, definining the Hiptmair smoother, and setting the coarse
 * solver. More information about MultiLevelPreconditioner can be found of the ML
 * user's guide.
 *****************************************************************************/

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"

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

#ifdef ML_MPI
#define COMMUNICATOR MPI_COMM_WORLD
#else
#define COMMUNICATOR AZ_NOT_MPI
#endif

/*****************************************************************************/
/* Function definitions.                                                     */
/*****************************************************************************/
extern void user_partition_edges(struct user_partition *,
				 struct user_partition *);
extern void user_partition_nodes(struct user_partition *Partition);
extern AZ_MATRIX   *user_Ke_build(struct user_partition *);
extern AZ_MATRIX   *user_Kn_build(struct user_partition *);
extern ML_Operator *user_T_build (struct user_partition *, 
                                  struct user_partition *, ML_Operator *, ML_Comm *);
int main(int argc, char *argv[])
{
  int    Nnodes=32*32;              /* Total number of nodes in the problem.*/
                                    /* 'Nnodes' must be a perfect square.   */

  struct       user_partition Edge_Partition = {NULL, NULL,0,0,NULL,0,0,0}, 
                                Node_Partition = {NULL, NULL,0,0,NULL,0,0,0};

  int          proc_config[AZ_PROC_SIZE];

  /* get processor information (id & # of procs) and set ML's printlevel. */

  // create communicators. In this example Epetra communicators are not
  // required (they are automatically defined in
  // ML_Operator2EpetraCrsMatrix()), but other codes may need to define
  // Epetra_MpiComm or Epetra_SerialComm.

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif
  AZ_set_proc_config(proc_config, COMMUNICATOR);
  ML_Comm * comm;
  ML_Comm_Create( &comm );

  // ================================= //
  // C R E A T E   P A R T I T I O N S //
  // ================================= //

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

  // ================================================= //
  // M A T R I X   C O N S T R U C T I O N   P H A S E //
  // ================================================= //
  
  // Matrix creation is here done as follows:
  // - as this example is generated from ml_simple_max.c, matrices are
  // first created as Aztec matrices. 
  // - then, they are converted to ML_Operator matrices
  // - finally, ML_Operator matrices are converted to Epetra_CrsMatrix
  // (derived from Epetra_RowMatrix).
  //
  // The first conversion simply defines some wrappers, so only an
  // additional amount of memory is really allocated for the new
  // ML_Operator's. Instead, Epetra_CrsMatrix's are created using
  // function ML_Operator2EpetraCrsMatrix(). After a call to
  // ML_Operator2EpetraCrsMatrix(), two matrices exists: the original
  // ML_Operator, and the Epetra matrix. This function supports
  // rectangular matrix as well (as required by Epetra_T).
  // build the matrices as Aztec matrices (as done in example
  // ml_simple_max.c).
  //
  // As the goal of this example is to show how
  // ML_Epetra::MultiLevelPreconditioner works with Epetra matrices, it
  // is here not important that some more memory is allocated in the
  // conversion. NOTE THAT ML_Epetra::MultiLevelPreconditioner DOES NOT
  // REQUIRE THE INPUT MATRICES TO BE DERIVED FROM ANY ML_OPERATOR OR
  // SIMILAR. ANY Epetra_RowMatrix CAN BE USED IN INPUT.
  
  AZ_MATRIX * AZ_Ke = user_Ke_build(&Edge_Partition);
  AZ_MATRIX * AZ_Kn = user_Kn_build(&Node_Partition);

  // convert (put wrappers) from Aztec matrices to ML_Operator's

  ML_Operator * ML_Ke, * ML_Kn, * ML_Tmat;

  ML_Ke = ML_Operator_Create( comm );
  ML_Kn = ML_Operator_Create( comm );

  AZ_convert_aztec_matrix_2ml_matrix(AZ_Ke,ML_Ke,proc_config);
  AZ_convert_aztec_matrix_2ml_matrix(AZ_Kn,ML_Kn,proc_config);

  ML_Tmat = user_T_build(&Edge_Partition, &Node_Partition, 
		      ML_Kn, comm);


  // convert to Epetra_CrsMatrix

  Epetra_CrsMatrix * Epetra_Kn, * Epetra_Ke, * Epetra_T;
  
  int MaxNumNonzeros;
  double CPUTime;

  ML_Operator2EpetraCrsMatrix(ML_Ke,Epetra_Ke,
			      MaxNumNonzeros,
			      true,CPUTime);

  ML_Operator2EpetraCrsMatrix(ML_Kn,
			      Epetra_Kn,MaxNumNonzeros,
			      true,CPUTime);

  ML_Operator2EpetraCrsMatrix(ML_Tmat,Epetra_T,MaxNumNonzeros,
			      true,CPUTime);  


  // ==================================================== //
  // S E T U P   O F    M L   P R E C O N D I T I O N E R //
  // ==================================================== //

  Teuchos::ParameterList MLList;
  ML_Epetra::SetDefaults("maxwell", MLList);
  
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("coarse: max size", 30);
  MLList.set("aggregation: threshold", 0.0);

  MLList.set("coarse: type", "Amesos_KLU");

  // this creates the multilevel hierarchy, set the smoother as
  // Hiptmair, prepare the coarse solver, etc...

  ML_Epetra::MultiLevelPreconditioner * MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(*Epetra_Ke, *Epetra_T, *Epetra_Kn,
					    MLList);

  // ========================================================= //
  // D E F I N I T I O N   O F   A Z T E C O O   P R O B L E M //
  // ========================================================= //

  // create left-hand side and right-hand side, and populate them with
  // random and constant values. Both vectors are defined on the domain
  // map of the edge matrix.
  // Epetra_Vectors can be created in View mode, to accept pointers to
  // double vectors.

  Epetra_Vector LHS(Epetra_Ke->DomainMap()); LHS.Random();
  Epetra_Vector RHS(Epetra_Ke->DomainMap()); RHS.PutScalar(1.0);
  
  // for AztecOO, we need an Epetra_LinearProblem
  Epetra_LinearProblem Problem(Epetra_Ke,&LHS,&RHS);
  // AztecOO Linear problem
  AztecOO solver(Problem);
  // set MLPrec as precondititoning operator for AztecOO linear problem
  solver.SetPrecOperator(MLPrec);

  // a few options for AztecOO
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(500, 1e-8);

  // =============== //
  // C L E A N   U P //
  // =============== //
  
  delete MLPrec;    // destroy phase prints out some information
  delete Epetra_Kn;
  delete Epetra_Ke;
  delete Epetra_T;
  
  ML_Operator_Destroy( &ML_Ke );
  ML_Operator_Destroy( &ML_Kn );
  ML_Comm_Destroy( &comm );

  if (Edge_Partition.my_local_ids != NULL) free(Edge_Partition.my_local_ids);
  if (Node_Partition.my_local_ids != NULL) free(Node_Partition.my_local_ids);
  if (Node_Partition.my_global_ids != NULL) free(Node_Partition.my_global_ids);
  if (Edge_Partition.my_global_ids != NULL) free(Edge_Partition.my_global_ids);
  if (Node_Partition.needed_external_ids != NULL) 
    free(Node_Partition.needed_external_ids);
  if (Edge_Partition.needed_external_ids != NULL) 
    free(Edge_Partition.needed_external_ids);

  if (AZ_Ke!= NULL) {
    AZ_free(AZ_Ke->bindx);
    AZ_free(AZ_Ke->val);
    AZ_free(AZ_Ke->data_org);
    AZ_matrix_destroy(&AZ_Ke);
  }
  if (AZ_Kn!= NULL) {
    AZ_free(AZ_Kn->bindx);
    AZ_free(AZ_Kn->val);
    AZ_free(AZ_Kn->data_org);
    AZ_matrix_destroy(&AZ_Kn);
  }

  ML_Operator_Destroy(&ML_Tmat);
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
  int    proc_config[AZ_PROC_SIZE];

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_input_update(NULL,&(Partition->Nlocal), &(Partition->my_global_ids),
		  proc_config, Partition->Nglobal, 1, AZ_linear);
  Partition->Nghost = 0;
}
/* Assign unknowns to processors */
void user_partition_edges(struct user_partition *Partition,
			  struct user_partition *NodePartition)
{
  int i;

  /* Edge partition is derived from node partition */

  Partition->Nlocal = NodePartition->Nlocal*2;
  Partition->my_global_ids = (int *) malloc(sizeof(int)*Partition->Nlocal);

  for (i = 0; i < NodePartition->Nlocal; i++)
    (Partition->my_global_ids)[i] = (NodePartition->my_global_ids)[i];
  for (i = 0; i < NodePartition->Nlocal; i++)
    (Partition->my_global_ids)[i+NodePartition->Nlocal] = 
               (NodePartition->my_global_ids)[i] + NodePartition->Nglobal;

}


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
