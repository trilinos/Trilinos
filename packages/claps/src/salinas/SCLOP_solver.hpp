#include <mpi.h>
#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "CLOP_interface.H"
#include "MpcLocal.H"
#include "CLOP_solver.hpp"

class SCLOP_solver
{
 public:
  SCLOP_solver(
  int nnode_,               // number of nodes in superstructure
  int nelem_,               // number of elements in superstructure
  int neq_,                 // number of of active dofs in superstructure
  int nadj_proc_,           // number of adjacent processors
  int ndofs_deleted_,       // number of deleted degrees of freedom
  int NumberMpc_,           // number of constraint equations
  const int E1_[],          // nodes (local numbering) for elements
  const int E2_[],          // pointer array for E1_
  const int adj_proc_[],    // adjacent processors
  const int H1_global_[],   // nodes on boundary (global numbering)
  const int H2_[],          // pointer array for H1_global
  const MPI_Comm *mpicomm_, // pointer to MPI communicator
  const int gnn_[],         // global node numbers in superstructure
  const unsigned short *dofmap_on_node_, // codes for active nodal dofs
  const double x_[],        // x-coordinates of nodes
  const double y_[],        // y-coordinates of nodes
  const double z_[],        // z-coordinates of nodes
  const int rowbeg_KU_[],   // beginning of rows      in upper t sparse matrix
  const int colidx_KU_[],   // nonzero column numbers in upper t sparse matrix
  const double KU_[],       // nonzero entries        in upper t sparse matrix
  const int mapdof2row_[],  // a degree of freedom map
  const int bc_node_ids_[], // array of deleted dof node numbers
  const int bc_node_dof_[], // array of deleted dof local dof numbers
  const MpcLocal* MpcVector_);

  ~SCLOP_solver();
  void construct_H1();
  void determine_ownership();
  void construct_K_base();
  void CLOP_solver_init(int overlap, double solver_tol, int maxiter, 
			int atype, int ndim, int local_solver, int max_orthog,
			int prt_debug, int prt_summary);
  void solve(double f[], double u[], int & number_iterations, 
	     int & SCLOP_status);
  void MpcForces( double *cvals);
  void EPmat_datfile(Epetra_CrsMatrix* A, char fname[]);
 private: // variables
  int nnode, nelem, neq, nadj_proc, ndofs_deleted, NumberMpc;
  const int *E1, *E2, *adj_proc, *H1_global, *H2, *gnn, *rowbeg_KU, *colidx_KU;
  const int *mapdof2row, *bc_node_ids, *bc_node_dof;
  const unsigned short *dofmap_on_node;
  const double *x, *y, *z, *KU;
  const MpcLocal *MpcVector;
  MPI_Comm mpicomm;

  CLOP_solver *CS;
  Epetra_MpiComm *Comm;
  Epetra_CrsMatrix *AStandard, *ConStandard;
  Epetra_MultiVector *CStandard;
  Epetra_Map *SubMap, *StandardMap, *MpcLocalMap;
  Epetra_IntVector *LDStandard, *GNStandard;
  Epetra_Vector *uStand, *fStand, *uSub, *uLocal;
  Epetra_Import *ImporterSt2Sub, *ImporterStLam2Loc;
  Epetra_Export *ExporterSub2St;
  int maxdofnode, ndof, *H1, ndof_mine, nmpc_loc, MyPID;
  double *uvec, *fvec, *subvec, *locvec;
};
