//@HEADER
// ************************************************************************
// 
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

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
#include "CLIP_interface.H"
#include "MpcLocal.H"
#include "CLIP_solver.hpp"

class SCLIP_solver
{
 public:
  SCLIP_solver(
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

  ~SCLIP_solver();
  void construct_H1();
  void determine_ownership();
  void construct_K_base();
  void CLIP_solver_init(int cdof_option, double solver_tol, int maxiter, 
			int atype, int ndim, int local_solver, int max_orthog,
			int prt_debug, int prt_summary, 
			int chk_sub_singularity, int krylov_method,
			int scale_option, int num_rigid_mode);
  void solve(double f[], double u[], int & number_iterations, 
	     int & SCLIP_status, int & max_added_corner);
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

  CLIP_solver *CS;
  Epetra_MpiComm *Comm;
  Epetra_CrsMatrix *ASub, *ConStandard;
  Epetra_MultiVector *Coords;
  Epetra_Map *SubMap, *StandardMap, *MpcLocalMap;
  Epetra_IntVector *NodalDofs;
  Epetra_Vector *uStand, *fStand, *uSub, *uLocal;
  Epetra_Import *ImporterSt2Sub, *ImporterStLam2Loc;
  Epetra_Export *ExporterSub2St;
  int maxdofnode, ndof, *H1, ndof_mine, nmpc_loc, MyPID;
  double *uvec, *fvec, *subvec, *locvec;
};
