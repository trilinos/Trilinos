/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 

#ifndef ML_NOX_PRECONDITIONER_H
#define ML_NOX_PRECONDITIONER_H

#include <vector>

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"

// NOX-headers
#include "NOX_Epetra_Interface.H"

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_agg_VBMETIS.h"
// a matrixfree level
#include "ml_nox_matrixfreelevel.H"
// a nonlinear level
#include "ml_nox_nonlinearlevel.H"
// this is the wrapper for the Epetra_CrsGraph used by ml in matrix-free method
// it is derived from Epetra_RowMatrix so the ml - Epetra_RowMatrix - wrapper will take it
#include "ml_Epetra_wrap_CrsGraph_as_RowMatrix.H"



// interface to user code on fine level
#include "ml_nox_fineinterface.H"



// ---------- Forward Declarations ----------
class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;
class Nox_Problem_Interface;

namespace ML_NOX
{
class  ML_Nox_Preconditioner : public virtual Epetra_Operator, 
                               public NOX::EpetraNew::Interface::Preconditioner
{
public:
  // constructor
  explicit ML_Nox_Preconditioner(ML_NOX::Ml_Nox_Fineinterface&    interface,
                                 bool                             matrixfree,
                                 bool                             matfreelev0,
                                 double                           fd_alpha,
                                 double                           fd_beta,
                                 bool                             fd_centered,
                                 bool                             linearPrec,
                                 bool                             isnlnCG,
                                 int                              nitersCG,
                                 double                           FAS_normF,      
                                 double                           FAS_nupdate,
                                 int                              FAS_prefinesmooth,    
                                 int                              FAS_presmooth,  
                                 int                              FAS_postsmooth,
                                 int                              FAS_postfinesmooth, 
                                 int                              FAS_maxcycle,   
                                 int                              N_levels,
                                 int                              ml_printlevel,
                                 int                              ml_numPDE,
                                 int                              ml_dim_nullsp,
                                 string                           coarsentype,
                                 int                              nnodeperagg,
                                 string                           smoothertype,
                                 string                           finesmoothertype,
                                 int                              nsmooth[],
                                 string                           coarsesolve,
                                 int                              maxcoarsesize,
                                 int                              offset,
                                 Epetra_Map&                      dm, 
                                 Epetra_Map&                      rm, 
                                 Epetra_Comm&                     comm);

  // destructor
  ~ML_Nox_Preconditioner();

  // derived from Epetra_Operator
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  // derived from Epetra_Operator
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { cout << "ML_Nox_Preconditioner::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";throw -1;return(-1);}
    
    
  int SetUseTranspose(bool UseTranspose);

  double NormInf() const {return(-1.0);}

  char* Label() const {return(&(label_[0]));}

  bool UseTranspose() const {return(usetranspose_);}

  bool HasNormInf() const {return(false);}

  const Epetra_Comm & Comm() const {return(comm_);}

  const Epetra_Map & OperatorDomainMap() const {return(DomainMap_);}

  const Epetra_Map & OperatorRangeMap() const {return(RangeMap_);}
  
  const bool isinit() const {return(isinit_);}
  
  const bool setinit(bool flag) {isinit_=flag; return(true);}

  int getoffset() {return(offset_newPrec_);} 
  
  bool set_nox_solver(NOX::Solver::Manager* noxsolver) {noxsolver_ = noxsolver; return true;} 
  
  // derived from NOX::EpetraNew::Interface::Preconditioner
  virtual bool computePreconditioner(const Epetra_Vector& x,
				     NOX::Parameter::List* precParams = 0);
  
  // nonlinear multigrid FAS-solver
  int solve() ;

private:
  // don't want copy-ctor and = operator
  ML_Nox_Preconditioner(ML_Nox_Preconditioner& old);
  ML_Nox_Preconditioner operator = (const ML_Nox_Preconditioner& old);
  
  // Compute Preconditioner functions, 
  // depending on type of preconditioner and whether the Jacobian was supplied or not
  bool compPrec(const Epetra_Vector& x);
  bool ML_Nox_compute_Jacobian_Linearpreconditioner(const Epetra_Vector& xfine);      
  bool ML_Nox_compute_Matrixfree_Linearpreconditioner(const Epetra_Vector& xfine);    
  bool ML_Nox_compute_Jacobian_Nonlinearpreconditioner(const Epetra_Vector& xfine);   
  bool ML_Nox_compute_Matrixfree_Nonlinearpreconditioner(const Epetra_Vector& xfine); 
  Epetra_CrsMatrix* ML_Nox_computeFineLevelJacobian(const Epetra_Vector& x);

  // utils
  Epetra_CrsGraph* deepcopy_graph(const Epetra_CrsGraph* oldgraph);
  
  // ApplyInverse functions, 
  // depending on type of preconditioner
  int ML_Nox_ApplyInverse_Linear(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  int ML_Nox_ApplyInverse_NonLinear(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  // nonlinear multigrid functions
  //the FAS-preconditioner
  bool ML_Nox_FAS_cycle(Epetra_Vector* f, Epetra_Vector* x, int level, bool* converged) const;
  // the FAS-solver
  bool ML_Nox_FAS_cycle1(Epetra_Vector* f, Epetra_Vector* x, int level, bool* converged);

  // destructor of vector of matfreelevels
  bool destroy_matfreelevels(ML_NOX::ML_Nox_MatrixfreeLevel*** m, int nlevel);
  // destructor of vector of nonlinearlevels
  bool destroy_nonlinearlevels(ML_NOX::ML_Nox_NonlinearLevel*** nln, int nlevel);

private:
  // ---general data---
  char*                                 label_;             // the label of this class
  bool                                  usetranspose_;     // use transpose of operator (not supported)
  bool                                  isinit_;           // flag whether preconditioner is valid
  bool                                  ismatrixfree_;     // use Graph or Jacobian in building the preconditioner
  bool                                  matfreelev0_;
  bool                                  islinearPrec_;     // do standard linear ML-SA-Preconditioner (not matrixfree!)) 
  bool                                  isnlnCG_;          // switch between nlnCG and Newton's method

  Epetra_Map                            DomainMap_;        // the operators domain map
  Epetra_Map                            RangeMap_;         // the operators range map
  Epetra_Comm&                          comm_;             // the communicator;
  ML_NOX::Ml_Nox_Fineinterface&         interface_;        // the user interface
  bool                                  destroyfineJac_;   // flag to indicate this class destroys the fineJac_
  Epetra_CrsMatrix*                     fineJac_;          // fine grid matrix (only used when ismatrixfree_=false)
  int                                   offset_newPrec_;   // how often to recompute the precond from scratch
  int                                   ncalls_NewPrec_;   // how often computePreconditioner is called

  // ---data related to the ML-hierarchy---
  ML*                                   ml_;               // the ml handle
  ML_Aggregate*                         ag_;               // the aggregate handle
  int                                   ml_N_levels_;      // max. number of levels
  int                                   ml_printlevel_;    // ML's printlevel (0 - 10)
  int                                   ml_numPDE_;        // number of pdes per node
  int                                   ml_dim_nullsp_;    // dimension of nullspace
  string                                ml_coarsentype_;   // type of coarsening algorithm: Uncoupled MIS METIS VBMETIS
  int                                   ml_nnodeperagg_;   // objective number of nodes per aggregate with coarsentypes METIS & VBMETIS
  string                                ml_smoothertype_;  // type of smoothing algorithm: SGS BSGS Jacobi
  string                                ml_fsmoothertype_; // type of smoothing algorithm: SGS BSGS Jacobi on level 0
  int*                                  ml_nsmooth_;       // number of smoothing steps for each level
  string                                ml_coarsesolve_;   // type of coarse solver: SGS BSGS AmesosKLU Jacobi
  int                                   ml_maxcoarsesize_; // max. coarse grid size
  int                                   ml_nlevel_;        // number of level generated by ML_Gen_MGHierarchy_UsingAggregation
  int                                   ml_coarsestlev_;   // index of the coarsest level;
  int                                   ml_nblocks_;
  int*                                  ml_blocks_;        // block information for Block-smoothers and aggregation scheme VBMETIS
  int*                                  ml_block_pde_;     // information which pde each row belongs to

  // ---the linear preconditioner---
  ML_Epetra::MultiLevelOperator*        ml_linPrec_;       // the linear preconditioner class

  // ---matrix-free method--- 
  const Epetra_CrsGraph*                fineGraph_;        // the fine grid graph (only used when ismatrixfree_=true)
  ML_Epetra::CrsGraphWrapper*           ml_graphwrap_;     // in the matrix-free case this wraps the input fine level graph as an Epetra_RowMatrix
  int                                   nmatfreelevel_;    // number of matfreelevels, should be equal to ml_nlevel_
  ML_NOX::ML_Nox_MatrixfreeLevel**      ml_matfreelevel_;  // vector of matrix-free levels (including fine level)
  double                                fd_alpha_;         // FiniteDifferencing parameters (for matrixfree only)
  double                                fd_beta_;
  bool                                  fd_centered_;      // use centered finite differencing

  // --- data special to the nonlinear preconditioner ---
  int                                   n_nlnlevel_;
  ML_NOX::ML_Nox_NonlinearLevel**       nlnLevel_;         // ptr to vector of coarse interfaces
  double                                FAS_normF_;       // coarse NOX's convergence criteria
  double                                FAS_nupdate_;
  int                                   FAS_maxiter_;
  int                                   FAS_prefinesmooth_;
  int                                   FAS_presmooth_;
  int                                   FAS_postsmooth_;
  int                                   FAS_postfinesmooth_;
  int                                   FAS_maxcycle_; 
  NOX::Solver::Manager*                 noxsolver_;
  int                                   nitersCG_;
}; // class  ML_Nox_Preconditioner : public virtual Epetra_Operator
}  // namespace ML_NOX


#endif // ML_NOX_PRECONDITIONER_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
