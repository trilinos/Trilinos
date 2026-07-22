// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_MATRIXFREELEVEL_H
#define ML_NOX_MATRIXFREELEVEL_H

#include <vector>

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
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
#include "Epetra_MapColoring.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX-headers
#include <NOX.H>
#include <NOX_Epetra.H> 
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"
#include "NOX_EpetraNew_FiniteDifference.H"
#include "NOX_EpetraNew_FiniteDifferenceColoring.H"

// ML-headers
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_nox_fineinterface.H"    // interface to user code on fine level
#include "ml_nox_coarseinterface.H"  // interface to user code on coarse levels
#include "ml_nox_preconditioner_utils.H"  // interface to user code on coarse levels

// the graph coloring class 
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

// ---------- Forward Declarations ----------
class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;
class Nox_Problem_Interface;
class Nox_CoarseProblem_Interface;
class Epetra_MapColoring;
class Epetra_IntVector;

namespace ML_NOX
{
class  ML_Nox_MatrixfreeLevel 
{

public:
  // constructor
  explicit ML_Nox_MatrixfreeLevel(int level, int nlevel, int plevel, ML* ml, ML_Aggregate* ag, 
                                  Epetra_CrsMatrix** P, ML_NOX::Ml_Nox_Fineinterface& interface, 
                                  const Epetra_Comm& comm, const Epetra_Vector& xfine,
                                  double fd_alpha, double fd_beta, bool fd_centered, 
                                  bool isDiagonalOnly, int bsize);

  // a reconstructor
  bool recreateLevel(int level, int nlevel, int plevel, ML* ml, ML_Aggregate* ag,
                     Epetra_CrsMatrix** P, ML_NOX::Ml_Nox_Fineinterface& interface,
                     const Epetra_Comm& comm, const Epetra_Vector& xfine);
  // destructor
  ~ML_Nox_MatrixfreeLevel();
  
  // return ptr to underlying Epetra_CrsMatrix
  Epetra_CrsMatrix* getunderlyingMatrix() { return A_; }
  
  // return ptr to coarse interface
  ML_NOX::Nox_CoarseProblem_Interface* getCoarseinterface() { return coarseinterface_; }
  
  // set the ptr to the P-hierarchy in the coarseinterface_
  bool setP(Epetra_CrsMatrix** P) { if (coarseinterface_) return coarseinterface_->setP(P); 
                                    else return false; }

  // destroy the P-hierarchy in the coarse interface
  bool destroyP() { if (coarseinterface_) return coarseinterface_->destroyP(); 
                    else                  return false; }

private:
  // don't want copy-ctor and = operator
  ML_Nox_MatrixfreeLevel(ML_NOX::ML_Nox_MatrixfreeLevel& old);
  ML_Nox_MatrixfreeLevel operator = (const ML_NOX::ML_Nox_MatrixfreeLevel& old);
  
   // destroy the FD object
  bool destroyFD() { if (FD_) delete FD_; FD_ = 0; A_ = 0; graph_ = 0; return true;}
  // compare pointwise the global Indices in 2 graphs, returns false if not equal
  bool compare_graphs(const Epetra_CrsGraph* newgraph, const Epetra_CrsGraph* oldgraph);
  
  int                                       level_;
  int                                       nlevel_;
  int                                       ml_printlevel_;
  ML*                                       ml_;
  ML_Aggregate*                             ag_;
  ML_NOX::Ml_Nox_Fineinterface&             fineinterface_;
  const Epetra_Comm&                        comm_;
  ML_NOX::Nox_CoarseProblem_Interface*      coarseinterface_;
  Epetra_CrsGraph*                          graph_;
  Epetra_CrsMatrix*                         A_;
  //EpetraExt::CrsGraph_MapColoring*          MapColoring_;
  Epetra_MapColoring*                       colorMap_;
  EpetraExt::CrsGraph_MapColoringIndex*     colorMapIndex_;
  vector<Epetra_IntVector>*                 colorcolumns_;
  double                                    fd_alpha_;        // FiniteDifferencing parameters (for matrixfree only)
  double                                    fd_beta_;
  bool                                      fd_centered_;
  bool                                      isDiagonalOnly_;
  int                                       bsize_;
#if 1
  NOX::EpetraNew::FiniteDifferenceColoring* FD_;
#else
  NOX::EpetraNew::FiniteDifference*         FD_;
#endif
  
}; // class  ML_Nox_MatrixfreeLevel 
}  // namespace ML_NOX


#endif // ML_NOX_MATRIXFREELEVEL_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
