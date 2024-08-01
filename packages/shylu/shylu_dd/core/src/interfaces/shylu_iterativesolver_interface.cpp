// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "shylu_iterativesolver_interface_decl.hpp"
#include "shylu_iterativesolver_interface_def.hpp"

#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include "Tpetra_ConfigDefs.hpp"
#endif

namespace ShyLU {
  
  template class IterativeSolverInterface<Epetra_CrsMatrix, Epetra_MulitVector>;
  
#if defined(HAVE_SHYLU_TPETRA)

#if defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)
  template class IterativeSolverInterface<
    Tpetra::CrsMatrix<float, int, int>,
    Tpetra::MultiVector<float, int, int> >;
#endif // defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
  template class IterativeSolverInterface<
    Tpetra::CrsMatrix<double, int, int> ,
    Tpetra::MultiVector<double, int, int> >;
#endif // defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)

#endif // defined(HAVE_SHYLU_DDCORE_TPETRA)

}//end namespace Shylu


