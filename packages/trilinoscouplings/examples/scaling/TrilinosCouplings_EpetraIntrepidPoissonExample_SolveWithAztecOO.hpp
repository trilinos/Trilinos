// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_IntrepidPoissonExample_SolveWithAztecOO_hpp
#define __TrilinosCouplings_IntrepidPoissonExample_SolveWithAztecOO_hpp

/// \file TrilinosCouplings_IntrepidPoissonExample_SolveWithAztecOO.hpp
/// \brief AztecOO solver for the Intrepid Poisson test problem example.

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

namespace TrilinosCouplings {
namespace EpetraIntrepidPoissonExample {

/// \brief Solve the linear system(s) AX=B with AztecOO.
///
/// This interface will change in the future to accept the name of the
/// Aztec solver to use.  For now, the solver is hard-coded to
/// CG.
///
/// \param converged [out] Whether Aztec reported that the iterative
///   method converged (solved the linear system to the desired
///   tolerance).
///
/// \param numItersPerformed [out] Number of iterations that the Aztec
///   solver performed.
///
/// \param tol [in] Convergence tolerance for the iterative method.
///   The meaning of this depends on the particular iterative method.
///
/// \param maxNumIters [in] Maximum number of iterations that the
///   iterative method should perform, regardless of whether it
///   converged.
///
/// \param num_steps [in] Number of "time steps", i.e., the number of
//    times the solver is called in a fake time-step loop.
///
/// \param X [in/out] On input: the initial guess(es) for the iterative
///   method.  On output: the computed approximate solution.
///
/// \param A [in] The matrix in the linear system(s) AX=B to solve.
///
/// \param B [in] The right-hand side in the linear system AX=B to solve.
///
/// \param M_right [in] If nonnull, a right preconditioner that the
///   iterative method may use.  If null, the iterative method will
///   not use a right preconditioner.
void
solveWithAztecOO (bool& converged,
		  int& numItersPerformed,
		  const double tol,
		  const int maxNumIters,
		  const int num_steps,
		  const Teuchos::RCP<Epetra_Vector>& X,
		  const Teuchos::RCP<Epetra_Operator>& A,
		  const Teuchos::RCP<Epetra_Vector>& B,
		  const Teuchos::RCP<Epetra_Operator>& M_right);

} // namespace EpetraIntrepidPoissonExample
} // namespace TrilinosCouplings

#endif // __TrilinosCouplings_IntrepidPoissonExample_SolveWithAztecOO_hpp
