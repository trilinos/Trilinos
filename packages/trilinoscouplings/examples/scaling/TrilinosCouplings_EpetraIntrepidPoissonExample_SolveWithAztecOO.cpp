// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TrilinosCouplings_EpetraIntrepidPoissonExample_SolveWithAztecOO.hpp"

// AztecOO solver
#include "AztecOO.h"

namespace TrilinosCouplings {
namespace EpetraIntrepidPoissonExample {

void
solveWithAztecOO (bool& converged,
		  int& numItersPerformed,
		  const double tol,
		  const int maxNumIters,
		  const int num_steps,
		  const Teuchos::RCP<Epetra_Vector>& X,
		  const Teuchos::RCP<Epetra_Operator>& A,
		  const Teuchos::RCP<Epetra_Vector>& B,
		  const Teuchos::RCP<Epetra_Operator>& M_right)
{
  // Set these in advance, so that if the AztecOO solver throws an
  // exception for some reason, these will have sensible values.
  converged = false;
  numItersPerformed = 0;

  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null () || X.is_null () || B.is_null (),
    std::invalid_argument, "solveWithBelos: The A, X, and B arguments must all "
    "be nonnull.");

  AztecOO solver;
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_kspace, 200);
  solver.SetAztecOption(AZ_conv, AZ_r0);
  solver.SetAztecOption(AZ_output, 10);
  solver.SetUserOperator(A.get());
  if (M_right != Teuchos::null)
    solver.SetPrecOperator(M_right.get());
  solver.SetRHS(B.get());

  // Enter "time step" loop -- we're really solving the same system repeatedly
  converged = true;
  numItersPerformed = 0;
  for (int step=0; step<num_steps; ++step) {

    // Set x
    X->PutScalar(0.0);

    // Reset problem
    solver.SetLHS(X.get());

    int result = solver.Iterate(maxNumIters, tol);

    converged = converged && (result == 0);
    numItersPerformed += solver.NumIters ();

  }
}

} // namespace EpetraIntrepidPoissonExample
} // namespace TrilinosCouplings
