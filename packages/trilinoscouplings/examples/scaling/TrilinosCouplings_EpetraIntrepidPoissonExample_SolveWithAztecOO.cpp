// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
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
