// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TrilinosCouplings_IntrepidPoissonExample_SolveWithBelos.hpp"
#include "TrilinosCouplings_EpetraIntrepidPoissonExample.hpp"
#include <BelosEpetraAdapter.hpp>


namespace TrilinosCouplings {
namespace EpetraIntrepidPoissonExample {

void
solveWithBelos (bool& converged,
                int& numItersPerformed,
                const std::string& solverName,
                const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
                const int maxNumIters,
                const int num_steps,
                const Teuchos::RCP<multivector_type>& X,
                const Teuchos::RCP<const sparse_matrix_type>& A,
                const Teuchos::RCP<const multivector_type>& B,
                const Teuchos::RCP<operator_type>& M_left,
                const Teuchos::RCP<operator_type>& M_right)
{
  typedef multivector_type MV;
  typedef operator_type OP;

  // Create prec operator out of M (Apply->ApplyInverse)
  Teuchos::RCP<const operator_type> Mp_left = M_left;
  Teuchos::RCP<const operator_type> Mp_right = M_right;
  if (! M_left.is_null ()) {
    Mp_left = Teuchos::rcp (new Belos::EpetraPrecOp (M_left));
  }
  if (! M_right.is_null ()) {
    Mp_right = Teuchos::rcp (new Belos::EpetraPrecOp (M_right));
  }

  // Invoke the generic solve routine.
  IntrepidPoissonExample::solveWithBelos<ST, MV, OP> (converged, numItersPerformed,
                                                      solverName, tol, maxNumIters,
                                                      num_steps,
                                                      X, A, B, Mp_left, Mp_right);
}

} // namespace EpetraIntrepidPoissonExample
} // namespace TrilinosCouplings
