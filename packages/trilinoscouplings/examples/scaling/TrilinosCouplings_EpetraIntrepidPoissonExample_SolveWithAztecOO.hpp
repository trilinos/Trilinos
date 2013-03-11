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
