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

#include "TrilinosCouplings_IntrepidPoissonExample_SolveWithBelos.hpp"
#include "TrilinosCouplings_EpetraIntrepidPoissonExample.hpp"
#include <BelosEpetraAdapter.hpp>


namespace TrilinosCouplings {
namespace EpetraIntrepidPoissonExample {

void
solveWithBelos (bool& converged,
                int& numItersPerformed,
                const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
                const int maxNumIters,
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
  if (M_left != Teuchos::null)
    Mp_left = Teuchos::rcp(new Belos::EpetraPrecOp(M_left));
  if (M_right != Teuchos::null)
    Mp_right = Teuchos::rcp(new Belos::EpetraPrecOp(M_right));

  // Invoke the generic solve routine.
  IntrepidPoissonExample::solveWithBelos<ST, MV, OP> (converged, numItersPerformed, tol, maxNumIters, X, A, B, Mp_left, Mp_right);
}

} // namespace EpetraIntrepidPoissonExample
} // namespace TrilinosCouplings
