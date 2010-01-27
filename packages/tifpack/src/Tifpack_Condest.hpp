/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_CONDEST_HPP
#define TIFPACK_CONDEST_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_CondestType.hpp"
#include "Tifpack_Preconditioner.hpp"
#include <Teuchos_Ptr.hpp>

namespace Tifpack {

template<class Scalar,class LocalOrdinal,class GlobalOrdinal, class Node>
Scalar Condest(const Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>& TIFP,
		           const Tifpack::CondestType CT,
               const int MaxIters = 1550,
               const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& Tol = 1e-9,
               const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix_in = Teuchos::null)
{
  Scalar ConditionNumberEstimate = -1.0;
  Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > matrix = matrix_in;
  if (matrix_in == Teuchos::null) {
    matrix = TIFP.getMatrix().ptr();
  }

  if (CT == Tifpack::Cheap) {

    // Create a vector with all values equal to one
    Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Ones(TIFP.getDomainMap());
    Ones.putScalar(1.0);
    // Create the vector of results
    Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> OnesResult(TIFP.getRangeMap());
    // Compute the effect of the solve on the vector of ones
    TIFP.apply(Ones, OnesResult);
    // Make all values non-negative
    OnesResult.abs(OnesResult);
    // Get the maximum value across all processors
    Teuchos::Array<Scalar> cne(1,0);
    OnesResult.normInf(cne());
    ConditionNumberEstimate = cne[0];
  }
  else if (CT == Tifpack::CG) {

#ifdef HAVE_TIFPACK_AZTECOO
this code not yet converted!!!
    Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> LHS(IFP.getDomainMap());
    LHS.PutScalar(0.0);
    Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> RHS(IFP.getRangeMap());
    RHS.Random();
    Tpetra_LinearProblem Problem;
    Problem.SetOperator(matrix);
    Problem.SetLHS(&LHS);
    Problem.SetRHS(&RHS);

    AztecOO Solver(Problem);
    Solver.SetAztecOption(AZ_output,AZ_none);
    Solver.SetAztecOption(AZ_solver,AZ_cg_condnum);
    Solver.Iterate(MaxIters,Tol);

    const double* status = Solver.GetAztecStatus();
    ConditionNumberEstimate = status[AZ_condnum];
#endif

  } else if (CT == Tifpack::GMRES) {

#ifdef HAVE_TIFPACK_AZTECOO
this code not yet converted!!!
    Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> LHS(IFP.getDomainMap());
    LHS.PutScalar(0.0);
    Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> RHS(IFP.getRangeMap());
    RHS.Random();
    Tpetra_LinearProblem Problem;
    Problem.SetOperator(matrix);
    Problem.SetLHS(&LHS);
    Problem.SetRHS(&RHS);

    AztecOO Solver(Problem);
    Solver.SetAztecOption(AZ_solver,AZ_gmres_condnum);
    Solver.SetAztecOption(AZ_output,AZ_none);
    // the following can be problematic for large problems,
    // but any restart would destroy useful information about
    // the condition number.
    Solver.SetAztecOption(AZ_kspace,MaxIters);
    Solver.Iterate(MaxIters,Tol);

    const double* status = Solver.GetAztecStatus();
    ConditionNumberEstimate = status[AZ_condnum];
#endif
  }

  return(ConditionNumberEstimate);
}

}//namespace Tifpack

#endif // TIFPACK_CONDEST_HPP

