/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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

#include "Thyra_AztecOOLinearOpWithSolve.hpp"

namespace Thyra {

// Constructors/initializers/accessors

AztecOOLinearOpWithSolve::AztecOOLinearOpWithSolve()
{}

// Overridden from LinearOpBase

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AztecOOLinearOpWithSolve::range() const
{
  return epetraOp_.range();
}

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
AztecOOLinearOpWithSolve::domain() const
{
  return epetraOp_.domain();
}

Teuchos::RefCountPtr<const LinearOpBase<double> >
AztecOOLinearOpWithSolve::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// protected

// Overridden from SingleScalarLinearOpBase

bool AztecOOLinearOpWithSolve::opSupported(ETransp M_trans) const
{
  return true; // ToDo: Determine if the Epetra_Operator supports adjoints or not!
}

// Overridden from SingleRhsLinearOpBase

void AztecOOLinearOpWithSolve::apply(
  const ETransp                M_trans
  ,const VectorBase<double>    &x
  ,VectorBase<double>          *y
  ,const double                alpha
  ,const double                beta
  ) const
{
  Thyra::apply( epetraOp_, M_trans, x, y, alpha, beta );
}

// Overridden from SingleScalarLinearOpWithSolveBase

bool AztecOOLinearOpWithSolve::solveSupportsTrans(ETransp M_trans) const
{
  return true; // ToDo: Determine if the solver supports adjoints or not!
}

bool AztecOOLinearOpWithSolve::solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const
{
  return true; // I am a direct solver so I should be able to do it all!
}

// Overridden from SingleRhsLinearOpWithSolveBase

SolveStatus<double>
AztecOOLinearOpWithSolve::solve(
  const ETransp                         M_trans
  ,const VectorBase<double>             &b
  ,VectorBase<double>                   *x
  ,const SolveCriteria<double>          *solveCriteria
  ) const
{
  TEST_FOR_EXCEPT(true); // ToDo: Fill this in!
  SolveStatus<double> solveStatus;
  solveStatus.solveStatus = SOLVE_STATUS_UNKNOWN;
  solveStatus.achievedTol = SolveStatus<double>::unknownTolerance();
  solveStatus.iterations = 1;
  return solveStatus;
}

}	// end namespace Thyra
