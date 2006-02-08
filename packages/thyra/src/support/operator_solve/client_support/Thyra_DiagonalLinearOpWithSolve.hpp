// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_DiagonalLinearOpWithSolveDecl.hpp"
#include "Thyra_DiagonalLinearOp.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_TestingTools.hpp" // ToDo: I need to have a better way to get eps()!

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DiagonalLinearOpWithSolve<Scalar>::DiagonalLinearOpWithSolve()
{}

template<class Scalar>
DiagonalLinearOpWithSolve<Scalar>::DiagonalLinearOpWithSolve(
  const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
  ,const Scalar                                           &gamma
  )
{
  initialize(diag,gamma);
}

// protected

// Overridden from SingleScalarLinearOpWithSolveBase

template<class Scalar>
bool DiagonalLinearOpWithSolve<Scalar>::solveSupportsTrans(ETransp M_trans) const
{
  return true; // ToDo: Update this!
}

template<class Scalar>
bool DiagonalLinearOpWithSolve<Scalar>::solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const
{
  return true;
}

// Overridden from SingleRhsLinearOpWithSolveBase

template<class Scalar>
SolveStatus<Scalar> DiagonalLinearOpWithSolve<Scalar>::solve(
  const ETransp                         M_trans
  ,const VectorBase<Scalar>             &b
  ,VectorBase<Scalar>                   *x
  ,const SolveCriteria<Scalar>          *solveCriteria
  ) const
{
  // RAB: 4/16/2005: Warning! this does not work if Scalar is a complex type
  // and M_trans==CONJTRANS!
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  typedef SolveCriteria<Scalar> SC;
  typedef SolveStatus<Scalar>   SS;
  assign(x,ST::zero());
  ele_wise_divide( Scalar(ST::one()/this->gamma()), b, *this->diag(), x );
  SS solveStatus;
  solveStatus.solveStatus
    = (solveCriteria && solveCriteria->solveTolType!=SOLVE_TOL_DEFAULT
       ? SOLVE_STATUS_CONVERGED : SOLVE_STATUS_UNKNOWN );
  solveStatus.achievedTol = SolveStatus<Scalar>::unknownTolerance();
  solveStatus.iterations = 1;
  return solveStatus;
}

}	// end namespace Thyra

#endif	// THYRA_DIAGONAL_LINEAR_OP_WITH_SOLVE_HPP
