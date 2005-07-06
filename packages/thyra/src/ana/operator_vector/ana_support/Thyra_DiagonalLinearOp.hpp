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

#ifndef THYRA_DIAGONAL_LINEAR_OP_HPP
#define THYRA_DIAGONAL_LINEAR_OP_HPP

#include "Thyra_DiagonalLinearOpDecl.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_TestingTools.hpp" // ToDo: I need to have a better way to get eps()!

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DiagonalLinearOp<Scalar>::DiagonalLinearOp()
{}

template<class Scalar>
DiagonalLinearOp<Scalar>::DiagonalLinearOp(
  const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
  ,const Scalar                                           &gamma
  )
{
  initialize(diag,gamma);
}

template<class Scalar>
void DiagonalLinearOp<Scalar>::initialize(
  const Teuchos::RefCountPtr<const VectorBase<Scalar> >   &diag
  ,const Scalar                                           &gamma
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(diag.get()==NULL);
#endif
  diag_ = diag;
  gamma_ = gamma;
}

template<class Scalar>
void DiagonalLinearOp<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const VectorBase<Scalar> >  *diag
  ,Scalar                                          *gamma
  )
{
  if(diag) *diag = diag_;
  if(gamma) *gamma = gamma_;
  diag_ = Teuchos::null;
}

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DiagonalLinearOp<Scalar>::range() const
{
  return diag_->space();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DiagonalLinearOp<Scalar>::domain() const
{
  return diag_->space();
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DiagonalLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool DiagonalLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
  return true; // ToDo: Update this!
}

// Overridden from SingleRhsLinearOpBase

template<class Scalar>
void DiagonalLinearOp<Scalar>::apply(
  const ETransp                M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  ,const Scalar                beta
  ) const
{
  // RAB: 4/16/2005: Warning! this does not work if Scalar is a complex type
  // and M_trans==CONJTRANS!
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if( beta != ST::one() ) Vt_S( y, beta );
  ele_wise_prod( Scalar(alpha*gamma_), x, *diag_, y );
}

// Overridden from SingleScalarLinearOpWithSolveBase

template<class Scalar>
bool DiagonalLinearOp<Scalar>::solveSupported(ETransp M_trans) const
{
  return true; // ToDo: Update this!
}

// Overridden from SingleRhsLinearOpWithSolveBase

template<class Scalar>
SolveStatus<Scalar> DiagonalLinearOp<Scalar>::solve(
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
  assign(x,ST::zero());
  ele_wise_divide( Scalar(ST::one()/gamma_), b, *diag_, x );
  SolveStatus<Scalar> solveStatus;
  if( solveCriteria && solveCriteria->requestedTol!=SolveCriteria<Scalar>::defaultTolerance() ) {
    const ScalarMag eps = relErrSmallNumber<SMT::hasMachineParameters,ScalarMag>::smallNumber();
    if( solveCriteria->requestedTol <= eps ) {
      solveStatus.solveStatus = SOLVE_STATUS_CONVERGED;
    }
    else {
      solveStatus.solveStatus = SOLVE_STATUS_UNCONVERGED;
    }
    solveStatus.achievedTol = ScalarMag(2)*relErrSmallNumber<SMT::hasMachineParameters,ScalarMag>::smallNumber();
  }
  else {
    solveStatus.solveStatus = SOLVE_STATUS_UNKNOWN;
    solveStatus.achievedTol = SolveStatus<Scalar>::unknownTolerance();
  }
  solveStatus.numIterations   = 1;
  return solveStatus;
}

}	// end namespace Thyra

#endif	// THYRA_DIAGONAL_LINEAR_OP_HPP
