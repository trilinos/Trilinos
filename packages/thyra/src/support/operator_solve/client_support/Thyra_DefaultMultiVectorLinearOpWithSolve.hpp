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

#ifndef THYRA_MULTI_VECTOR_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_MULTI_VECTOR_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_DefaultMultiVectorLinearOpWithSolveDecl.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_dyn_cast.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultMultiVectorLinearOpWithSolve<Scalar>::DefaultMultiVectorLinearOpWithSolve()
{}


template<class Scalar>
void DefaultMultiVectorLinearOpWithSolve<Scalar>::initialize(
  const Teuchos::RCP<LinearOpWithSolveBase<Scalar> > &lows,
  const Teuchos::RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
  const Teuchos::RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
  )
{
  validateInitialize(lows,multiVecRange,multiVecDomain);
  lows_ = lows;
  multiVecRange_ = multiVecRange;
  multiVecDomain_ = multiVecDomain;
}


template<class Scalar>
void DefaultMultiVectorLinearOpWithSolve<Scalar>::initialize(
  const Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > &lows,
  const Teuchos::RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
  const Teuchos::RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
  )
{
  validateInitialize(lows,multiVecRange,multiVecDomain);
  lows_ = lows;
  multiVecRange_ = multiVecRange;
  multiVecDomain_ = multiVecDomain;
}


template<class Scalar>
Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
DefaultMultiVectorLinearOpWithSolve<Scalar>::getNonconstLinearOpWithSolve()
{
  return lows_.getNonconstObj();
}


template<class Scalar>
Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
DefaultMultiVectorLinearOpWithSolve<Scalar>::getLinearOpWithSolve() const
{
  return lows_.getConstObj();
}


template<class Scalar>
void DefaultMultiVectorLinearOpWithSolve<Scalar>::uninitialize()
{
  lows_.uninitialize();
  multiVecRange_ = Teuchos::null;
  multiVecDomain_ = Teuchos::null;
}


// Overridden from LinearOpBase


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultMultiVectorLinearOpWithSolve<Scalar>::range() const
{
  return multiVecRange_;
}


template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultMultiVectorLinearOpWithSolve<Scalar>::domain() const
{
  return multiVecDomain_;
}


template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultMultiVectorLinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null; // ToDo: Implement if needed ???
}


// protected


// Overridden from SingleScalarLinearOpBase


template<class Scalar>
bool DefaultMultiVectorLinearOpWithSolve<Scalar>::opSupported(
  EOpTransp M_trans
  ) const
{
  return Thyra::opSupported(*lows_.getConstObj(),M_trans);
}


// Overridden from SingleRhsLinearOpBase


template<class Scalar>
void DefaultMultiVectorLinearOpWithSolve<Scalar>::apply(
  const EOpTransp M_trans,
  const VectorBase<Scalar> &x,
  VectorBase<Scalar> *y,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  typedef DefaultMultiVectorProductVector<Scalar> MVPV;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==y);
#endif
 
  RCP<const MultiVectorBase<Scalar> >
    X = dyn_cast<const MVPV>(x).getMultiVector().assert_not_null();
  RCP<MultiVectorBase<Scalar> >
    Y = dyn_cast<MVPV>(*y).getNonconstMultiVector().assert_not_null();

  Thyra::apply( *lows_.getConstObj(), M_trans, *X, &*Y, alpha, beta );

}


// Overridden from SingleScalarLinearOpWithSolveBase


template<class Scalar>
bool DefaultMultiVectorLinearOpWithSolve<Scalar>::solveSupportsTrans(
  EOpTransp M_trans
  ) const
{
  return Thyra::solveSupports(*lows_.getConstObj(),M_trans);
}


template<class Scalar>
bool DefaultMultiVectorLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureType(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  return Thyra::solveSupportsSolveMeasureType(
    *lows_.getConstObj(),M_trans,solveMeasureType);
}


// Overridden from SingleRhsLinearOpWithSolveBase


template<class Scalar>
SolveStatus<Scalar>
DefaultMultiVectorLinearOpWithSolve<Scalar>::solve(
  const EOpTransp M_trans,
  const VectorBase<Scalar> &b,
  VectorBase<Scalar> *x,
  const SolveCriteria<Scalar> *solveCriteria
  ) const
{

  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  typedef DefaultMultiVectorProductVector<Scalar> MVPV;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==x);
#endif
 
  RCP<const MultiVectorBase<Scalar> >
    B = dyn_cast<const MVPV>(b).getMultiVector().assert_not_null();
  RCP<MultiVectorBase<Scalar> >
    X = dyn_cast<MVPV>(*x).getNonconstMultiVector().assert_not_null();

  return Thyra::solve(
    *lows_.getConstObj(), M_trans,
    *B, &*X, solveCriteria
    );

}


// private


template<class Scalar>
void DefaultMultiVectorLinearOpWithSolve<Scalar>::validateInitialize(
  const Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > &lows,
  const Teuchos::RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
  const Teuchos::RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(lows));
  TEST_FOR_EXCEPT(is_null(multiVecRange));
  TEST_FOR_EXCEPT(is_null(multiVecDomain));
  TEST_FOR_EXCEPT( multiVecRange->numBlocks() != multiVecDomain->numBlocks() );
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorLinearOpWithSolve<Scalar>::initialize(lows,multiVecRange,multiVecDomain)",
    *lows->range(), *multiVecRange->getBlock(0) );
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorLinearOpWithSolve<Scalar>::initialize(lows,multiVecRange,multiVecDomain)",
    *lows->domain(), *multiVecDomain->getBlock(0) );
#endif
}


}	// end namespace Thyra


#endif	// THYRA_MULTI_VECTOR_LINEAR_OP_WITH_SOLVE_HPP
