// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_DEF_HPP
#define THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_DEF_HPP


#include "Thyra_DefaultAdjointLinearOpWithSolve_decl.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultAdjointLinearOpWithSolve<Scalar>::DefaultAdjointLinearOpWithSolve()
  : transp_(NOTRANS)
{}


template<class Scalar>
void DefaultAdjointLinearOpWithSolve<Scalar>::initialize(
  const RCP<LinearOpWithSolveBase<Scalar> > &lows,
  const EOpTransp transp )
{
  lows_ = lows;
  transp_ = transp;
}


template<class Scalar>
void DefaultAdjointLinearOpWithSolve<Scalar>::initialize(
  const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
  const EOpTransp transp )
{
  lows_ = lows;
  transp_ = transp;
}


template<class Scalar>
const RCP<LinearOpWithSolveBase<Scalar> >
DefaultAdjointLinearOpWithSolve<Scalar>::getNonconstOp()
{
  return lows_.getNonconstObj();
}


template<class Scalar>
const RCP<const LinearOpWithSolveBase<Scalar> >
DefaultAdjointLinearOpWithSolve<Scalar>::getOp() const
{
  return lows_.getConstObj();
}


// Overridden from LinearOpBase */


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultAdjointLinearOpWithSolve<Scalar>::range() const
{
  return ( real_trans(transp_) == NOTRANS
    ? lows_()->range() : lows_()->domain() );
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultAdjointLinearOpWithSolve<Scalar>::domain() const
{
  return ( real_trans(transp_) == NOTRANS
    ? lows_()->domain() : lows_()->range() );
}


// protected

  
// Overridden from LinearOpBase


template<class Scalar>
bool DefaultAdjointLinearOpWithSolve<Scalar>::opSupportedImpl(
  EOpTransp M_trans) const
{
  return Thyra::opSupported(*lows_(), trans_trans(transp_, M_trans));
}


template<class Scalar>
void DefaultAdjointLinearOpWithSolve<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply( *lows_(), trans_trans(transp_, M_trans),
    X, Y, alpha, beta );
}


// Overridden from LinearOpWithSolveBase


template<class Scalar>
bool DefaultAdjointLinearOpWithSolve<Scalar>::solveSupportsImpl(EOpTransp M_trans) const
{
  return Thyra::solveSupports(*lows_(), trans_trans(transp_, M_trans));
}


template<class Scalar>
bool DefaultAdjointLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const
{
  return Thyra::solveSupportsSolveMeasureType(*lows_(),
    trans_trans(transp_, M_trans), solveMeasureType );
}



template<class Scalar>
SolveStatus<Scalar>
DefaultAdjointLinearOpWithSolve<Scalar>::solveImpl(
  const EOpTransp transp,
  const MultiVectorBase<Scalar> &B,
  const Ptr<MultiVectorBase<Scalar> > &X,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  ) const
{
  return Thyra::solve( *lows_(), trans_trans(transp_, transp),
    B, X, solveCriteria );
}


} // namespace Thyra


#endif // THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_DEF_HPP
