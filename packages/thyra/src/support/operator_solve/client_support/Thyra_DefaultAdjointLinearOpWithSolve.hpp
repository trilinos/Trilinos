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


#ifndef THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_HPP


#include "Thyra_DefaultAdjointLinearOpWithSolveDecl.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_SingleScalarLinearOpWithSolveBase.hpp"


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

  
// Overridden from SingleScalarLinearOpBase


template<class Scalar>
bool DefaultAdjointLinearOpWithSolve<Scalar>::opSupported(EOpTransp M_trans) const
{
  return Thyra::opSupported(*lows_(), trans_trans(transp_, M_trans));
}


template<class Scalar>
void DefaultAdjointLinearOpWithSolve<Scalar>::apply(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply( *lows_(), trans_trans(transp_, M_trans),
    X, Y, alpha, beta );
}


// Overridden from SingleScalarLinearOpWithSolveBase


template<class Scalar>
bool DefaultAdjointLinearOpWithSolve<Scalar>::solveSupportsTrans(EOpTransp M_trans) const
{
  return Thyra::solveSupports(*lows_(), trans_trans(transp_, M_trans));
}


template<class Scalar>
bool DefaultAdjointLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureType(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const
{
  return Thyra::solveSupportsSolveMeasureType(*lows_(),
    trans_trans(transp_, M_trans), solveMeasureType );
}


// Overridden from SingleScalarLinearOpWithSolveBase


template<class Scalar>
void DefaultAdjointLinearOpWithSolve<Scalar>::solve(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &B,
  MultiVectorBase<Scalar> *X,
  const int numBlocks,
  const BlockSolveCriteria<Scalar> blockSolveCriteria[],
  SolveStatus<Scalar> blockSolveStatus[]
  ) const
{
  return Thyra::solve( *lows_(), trans_trans(transp_, M_trans),
    B, X, numBlocks, blockSolveCriteria, blockSolveStatus );
}


} // namespace Thyra


#endif // THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_HPP
