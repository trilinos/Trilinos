// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
