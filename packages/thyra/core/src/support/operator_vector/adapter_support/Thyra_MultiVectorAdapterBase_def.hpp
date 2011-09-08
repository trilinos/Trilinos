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

#ifndef THYRA_MULTI_VECTOR_ADAPTER_BASE_DEF_HPP
#define THYRA_MULTI_VECTOR_ADAPTER_BASE_DEF_HPP

#include "Thyra_MultiVectorAdapterBase_decl.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_ScalarProdBase.hpp"


namespace Thyra {


// Overridden functions from LinearOp


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
MultiVectorAdapterBase<Scalar>::range() const
{
  return rangeScalarProdVecSpc();
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
MultiVectorAdapterBase<Scalar>::domain() const
{
  return domainScalarProdVecSpc();
}


// Overridden protected functions from LinearOpBase


template<class Scalar>
bool MultiVectorAdapterBase<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  if (ScalarTraits<Scalar>::isComplex)
    return (M_trans == NOTRANS || M_trans == CONJTRANS);
  return true;
}


template<class Scalar>
void MultiVectorAdapterBase<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  //
  // Perform:
  //
  //   NOTRANS:    Y = beta*Y + alpha * M * Q_D * X
  //
  //   CONJTRANS:  Y = beta*Y + alpha * M^H * Q_R * X
  //
  // where T = Q_D * X or Q_R * X
  //
  RCP<const ScalarProdVectorSpaceBase<Scalar> > scalarProdVecSpc =
    ( real_trans(M_trans) == NOTRANS
      ? domainScalarProdVecSpc()
      : rangeScalarProdVecSpc() );
  RCP<const ScalarProdBase<Scalar> > scalarProd = scalarProdVecSpc->getScalarProd();
  if (scalarProd->isEuclidean()) {
    // Y = beta*Y + alpha * op(M) * X
    this->euclideanApply(M_trans, X, Y, alpha, beta);
  }
  else {
    // T = Q * X
    RCP<MultiVectorBase<Scalar> > T = createMembers(X.range(), X.domain());
    ::Thyra::apply(*scalarProd->getLinearOp(), NOTRANS, X, T.ptr());
    // Y = beta*Y + alpha * op(M) * T
    this->euclideanApply(M_trans, *T, Y, alpha, beta);
  }
}


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_ADAPTER_BASE_DEF_HPP
