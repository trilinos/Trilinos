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
