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

#ifndef THYRA_EUCLIDEAN_LINEAR_OP_HPP
#define THYRA_EUCLIDEAN_LINEAR_OP_HPP

#include "Thyra_EuclideanLinearOpBaseDecl.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_ScalarProdBase.hpp"
#include "Thyra_AssertOp.hpp"

namespace Thyra {

// Virtual functions with default implementations

template<class RangeScalar, class DomainScalar>
void EuclideanLinearOpBase<RangeScalar,DomainScalar>::euclideanApplyTranspose(
  const EConj                            conj
  ,const MultiVectorBase<RangeScalar>    &X
  ,MultiVectorBase<DomainScalar>         *Y
  ,const DomainScalar                    alpha
  ,const DomainScalar                    beta
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"EuclideanLinearOpBase<"<<Teuchos::ScalarTraits<RangeScalar>::name()<<","<<Teuchos::ScalarTraits<RangeScalar>::name()<<">::applyTranspose(...): "
    "Error, the concrete subclass described as { " << this->description() << " } "
    " with this->applyTransposeSupports("<<toString(conj)<<")="<<this->applyTransposeSupports(conj)
    << " did not override this function and does not support transposes."
    );
}

// Overridden functions from OpBase */

template<class RangeScalar, class DomainScalar>
Teuchos::RCP<const VectorSpaceBase<RangeScalar> >
EuclideanLinearOpBase<RangeScalar,DomainScalar>::range() const
{
  return this->rangeScalarProdVecSpc();
}

template<class RangeScalar, class DomainScalar>
Teuchos::RCP<const VectorSpaceBase<DomainScalar> >
EuclideanLinearOpBase<RangeScalar,DomainScalar>::domain() const
{
  return this->domainScalarProdVecSpc();
}

// Overridden functions from LinearOpBase

template<class RangeScalar, class DomainScalar>
void EuclideanLinearOpBase<RangeScalar,DomainScalar>::apply(
  const EConj                            conj
  ,const MultiVectorBase<DomainScalar>   &X
  ,MultiVectorBase<RangeScalar>          *Y
  ,const RangeScalar                     alpha
  ,const RangeScalar                     beta
  ) const
{
  euclidean_apply_impl(conj,X,Y,alpha,beta);
}

template<class RangeScalar, class DomainScalar>
void EuclideanLinearOpBase<RangeScalar,DomainScalar>::applyTranspose(
  const EConj                            conj
  ,const MultiVectorBase<RangeScalar>    &X
  ,MultiVectorBase<DomainScalar>         *Y
  ,const DomainScalar                    alpha
  ,const DomainScalar                    beta
  ) const
{
  euclidean_applyTranspose_impl(conj,X,Y,alpha,beta);
}

// protected

template<class RangeScalar, class DomainScalar>
void EuclideanLinearOpBase<RangeScalar,DomainScalar>::euclidean_apply_impl(
  const EConj                            conj
  ,const MultiVectorBase<DomainScalar>   &X
  ,MultiVectorBase<RangeScalar>          *Y
  ,const RangeScalar                     alpha
  ,const RangeScalar                     beta
  ) const
{
  this->domainScalarProdVecSpc()->getScalarProd()->euclideanApply(
    *this, applyConjToTrans(conj), X, Teuchos::ptr(Y), alpha, beta);
}

template<class RangeScalar, class DomainScalar>
void EuclideanLinearOpBase<RangeScalar,DomainScalar>::euclidean_applyTranspose_impl(
  const EConj                            conj
  ,const MultiVectorBase<RangeScalar>    &X
  ,MultiVectorBase<DomainScalar>         *Y
  ,const DomainScalar                    alpha
  ,const DomainScalar                    beta
  ) const
{
  this->rangeScalarProdVecSpc()->getScalarProd()->euclideanApply(
    *this, applyTransposeConjToTrans(conj), X, Teuchos::ptr(Y), alpha, beta);
}

} // namespace Thyra

#endif // THYRA_EUCLIDEAN_LINEAR_OP_HPP
