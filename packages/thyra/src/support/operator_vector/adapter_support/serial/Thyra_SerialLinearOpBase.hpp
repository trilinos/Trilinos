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

#ifndef THYRA_SERIAL_LINEAR_OP_BASE_HPP
#define THYRA_SERIAL_LINEAR_OP_BASE_HPP

#include "Thyra_SerialLinearOpBaseDecl.hpp"
#include "Thyra_EuclideanLinearOpBase.hpp"
#include "Thyra_DefaultSerialVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

namespace Thyra {

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RCP< const ScalarProdVectorSpaceBase<Scalar> >
SerialLinearOpBase<Scalar>::rangeScalarProdVecSpc() const
{
  return range_;
}

template<class Scalar>
Teuchos::RCP< const ScalarProdVectorSpaceBase<Scalar> >
SerialLinearOpBase<Scalar>::domainScalarProdVecSpc() const
{
  return domain_;
}

template<class Scalar>
void SerialLinearOpBase<Scalar>::euclideanApply(
  const EOpTransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES("SerialLinearOpBase<Scalar>::apply()",*this,M_trans,X,Y);
#endif
  const ConstDetachedMultiVectorView<Scalar>         X_ev(X);
  const DetachedMultiVectorView<Scalar>  Y_ev(*Y);
  this->euclideanApply(M_trans,X_ev.smv(),&Y_ev.smv(),alpha,beta);
}

// protected

// Constructors/initializers/accessors

template<class Scalar>
SerialLinearOpBase<Scalar>::SerialLinearOpBase()
  :forceUnitStride_(true)
{}

template<class Scalar>
void SerialLinearOpBase<Scalar>::setSpaces(
  const Teuchos::RCP<const ScalarProdVectorSpaceBase<Scalar> >     &range
  ,const Teuchos::RCP<const ScalarProdVectorSpaceBase<Scalar> >    &domain
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
#endif
  range_  = range;
  domain_ = domain;
}

template<class Scalar>
void SerialLinearOpBase<Scalar>::setDimensions(
  const Index                                                dimRange
  ,const Index                                               dimDomain
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( dimRange <= 0 );
  TEST_FOR_EXCEPT( dimDomain <= 0 );
#endif
  range_  = Teuchos::rcp(new DefaultSerialVectorSpace<Scalar>(dimRange));
  domain_ = Teuchos::rcp(new DefaultSerialVectorSpace<Scalar>(dimDomain));
}

// Virtual functions to be overridden by subclasses

template<class Scalar>
void SerialLinearOpBase<Scalar>::euclideanApply(
  const EOpTransp                                     M_trans
  ,const RTOpPack::ConstSubMultiVectorView<Scalar>          &X
  ,const RTOpPack::SubMultiVectorView<Scalar>   *Y
  ,const Scalar                                     alpha
  ,const Scalar                                     beta
  ) const
{
  for(Index j = 0; j < X.numSubCols(); ++j ) {
    const RTOpPack::ConstSubVectorView<Scalar>         X_j = X.col(j);
    const RTOpPack::SubVectorView<Scalar>  Y_j = Y->col(j);
    this->euclideanApply(M_trans,X_j,&Y_j,alpha,beta);
  }
}

}	// end namespace Thyra

#endif	// THYRA_SERIAL_LINEAR_OP_BASE_HPP
