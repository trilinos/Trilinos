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
#include "Thyra_SerialVectorSpaceStd.hpp"
#include "Thyra_ExplicitVectorView.hpp"
#include "Thyra_ExplicitMultiVectorView.hpp"

namespace Thyra {

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
SerialLinearOpBase<Scalar>::rangeScalarProdVecSpc() const
{
  return range_;
}

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
SerialLinearOpBase<Scalar>::domainScalarProdVecSpc() const
{
  return domain_;
}

// Overridden from LinearOpBase

template<class Scalar>
void SerialLinearOpBase<Scalar>::euclideanApply(
  const ETransp            M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar            alpha
  ,const Scalar            beta
  ) const
{
#ifdef _DEBUG
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("SerialLinearOpBase<Scalar>::apply()",*this,M_trans,x,y);
#endif
  const ExplicitVectorView<Scalar>         x_ev(x,Range1D(),forceUnitStride());
  const ExplicitMutableVectorView<Scalar>  y_ev(*y,Range1D(),forceUnitStride());
  this->euclideanApply(M_trans,x_ev.sv(),&y_ev.sv(),alpha,beta);
}

template<class Scalar>
void SerialLinearOpBase<Scalar>::euclideanApply(
  const ETransp                 M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                 alpha
  ,const Scalar                 beta
  ) const
{
#ifdef _DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES("SerialLinearOpBase<Scalar>::apply()",*this,M_trans,X,Y);
#endif
  const ExplicitMultiVectorView<Scalar>         X_ev(X);
  const ExplicitMutableMultiVectorView<Scalar>  Y_ev(*Y);
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
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >     &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >    &domain
  )
{
#ifdef _DEBUG
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
#ifdef _DEBUG
  TEST_FOR_EXCEPT( dimRange <= 0 );
  TEST_FOR_EXCEPT( dimDomain <= 0 );
#endif
  range_  = Teuchos::rcp(new SerialVectorSpaceStd<Scalar>(dimRange));
  domain_ = Teuchos::rcp(new SerialVectorSpaceStd<Scalar>(dimDomain));
}

// Virtual functions to be overridden by subclasses

template<class Scalar>
void SerialLinearOpBase<Scalar>::euclideanApply(
  const ETransp                                     M_trans
  ,const RTOpPack::SubMultiVectorT<Scalar>          &X
  ,const RTOpPack::MutableSubMultiVectorT<Scalar>   *Y
  ,const Scalar                                     alpha
  ,const Scalar                                     beta
  ) const
{
  for(Index j = 1; j <= X.numSubCols(); ++j ) {
    const RTOpPack::SubVectorT<Scalar>         X_j = X.col(j);
    const RTOpPack::MutableSubVectorT<Scalar>  Y_j = Y->col(j);
    this->euclideanApply(M_trans,X_j,&Y_j,alpha,beta);
  }
}

}	// end namespace Thyra

#endif	// THYRA_SERIAL_LINEAR_OP_BASE_HPP
