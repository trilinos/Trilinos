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

#ifndef THYRA_SPMD_LINEAR_OP_BASE_HPP
#define THYRA_SPMD_LINEAR_OP_BASE_HPP

#include "Thyra_SpmdLinearOpBaseDecl.hpp"
#include "Thyra_SingleScalarEuclideanLinearOpBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

namespace Thyra {

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
SpmdLinearOpBase<Scalar>::rangeScalarProdVecSpc() const
{
  return sp_range_;
}

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
SpmdLinearOpBase<Scalar>::domainScalarProdVecSpc() const
{
  return sp_domain_;
}

template<class Scalar>
void SpmdLinearOpBase<Scalar>::euclideanApply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "SpmdLinearOpBase<Scalar>::euclideanApply()",*this,M_trans,X,Y
    );
#endif
  const Teuchos::RefCountPtr<const SpmdVectorSpaceBase<Scalar> >
    &Op_range  = ( M_trans == NOTRANS ? range_  : domain_ ),
    &Op_domain = ( M_trans == NOTRANS ? domain_ : range_  );
  const Index
    localOffsetRange  = Op_range->localOffset(), 
    localDimRange     = Op_range->localSubDim(),
    localOffsetDomain = Op_domain->localOffset(), 
    localDimDomain    = Op_domain->localSubDim(); 
  const ConstDetachedMultiVectorView<Scalar>
    local_X_ev(X,Range1D(localOffsetDomain,localOffsetDomain+localDimDomain-1));
  const DetachedMultiVectorView<Scalar>
    local_Y_ev(*Y,Range1D(localOffsetRange,localOffsetRange+localDimRange-1));
  this->euclideanApply(M_trans,local_X_ev.smv(),&local_Y_ev.smv(),alpha,beta);
}

// protected

// Constructors/initializers/accessors

template<class Scalar>
SpmdLinearOpBase<Scalar>::SpmdLinearOpBase()
  :forceUnitStride_(true)
{}

template<class Scalar>
void SpmdLinearOpBase<Scalar>::setSpaces(
  const Teuchos::RefCountPtr<const SpmdVectorSpaceBase<Scalar> >     &range
  ,const Teuchos::RefCountPtr<const SpmdVectorSpaceBase<Scalar> >    &domain
  )
{
  // Validate input
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
#endif
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >
    sp_range = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(range,true),
    sp_domain = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(domain,true);
  // Set state
  range_  = range;
  domain_ = domain;
  sp_range_  = sp_range;
  sp_domain_ = sp_domain;
}

template<class Scalar>
void SpmdLinearOpBase<Scalar>::setLocalDimensions(
  const Teuchos::RefCountPtr<const Teuchos::Comm<Index> >     &comm
  ,const Index                                                localDimRange
  ,const Index                                                localDimDomain
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( localDimRange <= 0 );
  TEST_FOR_EXCEPT( localDimDomain <= 0 );
#endif
  Teuchos::RefCountPtr<const DefaultSpmdVectorSpace<Scalar> >
    range  = Teuchos::rcp(new DefaultSpmdVectorSpace<Scalar>(comm,localDimRange,-1)),
    domain = Teuchos::rcp(new DefaultSpmdVectorSpace<Scalar>(comm,localDimDomain,-1));
  range_  = range;
  domain_ = domain;
  sp_range_  = range;
  sp_domain_ = domain;
}

// Virtual functions to be overridden by subclasses

template<class Scalar>
void SpmdLinearOpBase<Scalar>::euclideanApply(
  const ETransp                                     M_trans
  ,const RTOpPack::ConstSubMultiVectorView<Scalar>  &local_X
  ,const RTOpPack::SubMultiVectorView<Scalar>       *local_Y
  ,const Scalar                                     alpha
  ,const Scalar                                     beta
  ) const
{
  for(Index j = 0; j < local_X.numSubCols(); ++j ) {
    const RTOpPack::ConstSubVectorView<Scalar>
      local_X_j = local_X.col(j);
    const RTOpPack::SubVectorView<Scalar>
      local_Y_j = local_Y->col(j);
    this->euclideanApply(M_trans,local_X_j,&local_Y_j,alpha,beta);
  }
}

}	// end namespace Thyra

#endif	// THYRA_SPMD_LINEAR_OP_BASE_HPP
