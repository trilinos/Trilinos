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

#ifndef THYRA_MPI_LINEAR_OP_BASE_HPP
#define THYRA_MPI_LINEAR_OP_BASE_HPP

#include "Thyra_MPILinearOpBaseDecl.hpp"
#include "Thyra_EuclideanLinearOpBase.hpp"
#include "Thyra_MPIVectorSpaceStd.hpp"
#include "Thyra_ExplicitVectorView.hpp"
#include "Thyra_ExplicitMultiVectorView.hpp"

namespace Thyra {

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
MPILinearOpBase<Scalar>::rangeScalarProdVecSpc() const
{
  return range_;
}

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
MPILinearOpBase<Scalar>::domainScalarProdVecSpc() const
{
  return domain_;
}

template<class Scalar>
void MPILinearOpBase<Scalar>::euclideanApply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
#ifdef _DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES("MPILinearOpBase<Scalar>::euclideanApply()",*this,M_trans,X,Y);
#endif
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >
    &Op_range  = ( M_trans == NOTRANS ? range_  : domain_ ),
    &Op_domain = ( M_trans == NOTRANS ? domain_ : range_  );
  const Index
    localOffsetRange  = Op_range->localOffset(), 
    localDimRange     = Op_range->localSubDim(),
    localOffsetDomain = Op_domain->localOffset(), 
    localDimDomain    = Op_domain->localSubDim(); 
  const ExplicitMultiVectorView<Scalar>         local_X_ev(X,Range1D(localOffsetDomain+1,localOffsetDomain+localDimDomain));
  const ExplicitMutableMultiVectorView<Scalar>  local_Y_ev(*Y,Range1D(localOffsetRange+1,localOffsetRange+localDimRange));
  this->euclideanApply(M_trans,local_X_ev.smv(),&local_Y_ev.smv(),alpha,beta);
}

// protected

// Constructors/initializers/accessors

template<class Scalar>
MPILinearOpBase<Scalar>::MPILinearOpBase()
  :forceUnitStride_(true)
{}

template<class Scalar>
void MPILinearOpBase<Scalar>::setSpaces(
  const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >     &range
  ,const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    &domain
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
void MPILinearOpBase<Scalar>::setLocalDimensions(
  MPI_Comm                                                    mpiComm
  ,const Index                                                localDimRange
  ,const Index                                                localDimDomain
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( localDimRange <= 0 );
  TEST_FOR_EXCEPT( localDimDomain <= 0 );
#endif
  range_  = Teuchos::rcp(new MPIVectorSpaceStd<Scalar>(mpiComm,localDimRange,-1));
  domain_ = Teuchos::rcp(new MPIVectorSpaceStd<Scalar>(mpiComm,localDimDomain,-1));
}

// Virtual functions to be overridden by subclasses

template<class Scalar>
void MPILinearOpBase<Scalar>::euclideanApply(
  const ETransp                                     M_trans
  ,const RTOpPack::SubMultiVectorT<Scalar>          &local_X
  ,const RTOpPack::MutableSubMultiVectorT<Scalar>   *local_Y
  ,const Scalar                                     alpha
  ,const Scalar                                     beta
  ) const
{
  for(Index j = 1; j <= local_X.numSubCols(); ++j ) {
    const RTOpPack::SubVectorT<Scalar>         local_X_j = local_X.col(j);
    const RTOpPack::MutableSubVectorT<Scalar>  local_Y_j = local_Y->col(j);
    this->euclideanApply(M_trans,local_X_j,&local_Y_j,alpha,beta);
  }
}

}	// end namespace Thyra

#endif	// THYRA_MPI_LINEAR_OP_BASE_HPP
