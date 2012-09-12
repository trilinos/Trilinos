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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP

#include "Thyra_DefaultSpmdVectorSpace_decl.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


template<class Scalar>
RCP<DefaultSpmdVectorSpace<Scalar> >
DefaultSpmdVectorSpace<Scalar>::create()
{ 
  const RCP<DefaultSpmdVectorSpace<Scalar> > vs(new DefaultSpmdVectorSpace<Scalar>);
  vs->weakSelfPtr_ = vs.create_weak();
  return vs;
}


template<class Scalar>
void DefaultSpmdVectorSpace<Scalar>::initialize(
  const Ordinal dim_in
  )
{
  this->initialize(Teuchos::null, dim_in, dim_in);
}


template<class Scalar>
void DefaultSpmdVectorSpace<Scalar>::initialize(
  const RCP<const Teuchos::Comm<Ordinal> > &comm
  ,const Ordinal localSubDim_in, const Ordinal globalDim
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( localSubDim_in >= 0 ) );
#endif
  comm_ = comm;
  localSubDim_ = localSubDim_in;
  if (!is_null(comm)) {
    numProc_ = size(*comm);
    procRank_ = rank(*comm);
  }
  else {
    numProc_ = 1;
    procRank_ = 0;
  }
  this->updateState(globalDim);
}


template<class Scalar>
void DefaultSpmdVectorSpace<Scalar>::uninitialize()
{
  comm_         = Teuchos::null;
  localSubDim_  = 0;
}


// Overridden from VectorSpace


template<class Scalar>
RCP<VectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMember() const
{
  ArrayRCP<Scalar> values;
  if (localSubDim_)
    values = Teuchos::arcp<Scalar>(localSubDim_);
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      weakSelfPtr_.create_strong(),
      values,
      1 // stride
      )
    );
}


template<class Scalar>
RCP< MultiVectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      weakSelfPtr_.create_strong(),
      Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
        this->smallVecSpcFcty()->createVecSpc(numMembers),true
        )
      )
    );
}


template<class Scalar>
RCP<VectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMemberView(
  const RTOpPack::SubVectorView<Scalar> &raw_v
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( localSubDim_ != raw_v.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      weakSelfPtr_.create_strong(),
      Teuchos::arcp(raw_v.values().get(),0,raw_v.subDim(),false),
      raw_v.stride()
      )
    );
}


template<class Scalar>
RCP<const VectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMemberView(
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( localSubDim_ != raw_v.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      weakSelfPtr_.create_strong(),
      Teuchos::arcp(const_cast<Scalar*>(raw_v.values().get()),0,raw_v.subDim(),false),
      raw_v.stride()
      )
    );
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMembersView(
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( localSubDim_ != raw_mv.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      weakSelfPtr_.create_strong(),
      Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
        this->smallVecSpcFcty()->createVecSpc(raw_mv.numSubCols()),true),
      Teuchos::arcp(raw_mv.values().get(),0,raw_mv.leadingDim()*raw_mv.numSubCols(),false),
      raw_mv.leadingDim()
      )
    );
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMembersView(
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( localSubDim_ != raw_mv.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      weakSelfPtr_.create_strong(),
      Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
        this->smallVecSpcFcty()->createVecSpc(raw_mv.numSubCols()),true),
      Teuchos::arcp(
        const_cast<Scalar*>(raw_mv.values().get()),
        0,raw_mv.leadingDim()*raw_mv.numSubCols(),false),
      raw_mv.leadingDim()
      )
    );
}


template<class Scalar>
bool DefaultSpmdVectorSpace<Scalar>::hasInCoreView(
  const Range1D& rng_in, const EViewType viewType, const EStrideType strideType
  ) const
{
  const Range1D rng = full_range(rng_in,0,this->dim()-1);
  const Ordinal l_localOffset = this->localOffset();
  return ( l_localOffset<=rng.lbound() && rng.ubound()<l_localOffset+localSubDim_ );
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::clone() const
{
  return defaultSpmdVectorSpace<Scalar>(comm_,localSubDim_,this->dim());
}


// Overridden from SpmdVectorSpaceDefaultBase


template<class Scalar>
RCP<const Teuchos::Comm<Ordinal> >
DefaultSpmdVectorSpace<Scalar>::getComm() const
{
  return comm_;
}


template<class Scalar>
Ordinal DefaultSpmdVectorSpace<Scalar>::localSubDim() const
{
  return localSubDim_;
}


// private


template<class Scalar>
DefaultSpmdVectorSpace<Scalar>::DefaultSpmdVectorSpace()
  :localSubDim_(-1), numProc_(-1), procRank_(-1)
{
  // The base classes should automatically default initialize to a safe
  // uninitialized state.
}


#ifndef THYRA_HIDE_DEPRECATED_CODE
// Deprecated

template<class Scalar>
DefaultSpmdVectorSpace<Scalar>::DefaultSpmdVectorSpace(
  const Ordinal dim_in
  )
  :localSubDim_(-1), numProc_(-1), procRank_(-1)
{
  initialize(dim_in);
  weakSelfPtr_ = Teuchos::rcpFromRef(*this);
}

template<class Scalar>
DefaultSpmdVectorSpace<Scalar>::DefaultSpmdVectorSpace(
  const RCP<const Teuchos::Comm<Ordinal> > &comm,
  const Ordinal my_localSubDim, const Ordinal globalDim
  )
  :localSubDim_(-1), numProc_(-1), procRank_(-1)
{
  initialize(comm, my_localSubDim, globalDim);
  weakSelfPtr_ = Teuchos::rcpFromRef(*this);
}


#endif // THYRA_HIDE_DEPRECATED_CODE
} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP
