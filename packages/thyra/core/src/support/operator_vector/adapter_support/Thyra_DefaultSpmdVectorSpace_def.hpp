// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const RCP<const Teuchos::Comm<Ordinal> > &comm,
  const Ordinal localSubDim_in, const Ordinal globalDim,
  const bool isLocallyReplicated_in
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( localSubDim_in >= 0 ) );
#endif
  comm_ = comm;
  localSubDim_ = localSubDim_in;
  if (! comm.is_null ()) {
    numProc_ = comm->getSize ();
    procRank_ = comm->getRank ();
  }
  else {
    numProc_ = 1;
    procRank_ = 0;
  }
  this->updateState(globalDim, isLocallyReplicated_in);
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
  const Range1D& rng_in, const EViewType /* viewType */, const EStrideType /* strideType */
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
  return defaultSpmdVectorSpace<Scalar>(comm_, localSubDim_, this->dim(),
    this->isLocallyReplicated());
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


} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP
