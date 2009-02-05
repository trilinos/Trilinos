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

#ifndef THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP

#include "Thyra_DefaultSpmdVectorSpace_decl.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


template<class Scalar>
DefaultSpmdVectorSpace<Scalar>::DefaultSpmdVectorSpace()
  :localSubDim_(0), numProc_(0), procRank_(0)
{
  // The base classes should automatically default initialize to a safe
  // uninitialized state.
}


template<class Scalar>
DefaultSpmdVectorSpace<Scalar>::DefaultSpmdVectorSpace(
  const Index dim
  )
  :localSubDim_(0), numProc_(0), procRank_(0)
{
  initialize(dim);
}


template<class Scalar>
DefaultSpmdVectorSpace<Scalar>::DefaultSpmdVectorSpace(
  const RCP<const Teuchos::Comm<Index> > &comm_in
  ,const Index localSubDim_in, const Index globalDim_in
  )
  :localSubDim_(0), numProc_(0), procRank_(0)
{
  initialize(comm_in, localSubDim_in, globalDim_in);
}


template<class Scalar>
void DefaultSpmdVectorSpace<Scalar>::initialize(
  const Index dim
  )
{
  this->initialize(Teuchos::null, dim, dim);
}


template<class Scalar>
void DefaultSpmdVectorSpace<Scalar>::initialize(
  const RCP<const Teuchos::Comm<Index> > &comm
  ,const Index localSubDim_in, const Index globalDim
  )
{
#ifdef TEUCHOS_DEBUG
  //TEST_FOR_EXCEPT( !( localSubDim_in > 0 ) );
  TEST_FOR_EXCEPT( !( localSubDim_in >= 0 ) );
#endif
  comm_        = comm;
  localSubDim_ = localSubDim_in;
  if( comm.get() ) {
    numProc_  = size(*comm);
    procRank_ = rank(*comm);
  }
  else {
    numProc_  = 1;
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
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      Teuchos::rcp(this,false),
      Teuchos::arcp<Scalar>(localSubDim_),
      1
      )
    );
}


template<class Scalar>
RCP< MultiVectorBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::createMembers(int numMembers) const
{
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
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
  TEST_FOR_EXCEPT( localSubDim_ != raw_v.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      Teuchos::rcp(this,false),
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
  TEST_FOR_EXCEPT( localSubDim_ != raw_v.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdVector<Scalar>(
      Teuchos::rcp(this,false),
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
  TEST_FOR_EXCEPT( localSubDim_ != raw_mv.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      Teuchos::rcp(this,false),
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
  TEST_FOR_EXCEPT( localSubDim_ != raw_mv.subDim() );
#endif
  return Teuchos::rcp(
    new DefaultSpmdMultiVector<Scalar>(
      Teuchos::rcp(this,false),
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
  const Index localOffset = this->localOffset();
  return ( localOffset<=rng.lbound() && rng.ubound()<localOffset+localSubDim_ );
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultSpmdVectorSpace<Scalar>::clone() const
{
  return Teuchos::rcp(
    new DefaultSpmdVectorSpace<Scalar>(comm_,localSubDim_,this->dim())
    );
}


// Overridden from SpmdVectorSpaceDefaultBase


template<class Scalar>
RCP<const Teuchos::Comm<Index> >
DefaultSpmdVectorSpace<Scalar>::getComm() const
{
  return comm_;
}


template<class Scalar>
Index DefaultSpmdVectorSpace<Scalar>::localSubDim() const
{
  return localSubDim_;
}


} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_DEF_HPP
