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


#ifndef THYRA_TPETRA_VECTOR_SPACE_HPP
#define THYRA_TPETRA_VECTOR_SPACE_HPP


#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Thyra_TpetraMultiVector.hpp"


namespace Thyra {


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::create()
{
  const RCP<this_t> vs(new this_t);
  vs->weakSelfPtr_ = vs.create_weak();
  return vs;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize(
  const RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &tpetraMap
  )
{
  comm_ = convertTpetraToThyraComm(tpetraMap->getComm());
  tpetraMap_ = tpetraMap;
  localSubDim_ = tpetraMap->getNodeNumElements();
  numProc_ = comm_->getSize();
  procRank_ = comm_->getRank();
  this->updateState(tpetraMap->getGlobalNumElements());
}


// Overridden from VectorSpace


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<VectorBase<Scalar> >
TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createMember() const
{
  return tpetraVector<Scalar>(
    weakSelfPtr_.create_strong().getConst(),
    Teuchos::rcp(
      new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tpetraMap_, false)
      )
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< MultiVectorBase<Scalar> >
TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createMembers(int numMembers) const
{
  return tpetraMultiVector<Scalar>(
    weakSelfPtr_.create_strong().getConst(),
    tpetraVectorSpace<Scalar>(
      Tpetra::createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(
        numMembers, tpetraMap_->getComm(), tpetraMap_->getNode()
        )
      ),
    Teuchos::rcp(
      new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
        tpetraMap_, numMembers, false)
      )
    );
  // ToDo: Create wrapper function to create locally replicated vector space
  // and use it.
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasInCoreView(
  const Range1D& rng_in, const EViewType viewType, const EStrideType strideType
  ) const
{
  const Range1D rng = full_range(rng_in,0,this->dim()-1);
  const Ordinal l_localOffset = this->localOffset();
  return ( l_localOffset<=rng.lbound() && rng.ubound()<l_localOffset+localSubDim_ );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const VectorSpaceBase<Scalar> >
TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::clone() const
{
  return tpetraVectorSpace<Scalar>(tpetraMap_);
}


// Overridden from SpmdVectorSpaceDefaultBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Teuchos::Comm<Ordinal> >
TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getComm() const
{
  return comm_;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Ordinal TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::localSubDim() const
{
  return localSubDim_;
}


// private


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraVectorSpace()
  :localSubDim_(-1), numProc_(-1), procRank_(-1)
{
  // The base classes should automatically default initialize to a safe
  // uninitialized state.
}


} // end namespace Thyra


#endif // THYRA_TPETRA_VECTOR_SPACE_HPP
