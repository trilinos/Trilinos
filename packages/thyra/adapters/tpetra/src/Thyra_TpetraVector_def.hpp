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

#ifndef THYRA_TPETRA_VECTOR_HPP
#define THYRA_TPETRA_VECTOR_HPP


#include "Thyra_TpetraVector_decl.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraVector()
{}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
  )
{
  initializeImpl(tpetraVectorSpace, tpetraVector);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::constInitialize(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
  )
{
  initializeImpl(tpetraVectorSpace, tpetraVector);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetraVector()
{
  return tpetraVector_.getNonconstObj();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraVector() const
{
  return tpetraVector_;
}


// Overridden from SpmdMultiVectorBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const SpmdVectorSpaceBase<Scalar> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::spmdSpaceImpl() const
{
  return tpetraVectorSpace_;
}


// Overridden from SpmdVectorBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstLocalVectorDataImpl(
  const Ptr<ArrayRCP<Scalar> > &localValues )
{
  *localValues = tpetraVector_.getNonconstObj()->get1dViewNonConst();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalVectorDataImpl(
  const Ptr<ArrayRCP<const Scalar> > &localValues ) const
{
  *localValues = tpetraVector_->get1dView();
}


// Overridden from VectorBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::assignVecImpl(
  const VectorBase<Scalar>& x
  )
{
  // Try to cast x to this type.
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  RCP<const TV> tx = Teuchos::rcp_dynamic_cast<const TV>(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    tpetraVector_.getNonconstObj()->assign(*tx->getConstTpetraVector());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::assignVecImpl(x);
  }  
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::randomizeImpl(
  Scalar l,
  Scalar u
  )
{
  tpetraVector_.getNonconstObj()->randomize(l, u);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::absImpl(
  const VectorBase<Scalar>& x
  )
{
  // Try to cast x to this type.
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  RCP<const TV> tx = Teuchos::rcp_dynamic_cast<const TV>(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    tpetraVector_.getNonconstObj()->abs(*tx->getConstTpetraVector());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::absImpl(x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::reciprocalImpl(
  const VectorBase<Scalar>& x
  )
{
  // Try to cast x to this type.
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  RCP<const TV> tx = Teuchos::rcp_dynamic_cast<const TV>(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    tpetraVector_.getNonconstObj()->reciprocal(*tx->getConstTpetraVector());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::reciprocalImpl(x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::eleWiseScaleImpl(
  const VectorBase<Scalar>& x
  )
{
  // Try to cast x to this type.
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  RCP<const TV> tx = Teuchos::rcp_dynamic_cast<const TV>(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    tpetraVector_.getNonconstObj()->elementWiseMultiply(
      ST::one(), *tx->getConstTpetraVector(), *tpetraVector_.getConstObj(), ST::zero());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::eleWiseScaleImpl(x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateVecImpl(
  Scalar alpha,
  const VectorBase<Scalar>& x
  )
{
  // Try to cast x to this type.
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  RCP<const TV> tx = Teuchos::rcp_dynamic_cast<const TV>(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    tpetraVector_.getNonconstObj()->update(alpha, *tx->getConstTpetraVector(), ST::one());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::updateVecImpl(alpha, x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::linearCombinationVecImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const VectorBase<Scalar> > >& x,
  const Scalar& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), x.size());
#endif

  // Try to cast x to an array of this type
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  Teuchos::Array<Ptr<const TV> > tx(x.size());
  bool allCastsSuccessful = true;
  {
    auto xIter = x.begin();
    auto txIter = tx.begin();
    for (; txIter != tx.end(); ++xIter, ++txIter) {
      *txIter = Teuchos::ptr_dynamic_cast<const TV>(*xIter);
      if (is_null(*txIter)) {
        allCastsSuccessful = false;
        break;
      }
    }
  }

  // If each cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (allCastsSuccessful) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    auto txIter = tx.begin();
    auto alphaIter = alpha.begin();
    // We add two vecs at a time, so only scale if even num vecs,
    // and additionally do the first addition if odd num vecs.
    if ((tx.size() % 2) == 0) {
      tpetraVector_.getNonconstObj()->scale(beta);
    } else {
      tpetraVector_.getNonconstObj()->update(
        *alphaIter, *(*txIter)->getConstTpetraVector(), beta);
      ++txIter;
      ++alphaIter;
    }
    for (; txIter != tx.end(); txIter+=2, alphaIter+=2) {
      tpetraVector_.getNonconstObj()->update(
        *alphaIter, *(*txIter)->getConstTpetraVector(),
        *(alphaIter+1), *(*(txIter+1))->getConstTpetraVector(), ST::one());
    }
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::linearCombinationVecImpl(alpha, x, beta);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1Impl() const
{
  return tpetraVector_.getConstObj()->norm1();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2Impl() const
{
  return tpetraVector_.getConstObj()->norm2();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2WeightedImpl(
  const VectorBase<Scalar>& x
  ) const
{
  // Try to cast weighting vector to this type.
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  RCP<const TV> tx = Teuchos::rcp_dynamic_cast<const TV>(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    // Weighted 2-norm function for Tpetra vector seems to be deprecated...
    typedef Teuchos::ScalarTraits<Scalar> ST;
    RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > temp
      = Tpetra::createVector<Scalar>(tx->getConstTpetraVector()->getMap());
    temp->elementWiseMultiply(
      ST::one(), *tx->getConstTpetraVector(), *tpetraVector_.getConstObj(), ST::zero());
    return ST::magnitude(ST::squareroot(tpetraVector_.getConstObj()->dot(*temp)));
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    return VectorDefaultBase<Scalar>::norm2WeightedImpl(x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInfImpl() const
{
  return tpetraVector_.getConstObj()->normInf();
}


// private


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template<class TpetraVector_t>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initializeImpl(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<TpetraVector_t> &tpetraVector
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(nonnull(tpetraVectorSpace));
  TEUCHOS_ASSERT(nonnull(tpetraVector));
#endif
  tpetraVectorSpace_ = tpetraVectorSpace;
  tpetraVector_.initialize(tpetraVector);
  this->updateSpmdSpace();
}


} // end namespace Thyra


#endif // THYRA_TPETRA_VECTOR_HPP
