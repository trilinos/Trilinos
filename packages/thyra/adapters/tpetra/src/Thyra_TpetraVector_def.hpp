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
#include "Thyra_TpetraMultiVector.hpp"

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


// Overridden from VectorDefaultBase
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const VectorSpaceBase<Scalar> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::domain() const
{
  if (domainSpace_.is_null()) {
    domainSpace_ = tpetraVectorSpace<Scalar>(
      Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal>(
        1,
        tpetraVector_.getConstObj()->getMap()->getComm(),
        tpetraVector_.getConstObj()->getMap()->getNode()
        )
      );
  }
  return domainSpace_;
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
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::randomizeImpl(
  Scalar l,
  Scalar u
  )
{
  // Tpetra randomizes with different seed for each proc, so need a global
  // reduction to get locally-replicated random vector same on each proc.
  if (!tpetraVector_.getNonconstObj()->isDistributed()) {
    auto comm = tpetraVector_.getNonconstObj()->getMap()->getComm();
    if (tpetraVector_.getConstObj()->getMap()->getComm()->getRank() == 0)
      tpetraVector_.getNonconstObj()->randomize(l, u);
    else
      tpetraVector_.getNonconstObj()->putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    tpetraVector_.getNonconstObj()->reduce();
  } else {
    tpetraVector_.getNonconstObj()->randomize(l, u);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::absImpl(
  const VectorBase<Scalar>& x
  )
{
  auto tx = this->getConstTpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    tpetraVector_.getNonconstObj()->abs(*tx);
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
  auto tx = this->getConstTpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    tpetraVector_.getNonconstObj()->reciprocal(*tx);
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
  auto tx = this->getConstTpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    tpetraVector_.getNonconstObj()->elementWiseMultiply(
      ST::one(), *tx, *tpetraVector_.getConstObj(), ST::zero());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::eleWiseScaleImpl(x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2WeightedImpl(
  const VectorBase<Scalar>& x
  ) const
{
  auto tx = this->getConstTpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    // Weighted 2-norm function for Tpetra vector seems to be deprecated...
    typedef Teuchos::ScalarTraits<Scalar> ST;
    RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > temp
      = Tpetra::createVector<Scalar>(tx->getMap());
    temp->elementWiseMultiply(
      ST::one(), *tx, *tpetraVector_.getConstObj(), ST::zero());
    return ST::magnitude(ST::squareroot(tpetraVector_.getConstObj()->dot(*temp)));
  } else {
    // This version will require the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    return VectorDefaultBase<Scalar>::norm2WeightedImpl(x);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset
  ) const
{
  // Sync any non-target Tpetra vecs to host space
  for (auto itr = vecs.begin(); itr != vecs.end(); ++itr) {
    auto tv = this->getConstTpetraVector(Teuchos::rcpFromPtr(*itr));
    if (nonnull(tv)) {
      typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
      Teuchos::rcp_const_cast<TV>(tv)->template sync<Kokkos::HostSpace>();
    }
  }

  // Sync any target Tpetra vecs and mark modified on host
  for (auto itr = targ_vecs.begin(); itr != targ_vecs.end(); ++itr) {
    auto tv = this->getTpetraVector(Teuchos::rcpFromPtr(*itr));
    if (nonnull(tv)) {
      tv->template sync<Kokkos::HostSpace>();
      tv->template modify<Kokkos::HostSpace>();
    }
  }

  SpmdVectorDefaultBase<Scalar>::applyOpImpl(op, vecs, targ_vecs, reduct_obj, global_offset);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
acquireDetachedVectorViewImpl(
  const Range1D& rng,
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  // Only viewing data, so just sync dual view to host space
  typedef typename Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  Teuchos::rcp_const_cast<TV>(
    tpetraVector_.getConstObj())->template sync<Kokkos::HostSpace>();

  SpmdVectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(rng, sub_vec);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng,
  RTOpPack::SubVectorView<Scalar>* sub_vec 
  )
{
  // Sync to host and mark as modified
  tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
  tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();

  SpmdVectorDefaultBase<Scalar>::acquireNonconstDetachedVectorViewImpl(rng, sub_vec);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  SpmdVectorDefaultBase<Scalar>::commitNonconstDetachedVectorViewImpl(sub_vec);

  // Sync changes from host view to execution space
  typedef typename Tpetra::Vector<
    Scalar,LocalOrdinal,GlobalOrdinal,Node>::execution_space execution_space;
  tpetraVector_.getNonconstObj()->template sync<execution_space>();
}


// Overridden protected functions from MultiVectorBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::assignImpl(Scalar alpha)
{
  tpetraVector_.getNonconstObj()->putScalar(alpha);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
assignMultiVecImpl(const MultiVectorBase<Scalar>& mv)
{
  auto tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tmv)) {
    tpetraVector_.getNonconstObj()->assign(*tmv);
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    MultiVectorDefaultBase<Scalar>::assignMultiVecImpl(mv);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scaleImpl(Scalar alpha)
{
  tpetraVector_.getNonconstObj()->scale(alpha);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateImpl(
  Scalar alpha,
  const MultiVectorBase<Scalar>& mv
  )
{
  auto tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if (nonnull(tmv)) {
    tpetraVector_.getNonconstObj()->update(alpha, *tmv, ST::one());
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    MultiVectorDefaultBase<Scalar>::updateImpl(alpha, mv);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::linearCombinationImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
  const Scalar& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
#endif

  // Try to cast mv to Tpetra objects
  Teuchos::Array<RCP<const TpetraMultiVector_t> > tmvs(mv.size());
  RCP<const TpetraMultiVector_t> tmv;
  bool allCastsSuccessful = true;
  {
    auto mvIter = mv.begin();
    auto tmvIter = tmvs.begin();
    for (; mvIter != mv.end(); ++mvIter, ++tmvIter) {
      tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromPtr(*mvIter));
      if (nonnull(tmv)) {
        *tmvIter = tmv;
      } else {
        allCastsSuccessful = false;
        break;
      }
    }
  }

  // If casts succeeded, or input arrays are size 0, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  auto len = mv.size();
  if (len == 0) {
    tpetraVector_.getNonconstObj()->scale(beta);
  } else if (len == 1 && allCastsSuccessful) {
    tpetraVector_.getNonconstObj()->update(alpha[0], *tmvs[0], beta);
  } else if (len == 2 && allCastsSuccessful) {
    tpetraVector_.getNonconstObj()->update(alpha[0], *tmvs[0], alpha[1], *tmvs[1], beta);
  } else if (allCastsSuccessful) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    auto tmvIter = tmvs.begin();
    auto alphaIter = alpha.begin();

    // Check if any entry of tmvs aliases this object's wrapped vector.
    // If so, replace that entry in the array with a copy.
    tmv = Teuchos::null;
    for (; tmvIter != tmvs.end(); ++tmvIter) {
      if (tmvIter->getRawPtr() == tpetraVector_.getConstObj().getRawPtr()) {
        if (tmv.is_null()) {
          tmv = Teuchos::rcp(new TpetraMultiVector_t(
            *tpetraVector_.getConstObj(), Teuchos::Copy));
        }
        *tmvIter = tmv;
      }
    }
    tmvIter = tmvs.begin();

    // We add two MVs at a time, so only scale if even num MVs,
    // and additionally do the first addition if odd num MVs.
    if ((tmvs.size() % 2) == 0) {
      tpetraVector_.getNonconstObj()->scale(beta);
    } else {
      tpetraVector_.getNonconstObj()->update(*alphaIter, *(*tmvIter), beta);
      ++tmvIter;
      ++alphaIter;
    }
    for (; tmvIter != tmvs.end(); tmvIter+=2, alphaIter+=2) {
      tpetraVector_.getNonconstObj()->update(
        *alphaIter, *(*tmvIter), *(alphaIter+1), *(*(tmvIter+1)), ST::one());
    }
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    MultiVectorDefaultBase<Scalar>::linearCombinationImpl(alpha, mv, beta);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dotsImpl(
    const MultiVectorBase<Scalar>& mv,
    const ArrayView<Scalar>& prods
    ) const
{
  auto tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tmv)) {
    tpetraVector_.getConstObj()->dot(*tmv, prods);
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    MultiVectorDefaultBase<Scalar>::dotsImpl(mv, prods);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norms1Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  tpetraVector_.getConstObj()->norm1(norms);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norms2Impl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const
{
  tpetraVector_.getConstObj()->norm2(norms);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normsInfImpl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  tpetraVector_.getConstObj()->normInf(norms);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  // Try to extract Tpetra objects from X and Y
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;
  Teuchos::RCP<const TMV> X_tpetra = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(X));
  Teuchos::RCP<TMV> Y_tpetra = this->getTpetraMultiVector(Teuchos::rcpFromPtr(Y));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the default implementation.
  if (nonnull(X_tpetra) && nonnull(Y_tpetra)) {
    // Sync everything to the execution space
    typedef typename TMV::execution_space execution_space;
    Teuchos::rcp_const_cast<TMV>(X_tpetra)->template sync<execution_space>();
    Y_tpetra->template sync<execution_space>();
    Teuchos::rcp_const_cast<TV>(tpetraVector_.getConstObj())->template sync<execution_space>();

    typedef Teuchos::ScalarTraits<Scalar> ST;
    TEUCHOS_TEST_FOR_EXCEPTION(ST::isComplex && (M_trans == CONJ),
      std::logic_error,
      "Error, conjugation without transposition is not allowed for complex scalar types!");

    Teuchos::ETransp trans = Teuchos::NO_TRANS;
    switch (M_trans) {
      case NOTRANS:
        trans = Teuchos::NO_TRANS;
        break;
      case CONJ:
        trans = Teuchos::NO_TRANS;
        break;
      case TRANS:
        trans = Teuchos::TRANS;
        break;
      case CONJTRANS:
        trans = Teuchos::CONJ_TRANS;
        break;
    }

    Y_tpetra->template modify<execution_space>();
    Y_tpetra->multiply(trans, Teuchos::NO_TRANS, alpha, *tpetraVector_.getConstObj(), *X_tpetra, beta);
  } else {
    Teuchos::rcp_const_cast<TV>(tpetraVector_.getConstObj())->template sync<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::applyImpl(M_trans, X, Y, alpha, beta);
  }

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


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetraMultiVector(const RCP<MultiVectorBase<Scalar> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

  RCP<TMV> tmv = rcp_dynamic_cast<TMV>(mv);
  if (nonnull(tmv)) {
    return tmv->getTpetraMultiVector();
  }

  RCP<TV> tv = rcp_dynamic_cast<TV>(mv);
  if (nonnull(tv)) {
    return tv->getTpetraVector();
  }

  return Teuchos::null;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

  RCP<const TMV> tmv = rcp_dynamic_cast<const TMV>(mv);
  if (nonnull(tmv)) {
    return tmv->getConstTpetraMultiVector();
  }

  RCP<const TV> tv = rcp_dynamic_cast<const TV>(mv);
  if (nonnull(tv)) {
    return tv->getConstTpetraVector();
  }

  return Teuchos::null;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetraVector(const RCP<VectorBase<Scalar> >& v) const
{
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TV;
  RCP<TV> tv = Teuchos::rcp_dynamic_cast<TV>(v);
  if (nonnull(tv)) {
    return tv->getTpetraVector();
  } else {
    return Teuchos::null;
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraVector(const RCP<const VectorBase<Scalar> >& v) const
{
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TV;
  RCP<const TV> tv = Teuchos::rcp_dynamic_cast<const TV>(v);
  if (nonnull(tv)) {
    return tv->getConstTpetraVector();
  } else {
    return Teuchos::null;
  }
}


} // end namespace Thyra


#endif // THYRA_TPETRA_VECTOR_HPP
