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
  this->assignMultiVecImpl(x);
  /*auto tx = this->getConstTpetraVector(Teuchos::rcpFromRef(x));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tx)) {
    tpetraVector_.getNonconstObj()->assign(*tx);
  } else {
    // This version will require/modify the host view of this vector.
    tpetraVector_.getNonconstObj()->template sync<Kokkos::HostSpace>();
    tpetraVector_.getNonconstObj()->template modify<Kokkos::HostSpace>();
    VectorDefaultBase<Scalar>::assignVecImpl(x);
  }*/
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
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateVecImpl(
  Scalar alpha,
  const VectorBase<Scalar>& x
  )
{
  this->updateImpl(alpha, x);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::linearCombinationVecImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const VectorBase<Scalar> > >& x,
  const Scalar& beta
  )
{
  Teuchos::Array<Ptr<const MultiVectorBase<Scalar> > > mv(x.size());
  for (typename Teuchos::Array<Ptr<const MultiVectorBase<Scalar> > >::size_type i = 0; i < mv.size(); ++i)
    mv[i] = Teuchos::ptr_static_cast<const MultiVectorBase<Scalar> >(x[i]);
  this->linearCombinationImpl(alpha, mv, beta);
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
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInfImpl() const
{
  return tpetraVector_.getConstObj()->normInf();
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
      // Tpetra frequently const_casts to sync, so I'm assuming this is alright
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

  // If each cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (allCastsSuccessful) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    auto tmvIter = tmvs.begin();
    auto alphaIter = alpha.begin();
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


/*template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

  RCP<const TpetraMultiVector_t> X_tpetra;
  if (nonnull(Teuchos::ptr_dynamic_cast<const TMV>(Teuchos::ptrFromRef(X)))) {
    X_tpetra = dynamic_cast<const TMV&>(X).getConstTpetraMultiVector();
  } else if (nonnull(Teuchos::ptr_dynamic_cast<const TV>(Teuchos::ptrFromRef(X)))) {
    X_tpetra = dynamic_cast<const TV&>(X).getConstTpetraVector();
  }

  RCP<TpetraMultiVector_t> Y_tpetra;
  if (nonnull(Teuchos::ptr_dynamic_cast<TMV>(Y))) {
    Y_tpetra = Teuchos::ptr_dynamic_cast<TMV>(Y)->getTpetraMultiVector();
  } else if (nonnull(Teuchos::ptr_dynamic_cast<TV>(Y))) {
    Y_tpetra = Teuchos::ptr_dynamic_cast<TV>(Y)->getTpetraVector();
  }

  // If we succeeded, call Tpetra directly,
  // otherwise fall back on the default implementation
  if (nonnull(X_tpetra) && nonnull(Y_tpetra)) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(X_tpetra->getNumVectors() == Y_tpetra->getNumVectors());
#endif
    if (M_trans == NOTRANS || (M_trans == CONJ && !ST::isComplex) ) {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(tpetraVector_.getConstObj()->getMap()->isCompatible(*Y_tpetra->getMap()));
      TEUCHOS_ASSERT(X_tpetra->getLocalLength() == 1);
      TEUCHOS_ASSERT(!(X_tpetra->isDistributed()));
#endif
      // Do this without copying X to host?
      for (std::size_t i = 0; i < Y_tpetra->getNumVectors(); ++i) {
        Y_tpetra->getVectorNonConst(i)->update(alpha * X_tpetra->getData(i)[0],
          *tpetraVector_.getConstObj(), beta);
      }
    } else if(M_trans == CONJTRANS || (M_trans == TRANS && !ST::isComplex) ) {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(tpetraVector_.getConstObj()->getMap()->isCompatible(*X_tpetra->getMap()));
      TEUCHOS_ASSERT(Y_tpetra->getLocalLength() == 1);
      TEUCHOS_ASSERT(!(Y_tpetra->isDistributed()));
#endif
      // Do this without copying Y and scalar prods to host?
      for (std::size_t i = 0; i < Y_tpetra->getNumVectors(); ++i) {
        Scalar& value = Y_tpetra->getDataNonConst(i)[0];
        // Argument of dot conjugated if complex.
        // For complex Scalar, getting an error adding Kokkos:complex and
        // std::complex, so need to explicitly cast the dot result.
        value = alpha * static_cast<Scalar>(
          X_tpetra->getVector(i)->dot(*tpetraVector_.getConstObj())) + beta * value;
      }
      Y_tpetra->template modify<Kokkos::HostSpace>();
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "TpetraVector<"<<ST::name()<<">::apply(M_trans,...): Error, M_trans="
        <<toString(M_trans)<<" not supported!" );
    }
  } else {
    VectorDefaultBase<Scalar>::applyImpl(M_trans, X, Y, alpha, beta);
  }
}*/


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
