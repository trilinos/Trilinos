// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_TPETRA_MULTIVECTOR_HPP
#define THYRA_TPETRA_MULTIVECTOR_HPP

#include "Thyra_TpetraMultiVector_decl.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TpetraMultiVector()
{}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
  )
{
  initializeImpl(tpetraVectorSpace, domainSpace, tpetraMultiVector);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::constInitialize(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
  )
{
  initializeImpl(tpetraVectorSpace, domainSpace, tpetraMultiVector);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getTpetraMultiVector()
{
  return tpetraMultiVector_.getNonconstObj();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraMultiVector() const
{
  return tpetraMultiVector_;
}


// Overridden public functions form MultiVectorAdapterBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP< const ScalarProdVectorSpaceBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::domainScalarProdVecSpc() const
{
  return domainSpace_;
}


// Overridden protected functions from MultiVectorBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::assignImpl(Scalar alpha)
{
  tpetraMultiVector_.getNonconstObj()->putScalar(alpha);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
assignMultiVecImpl(const MultiVectorBase<Scalar>& mv)
{
  auto tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tmv)) {
    tpetraMultiVector_.getNonconstObj()->assign(*tmv);
  } else {
    MultiVectorDefaultBase<Scalar>::assignMultiVecImpl(mv);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::scaleImpl(Scalar alpha)
{
  tpetraMultiVector_.getNonconstObj()->scale(alpha);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateImpl(
  Scalar alpha,
  const MultiVectorBase<Scalar>& mv
  )
{
  auto tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tmv)) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    tpetraMultiVector_.getNonconstObj()->update(alpha, *tmv, ST::one());
  } else {
    MultiVectorDefaultBase<Scalar>::updateImpl(alpha, mv);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::linearCombinationImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
  const Scalar& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
#endif

  // Try to cast mv to an array of this type
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  Teuchos::Array<RCP<const TMV> > tmvs(mv.size());
  RCP<const TMV> tmv;
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
  auto len = tmvs.size();
  if (len == 0) {
    tpetraMultiVector_.getNonconstObj()->scale(beta);
  } else if (len == 1 && allCastsSuccessful) {
    tpetraMultiVector_.getNonconstObj()->update(alpha[0], *tmvs[0], beta);
  } else if (len == 2 && allCastsSuccessful) {
    tpetraMultiVector_.getNonconstObj()->update(alpha[0], *tmvs[0], alpha[1], *tmvs[1], beta);
  } else if (allCastsSuccessful) {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    auto tmvIter = tmvs.begin();
    auto alphaIter = alpha.begin();

    // Check if any entry of tmvs aliases this object's wrapped vector.
    // If so, replace that entry in the array with a copy.
    tmv = Teuchos::null;
    for (; tmvIter != tmvs.end(); ++tmvIter) {
      if (tmvIter->getRawPtr() == tpetraMultiVector_.getConstObj().getRawPtr()) {
        if (tmv.is_null()) {
          tmv = Teuchos::rcp(new TMV(*tpetraMultiVector_.getConstObj(), Teuchos::Copy));
        }
        *tmvIter = tmv;
      }
    }
    tmvIter = tmvs.begin();

    // We add two MVs at a time, so only scale if even num MVs,
    // and additionally do the first addition if odd num MVs.
    if ((tmvs.size() % 2) == 0) {
      tpetraMultiVector_.getNonconstObj()->scale(beta);
    } else {
      tpetraMultiVector_.getNonconstObj()->update(*alphaIter, *(*tmvIter), beta);
      ++tmvIter;
      ++alphaIter;
    }
    for (; tmvIter != tmvs.end(); tmvIter+=2, alphaIter+=2) {
      tpetraMultiVector_.getNonconstObj()->update(
        *alphaIter, *(*tmvIter), *(alphaIter+1), *(*(tmvIter+1)), ST::one());
    }
  } else {
    MultiVectorDefaultBase<Scalar>::linearCombinationImpl(alpha, mv, beta);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dotsImpl(
    const MultiVectorBase<Scalar>& mv,
    const ArrayView<Scalar>& prods
    ) const
{
  auto tmv = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(mv));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the RTOp implementation.
  if (nonnull(tmv)) {
    tpetraMultiVector_.getConstObj()->dot(*tmv, prods);
  } else {
    MultiVectorDefaultBase<Scalar>::dotsImpl(mv, prods);
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norms1Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  tpetraMultiVector_.getConstObj()->norm1(norms);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norms2Impl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const
{
  tpetraMultiVector_.getConstObj()->norm2(norms);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normsInfImpl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  tpetraMultiVector_.getConstObj()->normInf(norms);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const VectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::colImpl(Ordinal j) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, this->domain()->dim());
#endif
  return constTpetraVector<Scalar>(
    tpetraVectorSpace_,
    tpetraMultiVector_->getVector(j)
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<VectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconstColImpl(Ordinal j)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, this->domain()->dim());
#endif
  return tpetraVector<Scalar>(
    tpetraVectorSpace_,
    tpetraMultiVector_.getNonconstObj()->getVectorNonConst(j)
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const MultiVectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::contigSubViewImpl(
  const Range1D& col_rng_in
  ) const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nTpetraMultiVector::subView(Range1D) const called!\n";
#endif
  const Range1D colRng = this->validateColRange(col_rng_in);

  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraView =
    this->getConstTpetraMultiVector()->subView(colRng);

  const RCP<const ScalarProdVectorSpaceBase<Scalar> > viewDomainSpace =
    tpetraVectorSpace<Scalar>(
        Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
          tpetraView->getNumVectors(),
          tpetraView->getMap()->getComm()
          )
        );

  return constTpetraMultiVector(
      tpetraVectorSpace_,
      viewDomainSpace,
      tpetraView
      );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MultiVectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconstContigSubViewImpl(
  const Range1D& col_rng_in
  )
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nTpetraMultiVector::subView(Range1D) called!\n";
#endif
  const Range1D colRng = this->validateColRange(col_rng_in);

  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraView =
    this->getTpetraMultiVector()->subViewNonConst(colRng);

  const RCP<const ScalarProdVectorSpaceBase<Scalar> > viewDomainSpace =
    tpetraVectorSpace<Scalar>(
        Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
          tpetraView->getNumVectors(),
          tpetraView->getMap()->getComm()
          )
        );

  return tpetraMultiVector(
      tpetraVectorSpace_,
      viewDomainSpace,
      tpetraView
      );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const MultiVectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonContigSubViewImpl(
  const ArrayView<const int>& cols_in
  ) const
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nTpetraMultiVector::subView(ArrayView) const called!\n";
#endif
  // Tpetra wants col indices as size_t
  Array<std::size_t> cols(cols_in.size());
  for (Array<std::size_t>::size_type i = 0; i < cols.size(); ++i)
    cols[i] = static_cast<std::size_t>(cols_in[i]);

  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraView =
    this->getConstTpetraMultiVector()->subView(cols());

  const RCP<const ScalarProdVectorSpaceBase<Scalar> > viewDomainSpace =
    tpetraVectorSpace<Scalar>(
        Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
          tpetraView->getNumVectors(),
          tpetraView->getMap()->getComm()
          )
        );

  return constTpetraMultiVector(
      tpetraVectorSpace_,
      viewDomainSpace,
      tpetraView
      );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MultiVectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconstNonContigSubViewImpl(
  const ArrayView<const int>& cols_in
  )
{
#ifdef THYRA_DEFAULT_SPMD_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nTpetraMultiVector::subView(ArrayView) called!\n";
#endif
  // Tpetra wants col indices as size_t
  Array<std::size_t> cols(cols_in.size());
  for (Array<std::size_t>::size_type i = 0; i < cols.size(); ++i)
    cols[i] = static_cast<std::size_t>(cols_in[i]);

  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraView =
    this->getTpetraMultiVector()->subViewNonConst(cols());

  const RCP<const ScalarProdVectorSpaceBase<Scalar> > viewDomainSpace =
    tpetraVectorSpace<Scalar>(
        Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(
          tpetraView->getNumVectors(),
          tpetraView->getMap()->getComm()
          )
        );

  return tpetraMultiVector(
      tpetraVectorSpace_,
      viewDomainSpace,
      tpetraView
      );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
mvMultiReductApplyOpImpl(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal primary_global_offset
  ) const
{

  MultiVectorAdapterBase<Scalar>::mvMultiReductApplyOpImpl(
    primary_op, multi_vecs, targ_multi_vecs, reduct_objs, primary_global_offset);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
acquireDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
  ) const
{
  SpmdMultiVectorDefaultBase<Scalar>::
    acquireDetachedMultiVectorViewImpl(rowRng, colRng, sub_mv);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
acquireNonconstDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::SubMultiVectorView<Scalar>* sub_mv
  )
{
  SpmdMultiVectorDefaultBase<Scalar>::
    acquireNonconstDetachedMultiVectorViewImpl(rowRng, colRng, sub_mv);
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
commitNonconstDetachedMultiVectorViewImpl(
  RTOpPack::SubMultiVectorView<Scalar>* sub_mv
  )
{
  SpmdMultiVectorDefaultBase<Scalar>::
    commitNonconstDetachedMultiVectorViewImpl(sub_mv);

}


/* ToDo: Implement these?


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const MultiVectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonContigSubViewImpl(
  const ArrayView<const int> &cols
  ) const
{
  THYRA_DEBUG_ASSERT_MV_COLS("nonContigSubViewImpl(cols)", cols);
  const int numCols = cols.size();
  const ArrayRCP<Scalar> localValuesView = createContiguousCopy(cols);
  return defaultSpmdMultiVector<Scalar>(
    spmdRangeSpace_,
    createSmallScalarProdVectorSpaceBase<Scalar>(spmdRangeSpace_, numCols),
    localValuesView
    );
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MultiVectorBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconstNonContigSubViewImpl(
  const ArrayView<const int> &cols )
{
  THYRA_DEBUG_ASSERT_MV_COLS("nonContigSubViewImpl(cols)", cols);
  const int numCols = cols.size();
  const ArrayRCP<Scalar> localValuesView = createContiguousCopy(cols);
  const Ordinal localSubDim = spmdRangeSpace_->localSubDim();
  RCP<CopyBackSpmdMultiVectorEntries<Scalar> > copyBackView =
    copyBackSpmdMultiVectorEntries<Scalar>(cols, localValuesView.getConst(),
      localSubDim, localValues_.create_weak(), leadingDim_);
  return Teuchos::rcpWithEmbeddedObjPreDestroy(
    new TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
      spmdRangeSpace_,
      createSmallScalarProdVectorSpaceBase<Scalar>(spmdRangeSpace_, numCols),
      localValuesView),
    copyBackView
    );
}

*/


// Overridden protected members from SpmdMultiVectorBase


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const SpmdVectorSpaceBase<Scalar> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::spmdSpaceImpl() const
{
  return tpetraVectorSpace_;
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstLocalMultiVectorDataImpl(
  const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
  )
{
  *localValues = tpetraMultiVector_.getNonconstObj()->get1dViewNonConst();
  *leadingDim = tpetraMultiVector_->getStride();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMultiVectorDataImpl(
  const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
  ) const
{
  *localValues = tpetraMultiVector_->get1dView();
  *leadingDim = tpetraMultiVector_->getStride();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::euclideanApply(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  // Try to extract Tpetra objects from X and Y
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  Teuchos::RCP<const TMV> X_tpetra = this->getConstTpetraMultiVector(Teuchos::rcpFromRef(X));
  Teuchos::RCP<TMV> Y_tpetra = this->getTpetraMultiVector(Teuchos::rcpFromPtr(Y));

  // If the cast succeeded, call Tpetra directly.
  // Otherwise, fall back to the default implementation.
  if (nonnull(X_tpetra) && nonnull(Y_tpetra)) {
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

    Y_tpetra->multiply(trans, Teuchos::NO_TRANS, alpha, *tpetraMultiVector_.getConstObj(), *X_tpetra, beta);
    Kokkos::fence();
  } else {
    SpmdMultiVectorDefaultBase<Scalar>::euclideanApply(M_trans, X, Y, alpha, beta);
  }

}

// private


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template<class TpetraMultiVector_t>
void TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initializeImpl(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<TpetraMultiVector_t> &tpetraMultiVector
  )
{
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(tpetraVectorSpace));
  TEUCHOS_ASSERT(nonnull(domainSpace));
  TEUCHOS_ASSERT(nonnull(tpetraMultiVector));
  // ToDo: Check to make sure that tpetraMultiVector is compatible with
  // tpetraVectorSpace.
#endif
  tpetraVectorSpace_ = tpetraVectorSpace;
  domainSpace_ = domainSpace;
  tpetraMultiVector_.initialize(tpetraMultiVector);
  this->updateSpmdSpace();
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetraMultiVector(const RCP<MultiVectorBase<Scalar> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

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
TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const
{
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TV;

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


} // end namespace Thyra


#endif // THYRA_TPETRA_MULTIVECTOR_HPP
