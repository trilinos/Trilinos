// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_HPP
#define THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_HPP


#include "Thyra_DefaultMultiVectorProductVector_decl.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultMultiVectorProductVector<Scalar>::DefaultMultiVectorProductVector()
{
  uninitialize();
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::initialize(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace_in,
  const RCP<MultiVectorBase<Scalar> > &multiVec
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(productSpace_in));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVec));
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorProductVector<Scalar>::initialize(productSpace,multiVec)",
    *multiVec->range(), *productSpace_in->getBlock(0)
    );
  TEUCHOS_ASSERT_EQUALITY( multiVec->domain()->dim(), productSpace_in->numBlocks());
#endif

  numBlocks_ = productSpace_in->numBlocks();

  productSpace_ = productSpace_in;

  multiVec_ = multiVec;

}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::initialize(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace_in,
  const RCP<const MultiVectorBase<Scalar> > &multiVec
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(productSpace_in));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVec));
  THYRA_ASSERT_VEC_SPACES(
    "DefaultMultiVectorProductVector<Scalar>::initialize(productSpace_in,multiVec)",
    *multiVec->range(), *productSpace_in->getBlock(0)
    );
  TEUCHOS_ASSERT_EQUALITY( multiVec->domain()->dim(), productSpace_in->numBlocks() );
#endif

  numBlocks_ = productSpace_in->numBlocks();

  productSpace_ = productSpace_in;

  multiVec_ = multiVec;

}


template <class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstMultiVector()
{
  return multiVec_.getNonconstObj();
}


template <class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getMultiVector() const
{
  return multiVec_.getConstObj();
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::uninitialize()
{
  numBlocks_ = 0;
  productSpace_ = Teuchos::null;
  multiVec_.uninitialize();
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultMultiVectorProductVector<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description()
    << "{"
    << "dim="<<this->space()->dim()
    << ",numColumns = "<<numBlocks_
    << "}";
  return oss.str();
}

template<class Scalar>
void DefaultMultiVectorProductVector<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::describe;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description() << "{"
        << "dim=" << this->space()->dim()
        << "}\n";
      OSTab tab2(out);
      *out <<  "multiVec = " << Teuchos::describe(*multiVec_.getConstObj(),verbLevel);
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// Overridden from ProductVectorBase


template <class Scalar>
RCP<VectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( k, 0, numBlocks_ );
#endif
  return multiVec_.getNonconstObj()->col(k);
}


template <class Scalar>
RCP<const VectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getVectorBlock(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( k, 0, numBlocks_ );
#endif
  return multiVec_.getConstObj()->col(k);
}


// Overridden from ProductMultiVectorBase


template <class Scalar>
RCP<const ProductVectorSpaceBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template <class Scalar>
bool DefaultMultiVectorProductVector<Scalar>::blockIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( k, 0, numBlocks_ );
#else
  (void)k;
#endif
  return multiVec_.isConst();
}


template <class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}


template <class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}


// Overridden public functions from VectorBase


template <class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultMultiVectorProductVector<Scalar>::space() const
{
  return productSpace_;
}


// protected


// Overridden protected functions from VectorBase


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::randomizeImpl(
  Scalar l,
  Scalar u
  )
{
  for(int k = 0; k < numBlocks_; ++k) {
    multiVec_.getNonconstObj()->col(k)->randomize(l, u);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::absImpl(
  const VectorBase<Scalar>& x
  )
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(productSpace_->isCompatible(*x.space()));
#endif
  // Apply operation block-by-block if mv is also a ProductVector
  const RCP<const ProductVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductVectorBase<Scalar> >(Teuchos::rcpFromRef(x));
  if (nonnull(pv)) {
    for(int k = 0; k < numBlocks_; ++k) {
      multiVec_.getNonconstObj()->col(k)->abs(*pv->getVectorBlock(k));
    }
  } else {
    VectorDefaultBase<Scalar>::absImpl(x);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::reciprocalImpl(
  const VectorBase<Scalar>& x
  )
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(productSpace_->isCompatible(*x.space()));
#endif
  // Apply operation block-by-block if mv is also a ProductVector
  const RCP<const ProductVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductVectorBase<Scalar> >(Teuchos::rcpFromRef(x));
  if (nonnull(pv)) {
    for(int k = 0; k < numBlocks_; ++k) {
      multiVec_.getNonconstObj()->col(k)->reciprocal(*pv->getVectorBlock(k));
    }
  } else {
    VectorDefaultBase<Scalar>::reciprocalImpl(x);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::eleWiseScaleImpl(
  const VectorBase<Scalar>& x
  )
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(productSpace_->isCompatible(*x.space()));
#endif
  // Apply operation block-by-block if mv is also a ProductVector
  const RCP<const ProductVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductVectorBase<Scalar> >(Teuchos::rcpFromRef(x));
  if (nonnull(pv)) {
    for(int k = 0; k < numBlocks_; ++k) {
      multiVec_.getNonconstObj()->col(k)->ele_wise_scale(*pv->getVectorBlock(k));
    }
  } else {
    VectorDefaultBase<Scalar>::eleWiseScaleImpl(x);
  }
}


template <class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
DefaultMultiVectorProductVector<Scalar>::norm2WeightedImpl(
  const VectorBase<Scalar>& x
  ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT(productSpace_->isCompatible(*x.space()));
#endif
  // Apply operation block-by-block if mv is also a ProductVector
  typedef ScalarTraits<Scalar> ST;
  const RCP<const ProductVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductVectorBase<Scalar> >(Teuchos::rcpFromRef(x));
  if (nonnull(pv)) {
    typename ST::magnitudeType norm = ST::magnitude(ST::zero());
    for(int k = 0; k < numBlocks_; ++k) {
      typename ST::magnitudeType subNorm
        = multiVec_.getConstObj()->col(k)->norm_2(*pv->getVectorBlock(k));
      norm += subNorm*subNorm;
    }
    return ST::magnitude(ST::squareroot(norm));
  } else {
    return VectorDefaultBase<Scalar>::norm2WeightedImpl(x);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset
  ) const
{
  this->getDefaultProductVector()->applyOp(
    op, vecs, targ_vecs, reduct_obj, global_offset );
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::acquireDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  this->getDefaultProductVector()->acquireDetachedView(rng_in,sub_vec);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  this->getDefaultProductVector()->releaseDetachedView(sub_vec);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::acquireNonconstDetachedVectorViewImpl(
  const Range1D& /* rng_in */, RTOpPack::SubVectorView<Scalar>* /* sub_vec */
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::acquireNonconstDetachedVectorViewImpl(...)!");
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* /* sub_vec */
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::commitNonconstDetachedVectorViewImpl(...)!");
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::setSubVectorImpl(
  const RTOpPack::SparseSubVectorT<Scalar>& /* sub_vec */
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement DefaultMultiVectorProductVector<Scalar>::setSubVector(...)!");
}


// Overridden protected functions from VectorBase


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::assignImpl(Scalar alpha)
{
  multiVec_.getNonconstObj()->assign(alpha);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::assignMultiVecImpl(
  const MultiVectorBase<Scalar>& mv
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(mv.domain()->dim(), 1);
  TEUCHOS_ASSERT(productSpace_->isCompatible(*mv.range()));
#endif
  const RCP<const ProductMultiVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(Teuchos::rcpFromRef(mv));
  if (nonnull(pv)) {
    for(int k = 0; k < numBlocks_; ++k) {
      multiVec_.getNonconstObj()->col(k)->assign(*pv->getMultiVectorBlock(k));
    }
  } else {
    MultiVectorDefaultBase<Scalar>::assignMultiVecImpl(mv);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::scaleImpl(Scalar alpha)
{
  multiVec_.getNonconstObj()->scale(alpha);
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::updateImpl(
  Scalar alpha,
  const MultiVectorBase<Scalar>& mv
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(mv.domain()->dim(), 1);
  TEUCHOS_ASSERT(productSpace_->isCompatible(*mv.range()));
#endif
  const RCP<const ProductMultiVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(Teuchos::rcpFromRef(mv));
  if (nonnull(pv)) {
    for(int k = 0; k < numBlocks_; ++k) {
      multiVec_.getNonconstObj()->col(k)->update(alpha, *pv->getMultiVectorBlock(k));
    }
  } else {
    MultiVectorDefaultBase<Scalar>::updateImpl(alpha, mv);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::linearCombinationImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
  const Scalar& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
  for (Ordinal i = 0; i < mv.size(); ++i) {
    TEUCHOS_ASSERT_EQUALITY(mv[i]->domain()->dim(), 1);
    TEUCHOS_ASSERT(productSpace_->isCompatible(*mv[i]->range()));
  }
#endif

  // Apply operation block-by-block if each element of mv is also a ProductMultiVector
  bool allCastsSuccessful = true;
  Array<Ptr<const ProductMultiVectorBase<Scalar> > > pvs(mv.size());
  for (Ordinal i = 0; i < mv.size(); ++i) {
    pvs[i] = Teuchos::ptr_dynamic_cast<const ProductMultiVectorBase<Scalar> >(mv[i]);
    if (pvs[i].is_null()) {
      allCastsSuccessful = false;
      break;
    }
  }

  if (allCastsSuccessful) {
    Array<RCP<const MultiVectorBase<Scalar> > > blocks_rcp(pvs.size());
    Array<Ptr<const MultiVectorBase<Scalar> > > blocks(pvs.size());
    for (int k = 0; k < numBlocks_; ++k) {
      for (Ordinal i = 0; i < pvs.size(); ++i) {
        blocks_rcp[i] = pvs[i]->getMultiVectorBlock(k);
        blocks[i] = blocks_rcp[i].ptr();
      }
      multiVec_.getNonconstObj()->col(k)->linear_combination(alpha, blocks(), beta);
    }
  } else {
    MultiVectorDefaultBase<Scalar>::linearCombinationImpl(alpha, mv, beta);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::dotsImpl(
  const MultiVectorBase<Scalar>& mv,
  const ArrayView<Scalar>& prods
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(mv.domain()->dim(), 1);
  TEUCHOS_ASSERT_EQUALITY(prods.size(), 1);
  TEUCHOS_ASSERT(productSpace_->isCompatible(*mv.range()));
#endif
  const RCP<const ProductMultiVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(Teuchos::rcpFromRef(mv));
  if (nonnull(pv)) {
    prods[0] = ScalarTraits<Scalar>::zero();
    for(int k = 0; k < numBlocks_; ++k) {
      Scalar prod;
      multiVec_.getConstObj()->col(k)->dots(*pv->getMultiVectorBlock(k), Teuchos::arrayView(&prod, 1));
      prods[0] += prod;
    }
  } else {
    MultiVectorDefaultBase<Scalar>::dotsImpl(mv, prods);
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::norms1Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_EQUALITY(norms.size(), 1);
#endif
  typedef ScalarTraits<Scalar> ST;
  norms[0] = ST::magnitude(ST::zero());
  for(int k = 0; k < numBlocks_; ++k) {
    norms[0] += multiVec_.getConstObj()->col(k)->norm_1();
  }
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::norms2Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_EQUALITY(norms.size(), 1);
#endif
  typedef ScalarTraits<Scalar> ST;
  norms[0] = ST::magnitude(ST::zero());
  for(int k = 0; k < numBlocks_; ++k) {
    typename ST::magnitudeType subNorm = multiVec_.getConstObj()->col(k)->norm_2();
    norms[0] += subNorm*subNorm;
  }
  norms[0] = ST::magnitude(ST::squareroot(norms[0]));
}


template <class Scalar>
void DefaultMultiVectorProductVector<Scalar>::normsInfImpl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_EQUALITY(norms.size(), 1);
#endif
  typedef ScalarTraits<Scalar> ST;
  norms[0] = ST::magnitude(ST::zero());
  for(int k = 0; k < numBlocks_; ++k) {
    typename ST::magnitudeType subNorm = multiVec_.getConstObj()->col(k)->norm_inf();
    if (subNorm > norms[0])
      norms[0] = subNorm;
  }
}


// private


template <class Scalar>
RCP<const DefaultProductVector<Scalar> >
DefaultMultiVectorProductVector<Scalar>::getDefaultProductVector() const
{

  // This function exists since in general we can not create views of a column
  // vectors and expect the changes to be mirrored in the mulit-vector
  // automatically.  Later, we might be able to change this once we have a
  // Thyra::MultiVectorBase::hasDirectColumnVectorView() function and it
  // returns true.  Until then, this is the safe way to do this ...

  Array<RCP<const VectorBase<Scalar> > > vecArray;
  for ( int k = 0; k < numBlocks_; ++k) {
    vecArray.push_back(multiVec_.getConstObj()->col(k));
  }

  return Thyra::defaultProductVector<Scalar>(
    productSpace_->getDefaultProductVectorSpace(),
    vecArray()
    );

}


} // namespace Thyra


#endif // THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_HPP
