// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DEF_HPP
#define THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DEF_HPP

#include "Thyra_DefaultProductMultiVector_decl.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultProductMultiVector<Scalar>::DefaultProductMultiVector()
  :numBlocks_(0)
{}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::initialize(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
  const int numMembers
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( is_null(productSpace_in) );
  TEUCHOS_TEST_FOR_EXCEPT( numMembers <= 0 );
#endif
  Array<RCP<MultiVectorBase<Scalar> > > multiVecs;
  const int numBlocks = productSpace_in->numBlocks();
  for ( int k = 0; k < numBlocks; ++k )
    multiVecs.push_back(createMembers(productSpace_in->getBlock(k),numMembers));
  initialize(productSpace_in,multiVecs);
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::initialize(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
  const ArrayView<const RCP<MultiVectorBase<Scalar> > > &multiVecs
  )
{
  initializeImpl(productSpace_in,multiVecs);
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::initialize(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
  const ArrayView<const RCP<const MultiVectorBase<Scalar> > > &multiVecs
  )
{
  initializeImpl(productSpace_in,multiVecs);
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::uninitialize()
{
  productSpace_ = Teuchos::null;
  multiVecs_.resize(0);
  numBlocks_ = 0;
}


// Overridden public functions from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultProductMultiVector<Scalar>::description() const
{
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description()
    << "{"
    << "rangeDim="<<this->range()->dim()
    << ",domainDim="<<this->domain()->dim()
    << ",numBlocks = "<<numBlocks_
    << "}";
  return oss.str();
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::describe(
  FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::describe;
  if (verbLevel == Teuchos::VERB_NONE)
    return;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_NONE:
      break;
    case Teuchos::VERB_DEFAULT: // fall through
    case Teuchos::VERB_LOW: // fall through
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM: // fall through
    case Teuchos::VERB_HIGH: // fall through
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description() << "{"
        << "rangeDim="<<this->range()->dim()
        << ",domainDim="<<this->domain()->dim()
        << "}\n";
      OSTab tab2(out);
      *out
        <<  "numBlocks="<< numBlocks_ << std::endl
        <<  "Constituent multi-vector objects V[0], V[1], ... V[numBlocks-1]:\n";
      OSTab tab3(out);
      for( int k = 0; k < numBlocks_; ++k ) {
        *out << "V["<<k<<"] = "
             << Teuchos::describe(*multiVecs_[k].getConstObj(),verbLevel);
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// Overridden public functions from ProductMultiVectorBase


template<class Scalar>
RCP<const ProductVectorSpaceBase<Scalar> >
DefaultProductMultiVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template<class Scalar>
bool DefaultProductMultiVector<Scalar>::blockIsConst(const int k) const
{
  return multiVecs_[k].isConst();
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return multiVecs_[k].getNonconstObj();
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return multiVecs_[k].getConstObj();
}


// Overridden public functions from MultiVectorBase


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::clone_mv() const
{
  assertInitialized();
  Array<RCP<MultiVectorBase<Scalar> > > blocks;
  for ( int k = 0; k < numBlocks_; ++k )
    blocks.push_back(multiVecs_[k].getConstObj()->clone_mv());
  return defaultProductMultiVector<Scalar>(productSpace_, blocks());
}


// Overriden public functions from LinearOpBase


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultProductMultiVector<Scalar>::range() const
{
  return productSpace_;
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultProductMultiVector<Scalar>::domain() const
{
  if (is_null(productSpace_))
    return Teuchos::null;
  return multiVecs_[0].getConstObj()->domain();
}


// protected


// Overriden protected functions from MultiVectorBase


template<class Scalar>
void DefaultProductMultiVector<Scalar>::assignImpl(Scalar alpha)
{
  for ( int k = 0; k < numBlocks_; ++k ) {
    multiVecs_[k].getNonconstObj()->assign(alpha);
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::assignMultiVecImpl(
  const MultiVectorBase<Scalar>& mv
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), mv.domain()->dim());
  TEUCHOS_ASSERT(productSpace_->isCompatible(*mv.range()));
#endif

  // Apply operation block-by-block if mv is also a ProductMultiVector
  const RCP<const ProductMultiVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(Teuchos::rcpFromRef(mv));
  if (nonnull(pv)) {
    for (int k = 0; k < numBlocks_; ++k) {
      multiVecs_[k].getNonconstObj()->assign(*pv->getMultiVectorBlock(k));
    }
  } else {
    MultiVectorDefaultBase<Scalar>::assignMultiVecImpl(mv);
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::scaleImpl(Scalar alpha)
{
  for (int k = 0; k < numBlocks_; ++k) {
    multiVecs_[k].getNonconstObj()->scale(alpha);
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::updateImpl(
  Scalar alpha,
  const MultiVectorBase<Scalar>& mv
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), mv.domain()->dim());
  TEUCHOS_ASSERT(productSpace_->isCompatible(*mv.range()));
#endif

  // Apply operation block-by-block if mv is also a ProductMultiVector
  const RCP<const ProductMultiVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(Teuchos::rcpFromRef(mv));
  if (nonnull(pv)) {
    for (int k = 0; k < numBlocks_; ++k) {
      multiVecs_[k].getNonconstObj()->update(alpha, *pv->getMultiVectorBlock(k));
    }
  } else {
    MultiVectorDefaultBase<Scalar>::updateImpl(alpha, mv);
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::linearCombinationImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
  const Scalar& beta
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
  for (Ordinal i = 0; i < mv.size(); ++i) {
    TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), mv[i]->domain()->dim());
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
      multiVecs_[k].getNonconstObj()->linear_combination(alpha, blocks(), beta);
    }
  } else {
    MultiVectorDefaultBase<Scalar>::linearCombinationImpl(alpha, mv, beta);
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::dotsImpl(
  const MultiVectorBase<Scalar>& mv,
  const ArrayView<Scalar>& prods
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), mv.domain()->dim());
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), prods.size());
  TEUCHOS_ASSERT(productSpace_->isCompatible(*mv.range()));
#endif
  // Apply operation block-by-block if mv is also a ProductMultiVector
  const RCP<const ProductMultiVectorBase<Scalar> > pv
    = Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(Teuchos::rcpFromRef(mv));
  if (nonnull(pv)) {
    for (Ordinal i = 0; i < prods.size(); ++i)
      prods[i] = ScalarTraits<Scalar>::zero();

    Array<Scalar> subProds(prods.size());
    for (int k = 0; k < numBlocks_; ++k) {
      multiVecs_[k].getConstObj()->dots(*pv->getMultiVectorBlock(k), subProds());
      for (Ordinal i = 0; i < prods.size(); ++i) {
        prods[i] += subProds[i];
      }
    }
  } else {
    MultiVectorDefaultBase<Scalar>::dotsImpl(mv, prods);
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::norms1Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), norms.size());
#endif
  typedef ScalarTraits<Scalar> ST;
  for (Ordinal i = 0; i < norms.size(); ++i)
    norms[i] = ST::magnitude(ST::zero());

  Array<typename ST::magnitudeType> subNorms(norms.size());
  for (int k = 0; k < numBlocks_; ++k) {
    multiVecs_[k].getConstObj()->norms_1(subNorms());
    for (Ordinal i = 0; i < norms.size(); ++i) {
      norms[i] += subNorms[i];
    }
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::norms2Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), norms.size());
#endif
  typedef ScalarTraits<Scalar> ST;
  for (Ordinal i = 0; i < norms.size(); ++i)
    norms[i] = ST::magnitude(ST::zero());

  Array<typename ST::magnitudeType> subNorms(norms.size());
  for (int k = 0; k < numBlocks_; ++k) {
    multiVecs_[k].getConstObj()->norms_2(subNorms());
    for (Ordinal i = 0; i < norms.size(); ++i) {
        norms[i] += subNorms[i]*subNorms[i];
    }
  }

  for (Ordinal i = 0; i < norms.size(); ++i) {
    norms[i] = ST::magnitude(ST::squareroot(norms[i]));
  }
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::normsInfImpl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(this->domain()->dim(), norms.size());
#endif
  typedef ScalarTraits<Scalar> ST;
  for (Ordinal i = 0; i < norms.size(); ++i)
    norms[i] = ST::magnitude(ST::zero());

  Array<typename ST::magnitudeType> subNorms(norms.size());
  for (int k = 0; k < numBlocks_; ++k) {
    multiVecs_[k].getConstObj()->norms_inf(subNorms());
    for (Ordinal i = 0; i < norms.size(); ++i) {
      if (subNorms[i] > norms[i])
        norms[i] = subNorms[i];
    }
  }
}


template<class Scalar>
RCP<const VectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::colImpl(Ordinal j) const
{
  validateColIndex(j);
  Array<RCP<const VectorBase<Scalar> > > cols_;
  for ( int k = 0; k < numBlocks_; ++k )
    cols_.push_back(multiVecs_[k].getConstObj()->col(j));
  return defaultProductVector<Scalar>(productSpace_, cols_());
}


template<class Scalar>
RCP<VectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::nonconstColImpl(Ordinal j)
{
  validateColIndex(j);
  Array<RCP<VectorBase<Scalar> > > cols_;
  for ( int k = 0; k < numBlocks_; ++k )
    cols_.push_back(multiVecs_[k].getNonconstObj()->col(j));
  return defaultProductVector<Scalar>(productSpace_, cols_());
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::contigSubViewImpl( const Range1D& colRng ) const
{
  assertInitialized();
  Array<RCP<const MultiVectorBase<Scalar> > > blocks;
  for ( int k = 0; k < numBlocks_; ++k )
    blocks.push_back(multiVecs_[k].getConstObj()->subView(colRng));
  return defaultProductMultiVector<Scalar>(productSpace_, blocks());
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::nonconstContigSubViewImpl( const Range1D& colRng )
{
  assertInitialized();
  Array<RCP<MultiVectorBase<Scalar> > > blocks;
  for ( int k = 0; k < numBlocks_; ++k )
    blocks.push_back(multiVecs_[k].getNonconstObj()->subView(colRng));
  return defaultProductMultiVector<Scalar>(productSpace_, blocks());
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::nonContigSubViewImpl(
  const ArrayView<const int> &cols
  ) const
{
  assertInitialized();
  Array<RCP<const MultiVectorBase<Scalar> > > blocks;
  for ( int k = 0; k < numBlocks_; ++k )
    blocks.push_back(multiVecs_[k].getConstObj()->subView(cols));
  return defaultProductMultiVector<Scalar>(productSpace_, blocks());
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultProductMultiVector<Scalar>::nonconstNonContigSubViewImpl(
  const ArrayView<const int> &cols
  )
{
  assertInitialized();
  Array<RCP<MultiVectorBase<Scalar> > > blocks;
  for ( int k = 0; k < numBlocks_; ++k )
    blocks.push_back(multiVecs_[k].getNonconstObj()->subView(cols));
  return defaultProductMultiVector<Scalar>(productSpace_, blocks());
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::mvMultiReductApplyOpImpl(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs_in,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs_inout,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal primary_global_offset_in
  ) const
{

  using Teuchos::ptr_dynamic_cast;
  using Thyra::applyOp;

  assertInitialized();

#ifdef TEUCHOS_DEBUG
  for ( int j = 0; j < multi_vecs_in.size(); ++j ) {
    THYRA_ASSERT_VEC_SPACES(
      "DefaultProductMultiVector<Scalar>::mvMultiReductApplyOpImpl(...)",
      *this->range(), *multi_vecs_in[j]->range()
      );
    THYRA_ASSERT_VEC_SPACES(
      "DefaultProductMultiVector<Scalar>::mvMultiReductApplyOpImpl(...)",
      *this->domain(), *multi_vecs_in[j]->domain()
      );
  }
  for ( int j = 0; j < targ_multi_vecs_inout.size(); ++j ) {
    THYRA_ASSERT_VEC_SPACES(
      "DefaultProductMultiVector<Scalar>::mvMultiReductApplyOpImpl(...)",
      *this->range(), *targ_multi_vecs_inout[j]->range()
      );
    THYRA_ASSERT_VEC_SPACES(
      "DefaultProductMultiVector<Scalar>::mvMultiReductApplyOpImpl(...)",
      *this->domain(), *targ_multi_vecs_inout[j]->domain()
      );
  }
#endif //  TEUCHOS_DEBUG

  //
  // Try to dynamic cast all of the multi-vector objects to the
  // ProductMultiVectoBase interface.
  //

  bool allProductMultiVectors = true;

  Array<Ptr<const ProductMultiVectorBase<Scalar> > > multi_vecs;
  for ( int j = 0; j < multi_vecs_in.size() && allProductMultiVectors; ++j ) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( is_null(multi_vecs_in[j]) );
#endif
    const Ptr<const ProductMultiVectorBase<Scalar> >
      multi_vecs_j = ptr_dynamic_cast<const ProductMultiVectorBase<Scalar> >(
        multi_vecs_in[j]
        );
    if ( !is_null(multi_vecs_j) ) {
      multi_vecs.push_back(multi_vecs_j);
    }
    else {
      allProductMultiVectors = false;
    }
  }

  Array<Ptr<ProductMultiVectorBase<Scalar> > > targ_multi_vecs;
  for ( int j = 0; j < targ_multi_vecs_inout.size() && allProductMultiVectors; ++j )
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( is_null(targ_multi_vecs_inout[j]) );
#endif
    Ptr<ProductMultiVectorBase<Scalar> >
      targ_multi_vecs_j = ptr_dynamic_cast<ProductMultiVectorBase<Scalar> >(
        targ_multi_vecs_inout[j]
        );
    if (!is_null(targ_multi_vecs_j)) {
      targ_multi_vecs.push_back(targ_multi_vecs_j);
    }
    else {
      allProductMultiVectors = false;
    }
  }

  //
  // Do the reduction operations
  //
  
  if ( allProductMultiVectors ) {

    // All of the multi-vector objects support the ProductMultiVectorBase
    // interface so we can do the reductions block by block.  Note, this is
    // not the most efficient implementation in an SPMD program but this is
    // easy to code up and use!

    // We must set up temporary arrays for the pointers to the MultiVectorBase
    // blocks for each block of objects!  What a mess!
    Array<RCP<const MultiVectorBase<Scalar> > >
      multi_vecs_rcp_block_k(multi_vecs_in.size());
    Array<Ptr<const MultiVectorBase<Scalar> > >
      multi_vecs_block_k(multi_vecs_in.size());
    Array<RCP<MultiVectorBase<Scalar> > >
      targ_multi_vecs_rcp_block_k(targ_multi_vecs_inout.size());
    Array<Ptr<MultiVectorBase<Scalar> > >
      targ_multi_vecs_block_k(targ_multi_vecs_inout.size());

    Ordinal g_off = primary_global_offset_in;

    for ( int k = 0; k < numBlocks_; ++k ) {

      const Ordinal dim_k = productSpace_->getBlock(k)->dim();

      // Fill the MultiVector array objects for this block

      for ( int j = 0; j < multi_vecs_in.size(); ++j ) {
        RCP<const MultiVectorBase<Scalar> > multi_vecs_rcp_block_k_j =
          multi_vecs[j]->getMultiVectorBlock(k);
        multi_vecs_rcp_block_k[j] = multi_vecs_rcp_block_k_j;
        multi_vecs_block_k[j] = multi_vecs_rcp_block_k_j.ptr();
      }

      for ( int j = 0; j < targ_multi_vecs_inout.size(); ++j ) {
        RCP<MultiVectorBase<Scalar> > targ_multi_vecs_rcp_block_k_j =
          targ_multi_vecs[j]->getNonconstMultiVectorBlock(k);
        targ_multi_vecs_rcp_block_k[j] = targ_multi_vecs_rcp_block_k_j;
        targ_multi_vecs_block_k[j] = targ_multi_vecs_rcp_block_k_j.ptr();
      }

      // Apply the RTOp object to the MultiVectors for this block

      Thyra::applyOp<Scalar>(
        primary_op, multi_vecs_block_k(), targ_multi_vecs_block_k(),
        reduct_objs,
        g_off
        );

      g_off += dim_k;
    }

  }
  else {

    // All of the multi-vector objects do not support the
    // ProductMultiVectorBase interface but if we got here (in debug mode)
    // then the spaces said that they are compatible so fall back on the
    // column-by-column implementation that will work correctly in serial.

    MultiVectorDefaultBase<Scalar>::mvMultiReductApplyOpImpl(
      primary_op, multi_vecs_in(), targ_multi_vecs_inout(),
      reduct_objs, primary_global_offset_in);

  }

}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::acquireDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
  ) const
{
  return MultiVectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(
    rowRng, colRng, sub_mv );
  // ToDo: Override this implementation if needed!
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::releaseDetachedMultiVectorViewImpl(
  RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
  ) const
{
  return MultiVectorDefaultBase<Scalar>::releaseDetachedMultiVectorViewImpl(
    sub_mv );
  // ToDo: Override this implementation if needed!
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::acquireNonconstDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::SubMultiVectorView<Scalar> *sub_mv
  )
{
  return MultiVectorDefaultBase<Scalar>::acquireNonconstDetachedMultiVectorViewImpl(
    rowRng,colRng,sub_mv );
  // ToDo: Override this implementation if needed!
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::commitNonconstDetachedMultiVectorViewImpl(
  RTOpPack::SubMultiVectorView<Scalar>* sub_mv
  )
{
  return MultiVectorDefaultBase<Scalar>::commitNonconstDetachedMultiVectorViewImpl(sub_mv);
  // ToDo: Override this implementation if needed!
}


// Overridden protected functions from LinearOpBase


template<class Scalar>
bool DefaultProductMultiVector<Scalar>::opSupportedImpl(EOpTransp /* M_trans */) const
{
  return true; // We can do it all!
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X_in,
  const Ptr<MultiVectorBase<Scalar> > &Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::dyn_cast;
  using Thyra::apply;

#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultProductMultiVector<Scalar>::apply(...)",
    *this, M_trans, X_in, &*Y_inout );
#endif

  if ( real_trans(M_trans) == NOTRANS ) {
    //
    // Y = b*Y + a*M*X
    //
    //    =>
    //
    // Y[k] = b*Y[k] + a*M[k]*X, k = 0...numBlocks-1
    //
    ProductMultiVectorBase<Scalar>
      &Y = dyn_cast<ProductMultiVectorBase<Scalar> >(*Y_inout);
    for ( int k = 0; k < numBlocks_; ++k ) {
      Thyra::apply(
        *multiVecs_[k].getConstObj(), M_trans,
        X_in, Y.getNonconstMultiVectorBlock(k).ptr(),
        alpha, beta ); 
    }
  }
  else {
    //
    // Y = b*Y + a*trans(M)*X
    //
    //    =>
    //
    // Y = b*Y + sum( a * trans(M[k]) * X[k], k=0...numBlocks-1 )
    //
    const ProductMultiVectorBase<Scalar>
      &X = dyn_cast<const ProductMultiVectorBase<Scalar> >(X_in);
    for ( int k = 0; k < numBlocks_; ++k ) {
      RCP<const MultiVectorBase<Scalar> >
        M_k = multiVecs_[k].getConstObj(),
        X_k = X.getMultiVectorBlock(k);
      if ( 0 == k ) {
        Thyra::apply( *M_k, M_trans, *X_k, Y_inout.ptr(), alpha, beta ); 
      }
      else {
        Thyra::apply( *M_k, M_trans, *X_k, Y_inout.ptr(), alpha, ST::one() ); 
      }
    }
  }
}


// private


template<class Scalar>
template<class MultiVectorType>
void DefaultProductMultiVector<Scalar>::initializeImpl(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
  const ArrayView<const RCP<MultiVectorType> > &multiVecs
  )
{
  // This function provides the "strong" guarantee (i.e. if an exception is
  // thrown, then *this will be left in the original state as before the
  // function was called)!
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT(nonnull(productSpace_in));
  TEUCHOS_ASSERT_EQUALITY(multiVecs.size(), productSpace_in->numBlocks());
#endif // TEUCHOS_DEBUG
  const RCP<const VectorSpaceBase<Scalar> >
    theDomain = multiVecs[0]->domain();
  const int numBlocks = productSpace_in->numBlocks();
#ifdef TEUCHOS_DEBUG
  for ( int k = 0; k < numBlocks; ++k ) {
    THYRA_ASSERT_VEC_SPACES(
      Teuchos::TypeNameTraits<DefaultProductMultiVector<Scalar> >::name(),
      *theDomain, *multiVecs[k]->domain()
      );
  }
#endif
  productSpace_ = productSpace_in;
  numBlocks_ = numBlocks;
  multiVecs_.assign(multiVecs.begin(),multiVecs.end());
}


#ifdef TEUCHOS_DEBUG


template<class Scalar>
void DefaultProductMultiVector<Scalar>::assertInitialized() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    is_null(productSpace_), std::logic_error,
    "Error, this DefaultProductMultiVector object is not intialized!"
    );
}


template<class Scalar>
void DefaultProductMultiVector<Scalar>::validateColIndex(const int j) const
{
  assertInitialized();
  const int domainDim = multiVecs_[0].getConstObj()->domain()->dim();
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! ( 0 <= j && j < domainDim ), std::logic_error,
    "Error, the column index j = " << j << " does not fall in the range [0,"<<domainDim<<"]!"
    );
}


#endif // TEUCHOS_DEBUG


} // namespace Thyra


template<class Scalar>
Teuchos::RCP<Thyra::DefaultProductMultiVector<Scalar> >
Thyra::defaultProductMultiVector()
{
  return Teuchos::rcp(new DefaultProductMultiVector<Scalar>);
}


template<class Scalar>
Teuchos::RCP<Thyra::DefaultProductMultiVector<Scalar> >
Thyra::defaultProductMultiVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const int numMembers
  )
{
  RCP<DefaultProductMultiVector<Scalar> > pmv = defaultProductMultiVector<Scalar>();
  pmv->initialize(productSpace, numMembers);
  return pmv;
}


template<class Scalar>
Teuchos::RCP<Thyra::DefaultProductMultiVector<Scalar> >
Thyra::defaultProductMultiVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const ArrayView<const RCP<MultiVectorBase<Scalar> > > &multiVecs
  )
{
  const RCP<DefaultProductMultiVector<Scalar> > pmv =
    defaultProductMultiVector<Scalar>();
  pmv->initialize(productSpace, multiVecs);
  return pmv;
}


template<class Scalar>
Teuchos::RCP<Thyra::DefaultProductMultiVector<Scalar> >
Thyra::defaultProductMultiVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const ArrayView<const RCP<const MultiVectorBase<Scalar> > > &multiVecs
  )
{
  const RCP<DefaultProductMultiVector<Scalar> > pmv =
    defaultProductMultiVector<Scalar>();
  pmv->initialize(productSpace, multiVecs);
  return pmv;
}


template<class Scalar>
Teuchos::RCP<const Thyra::ProductMultiVectorBase<Scalar> >
Thyra::castOrCreateSingleBlockProductMultiVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const RCP<const MultiVectorBase<Scalar> > &mv
  )
{
  const RCP<const ProductMultiVectorBase<Scalar> > pmv =
    Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(mv);
  if (nonnull(pmv))
    return pmv;
  return defaultProductMultiVector<Scalar>(productSpace, Teuchos::tuple(mv)());
}


template<class Scalar>
Teuchos::RCP<Thyra::ProductMultiVectorBase<Scalar> >
Thyra::nonconstCastOrCreateSingleBlockProductMultiVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const RCP<MultiVectorBase<Scalar> > &mv
  )
{
  const RCP<ProductMultiVectorBase<Scalar> > pmv =
    Teuchos::rcp_dynamic_cast<ProductMultiVectorBase<Scalar> >(mv);
  if (nonnull(pmv))
    return pmv;
  return defaultProductMultiVector<Scalar>(productSpace, Teuchos::tuple(mv)());
}


//
// Explicit instantiation macro
//
// Must be expanded from within the Thyra namespace!
//


#define THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_INSTANT(SCALAR) \
  \
  template class DefaultProductMultiVector<SCALAR >; \
  \
  template RCP<DefaultProductMultiVector<SCALAR > >  \
  defaultProductMultiVector<SCALAR >();  \
    \
    \
  template RCP<DefaultProductMultiVector<SCALAR > >  \
  defaultProductMultiVector(  \
    const RCP<const DefaultProductVectorSpace<SCALAR > > &productSpace,  \
    const int numMembers  \
    );  \
    \
    \
  template RCP<DefaultProductMultiVector<SCALAR > >  \
  defaultProductMultiVector(  \
    const RCP<const DefaultProductVectorSpace<SCALAR > > &productSpace,  \
    const ArrayView<const RCP<MultiVectorBase<SCALAR > > > &multiVecs  \
    );  \
    \
    \
  template RCP<DefaultProductMultiVector<SCALAR > >  \
  defaultProductMultiVector(  \
    const RCP<const DefaultProductVectorSpace<SCALAR > > &productSpace,  \
    const ArrayView<const RCP<const MultiVectorBase<SCALAR > > > &multiVecs  \
    );  \
    \
    \
  template RCP<const ProductMultiVectorBase<SCALAR > >  \
  castOrCreateSingleBlockProductMultiVector(  \
    const RCP<const DefaultProductVectorSpace<SCALAR > > &productSpace, \
    const RCP<const MultiVectorBase<SCALAR > > &mv  \
    );  \
    \
    \
  template RCP<ProductMultiVectorBase<SCALAR > >  \
  nonconstCastOrCreateSingleBlockProductMultiVector(  \
    const RCP<const DefaultProductVectorSpace<SCALAR > > &productSpace, \
    const RCP<MultiVectorBase<SCALAR > > &mv  \
    );


#endif // THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DEF_HPP
