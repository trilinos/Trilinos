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

#ifndef THYRA_DEFAULT_PRODUCT_VECTOR_DEF_HPP
#define THYRA_DEFAULT_PRODUCT_VECTOR_DEF_HPP


#include "Thyra_DefaultProductVector_decl.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Teuchos_Workspace.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultProductVector<Scalar>::DefaultProductVector()
  : numBlocks_(0)
{
  uninitialize();
}


template <class Scalar>
DefaultProductVector<Scalar>::DefaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace_in
  )
  : numBlocks_(0)
{
  initialize(productSpace_in);
}


template <class Scalar>
void DefaultProductVector<Scalar>::initialize(
  const RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace_in
  )
{
  // ToDo: Validate input!
  numBlocks_ = productSpace_in->numBlocks();
  productSpace_ = productSpace_in;
  vecs_.resize(numBlocks_);
  for( int k = 0; k < numBlocks_; ++k )
    vecs_[k].initialize(createMember(productSpace_in->getBlock(k)));
}


template <class Scalar>
void DefaultProductVector<Scalar>::initialize(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
  const ArrayView<const RCP<VectorBase<Scalar> > > &vecs
  )
{
  using Teuchos::as;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY( as<Ordinal>(productSpace_in->numBlocks()),
    as<Ordinal>(vecs.size()) );
#endif
  numBlocks_ = productSpace_in->numBlocks();
  productSpace_ = productSpace_in;
  vecs_.resize(numBlocks_);
  for( int k = 0; k < numBlocks_; ++k )
    vecs_[k].initialize(vecs[k]);
}


template <class Scalar>
void DefaultProductVector<Scalar>::initialize(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
  const ArrayView<const RCP<const VectorBase<Scalar> > > &vecs
  )
{
  using Teuchos::as;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY( as<Ordinal>(productSpace_in->numBlocks()),
    as<Ordinal>(vecs.size()) );
#endif
  numBlocks_ = productSpace_in->numBlocks();
  productSpace_ = productSpace_in;
  vecs_.resize(numBlocks_);
  for( int k = 0; k < numBlocks_; ++k )
    vecs_[k].initialize(vecs[k]);
}


template <class Scalar>
void DefaultProductVector<Scalar>::uninitialize()
{
  productSpace_ = Teuchos::null;
  vecs_.resize(0);
  numBlocks_ = 0;
}


// Overridden from Teuchos::Describable

                                                
template<class Scalar>
std::string DefaultProductVector<Scalar>::description() const
{
  const RCP<const VectorSpaceBase<Scalar> > vs = this->space();
  std::ostringstream oss;
  oss
    << Teuchos::Describable::description()
    << "{"
    << "dim="<<(nonnull(vs) ? vs->dim() : 0)
    << ", numBlocks = "<<numBlocks_
    << "}";
  return oss.str();
}


template<class Scalar>
void DefaultProductVector<Scalar>::describe(
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_NONE:
      break;
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
      *out
        <<  "numBlocks="<< numBlocks_ << std::endl
        <<  "Constituent vector objects v[0], v[1], ... v[numBlocks-1]:\n";
      OSTab tab3(out);
      for( int k = 0; k < numBlocks_; ++k ) {
        *out << "v["<<k<<"] = " << describe(*vecs_[k].getConstObj(),verbLevel);
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
}


// Extensions to ProductVectorBase suitable for physically-blocked vectors


template <class Scalar>
void DefaultProductVector<Scalar>::setBlock(
  int i, const RCP<const VectorBase<Scalar> >& b
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(i < 0 || i >= numBlocks_);
  TEUCHOS_TEST_FOR_EXCEPT(!productSpace_->getBlock(i)->isCompatible(*(b->space())));
#endif
  vecs_[i] = b;
}


template <class Scalar>
void DefaultProductVector<Scalar>::setNonconstBlock(
  int i, const RCP<VectorBase<Scalar> >& b
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(i < 0 || i >= numBlocks_);
  TEUCHOS_TEST_FOR_EXCEPT(!productSpace_->getBlock(i)->isCompatible(*(b->space())));
#endif
  vecs_[i] = b;
}


// Overridden from ProductVectorBase


template <class Scalar>
RCP<VectorBase<Scalar> >
DefaultProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
#endif
  return vecs_[k].getNonconstObj();
}


template <class Scalar>
RCP<const VectorBase<Scalar> >
DefaultProductVector<Scalar>::getVectorBlock(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
#endif
  return vecs_[k].getConstObj();
}


// Overridden from ProductMultiVectorBase


template <class Scalar>
RCP<const ProductVectorSpaceBase<Scalar> >
DefaultProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template <class Scalar>
bool DefaultProductVector<Scalar>::blockIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( k < 0 || numBlocks_-1 < k);
#endif
  return vecs_[k].isConst();
}


template <class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}


template <class Scalar>
RCP<const MultiVectorBase<Scalar> >
DefaultProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}


// Overridden from VectorBase


template <class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultProductVector<Scalar>::space() const
{
  return productSpace_;
}


template <class Scalar>
void DefaultProductVector<Scalar>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset_in
  ) const
{

  // 2008/02/20: rabartl: ToDo: Upgrade Teuchos::Workspace<T> to implicitly
  // convert to Teuchos::ArrayView<T>.  This will allow the calls to
  // applyOp(...) with sub_vecs and sub_targ_vecs to work without trouble!
  // For now, I just want to get this done.  It is likely that this function
  // is going to change in major ways soon anyway!

  //using Teuchos::Workspace;
  using Teuchos::ptr_dynamic_cast;
  using Teuchos::describe;
  using Teuchos::null;

  //Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const Ordinal	n = productSpace_->dim();
  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();

  // Validate the compatibility of the vectors!
#ifdef TEUCHOS_DEBUG
  bool test_failed;
  for(int k = 0; k < num_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*vecs[k]->space());
    TEUCHOS_TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultProductVector::applyOp(...): Error vecs["<<k<<"]->space() = "
      <<vecs[k]->space()->description()<<"\' is not compatible with this "
      <<"vector space = "<<this->space()->description()<<"!"
      );
  }
  for(int k = 0; k < num_targ_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*targ_vecs[k]->space());
    TEUCHOS_TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultProductVector::applyOp(...): Error targ_vecs["<<k<<"]->space() = "
      <<targ_vecs[k]->space()->description()<<"\' is not compatible with this "
      <<"vector space = "<<this->space()->description()<<"!"
      );
  }
#endif

  //
  // Dynamic cast each of the vector arguments to the ProductVectorBase interface
  //
  // NOTE: If the constituent vector is not a product vector, then a product
  // vector of one component is created.
  //

  Array<RCP<const ProductVectorBase<Scalar> > > vecs_args_store(num_vecs);
  Array<Ptr<const ProductVectorBase<Scalar> > > vecs_args(num_vecs);
  for(int k = 0; k < num_vecs; ++k) {
    vecs_args_store[k] =
      castOrCreateProductVectorBase<Scalar>(rcpFromPtr(vecs[k]));
    vecs_args[k] = vecs_args_store[k].ptr();
  }

  Array<RCP<ProductVectorBase<Scalar> > > targ_vecs_args_store(num_targ_vecs);
  Array<Ptr<ProductVectorBase<Scalar> > > targ_vecs_args(num_targ_vecs);
  for(int k = 0; k < num_targ_vecs; ++k) {
    targ_vecs_args_store[k] =
      castOrCreateNonconstProductVectorBase<Scalar>(rcpFromPtr(targ_vecs[k]));
    targ_vecs_args[k] = targ_vecs_args_store[k].ptr();
  }

  //
  // If we get here, then we will implement the applyOpImpl(...) one vector
  // block at a time.
  //
  const Ordinal dim = n;
  Ordinal num_elements_remaining = dim;
  const int numBlocks = productSpace_->numBlocks();
  Array<RCP<const VectorBase<Scalar> > >
    sub_vecs_rcps(num_vecs);
  Array<Ptr<const VectorBase<Scalar> > >
    sub_vecs(num_vecs);
  Array<RCP<VectorBase<Scalar> > >
    sub_targ_vecs_rcps(num_targ_vecs);
  Array<Ptr<VectorBase<Scalar> > >
    sub_targ_vecs(num_targ_vecs);
  Ordinal g_off = 0;
  for(int k = 0; k < numBlocks; ++k) {
    const Ordinal dim_k = productSpace_->getBlock(k)->dim();
    // Fill constituent vectors for block k
    for( int i = 0; i < num_vecs; ++i ) {
      sub_vecs_rcps[i] = vecs_args[i]->getVectorBlock(k);
      sub_vecs[i] = sub_vecs_rcps[i].ptr();
    }
    // Fill constituent target vectors for block k
    for( int j = 0; j < num_targ_vecs; ++j ) {
      sub_targ_vecs_rcps[j] = targ_vecs_args[j]->getNonconstVectorBlock(k);
      sub_targ_vecs[j] = sub_targ_vecs_rcps[j].ptr();
    }
    Thyra::applyOp<Scalar>(
      op, sub_vecs(), sub_targ_vecs(),
      reduct_obj,
      global_offset_in + g_off
      );
    g_off += dim_k;
    num_elements_remaining -= dim_k;
  }
  TEUCHOS_TEST_FOR_EXCEPT(!(num_elements_remaining==0));

}


// protected


// Overridden protected functions from VectorBase


template <class Scalar>
void DefaultProductVector<Scalar>::acquireDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  const Range1D
    rng = rng_in.full_range() ? Range1D(0,productSpace_->dim()-1) : rng_in;
  int    kth_vector_space  = -1;
  Ordinal  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    rng.lbound() + rng.size()
    <= kth_global_offset + vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This involves only one sub-vector so just return it.
    const_cast<const VectorBase<Scalar>*>(
      &*vecs_[kth_vector_space].getConstObj()
      )->acquireDetachedView( rng - kth_global_offset, sub_vec );
    sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
  }
  else {
    // Just let the default implementation handle this.  ToDo: In the future
    // we could manually construct an explicit sub-vector that spanned
    // two or more constituent vectors but this would be a lot of work.
    // However, this would require the use of temporary memory but
    // so what.
    VectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(rng_in,sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
  if( sub_vec->values().get() == NULL ) return;
  int    kth_vector_space  = -1;
  Ordinal  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(sub_vec->globalOffset(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    sub_vec->globalOffset() + sub_vec->subDim()
    <= kth_global_offset +  vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This sub_vec was extracted from a single constituent vector
    sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
    vecs_[kth_vector_space].getConstObj()->releaseDetachedView(sub_vec);
  }
  else {
    // This sub_vec was created by the default implementation!
    VectorDefaultBase<Scalar>::releaseDetachedVectorViewImpl(sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  const Range1D
    rng = rng_in.full_range() ? Range1D(0,productSpace_->dim()-1) : rng_in;
  int    kth_vector_space  = -1;
  Ordinal  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    rng.lbound() + rng.size()
    <= kth_global_offset + vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This involves only one sub-vector so just return it.
    vecs_[kth_vector_space].getConstObj()->acquireDetachedView(
      rng - kth_global_offset, sub_vec
      );
    sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
  }
  else {
    // Just let the default implementation handle this.  ToDo: In the future
    // we could manually construct an explicit sub-vector that spanned
    // two or more constituent vectors but this would be a lot of work.
    // However, this would require the use of temporary memory but
    // so what.
    VectorDefaultBase<Scalar>::acquireNonconstDetachedVectorViewImpl(rng_in,sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
  if( sub_vec->values().get() == NULL ) return;
  int    kth_vector_space  = -1;
  Ordinal  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(sub_vec->globalOffset(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    sub_vec->globalOffset() + sub_vec->subDim()
    <= kth_global_offset +  vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This sub_vec was extracted from a single constituent vector
    sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
    vecs_[kth_vector_space].getNonconstObj()->commitDetachedView(sub_vec);
  }
  else {
    // This sub_vec was created by the default implementation!
    VectorDefaultBase<Scalar>::commitNonconstDetachedVectorViewImpl(sub_vec);
  }
}


template <class Scalar>
void DefaultProductVector<Scalar>::setSubVectorImpl(
  const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
  )
{
  int    kth_vector_space  = -1;
  Ordinal  kth_global_offset = 0;
  productSpace_->getVecSpcPoss(sub_vec.globalOffset(),&kth_vector_space,&kth_global_offset);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( 0 <= kth_vector_space && kth_vector_space <= numBlocks_ ) );
#endif
  if(
    sub_vec.globalOffset() + sub_vec.subDim()
    <= kth_global_offset + vecs_[kth_vector_space].getConstObj()->space()->dim()
    )
  {
    // This sub-vector fits into a single constituent vector
    RTOpPack::SparseSubVectorT<Scalar> sub_vec_g = sub_vec;
    sub_vec_g.setGlobalOffset( sub_vec_g.globalOffset() - kth_global_offset );
    vecs_[kth_vector_space].getNonconstObj()->setSubVector(sub_vec_g);
  }
  else {
    // Let the default implementation take care of this.  ToDo: In the future
    // it would be possible to manually set the relevant constituent
    // vectors with no temp memory allocations.
    VectorDefaultBase<Scalar>::setSubVector(sub_vec);
  }
}


} // namespace Thyra


template<class Scalar>
Teuchos::RCP<Thyra::ProductVectorBase<Scalar> >
Thyra::castOrCreateNonconstProductVectorBase(const RCP<VectorBase<Scalar> > v)
{
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::tuple;
  const RCP<ProductVectorBase<Scalar> > prod_v =
    rcp_dynamic_cast<ProductVectorBase<Scalar> >(v);
  if (nonnull(prod_v)) {
    return prod_v;
  }
  return defaultProductVector<Scalar>(
    productVectorSpace<Scalar>(tuple(v->space())()),
    tuple(v)()
    );
}


template<class Scalar>
Teuchos::RCP<const Thyra::ProductVectorBase<Scalar> >
Thyra::castOrCreateProductVectorBase(const RCP<const VectorBase<Scalar> > v)
{
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::tuple;
  const RCP<const ProductVectorBase<Scalar> > prod_v =
    rcp_dynamic_cast<const ProductVectorBase<Scalar> >(v);
  if (nonnull(prod_v)) {
    return prod_v;
  }
  return defaultProductVector<Scalar>(
    productVectorSpace<Scalar>(tuple(v->space())()),
    tuple(v)()
    );
}


//
// Explicit instant macro
//

#define THYRA_DEFAULT_PRODUCT_VECTOR_INSTANT(SCALAR) \
  \
  template class DefaultProductVector<SCALAR >; \
  \
  template RCP<ProductVectorBase<SCALAR > >  \
  castOrCreateNonconstProductVectorBase(const RCP<VectorBase<SCALAR > > v);  \
  \
  template RCP<const ProductVectorBase<SCALAR > >  \
  castOrCreateProductVectorBase(const RCP<const VectorBase<SCALAR > > v);  \



#endif // THYRA_DEFAULT_PRODUCT_VECTOR_DEF_HPP
