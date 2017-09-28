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

#ifndef THYRA_MULTI_VECTOR_DEFAULT_BASE_DEF_HPP
#define THYRA_MULTI_VECTOR_DEFAULT_BASE_DEF_HPP


#include "Thyra_MultiVectorDefaultBase_decl.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_DefaultColumnwiseMultiVector.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_ROpNorm1.hpp"
#include "RTOpPack_ROpNorm2.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "RTOpPack_TOpAssignVectors.hpp"
#include "RTOpPack_TOpAXPY.hpp"
#include "RTOpPack_TOpLinearCombination.hpp"
#include "RTOpPack_TOpScaleVector.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace Thyra {


// Overridden public member functions from MultiVectorBase


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::clone_mv() const
{
  const VectorSpaceBase<Scalar>
    &l_domain = *this->domain(),
    &l_range = *this->range();
  RCP<MultiVectorBase<Scalar> >
    copy = createMembers(l_range,l_domain.dim());
  ::Thyra::assign<Scalar>(copy.ptr(), *this);
  return copy;
}


// protected


// Overridden protected member functions from MultiVectorBase

template<class Scalar>
void MultiVectorDefaultBase<Scalar>::assignImpl(Scalar alpha)
{
  using Teuchos::tuple; using Teuchos::null;
  RTOpPack::TOpAssignScalar<Scalar> assign_scalar_op(alpha);
  Thyra::applyOp<Scalar>(assign_scalar_op,
    ArrayView<Ptr<const MultiVectorBase<Scalar> > >(null),
    tuple<Ptr<MultiVectorBase<Scalar> > >(ptr(this)), null);
}

template<class Scalar>
void MultiVectorDefaultBase<Scalar>::assignMultiVecImpl(const MultiVectorBase<Scalar>& mv)
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpAssignVectors<Scalar> assign_vectors_op;
  Thyra::applyOp<Scalar>(assign_vectors_op, tuple(ptrInArg(mv)),
    tuple<Ptr<MultiVectorBase<Scalar> > >(ptr(this)), null);
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::scaleImpl(Scalar alpha)
{
  using Teuchos::tuple; using Teuchos::null;
  typedef ScalarTraits<Scalar> ST;
  if (alpha==ST::zero()) {
    this->assign(ST::zero());
    //assign(tuple<Ptr<MultiVectorBase<Scalar> > >(ptr(this)), ST::zero());
    return;
  }
  if (alpha==ST::one()) {
    return;
  }
  RTOpPack::TOpScaleVector<Scalar> scale_vector_op(alpha);
  Thyra::applyOp<Scalar>(scale_vector_op,
    ArrayView<Ptr<const MultiVectorBase<Scalar> > >(null),
    tuple<Ptr<MultiVectorBase<Scalar> > >(ptr(this)), null );
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::updateImpl(
  Scalar alpha,
  const MultiVectorBase<Scalar>& mv
  )
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  RTOpPack::TOpAXPY<Scalar> axpy_op(alpha);
  Thyra::applyOp<Scalar>( axpy_op, tuple(ptrInArg(mv)),
    tuple<Ptr<MultiVectorBase<Scalar> > >(ptr(this)), null );
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::linearCombinationImpl(
  const ArrayView<const Scalar>& alpha,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
  const Scalar& beta
  )
{
  using Teuchos::tuple; using Teuchos::null;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), mv.size());
#endif
  const int m = alpha.size();
  if ( beta == ScalarTraits<Scalar>::one() && m == 1 ) {
    this->update(alpha[0], *mv[0]);
    return;
  }
  else if (m == 0) {
    this->scale(beta);
    return;
  }
  RTOpPack::TOpLinearCombination<Scalar> lin_comb_op(alpha, beta);
  Thyra::applyOp<Scalar>(lin_comb_op, mv,
    tuple<Ptr<MultiVectorBase<Scalar> > >(ptr(this)), null );
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::norms1Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  const int m = this->domain()->dim();
  RTOpPack::ROpNorm1<Scalar> norm_op;
  Array<RCP<RTOpPack::ReductTarget> > rcp_op_targs(m);
  Array<Ptr<RTOpPack::ReductTarget> > op_targs(m);
  for(int kc = 0; kc < m; ++kc) {
    rcp_op_targs[kc] = norm_op.reduct_obj_create();
    op_targs[kc] = rcp_op_targs[kc].ptr();
  }
  Thyra::applyOp<Scalar>(norm_op,
    tuple<Ptr<const MultiVectorBase<Scalar> > >(ptrInArg(*this)),
    ArrayView<Ptr<MultiVectorBase<Scalar> > >(null),
    op_targs);
  for(int kc = 0; kc < m; ++kc) {
    norms[kc] = norm_op(*op_targs[kc]);
  }
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::norms2Impl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  const int m = this->domain()->dim();
  RTOpPack::ROpNorm2<Scalar> norm_op;
  Array<RCP<RTOpPack::ReductTarget> > rcp_op_targs(m);
  Array<Ptr<RTOpPack::ReductTarget> > op_targs(m);
  for(int kc = 0; kc < m; ++kc) {
    rcp_op_targs[kc] = norm_op.reduct_obj_create();
    op_targs[kc] = rcp_op_targs[kc].ptr();
  }
  Thyra::applyOp<Scalar>(norm_op,
    tuple<Ptr<const MultiVectorBase<Scalar> > >(ptrInArg(*this)),
    ArrayView<Ptr<MultiVectorBase<Scalar> > >(null),
    op_targs);
  for(int kc = 0; kc < m; ++kc) {
    norms[kc] = norm_op(*op_targs[kc]);
  }
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::dotsImpl(
  const MultiVectorBase<Scalar>& mv,
  const ArrayView<Scalar>& prods
  ) const
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  const int m = this->domain()->dim();
  RTOpPack::ROpDotProd<Scalar> dot_op;
  Array<RCP<RTOpPack::ReductTarget> > rcp_dot_targs(m);
  Array<Ptr<RTOpPack::ReductTarget> > dot_targs(m);
  for( int kc = 0; kc < m; ++kc ) {
    rcp_dot_targs[kc] = dot_op.reduct_obj_create();
    dot_targs[kc] = rcp_dot_targs[kc].ptr();
  }
  Thyra::applyOp<Scalar>( dot_op,
    tuple<Ptr<const MultiVectorBase<Scalar> > >(ptrInArg(mv), ptrInArg(*this)),
    ArrayView<Ptr<MultiVectorBase<Scalar> > >(null),
    dot_targs );
  for( int kc = 0; kc < m; ++kc ) {
    prods[kc] = dot_op(*dot_targs[kc]);
  }
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::normsInfImpl(
  const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
  ) const
{
  using Teuchos::tuple; using Teuchos::ptrInArg; using Teuchos::null;
  const int m = this->domain()->dim();
  RTOpPack::ROpNormInf<Scalar> norm_op;
  Array<RCP<RTOpPack::ReductTarget> > rcp_op_targs(m);
  Array<Ptr<RTOpPack::ReductTarget> > op_targs(m);
  for(int kc = 0; kc < m; ++kc) {
    rcp_op_targs[kc] = norm_op.reduct_obj_create();
    op_targs[kc] = rcp_op_targs[kc].ptr();
  }
  Thyra::applyOp<Scalar>(norm_op,
    tuple<Ptr<const MultiVectorBase<Scalar> > >(ptrInArg(*this)),
    ArrayView<Ptr<MultiVectorBase<Scalar> > >(null),
    op_targs);
  for(int kc = 0; kc < m; ++kc) {
    norms[kc] = norm_op(*op_targs[kc]);
  }
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::contigSubViewImpl( const Range1D& colRng_in ) const
{
  using Teuchos::Workspace;
  using Teuchos::as;
  Teuchos::WorkspaceStore *wss = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar> &l_domain = *this->domain();
  const VectorSpaceBase<Scalar> &l_range = *this->range();
  const Ordinal dimDomain = l_domain.dim();
  const Range1D colRng = Teuchos::full_range(colRng_in,0,dimDomain-1);
  if( colRng.lbound() == 0 && as<Ordinal>(colRng.ubound()) == dimDomain-1 )
    return Teuchos::rcp(this,false); // Takes all of the columns!
  if( colRng.size() ) {
    // We have to create a view of a subset of the columns
    Workspace< RCP< VectorBase<Scalar> > > col_vecs(wss,colRng.size());
    for( Ordinal j = colRng.lbound(); j <= colRng.ubound(); ++j )
      col_vecs[j-colRng.lbound()] = Teuchos::rcp_const_cast<VectorBase<Scalar> >(this->col(j));
    return Teuchos::rcp(
      new DefaultColumnwiseMultiVector<Scalar>(
        this->range(),l_range.smallVecSpcFcty()->createVecSpc(colRng.size()),col_vecs
        )
      );
  }
  return Teuchos::null; // There was an empty set in colRng_in!
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::nonconstContigSubViewImpl( const Range1D& colRng_in )
{
  using Teuchos::Workspace;
  using Teuchos::as;
  Teuchos::WorkspaceStore *wss = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar> &l_domain = *this->domain();
  const VectorSpaceBase<Scalar> &l_range = *this->range();
  const Ordinal dimDomain = l_domain.dim();
  const Range1D colRng = Teuchos::full_range(colRng_in,0,dimDomain-1);
  if( colRng.lbound() == 0 && as<Ordinal>(colRng.ubound()) == dimDomain-1 )
    return Teuchos::rcp(this,false); // Takes all of the columns!
  if( colRng.size() ) {
    // We have to create a view of a subset of the columns
    Workspace< RCP< VectorBase<Scalar> > > col_vecs(wss,colRng.size());
    for( Ordinal j = colRng.lbound(); j <= colRng.ubound(); ++j )
      col_vecs[j-colRng.lbound()] = this->col(j);
    return Teuchos::rcp(
      new DefaultColumnwiseMultiVector<Scalar>(
        this->range(),l_range.smallVecSpcFcty()->createVecSpc(colRng.size()),col_vecs
        )
      );
  }
  return Teuchos::null; // There was an empty set in colRng_in!
}


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::nonContigSubViewImpl(
  const ArrayView<const int> &cols
  ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore *wss = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar> &l_range = *this->range();
  const int numCols = cols.size();
#ifdef TEUCHOS_DEBUG
  const VectorSpaceBase<Scalar> &l_domain = *this->domain();
  const Ordinal dimDomain = l_domain.dim();
  const char msg_err[] = "MultiVectorDefaultBase<Scalar>::subView(numCols,cols[]): Error!";
  TEUCHOS_TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
  // We have to create a view of a subset of the columns
  Workspace< RCP< VectorBase<Scalar> > > col_vecs(wss,numCols);
  for( int k = 0; k < numCols; ++k ) {
    const int col_k = cols[k];
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      !( 0 <= col_k && col_k < dimDomain ), std::invalid_argument
      ,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [0,"<<(dimDomain-1)<<"]!"
      );
#endif
    col_vecs[k] = Teuchos::rcp_const_cast<VectorBase<Scalar> >(this->col(col_k));
  }
  return Teuchos::rcp(
    new DefaultColumnwiseMultiVector<Scalar>(
      this->range(), l_range.smallVecSpcFcty()->createVecSpc(numCols), col_vecs
      )
    );
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
MultiVectorDefaultBase<Scalar>::nonconstNonContigSubViewImpl(
  const ArrayView<const int> &cols
  )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore *wss = Teuchos::get_default_workspace_store().get();
  const VectorSpaceBase<Scalar> &l_range = *this->range();
  const int numCols = cols.size();
#ifdef TEUCHOS_DEBUG
  const VectorSpaceBase<Scalar> &l_domain = *this->domain();
  const Ordinal dimDomain = l_domain.dim();
  const char msg_err[] = "MultiVectorDefaultBase<Scalar>::subView(numCols,cols[]): Error!";
  TEUCHOS_TEST_FOR_EXCEPTION( numCols < 1 || dimDomain < numCols, std::invalid_argument, msg_err );
#endif
  // We have to create a view of a subset of the columns
  Workspace< RCP< VectorBase<Scalar> > > col_vecs(wss,numCols);
  for( int k = 0; k < numCols; ++k ) {
    const int col_k = cols[k];
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      !( 0 <= col_k && col_k < dimDomain ), std::invalid_argument
      ,msg_err << " col["<<k<<"] = " << col_k << " is not in the range [0,"<<(dimDomain-1)<<"]!"
      );
#endif
    col_vecs[k] = this->col(col_k);
  }
  return Teuchos::rcp(
    new DefaultColumnwiseMultiVector<Scalar>(
      this->range(), l_range.smallVecSpcFcty()->createVecSpc(numCols), col_vecs
      )
    );
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::mvMultiReductApplyOpImpl(
  const RTOpPack::RTOpT<Scalar> &prim_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal prim_global_offset_in
  ) const
{

  using Teuchos::Workspace;
  using Teuchos::as;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const int num_multi_vecs = multi_vecs.size();
  const int num_targ_multi_vecs = targ_multi_vecs.size();

  // ToDo: Validate the input!

  const VectorSpaceBase<Scalar> &l_domain = *this->domain();

  // Get the primary and secondary dimensions.

  const Ordinal sec_dim = l_domain.dim();

  //
  // Apply the reduction/transformation operator and transform the
  // target vectors and reduce each of the reduction objects.
  //

  Workspace<RCP<const VectorBase<Scalar> > > vecs_s(wss, num_multi_vecs);
  Workspace<Ptr<const VectorBase<Scalar> > > vecs(wss, num_multi_vecs);
  Workspace<RCP<VectorBase<Scalar> > > targ_vecs_s(wss, num_targ_multi_vecs);
  Workspace<Ptr<VectorBase<Scalar> > > targ_vecs(wss, num_targ_multi_vecs);

  for(Ordinal j = 0; j < sec_dim; ++j) {
    // Fill the arrays of vector arguments
    {for(Ordinal k = 0; k < as<Ordinal>(num_multi_vecs); ++k) {
        vecs_s[k] = multi_vecs[k]->col(j);
        vecs[k] = vecs_s[k].ptr();
      }}
    {for(Ordinal k = 0; k < as<Ordinal>(num_targ_multi_vecs); ++k) {
        targ_vecs_s[k] = targ_multi_vecs[k]->col(j);
        targ_vecs[k] = targ_vecs_s[k].ptr();
      }}
    // Apply the reduction/transformation operator
    Thyra::applyOp(
      prim_op,
      vecs().getConst(),
      targ_vecs().getConst(),
      reduct_objs.size() ? reduct_objs[j] : Ptr<RTOpPack::ReductTarget>(),
      prim_global_offset_in);
  }
  // At this point all of the designated targ vectors in the target multi-vectors have
  // been transformed and all the reduction objects in reduct_obj[] have accumulated
  // the reductions.
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::mvSingleReductApplyOpImpl(
  const RTOpPack::RTOpT<Scalar> &prim_op,
  const RTOpPack::RTOpT<Scalar> &sec_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal prim_global_offset_in
  ) const
{

  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  // ToDo: Validate the input!

  const VectorSpaceBase<Scalar> &l_domain = *this->domain();

  // Get the primary and secondary dimensions.
  const Ordinal sec_dim = l_domain.dim();

  // Create a temporary buffer for the reduction objects of the primary reduction
  // so that we can call the companion version of this method.
  const int reduct_objs_size = (!is_null(reduct_obj) ? sec_dim : 0);
  Workspace<RCP<RTOpPack::ReductTarget> > rcp_reduct_objs(wss, reduct_objs_size);
  Workspace<Ptr<RTOpPack::ReductTarget> > reduct_objs(wss, reduct_objs_size);
  if (!is_null(reduct_obj)) {
    for(Ordinal k = 0; k < sec_dim; ++k) {
      rcp_reduct_objs[k] = prim_op.reduct_obj_create();
      reduct_objs[k] = rcp_reduct_objs[k].ptr();
    }
  }
 
  // Call the companion version that accepts an array of reduction objects
  this->applyOp(
    prim_op, multi_vecs, targ_multi_vecs, reduct_objs,
    prim_global_offset_in);

  // Reduce all the reduction objects using the secondary reduction operator
  // into one reduction object and free the intermediate reduction objects.
  if (!is_null(reduct_obj)) {
    for (Ordinal k = 0; k < sec_dim; ++k) {
      sec_op.reduce_reduct_objs( *reduct_objs[k], reduct_obj );
    }
  }
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(
  const Range1D &rowRng_in,
  const Range1D &colRng_in,
  RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
  ) const
{
  const Ordinal
    rangeDim = this->range()->dim(),
    domainDim = this->domain()->dim();
  const Range1D
    rowRng = rowRng_in.full_range() ? Range1D(0,rangeDim-1) : rowRng_in,
    colRng = colRng_in.full_range() ? Range1D(0,domainDim-1) : colRng_in;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(rowRng.ubound() < rangeDim), std::out_of_range
    ,"MultiVectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(...): Error, rowRng = ["
    <<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not in the range = [0,"
    <<(rangeDim-1)<<"]!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(colRng.ubound() < domainDim), std::out_of_range
    ,"MultiVectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(...): Error, colRng = ["
    <<colRng.lbound()<<","<<colRng.ubound()<<"] is not in the range = [0,"
    <<(domainDim-1)<<"]!"
    );
#endif
  // Allocate storage for the multi-vector (stored column-major)
  const ArrayRCP<Scalar> values = Teuchos::arcp<Scalar>(rowRng.size() * colRng.size());
  // Extract multi-vector values column by column
  RTOpPack::ConstSubVectorView<Scalar> sv; // uninitialized by default
  for( int k = colRng.lbound(); k <= colRng.ubound(); ++k ) {
    RCP<const VectorBase<Scalar> > col_k = this->col(k);
    col_k->acquireDetachedView( rowRng, &sv );
    for( int i = 0; i < rowRng.size(); ++i )
      values[ i + k*rowRng.size() ] = sv[i];
    col_k->releaseDetachedView( &sv );
  }
  // Initialize the multi-vector view object
  sub_mv->initialize(
    rowRng.lbound(), // globalOffset
    rowRng.size(), // subDim
    colRng.lbound(), // colOffset
    colRng.size(), // numSubCols
    values, // values
    rowRng.size() // leadingDim
    );
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::releaseDetachedMultiVectorViewImpl(
  RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
  ) const
{
  // Here we just need to free the view and that is it!
  sub_mv->uninitialize();
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::acquireNonconstDetachedMultiVectorViewImpl(
  const Range1D &rowRng,
  const Range1D &colRng,
  RTOpPack::SubMultiVectorView<Scalar> *sub_mv
  )
{
  using Teuchos::as;
  // Use the non-const implementation since it does exactly the
  // correct thing in this case also!
  MultiVectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(
    rowRng, colRng,
    as<RTOpPack::ConstSubMultiVectorView<Scalar>*>(sub_mv)
    // This cast will work as long as SubMultiVectorView
    // maintains no extra state over ConstSubMultiVectorView (which it
    // currently does not) but this is something that I should
    // technically check for some how.
    );
}


template<class Scalar>
void MultiVectorDefaultBase<Scalar>::commitNonconstDetachedMultiVectorViewImpl(
  RTOpPack::SubMultiVectorView<Scalar>* sub_mv
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    sub_mv==NULL, std::logic_error,
    "MultiVectorDefaultBase<Scalar>::commitNonconstDetachedMultiVectorViewImpl(...): Error!"
    );
#endif
  // Set back the multi-vector values column by column
  const Range1D rowRng(sub_mv->globalOffset(),sub_mv->globalOffset()+sub_mv->subDim()-1);
  RTOpPack::SubVectorView<Scalar> msv; // uninitialized by default
  for( int k = sub_mv->colOffset(); k < sub_mv->numSubCols(); ++k ) {
    RCP<VectorBase<Scalar> > col_k = this->col(k);
    col_k->acquireDetachedView( rowRng, &msv );
    for( int i = 0; i < rowRng.size(); ++i )
      msv[i] = sub_mv->values()[ i + k*rowRng.size() ];
    col_k->commitDetachedView( &msv );
  }
  // Zero out the view
  sub_mv->uninitialize();
}


} // end namespace Thyra


#endif // THYRA_MULTI_VECTOR_DEFAULT_BASE_DEF_HPP
