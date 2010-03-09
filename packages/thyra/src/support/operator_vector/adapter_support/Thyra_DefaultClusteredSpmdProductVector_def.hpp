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

#ifndef THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_HPP
#define THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_HPP

#include "Thyra_DefaultClusteredSpmdProductVector_decl.hpp"
#include "Thyra_DefaultClusteredSpmdProductVectorSpace.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "RTOp_parallel_helpers.h"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultClusteredSpmdProductVector<Scalar>::DefaultClusteredSpmdProductVector()
{
  uninitialize();
}


template <class Scalar>
DefaultClusteredSpmdProductVector<Scalar>::DefaultClusteredSpmdProductVector(
  const Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  &productSpace_in
  ,const Teuchos::RCP<VectorBase<Scalar> >                                   vecs[]
  )
{
  initialize(productSpace_in,vecs);
}


template <class Scalar>
void DefaultClusteredSpmdProductVector<Scalar>::initialize(
  const Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  &productSpace_in
  ,const Teuchos::RCP<VectorBase<Scalar> >                                   vecs[]
  )
{
  // ToDo: Validate input!
  productSpace_ = productSpace_in;
  const int numBlocks = productSpace_->numBlocks();
  vecs_.resize(numBlocks);
  if(vecs) {
    std::copy( vecs, vecs + numBlocks, &vecs_[0] );
  }
  else {
    for( int k = 0; k < numBlocks; ++k )
      vecs_[k] = createMember(productSpace_->getBlock(k));
  }
}


template <class Scalar>
void DefaultClusteredSpmdProductVector<Scalar>::uninitialize(
  Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  *productSpace_in
  ,Teuchos::RCP<VectorBase<Scalar> >                                   vecs[]
  )
{
  const int numBlocks = vecs_.size();
  if(productSpace_in) *productSpace_in = productSpace_;
  if(vecs) std::copy( &vecs_[0], &vecs_[0]+numBlocks, vecs );
  productSpace_ = Teuchos::null;
  vecs_.resize(0);
}


// Overridden from ProductVectorBase


template <class Scalar>
Teuchos::RCP<VectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
  using Teuchos::implicit_cast;
  TEST_FOR_EXCEPT( ! ( 0 <= k && k < implicit_cast<int>(vecs_.size()) ) );
  return vecs_[k];
}


template <class Scalar>
Teuchos::RCP<const VectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getVectorBlock(const int k) const
{
  using Teuchos::implicit_cast;
  TEST_FOR_EXCEPT( ! ( 0 <= k && k < implicit_cast<int>(vecs_.size()) ) );
  return vecs_[k];
}


// Overridden from ProductVectorBase


template <class Scalar>
Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}


template <class Scalar>
bool DefaultClusteredSpmdProductVector<Scalar>::blockIsConst(const int k) const
{
  using Teuchos::implicit_cast;
  TEST_FOR_EXCEPT( ! ( 0 <= k && k < implicit_cast<int>(vecs_.size()) ) );
  return false;
}


template <class Scalar>
Teuchos::RCP<MultiVectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}


template <class Scalar>
Teuchos::RCP<const MultiVectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}


// Overridden from VectorBase


template <class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::space() const
{
  return productSpace_;
}


// Overridden protected members from VectorBase


template <class Scalar>
void DefaultClusteredSpmdProductVector<Scalar>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset_in
  ) const
{

  const Ordinal first_ele_offset_in = 0;
  const Ordinal sub_dim_in =  -1;

  using Teuchos::null;
  using Teuchos::ptr_dynamic_cast;
  using Teuchos::tuple;

  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();

  // Validate input
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    global_offset_in < 0, std::invalid_argument,
    "DefaultClusteredSpmdProductVector::applyOp(...): Error, "
    "global_offset_in = " << global_offset_in << " is not acceptable" );
  bool test_failed;
  for (int k = 0; k < num_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*vecs[k]->space());
    TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces,
      "DefaultClusteredSpmdProductVector::applyOp(...): Error vecs["<<k<<"]->space() "
      <<"of type \'"<<typeName(*vecs[k]->space())<<"\' is not compatible with this "
      <<"\'VectorSpaceBlocked\' vector space!"
      );
  }
  for (int k = 0; k < num_targ_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*targ_vecs[k]->space());
    TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultClusteredSpmdProductVector::applyOp(...): Error targ_vecs["<<k<<"]->space() "
      <<"of type \'"<<typeName(*vecs[k]->space())<<"\' is not compatible with this "
      <<"\'VectorSpaceBlocked\' vector space!"
      );
  }
#endif

  //
  // Cast all of the vector arguments to DefaultClusteredSpmdProductVector and
  // make sure that they are all compatible!
  //
  Array<Ptr<const DefaultClusteredSpmdProductVector<Scalar> > > cl_vecs(num_vecs);
  for ( int k = 0; k < num_vecs; ++k ) {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(vecs[k]==null);
#endif
    cl_vecs[k] = ptr_dynamic_cast<const DefaultClusteredSpmdProductVector<Scalar> >(vecs[k],true);
  }
  Array<Ptr<DefaultClusteredSpmdProductVector<Scalar> > > cl_targ_vecs(num_targ_vecs);
  for ( int k = 0; k < num_targ_vecs; ++k ) {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(targ_vecs[k]==null);
#endif
    cl_targ_vecs[k] = ptr_dynamic_cast<DefaultClusteredSpmdProductVector<Scalar> >(targ_vecs[k],true);
  }

  //
  // Get the overlap of the element for this cluster that will participate in
  // the RTOp operation.
  //
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    intraClusterComm = productSpace_->intraClusterComm(),
    interClusterComm = productSpace_->interClusterComm();
  const Ordinal
    clusterSubDim = productSpace_->clusterSubDim(),
    clusterOffset = productSpace_->clusterOffset(),
    globalDim = productSpace_->dim();
  Ordinal  overlap_first_cluster_ele_off  = 0;
  Ordinal  overlap_cluster_sub_dim        = 0;
  Ordinal  overlap_global_off             = 0;
  if (clusterSubDim) {
    RTOp_parallel_calc_overlap(
      globalDim,clusterSubDim,clusterOffset,first_ele_offset_in,sub_dim_in
      ,global_offset_in
      ,&overlap_first_cluster_ele_off,&overlap_cluster_sub_dim,&overlap_global_off
      );
  }

  //
  // Perform the RTOp for each set of block vectors just within this cluster
  // of processes.
  //
  Teuchos::RCP<RTOpPack::ReductTarget> i_reduct_obj;
  if (!is_null(reduct_obj))
    i_reduct_obj = op.reduct_obj_create();
  // Note: i_reduct_obj will accumulate the reduction within this cluster of
  // processes.
  const int numBlocks = vecs_.size();
  if (overlap_first_cluster_ele_off >=0) {

    //
    // There is overlap with at least one element in one block
    // vector for this cluster
    //
    Array<Ptr<const VectorBase<Scalar> > >  v_vecs(num_vecs);
    Array<Ptr<VectorBase<Scalar> > > v_targ_vecs(num_targ_vecs);
    Ordinal overall_global_offset = overlap_global_off;
    for( int j = 0; j < numBlocks; ++j ) {
      const VectorBase<Scalar>
        &v = *vecs_[j];
      const VectorSpaceBase<Scalar>
        &v_space = *v.space();
      // Load up the constutuent block SPMD vectors
      for( int k = 0; k < num_vecs ; ++k )
        v_vecs[k] = cl_vecs[k]->vecs_[j].ptr();
      for( int k = 0; k < num_targ_vecs ; ++k )
        v_targ_vecs[k] = cl_targ_vecs[k]->vecs_[j].ptr();
      TEST_FOR_EXCEPTION(
        numBlocks > 1, std::logic_error
        ,"Error, Have not implemented general support for numBlocks > 1!"
        ); // ToDo: Fix the below code for numBlocks_ > 1!
      Ordinal v_global_offset = overall_global_offset;
      // Apply RTOp on just this cluster
      Thyra::applyOp<Scalar>(
        op, v_vecs(), v_targ_vecs(), i_reduct_obj.ptr(),
        v_global_offset);
      //
      overall_global_offset += v_space.dim();
    }

  }

  //
  // Perform the global reduction across all of the root processes in each of
  // the clusters and then move the global reduction out to each of the
  // processes within the cluster.
  //
  if (!is_null(reduct_obj)) {
    Teuchos::RCP<RTOpPack::ReductTarget>
      icl_reduct_obj = op.reduct_obj_create();
    // First, accumulate the global reduction across all of the elements by
    // just performing the global reduction involving the root processes of
    // each cluster.
    if (interClusterComm.get()) {
      RTOpPack::SPMD_all_reduce(
        &*interClusterComm,
        op,
        1,
        tuple<const RTOpPack::ReductTarget*>(&*i_reduct_obj).getRawPtr(),
        tuple<RTOpPack::ReductTarget*>(&*icl_reduct_obj).getRawPtr()
        );
    }
    // Now the root processes in each cluster have the full global reduction
    // across all elements stored in *icl_reduct_obj and the other processes
    // in each cluster have empty reductions in *icl_reduct_obj.  The last
    // thing to do is to just perform the reduction within each cluster of
    // processes and to add into the in/out *reduct_obj.
    RTOpPack::SPMD_all_reduce(
      &*intraClusterComm,
      op,
      1,
      tuple<const RTOpPack::ReductTarget*>(&*icl_reduct_obj).getRawPtr(),
      tuple<RTOpPack::ReductTarget*>(&*reduct_obj).getRawPtr()
      );
    // ToDo: Replace the above operation with a reduction across clustere into
    // reduct_obj in the root processes and then broadcast the result to all
    // of the slave processes.
  }

}


} // namespace Thyra


#endif // THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_HPP
