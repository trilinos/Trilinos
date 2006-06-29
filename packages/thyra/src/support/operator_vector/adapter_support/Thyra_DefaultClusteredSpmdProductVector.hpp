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

#include "Thyra_DefaultClusteredSpmdProductVectorDecl.hpp"
#include "Thyra_DefaultClusteredSpmdProductVectorSpace.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "RTOp_parallel_helpers.h"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_arrayArg.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template <class Scalar>
DefaultClusteredSpmdProductVector<Scalar>::DefaultClusteredSpmdProductVector()
{
  uninitialize();
}

template <class Scalar>
DefaultClusteredSpmdProductVector<Scalar>::DefaultClusteredSpmdProductVector(
  const Teuchos::RefCountPtr<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  &productSpace
  ,const Teuchos::RefCountPtr<VectorBase<Scalar> >                                   vecs[]
  )
{
  initialize(productSpace,vecs);
}

template <class Scalar>
void DefaultClusteredSpmdProductVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  &productSpace
  ,const Teuchos::RefCountPtr<VectorBase<Scalar> >                                   vecs[]
  )
{
  // ToDo: Validate input!
  productSpace_ = productSpace;
  const int numBlocks = productSpace_->numBlocks();
  vecs_.resize(numBlocks);
  if(vecs) {
    std::copy( vecs, vecs + numBlocks, &vecs_[0] );
  }
  else {
    for( int k = 0; k < numBlocks; ++k )
      vecs_[k] = createMember(productSpace->getBlock(k));
  }
}

template <class Scalar>
void DefaultClusteredSpmdProductVector<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  *productSpace
  ,Teuchos::RefCountPtr<VectorBase<Scalar> >                                   vecs[]
  )
{
  const int numBlocks = vecs_.size();
  if(productSpace) *productSpace = productSpace_;
  if(vecs) std::copy( &vecs_[0], &vecs_[0]+numBlocks, vecs );
  productSpace_ = Teuchos::null;
  vecs_.resize(0);
}

// Overridden from ProductVectorBase

template <class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getNonconstVectorBlock(const int k)
{
  TEST_FOR_EXCEPT( ! ( 0 <= k && k < vecs_.size() ) );
  return vecs_[k];
}

template <class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getVectorBlock(const int k) const
{
  TEST_FOR_EXCEPT( ! ( 0 <= k && k < vecs_.size() ) );
  return vecs_[k];
}

// Overridden from ProductVectorBase

template <class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::productSpace() const
{
  return productSpace_;
}

template <class Scalar>
bool DefaultClusteredSpmdProductVector<Scalar>::blockIsConst(const int k) const
{
  TEST_FOR_EXCEPT( ! ( 0 <= k && k < vecs_.size() ) );
  return false;
}

template <class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getNonconstMultiVectorBlock(const int k)
{
  return getNonconstVectorBlock(k);
}

template <class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::getMultiVectorBlock(const int k) const
{
  return getVectorBlock(k);
}

// Overridden from VectorBase

template <class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultClusteredSpmdProductVector<Scalar>::space() const
{
  return productSpace_;
}

template <class Scalar>
void DefaultClusteredSpmdProductVector<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>    &op
  ,const int                       num_vecs
  ,const VectorBase<Scalar>*const  vecs[]
  ,const int                       num_targ_vecs
  ,VectorBase<Scalar>*const        targ_vecs[]
  ,RTOpPack::ReductTarget          *reduct_obj
  ,const Index                     first_ele_offset_in
  ,const Index                     sub_dim_in
  ,const Index                     global_offset_in
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  //
  const Index	n = productSpace_->dim();
  // Validate input
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !(0 <= first_ele_offset_in && first_ele_offset_in < n), std::out_of_range
    ,"DefaultClusteredSpmdProductVector::applyOp(...): Error, "
    "first_ele_offset_in = " << first_ele_offset_in << " is not in range [0,"<<(n-1)<<"]" );
  TEST_FOR_EXCEPTION(
    global_offset_in < 0, std::invalid_argument
    ,"DefaultClusteredSpmdProductVector::applyOp(...): Error, "
    "global_offset_in = " << global_offset_in << " is not acceptable" );
  TEST_FOR_EXCEPTION(
    sub_dim_in > 0 && sub_dim_in - first_ele_offset_in > n, std::length_error
    ,"DefaultClusteredSpmdProductVector::applyOp(...): Error, "
    "global_offset_in = " << global_offset_in << ", sub_dim_in = " << sub_dim_in
    << "first_ele_offset_in = " << first_ele_offset_in << " and n = " << n
    << " are not compatible" );
  bool test_failed;
  for(int k = 0; k < num_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*vecs[k]->space());
    TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultClusteredSpmdProductVector::applyOp(...): Error vecs["<<k<<"]->space() "
      <<"of type \'"<<typeid(*vecs[k]->space()).name()<<"\' is not compatible with this "
      <<"\'VectorSpaceBlocked\' vector space!"
      );
  }
  for(int k = 0; k < num_targ_vecs; ++k) {
    test_failed = !this->space()->isCompatible(*targ_vecs[k]->space());
    TEST_FOR_EXCEPTION(
      test_failed, Exceptions::IncompatibleVectorSpaces
      ,"DefaultClusteredSpmdProductVector::applyOp(...): Error targ_vecs["<<k<<"]->space() "
      <<"of type \'"<<typeid(*vecs[k]->space()).name()<<"\' is not compatible with this "
      <<"\'VectorSpaceBlocked\' vector space!"
      );
  }
#endif
  //
  // Cast all of the vector arguments to DefaultClusteredSpmdProductVector and
  // make sure that they are all compatible!
  //
  Workspace<const DefaultClusteredSpmdProductVector<Scalar>*> cl_vecs(wss,num_vecs);
  for( int k = 0; k < num_vecs; ++k ) {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(vecs[k]==NULL);
#endif
    cl_vecs[k] = &dyn_cast<const DefaultClusteredSpmdProductVector<Scalar> >(*vecs[k]);
  }
  Workspace<DefaultClusteredSpmdProductVector<Scalar>*> cl_targ_vecs(wss,num_targ_vecs);
  for( int k = 0; k < num_targ_vecs; ++k ) {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPT(targ_vecs[k]==NULL);
#endif
    cl_targ_vecs[k] = &dyn_cast<DefaultClusteredSpmdProductVector<Scalar> >(*targ_vecs[k]);
  }
  //
  // Get the overlap of the element for this cluster that will participate in
  // the RTOp operation.
  //
  const Teuchos::RefCountPtr<const Teuchos::Comm<Index> >
    intraClusterComm = productSpace_->intraClusterComm(),
    interClusterComm = productSpace_->interClusterComm();
  const Index
    clusterSubDim = productSpace_->clusterSubDim(),
    clusterOffset = productSpace_->clusterOffset(),
    globalDim = productSpace_->dim();
  Index  overlap_first_cluster_ele_off  = 0;
  Index  overlap_cluster_sub_dim        = 0;
  Index  overlap_global_off             = 0;
  if(clusterSubDim) {
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
  Teuchos::RefCountPtr<RTOpPack::ReductTarget> i_reduct_obj;
  if(reduct_obj) i_reduct_obj = op.reduct_obj_create();
  // Note: i_reduct_obj will accumulate the reduction within this cluster of
  // processes.
  const int numBlocks = vecs_.size();
  if( overlap_first_cluster_ele_off >=0 ) {
    //
    // There is overlap with at least one element in one block
    // vector for this cluster
    //
    Workspace<const VectorBase<Scalar>*>  v_vecs(wss,num_vecs);
    Workspace<VectorBase<Scalar>*>        v_targ_vecs(wss,num_targ_vecs);
    Index overall_global_offset = overlap_global_off;
    for( int j = 0; j < numBlocks; ++j ) {
      const VectorBase<Scalar>
        &v = *vecs_[j];
      const VectorSpaceBase<Scalar>
        &v_space = *v.space();
      // Load up the constutuent block SPMD vectors
      for( int k = 0; k < num_vecs ; ++k )
        v_vecs[k] = &*cl_vecs[k]->vecs_[j];
      for( int k = 0; k < num_targ_vecs ; ++k )
        v_targ_vecs[k] = &*cl_targ_vecs[k]->vecs_[j];
      TEST_FOR_EXCEPTION(
        numBlocks > 1, std::logic_error
        ,"Error, Have not implemented general support for numBlocks > 1!"
        ); // ToDo: Fix the below code for numBlocks_ > 1!
      Index
        v_global_offset     = overall_global_offset,
        v_first_ele_offset  = overlap_first_cluster_ele_off,
        v_sub_dim           = overlap_cluster_sub_dim;
      // Apply RTOp on just this cluster
      Thyra::applyOp(
        op
        ,num_vecs, num_vecs ? &v_vecs[0] : NULL
        ,num_targ_vecs, num_targ_vecs ? &v_targ_vecs[0] : NULL
        ,i_reduct_obj.get() ? &*i_reduct_obj : NULL
        ,v_first_ele_offset,v_sub_dim,v_global_offset
        );
      //
      overall_global_offset += v_space.dim();
    }
  }
  //
  // Perform the global reduction across all of the root processes in each of
  // the clusters and then move the global reduction out to each of the
  // processes within the cluster.
  //
  if( reduct_obj ) {
    Teuchos::RefCountPtr<RTOpPack::ReductTarget>
      icl_reduct_obj = op.reduct_obj_create();
    // First, accumulate the global reduction across all of the elements by
    // just performing the global reduction involving the root processes of
    // each cluster.
    if(interClusterComm.get()) {
      RTOpPack::SPMD_all_reduce(
        &*interClusterComm
        ,op
        ,1
        ,Teuchos::arrayArg<const RTOpPack::ReductTarget*>(&*i_reduct_obj)()
        ,Teuchos::arrayArg<RTOpPack::ReductTarget*>(&*icl_reduct_obj)()
        );
    }
    // Now the root processes in each cluster have the full global reduction
    // across all elements stored in *icl_reduct_obj and the other processes
    // in each cluster have empty reductions in *icl_reduct_obj.  The last
    // thing to do is to just perform the reduction within each cluster of
    // processes and to add into the in/out *reduct_obj.
    RTOpPack::SPMD_all_reduce(
      &*intraClusterComm
      ,op
      ,1
      ,Teuchos::arrayArg<const RTOpPack::ReductTarget*>(&*icl_reduct_obj)()
      ,Teuchos::arrayArg<RTOpPack::ReductTarget*>(reduct_obj)()
      );
    // ToDo: Replace the above operation with a reduction across clustes into
    // reduct_obj in the root processes and then broadcast the result to all
    // of the slave processes.
  }
}

} // namespace Thyra

#endif // THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_HPP
