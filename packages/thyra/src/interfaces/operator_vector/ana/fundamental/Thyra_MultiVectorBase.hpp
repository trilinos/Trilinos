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

#ifndef THYRA_MULTI_VECTOR_BASE_HPP
#define THYRA_MULTI_VECTOR_BASE_HPP

#include "Thyra_MultiVectorBaseDecl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"

namespace Thyra {

// Provide access to the columns as VectorBase objects

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
MultiVectorBase<Scalar>::col(Index j) const
{
  return const_cast<MultiVectorBase*>(this)->col(j);
}

// Cloning

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
MultiVectorBase<Scalar>::clone_mv() const
{
  const VectorSpaceBase<Scalar>
    &domain = *this->domain(),
    &range  = *this->range();
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    copy = createMembers(range,domain.dim());
  assign( &*copy, *this );
  return copy;
}

// Collective applyOp() methods

template<class Scalar>
void MultiVectorBase<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>   &prim_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget*        reduct_objs[]
  ,const Index                    prim_first_ele_in
  ,const Index                    prim_sub_dim_in
  ,const Index                    prim_global_offset_in
  ,const Index                    sec_first_ele_in
  ,const Index                    sec_sub_dim_in
  ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  // ToDo: Validate the input!

  const VectorSpaceBase<Scalar>  &domain = *this->domain();

  // Get the primary and secondary dimensions.

  const Index sec_dim = domain.dim();
  const Index sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  -  sec_first_ele_in + 1  );
#ifdef _DEBUG
  const VectorSpaceBase<Scalar>  &range  = *this->range();
  const Index	prim_dim = range.dim();
  const Index prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim - prim_first_ele_in + 1 );
  const char err_msg[] = "MultiVectorBase<Scalar>::applyOp(...): Error!";
  TEST_FOR_EXCEPTION( !(0 < prim_sub_dim && prim_sub_dim <= prim_dim), std::invalid_argument, err_msg );
  TEST_FOR_EXCEPTION( !(0 < sec_sub_dim  && sec_sub_dim  <= sec_dim),  std::invalid_argument, err_msg );
#endif

  //
  // Apply the reduction/transformation operator and transform the
  // target vectors and reduce each of the reduction objects.
  //

  Workspace< Teuchos::RefCountPtr<const VectorBase<Scalar> > >   vecs_s(wss,num_multi_vecs);
  Workspace<const VectorBase<Scalar>*>                           vecs(wss,num_multi_vecs,false);
  Workspace< Teuchos::RefCountPtr<VectorBase<Scalar> > >         targ_vecs_s(wss,num_targ_multi_vecs);
  Workspace<VectorBase<Scalar>*>                                 targ_vecs(wss,num_targ_multi_vecs,false);

  for(Index j = sec_first_ele_in; j <= sec_first_ele_in - 1 + sec_sub_dim; ++j) {
    // Fill the arrays of vector arguments
    {for(Index k = 0; k < static_cast<Index>(num_multi_vecs); ++k) {
      vecs_s[k] = multi_vecs[k]->col(j);
      vecs[k] = vecs_s[k].get();
    }}
    {for(Index k = 0; k < static_cast<Index>(num_targ_multi_vecs); ++k) {
      targ_vecs_s[k] = targ_multi_vecs[k]->col(j);
      targ_vecs[k] = targ_vecs_s[k].get();
    }}
    // Apply the reduction/transformation operator
    Thyra::applyOp(
      prim_op
      ,num_multi_vecs,      (num_multi_vecs      ? &vecs[0]      : NULL)
      ,num_targ_multi_vecs, (num_targ_multi_vecs ? &targ_vecs[0] : NULL)
      ,reduct_objs ? reduct_objs[j-1] : NULL
      ,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
      );
  }
  // At this point all of the designated targ vectors in the target multi-vectors have
  // been transformed and all the reduction objects in reduct_obj[] have accumulated
  // the reductions.
}

template<class Scalar>
void MultiVectorBase<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>   &prim_op
  ,const RTOpPack::RTOpT<Scalar>  &sec_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const Index                    prim_first_ele_in
  ,const Index                    prim_sub_dim_in
  ,const Index                    prim_global_offset_in
  ,const Index                    sec_first_ele_in
  ,const Index                    sec_sub_dim_in
  ) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  // ToDo: Validate the input!

  const VectorSpaceBase<Scalar> &domain = *this->domain();

  // Get the primary and secondary dimensions.
  const Index sec_dim = domain.dim();
  const Index sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  -  sec_first_ele_in + 1  );
#ifdef _DEBUG
  const VectorSpaceBase<Scalar> &range = *this->range();
  const Index prim_dim = range.dim();
  const Index prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim - prim_first_ele_in + 1 );
  const char err_msg[] = "MultiVectorBase<Scalar>::applyOp(...): Error!";
  TEST_FOR_EXCEPTION( !(0 < prim_sub_dim && prim_sub_dim <= prim_dim), std::invalid_argument, err_msg );
  TEST_FOR_EXCEPTION( !(0 < sec_sub_dim  && sec_sub_dim  <= sec_dim),  std::invalid_argument, err_msg );
#endif

  // Create a temporary buffer for the reduction objects of the primary reduction
  // so that we can call the companion version of this method.
  Workspace<Teuchos::RefCountPtr<RTOpPack::ReductTarget> >
    rcp_reduct_objs(wss,reduct_obj!=NULL?sec_sub_dim:0);
  Workspace<RTOpPack::ReductTarget*>
    reduct_objs(wss,reduct_obj!=NULL?sec_sub_dim:0,false);
  if(reduct_obj) {
    for(Index k = 0; k < sec_sub_dim; ++k) {
      rcp_reduct_objs[k] = prim_op.reduct_obj_create();
      reduct_objs[k] = &*rcp_reduct_objs[k];
    }
  }
  
  // Call the companion version that accepts an array of reduction objects
  this->applyOp(
    prim_op
    ,num_multi_vecs,       multi_vecs
    ,num_targ_multi_vecs,  targ_multi_vecs
    ,reduct_obj ? &reduct_objs[0] : NULL
    ,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
    ,sec_first_ele_in,  sec_sub_dim_in
    );

  // Reduce all the reduction objects using the secondary reduction operator
  // into one reduction object and free the intermediate reduction objects.
  if(reduct_obj) {
    for(Index k = 0; k < sec_sub_dim; ++k) {
      sec_op.reduce_reduct_objs( *reduct_objs[k], reduct_obj );
    }
  }
}

// Explicit sub-multi-vector access

template<class Scalar>
void MultiVectorBase<Scalar>::getSubMultiVector(
  const Range1D                       &rowRng_in
  ,const Range1D                      &colRng_in
  ,RTOpPack::SubMultiVectorT<Scalar>  *sub_mv
  ) const
{
  const Index
    rangeDim  = this->range()->dim(),
    domainDim = this->domain()->dim();
  const Range1D
    rowRng = rowRng_in.full_range() ? Range1D(1,rangeDim)  : rowRng_in,
    colRng = colRng_in.full_range() ? Range1D(1,domainDim) : colRng_in;
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(
    rowRng.ubound() > rangeDim, std::out_of_range
    ,"MultiVectorBase<Scalar>::getSubMultiVector(...): Error, rowRng = ["
    <<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not in the range = [1,"
    <<rangeDim<<"]!"
    );
  TEST_FOR_EXCEPTION(
    colRng.ubound() > domainDim, std::out_of_range
    ,"MultiVectorBase<Scalar>::getSubMultiVector(...): Error, colRng = ["
    <<colRng.lbound()<<","<<colRng.ubound()<<"] is not in the range = [1,"
    <<domainDim<<"]!"
    );
#endif
  // Allocate storage for the multi-vector (stored column-major)
  Scalar *values = new Scalar[ rowRng.size() * colRng.size() ];
  // Extract multi-vector values column by column
  RTOpPack::SubVectorT<Scalar> sv; // uninitialized by default
  for( int k = colRng.lbound(); k <= colRng.ubound(); ++k ) {
    Teuchos::RefCountPtr<const VectorBase<Scalar> > col_k = this->col(k);
    col_k->getSubVector( rowRng, &sv );
    for( int i = 0; i < rowRng.size(); ++i )
      values[ i + (k-1)*rowRng.size() ] = sv[i];
    col_k->freeSubVector( &sv );
  }
  // Initialize the multi-vector view object
  sub_mv->initialize(
    rowRng.lbound()-1            // globalOffset
    ,rowRng.size()               // subDim
    ,colRng.lbound()-1           // colOffset
    ,colRng.size()               // numSubCols
    ,values                      // values
    ,rowRng.size()               // leadingDim
    );
}

template<class Scalar>
void MultiVectorBase<Scalar>::freeSubMultiVector(
  RTOpPack::SubMultiVectorT<Scalar>* sub_mv
  ) const
{
  // Here we just need to free the view and that is it!
  delete [] const_cast<Scalar*>(sub_mv->values());
  sub_mv->set_uninitialized();
}

template<class Scalar>
void MultiVectorBase<Scalar>::getSubMultiVector(
  const Range1D                                &rowRng
  ,const Range1D                               &colRng
  ,RTOpPack::MutableSubMultiVectorT<Scalar>    *sub_mv
  )
{
  // Use the non-const implementation since it does exactly the
  // correct thing in this case also!
  MultiVectorBase<Scalar>::getSubMultiVector(
    rowRng, colRng
    ,static_cast<RTOpPack::SubMultiVectorT<Scalar>*>(sub_mv)
    // This cast will work as long as MutableSubMultiVectorT
    // maintains no extra state over SubMultiVectorT (which it
    // currently does not) but this is something that I should
    // technically check for some how.
    );
}

template<class Scalar>
void MultiVectorBase<Scalar>::commitSubMultiVector(
  RTOpPack::MutableSubMultiVectorT<Scalar>* sub_mv
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(
    sub_mv==NULL, std::logic_error, "MultiVectorBase<Scalar>::commitSubMultiVector(...): Error!"
    );
#endif
  // Set back the multi-vector values column by column
  const Range1D rowRng(sub_mv->globalOffset()+1,sub_mv->globalOffset()+sub_mv->subDim());
  RTOpPack::MutableSubVectorT<Scalar> msv; // uninitialized by default
  for( int k = sub_mv->colOffset()+1; k <= sub_mv->numSubCols(); ++k ) {
    Teuchos::RefCountPtr<VectorBase<Scalar> > col_k = this->col(k);
    col_k->getSubVector( rowRng, &msv );
    for( int i = 0; i < rowRng.size(); ++i )
      msv[i] = sub_mv->values()[ i + (k-1)*rowRng.size() ];
    col_k->commitSubVector( &msv );
  }
  // Free the memory
  delete [] const_cast<Scalar*>(sub_mv->values());
  // Zero out the view
  sub_mv->set_uninitialized();
}

// Overridden methods from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
MultiVectorBase<Scalar>::clone() const
{
  return this->clone_mv();
}

} // end namespace Thyra

#endif // THYRA_MULTI_VECTOR_BASE_HPP
