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

#ifndef THYRA_SPMD_VECTOR_BASE_HPP
#define THYRA_SPMD_VECTOR_BASE_HPP

#include "Thyra_SpmdVectorBaseDecl.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_apply_op_helper.hpp"
#include "RTOp_parallel_helpers.h"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif // THYRA_SPMD_VECTOR_BASE_DUMP

namespace Thyra {

template<class Scalar>
SpmdVectorBase<Scalar>::SpmdVectorBase()
  :in_applyOp_(false)
  ,globalDim_(0)
  ,localOffset_(-1)
  ,localSubDim_(0)
{}

// Virtual methods with default implementations

template<class Scalar>
void SpmdVectorBase<Scalar>::getLocalData( const Scalar** values, Index* stride ) const
{
  const_cast<SpmdVectorBase<Scalar>*>(this)->getLocalData(const_cast<Scalar**>(values),stride);
}

template<class Scalar>
void SpmdVectorBase<Scalar>::freeLocalData( const Scalar* values ) const
{
  const_cast<SpmdVectorBase<Scalar>*>(this)->commitLocalData(const_cast<Scalar*>(values));
}

template<class Scalar>
void SpmdVectorBase<Scalar>::applyOp(
  const Teuchos::Comm<Index>      *comm_in
  ,const RTOpPack::RTOpT<Scalar>  &op
  ,const int                      num_vecs
  ,const VectorBase<Scalar>*const vecs[]
  ,const int                      num_targ_vecs
  ,VectorBase<Scalar>*const       targ_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const Index                    first_ele_offset_in
  ,const Index                    sub_dim_in
  ,const Index                    global_offset_in
  ) const
{
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(show_dump) {
    *out << "\nEntering SpmdVectorBase<Scalar>::applyOp(...) ...\n";
    *out
      << "\nop = " << typeName(op)
      << "\nnum_vecs = " << num_vecs
      << "\nnum_targ_vecs = " << num_targ_vecs
      << "\nreduct_obj = " << reduct_obj
      << "\nfirst_ele_offset_in = " << first_ele_offset_in
      << "\nsub_dim_in = " << sub_dim_in
      << "\nglobal_offset_in = " << global_offset_in
      << "\n"
      ;
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  const SpmdVectorSpaceBase<Scalar> &spmdSpc = *spmdSpace();
#ifdef TEUCHOS_DEBUG
  // ToDo: Validate input!
  TEST_FOR_EXCEPTION(
    in_applyOp_, std::invalid_argument
    ,"SpmdVectorBase<>::applyOp(...): Error, this method is being entered recursively which is a "
    "clear sign that one of the methods acquireDetachedView(...), releaseDetachedView(...) or commitDetachedView(...) "
    "was not implemented properly!"
    );
  Thyra::apply_op_validate_input(
    "SpmdVectorBase<>::applyOp(...)",*space()
    ,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele_offset_in,sub_dim_in,global_offset_in
    );
#endif
  Teuchos::RefCountPtr<const Teuchos::Comm<Index> > comm;
  if( comm_in != NULL )
    comm = Teuchos::rcp(comm_in,false);
  else
    comm = spmdSpc.getComm();
  // Flag that we are in applyOp()
  in_applyOp_ = true;
  // First see if this is a locally replicated vector in which case
  // we treat this as a local operation only.
  const bool locallyReplicated = ( comm_in == NULL && localSubDim_ == globalDim_ );
  // Get the overlap in the current process with the input logical sub-vector
  // from (first_ele_offset_in,sub_dim_in,global_offset_in)
  Index  overlap_first_local_ele_off  = 0;
  Index  overlap_local_sub_dim        = 0;
  Index  overlap_global_off           = 0;
  if(localSubDim_) {
    RTOp_parallel_calc_overlap(
      globalDim_, localSubDim_, localOffset_, first_ele_offset_in, sub_dim_in, global_offset_in
      ,&overlap_first_local_ele_off, &overlap_local_sub_dim, &overlap_global_off
      );
  }
  const Range1D local_rng = (
    overlap_first_local_ele_off>=0
    ? Range1D( localOffset_ + overlap_first_local_ele_off, localOffset_ + overlap_first_local_ele_off + overlap_local_sub_dim - 1 )
    : Range1D::Invalid
    );
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  if(show_dump) {
    *out
      << "\noverlap_first_local_ele_off = " << overlap_first_local_ele_off
      << "\noverlap_local_sub_dim = " << overlap_local_sub_dim
      << "\noverlap_global_off = " << overlap_global_off
      << "\nlocal_rng = ["<<local_rng.lbound()<<","<<local_rng.ubound()<<"]"
      << "\n"
      ;
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP
  // Create sub-vector views of all of the *participating* local data
  Workspace<RTOpPack::ConstSubVectorView<Scalar> > sub_vecs(wss,num_vecs);
  Workspace<RTOpPack::SubVectorView<Scalar> > sub_targ_vecs(wss,num_targ_vecs);
  if( overlap_first_local_ele_off >= 0 ) {
    {for(int k = 0; k < num_vecs; ++k ) {
      vecs[k]->acquireDetachedView( local_rng, &sub_vecs[k] );
      sub_vecs[k].setGlobalOffset( overlap_global_off );
    }}
    {for(int k = 0; k < num_targ_vecs; ++k ) {
      targ_vecs[k]->acquireDetachedView( local_rng, &sub_targ_vecs[k] );
      sub_targ_vecs[k].setGlobalOffset( overlap_global_off );
    }}
  }
  // Apply the RTOp operator object (all processors must participate)
  const bool validSubVecs = ( localSubDim_ > 0 && overlap_first_local_ele_off >= 0 );
  RTOpPack::SPMD_apply_op(
    locallyReplicated ? NULL : comm.get()                       // comm
    ,op                                                         // op
    ,num_vecs                                                   // num_vecs
    ,num_vecs && validSubVecs ? &sub_vecs[0] : NULL             // sub_vecs
    ,num_targ_vecs                                              // num_targ_vecs
    ,num_targ_vecs && validSubVecs ? &sub_targ_vecs[0] : NULL   // targ_sub_vecs
    ,reduct_obj                                                 // reduct_obj
    );
  // Free and commit the local data
  if( overlap_first_local_ele_off >= 0 ) {
    {for(int k = 0; k < num_vecs; ++k ) {
      sub_vecs[k].setGlobalOffset(local_rng.lbound());
      vecs[k]->releaseDetachedView( &sub_vecs[k] );
    }}
    {for(int k = 0; k < num_targ_vecs; ++k ) {
      sub_targ_vecs[k].setGlobalOffset(local_rng.lbound());
      targ_vecs[k]->commitDetachedView( &sub_targ_vecs[k] );
    }}
  }
  // Flag that we are leaving applyOp()
  in_applyOp_ = false;
#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
  if(show_dump) {
    *out << "\nLeaving SpmdVectorBase<Scalar>::applyOp(...) ...\n";
  }
#endif // THYRA_SPMD_VECTOR_BASE_DUMP
}

// Overridden from VectorBase

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
SpmdVectorBase<Scalar>::space() const
{
  return spmdSpace();
}

template<class Scalar>
void SpmdVectorBase<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>   &op
  ,const int                      num_vecs
  ,const VectorBase<Scalar>*const vecs[]
  ,const int                      num_targ_vecs
  ,VectorBase<Scalar>*const       targ_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const Index                    first_ele_offset
  ,const Index                    sub_dim
  ,const Index                    global_offset
  ) const
{
  applyOp(
    NULL
    ,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
    ,first_ele_offset,sub_dim,global_offset
    );
}

template<class Scalar>
void SpmdVectorBase<Scalar>::acquireDetachedView( const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const
{
  if( rng_in == Range1D::Invalid ) {
    // Just return an null view
    sub_vec->initialize(
      rng_in.lbound()    // globalOffset
      ,0                 // subDim
      ,0                 // values
      ,1                 // stride
      );
    return;
  }
  const Range1D rng = validateRange(rng_in);
  if( rng.lbound() < localOffset_ || localOffset_+localSubDim_-1 < rng.ubound() ) {
    // rng consists of off-processor elements so use the default implementation!
    VectorDefaultBase<Scalar>::acquireDetachedView(rng_in,sub_vec);
    return;
  }
  // rng consists of all local data so get it!
  const Scalar *localValues = NULL;
  Index stride = 0;
  this->getLocalData(&localValues,&stride);
  sub_vec->initialize(
    rng.lbound()                               // globalOffset
    ,rng.size()                                // subDim
    ,localValues+(rng.lbound()-localOffset_)   // values
    ,stride                                    // stride
    );
}

template<class Scalar>
void SpmdVectorBase<Scalar>::releaseDetachedView( RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    sub_vec==NULL || sub_vec->globalOffset() < 0 || sub_vec->globalOffset() + sub_vec->subDim() > globalDim_
    ,std::logic_error
    ,"SpmdVectorBase<Scalar>::releaseDetachedView(...) : Error, this sub vector was not gotten from acquireDetachedView(...)!"
    );
#endif
  if( sub_vec->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim() ) {
    // Let the default implementation handle it!
    VectorDefaultBase<Scalar>::releaseDetachedView(sub_vec);
    return;
  }
  // Nothing to deallocate!
  sub_vec->set_uninitialized();
}

template<class Scalar>
void SpmdVectorBase<Scalar>::acquireDetachedView( const Range1D& rng_in, RTOpPack::SubVectorView<Scalar>* sub_vec )
{
  if( rng_in == Range1D::Invalid ) {
    // Just return an null view
    sub_vec->initialize(
      rng_in.lbound()    // globalOffset
      ,0                 // subDim
      ,0                 // values
      ,1                 // stride
      );
    return;
  }
  const Range1D rng = validateRange(rng_in);
  if( rng.lbound() < localOffset_ || localOffset_+localSubDim_-1 < rng.ubound() ) {
    // rng consists of off-processor elements so use the default implementation!
    VectorDefaultBase<Scalar>::acquireDetachedView(rng_in,sub_vec);
    return;
  }
  // rng consists of all local data so get it!
  Scalar *localValues = NULL;
  Index stride = 0;
  this->getLocalData(&localValues,&stride);
  sub_vec->initialize(
    rng.lbound()                               // globalOffset
    ,rng.size()                                // subDim
    ,localValues+(rng.lbound()-localOffset_)   // values
    ,stride                                    // stride
    );
}

template<class Scalar>
void SpmdVectorBase<Scalar>::commitDetachedView( RTOpPack::SubVectorView<Scalar>* sub_vec )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    sub_vec==NULL || sub_vec->globalOffset() < 0 || sub_vec->globalOffset() + sub_vec->subDim() > globalDim_
    ,std::logic_error
    ,"SpmdVectorBase<Scalar>::commitDetachedView(...) : Error, this sub vector was not gotten from acquireDetachedView(...)!"
    );
#endif
  if( sub_vec->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim() ) {
    // Let the default implementation handle it!
    VectorDefaultBase<Scalar>::commitDetachedView(sub_vec);
    return;
  }
  sub_vec->set_uninitialized();  // Nothing to deallocate!
}

// protected

template<class Scalar>
void SpmdVectorBase<Scalar>::updateSpmdSpace()
{
  if(globalDim_ == 0) {
    const SpmdVectorSpaceBase<Scalar> *spmdSpace = this->spmdSpace().get();
    if(spmdSpace) {
      globalDim_    = spmdSpace->dim();
      localOffset_  = spmdSpace->localOffset();
      localSubDim_  = spmdSpace->localSubDim();
    }
    else {
      globalDim_    = 0;
      localOffset_  = -1;
      localSubDim_  = 0;
    }
  }
}

// private

template<class Scalar>
Range1D SpmdVectorBase<Scalar>::validateRange( const Range1D &rng_in ) const
{
  const Range1D rng = Teuchos::full_range(rng_in,0,globalDim_-1);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !(0 <= rng.lbound() && rng.ubound() < globalDim_), std::invalid_argument
    ,"SpmdVectorBase<Scalar>::validateRange(...): Error, the range ["
    <<rng.lbound()<<","<<rng.ubound()<<"] is not "
    "in the range [0,"<<(globalDim_-1)<<"]!"
    );
#endif
  return rng;
}

#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
template<class Scalar>
bool SpmdVectorBase<Scalar>::show_dump = false;
#endif // THYRA_SPMD_VECTOR_BASE_DUMP

} // end namespace Thyra

#endif // THYRA_SPMD_VECTOR_BASE_HPP
