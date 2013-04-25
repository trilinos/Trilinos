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

#ifndef THYRA_SPMD_VECTOR_DEFAULT_BASE_DEF_HPP
#define THYRA_SPMD_VECTOR_DEFAULT_BASE_DEF_HPP


#include "Thyra_SpmdVectorDefaultBase_decl.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_apply_op_helper.hpp"
#include "Thyra_SpmdLocalDataAccess.hpp"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


// Public interface functions


template<class Scalar>
SpmdVectorDefaultBase<Scalar>::SpmdVectorDefaultBase()
  :in_applyOpImpl_(false)
  ,globalDim_(0)
  ,localOffset_(-1)
  ,localSubDim_(0)
{}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::applyOpImplWithComm(
  const Ptr<const Teuchos::Comm<Ordinal> > &comm_in,
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset_in
  ) const
{

  using Teuchos::null;
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  using Teuchos::rcpFromPtr;

  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();

  Ptr<Teuchos::WorkspaceStore> wss = Teuchos::get_default_workspace_store().ptr();
  const SpmdVectorSpaceBase<Scalar> &spmdSpc = *this->spmdSpace();

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    in_applyOpImpl_, std::invalid_argument,
    "SpmdVectorDefaultBase<>::applyOp(...): Error, this method is being entered recursively"
    " which is a clear sign that one of the methods acquireDetachedView(...),"
    " releaseDetachedView(...) or commitDetachedView(...) was not implemented properly!"
    );
  Thyra::apply_op_validate_input(
    "SpmdVectorDefaultBase<>::applyOp(...)",*space(),
    op, vecs, targ_vecs, reduct_obj, global_offset_in);
#endif

  Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm;
  if (nonnull(comm_in))
    comm = Teuchos::rcpFromPtr(comm_in);
  else
    comm = spmdSpc.getComm();

  // Flag that we are in applyOp()
  in_applyOpImpl_ = true;

  // First see if this is a locally replicated vector in which case
  // we treat this as a local operation only.
  const bool locallyReplicated = spmdSpc.isLocallyReplicated();

  const Range1D local_rng(localOffset_, localOffset_+localSubDim_-1);

  // Create sub-vector views of all of the *participating* local data
  Workspace<RTOpPack::ConstSubVectorView<Scalar> > sub_vecs(wss.get(), num_vecs);
  Workspace<RTOpPack::SubVectorView<Scalar> > sub_targ_vecs(wss.get(), num_targ_vecs);
  for(int k = 0; k < num_vecs; ++k ) {
    sub_vecs[k] = getLocalSubVectorView<Scalar>(rcpFromPtr(vecs[k]));
    //vecs[k]->acquireDetachedView( local_rng, &sub_vecs[k] );
    sub_vecs[k].setGlobalOffset(localOffset_+global_offset_in);
  }
  for(int k = 0; k < num_targ_vecs; ++k ) {
    sub_targ_vecs[k] = getNonconstLocalSubVectorView<Scalar>(rcpFromPtr(targ_vecs[k]));
    //targ_vecs[k]->acquireDetachedView( local_rng, &sub_targ_vecs[k] );
    sub_targ_vecs[k].setGlobalOffset(localOffset_+global_offset_in);
  }

  // Apply the RTOp operator object (all processors must participate)
  RTOpPack::SPMD_apply_op(
    locallyReplicated ? NULL : &*comm,     // comm
    op,                                    // op
    num_vecs,                              // num_vecs
    sub_vecs.getRawPtr(),                  // sub_vecs
    num_targ_vecs,                         // num_targ_vecs
    sub_targ_vecs.getRawPtr(),             // targ_sub_vecs
    reduct_obj.get()                       // reduct_obj
    );

  // Free and commit the local data
  for (int k = 0; k < num_vecs; ++k ) {
    sub_vecs[k] = RTOpPack::ConstSubVectorView<Scalar>();
    //sub_vecs[k].setGlobalOffset(local_rng.lbound());
    //vecs[k]->releaseDetachedView(&sub_vecs[k]);
  }
  for (int k = 0; k < num_targ_vecs; ++k ) {
    sub_targ_vecs[k] = RTOpPack::SubVectorView<Scalar>();
    //sub_targ_vecs[k].setGlobalOffset(local_rng.lbound());
    //targ_vecs[k]->commitDetachedView(&sub_targ_vecs[k]);
  }

  // Flag that we are leaving applyOp()
  in_applyOpImpl_ = false;

}


// Overridden from Teuchos::Describable


template<class Scalar>
std::string SpmdVectorDefaultBase<Scalar>::description() const
{
  using Teuchos::RCP; using Teuchos::Comm; using Teuchos::null;
  using Teuchos::typeName;
  std::ostringstream ostr;
  ostr<<typeName(*this)<<"{spmdSpace="<<this->spmdSpace()->description()<<"}";
  return ostr.str();
}


// Overridden public functions from VectorBase


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
SpmdVectorDefaultBase<Scalar>::space() const
{
  return this->spmdSpace();
}


// protected


// Overridden protected functions from VectorBase


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::applyOpImpl(
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset
  ) const
{
  applyOpImplWithComm( Teuchos::null, op, vecs, targ_vecs, reduct_obj,
    global_offset);
}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(sub_vec);
#endif
  if( rng_in == Range1D::Invalid ) {
    // Just return an null view
    *sub_vec = RTOpPack::ConstSubVectorView<Scalar>();
    return;
  }
  const bool isFullRng = rng_in.full_range(); 
  const Range1D rng = validateRange(rng_in);
  const bool isLocallyReplicated = this->spmdSpace()->isLocallyReplicated();
  if (
    (
      rng.lbound() < localOffset_ 
      ||
      localOffset_+localSubDim_-1 < rng.ubound()
      ||
      isFullRng
      )
    &&
    !isLocallyReplicated
    )
  {
    // rng consists of off-processor elements so use the default implementation!
    VectorDefaultBase<Scalar>::acquireDetachedVectorViewImpl(rng_in, sub_vec);
    return;
  }
  // rng consists of all local data so get it!
  ArrayRCP<const Scalar>localValues;
  this->getLocalData(Teuchos::outArg(localValues));
  sub_vec->initialize(
    rng.lbound(),  // globalOffset
    rng.size(),  // subDim
    localValues.persistingView(rng.lbound()-localOffset_, rng.size()),
    1  // stride
    );
}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::releaseDetachedVectorViewImpl(
  RTOpPack::ConstSubVectorView<Scalar>* sub_vec
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    sub_vec==NULL || sub_vec->globalOffset() < 0 || sub_vec->globalOffset() + sub_vec->subDim() > globalDim_
    ,std::logic_error
    ,"SpmdVectorDefaultBase<Scalar>::releaseDetachedVectorViewImpl(...) : Error, this sub vector was not gotten from acquireDetachedView(...)!"
    );
#endif
  if(
    sub_vec->globalOffset() < localOffset_ 
    || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim()
    )
  {
    // Let the default implementation handle it!
    VectorDefaultBase<Scalar>::releaseDetachedVectorViewImpl(sub_vec);
    return;
  }
  // Nothing to deallocate!
  sub_vec->uninitialize();
}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::acquireNonconstDetachedVectorViewImpl(
  const Range1D& rng_in, RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(sub_vec);
#endif
  if( rng_in == Range1D::Invalid ) {
    // Just return an null view
    *sub_vec = RTOpPack::SubVectorView<Scalar>();
    return;
  }
  const bool isFullRng = rng_in.full_range(); 
  const Range1D rng = validateRange(rng_in);
  const bool isLocallyReplicated = this->spmdSpace()->isLocallyReplicated();
  if (
    (
      rng.lbound() < localOffset_ 
      ||
      localOffset_+localSubDim_-1 < rng.ubound()
      ||
      isFullRng
      )
    &&
    !isLocallyReplicated
    )
  {
    // rng consists of off-processor elements so use the default implementation!
    VectorDefaultBase<Scalar>::acquireNonconstDetachedVectorViewImpl(rng_in,sub_vec);
    return;
  }
  // rng consists of all local data so get it!
  ArrayRCP<Scalar> localValues;
  this->getNonconstLocalData(Teuchos::outArg(localValues));
  sub_vec->initialize(
    rng.lbound(),  // globalOffset
    rng.size(),  // subDim
    localValues.persistingView(rng.lbound()-localOffset_, rng.size()),
    1  // stride
    );
}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::commitNonconstDetachedVectorViewImpl(
  RTOpPack::SubVectorView<Scalar>* sub_vec
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    sub_vec==NULL || sub_vec->globalOffset() < 0 || sub_vec->globalOffset() + sub_vec->subDim() > globalDim_
    ,std::logic_error
    ,"SpmdVectorDefaultBase<Scalar>::commitDetachedView(...) : Error, this sub vector was not gotten from acquireDetachedView(...)!"
    );
#endif
  if(
    sub_vec->globalOffset() < localOffset_
    ||
    localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim()
    )
  {
    // Let the default implementation handle it!
    VectorDefaultBase<Scalar>::commitNonconstDetachedVectorViewImpl(sub_vec);
    return;
  }
  sub_vec->uninitialize();  // Nothing to deallocate!
}


// Overridden Protected functions from SpmdMultiVectorBase


template<class Scalar>
RTOpPack::SubMultiVectorView<Scalar>
SpmdVectorDefaultBase<Scalar>::getNonconstLocalSubMultiVectorImpl()
{
  ArrayRCP<Scalar> localValues;
  this->getNonconstLocalData(Teuchos::outArg(localValues));
  return RTOpPack::SubMultiVectorView<Scalar>(
    localOffset_, // globalOffset
    localSubDim_,
    0, // colOffset
    1, // numCols
    localValues,
    localSubDim_ // leadingDim
    );
}


template<class Scalar>
RTOpPack::ConstSubMultiVectorView<Scalar>
SpmdVectorDefaultBase<Scalar>::getLocalSubMultiVectorImpl() const
{
  using Teuchos::outArg;
  ArrayRCP<const Scalar> localValues;
  this->getLocalData(outArg(localValues));
  return RTOpPack::ConstSubMultiVectorView<Scalar>(
    localOffset_, // globalOffset
    localSubDim_,
    0, // colOffset
    1, // numCols
    localValues,
    localSubDim_ // leadingDim
    );
}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::getNonconstLocalMultiVectorDataImpl(
  const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim)
{
  this->getNonconstLocalData(localValues);
  *leadingDim = localValues->size();
}


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::getLocalMultiVectorDataImpl(
  const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim) const
{
  this->getLocalData(localValues);
  *leadingDim = localValues->size();
}


// Overridden Protected functions from SpmdVectorBase


template<class Scalar>
RTOpPack::SubVectorView<Scalar>
SpmdVectorDefaultBase<Scalar>::getNonconstLocalSubVectorImpl()
{
  ArrayRCP<Scalar> localValues;
  this->getNonconstLocalData(Teuchos::outArg(localValues));
  return RTOpPack::SubVectorView<Scalar>(
    localOffset_,
    localSubDim_,
    localValues,
    1 // stride
    );
}


template<class Scalar>
RTOpPack::ConstSubVectorView<Scalar>
SpmdVectorDefaultBase<Scalar>::getLocalSubVectorImpl() const
{
  ArrayRCP<const Scalar> localValues;
  this->getLocalData(Teuchos::outArg(localValues));
  return RTOpPack::ConstSubVectorView<Scalar>(
    localOffset_, // globalOffset?
    localSubDim_,
    localValues,
    1 // stride
    );
}


// Protected functions to be used by subclasses


template<class Scalar>
void SpmdVectorDefaultBase<Scalar>::updateSpmdSpace()
{
  if(globalDim_ == 0) {
    const SpmdVectorSpaceBase<Scalar> *l_spmdSpace = this->spmdSpace().get();
    if(l_spmdSpace) {
      globalDim_    = l_spmdSpace->dim();
      localOffset_  = l_spmdSpace->localOffset();
      localSubDim_  = l_spmdSpace->localSubDim();
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
Range1D SpmdVectorDefaultBase<Scalar>::validateRange( const Range1D &rng_in ) const
{
  const Range1D rng = Teuchos::full_range(rng_in,0,globalDim_-1);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= rng.lbound() && rng.ubound() < globalDim_), std::invalid_argument
    ,"SpmdVectorDefaultBase<Scalar>::validateRange(...): Error, the range ["
    <<rng.lbound()<<","<<rng.ubound()<<"] is not "
    "in the range [0,"<<(globalDim_-1)<<"]!"
    );
#endif
  return rng;
}


#ifdef THYRA_SPMD_VECTOR_BASE_DUMP
template<class Scalar>
bool SpmdVectorDefaultBase<Scalar>::show_dump = false;
#endif // THYRA_SPMD_VECTOR_BASE_DUMP


} // end namespace Thyra


#endif // THYRA_SPMD_VECTOR_DEFAULT_BASE_DEF_HPP
