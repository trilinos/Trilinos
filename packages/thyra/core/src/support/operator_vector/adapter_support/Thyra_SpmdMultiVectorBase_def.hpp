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

#ifndef THYRA_SPMD_MULTI_VECTOR_BASE_DEF_HPP
#define THYRA_SPMD_MULTI_VECTOR_BASE_DEF_HPP


#include "Thyra_SpmdMultiVectorBase_decl.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Thyra_MultiVectorAdapterBase.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_apply_op_helper.hpp"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "RTOp_parallel_helpers.h"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_CommHelpers.hpp"


// Define to see some timing output!
//#define THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES


namespace Thyra {


template<class Scalar>
SpmdMultiVectorBase<Scalar>::SpmdMultiVectorBase()
  :in_applyOp_(false),
   globalDim_(0),
   localOffset_(-1),
   localSubDim_(0),
   numCols_(0)
{}


// Overridden public functions from MultiVectorAdapterBase


template<class Scalar>
RCP< const ScalarProdVectorSpaceBase<Scalar> >
SpmdMultiVectorBase<Scalar>::rangeScalarProdVecSpc() const
{
  return Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(
    spmdSpace(), true
    );
}


// Protected functions overridden from MultiVectorBase


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::mvMultiReductApplyOpImpl(
  const RTOpPack::RTOpT<Scalar> &pri_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal pri_global_offset_in
  ) const
{

  using Teuchos::dyn_cast;
  using Teuchos::Workspace;

  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const Ordinal numCols = this->domain()->dim();
  const SpmdVectorSpaceBase<Scalar> &spmdSpc = *spmdSpace();

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    in_applyOp_, std::invalid_argument,
    "SpmdMultiVectorBase<>::mvMultiReductApplyOpImpl(...): Error, this method is"
    " being entered recursively which is a clear sign that one of the methods"
    " acquireDetachedView(...), releaseDetachedView(...) or commitDetachedView(...)"
    " was not implemented properly!"
    );
  apply_op_validate_input(
    "SpmdMultiVectorBase<>::mvMultiReductApplyOpImpl(...)", *this->domain(),
    *this->range(), pri_op, multi_vecs, targ_multi_vecs, reduct_objs,
    pri_global_offset_in);
#endif

  // Flag that we are in applyOp()
  in_applyOp_ = true;

  // First see if this is a locally replicated vector in which case
  // we treat this as a local operation only.
  const bool locallyReplicated = spmdSpc.isLocallyReplicated();

  const Range1D local_rng(localOffset_, localOffset_+localSubDim_-1);
  const Range1D col_rng(0, numCols-1);

  // Create sub-vector views of all of the *participating* local data
  Workspace<RTOpPack::ConstSubMultiVectorView<Scalar> >
    sub_multi_vecs(wss,multi_vecs.size());
  Workspace<RTOpPack::SubMultiVectorView<Scalar> >
    targ_sub_multi_vecs(wss,targ_multi_vecs.size());
  for(int k = 0; k < multi_vecs.size(); ++k ) {
    multi_vecs[k]->acquireDetachedView(local_rng, col_rng, &sub_multi_vecs[k]);
    sub_multi_vecs[k].setGlobalOffset(localOffset_+pri_global_offset_in);
  }
  for(int k = 0; k < targ_multi_vecs.size(); ++k ) {
    targ_multi_vecs[k]->acquireDetachedView(local_rng, col_rng, &targ_sub_multi_vecs[k]);
    targ_sub_multi_vecs[k].setGlobalOffset(localOffset_+pri_global_offset_in);
  }
  Workspace<RTOpPack::ReductTarget*> reduct_objs_ptr(wss, reduct_objs.size());
  for (int k = 0; k < reduct_objs.size(); ++k) {
    reduct_objs_ptr[k] = &*reduct_objs[k];
  }

  // Apply the RTOp operator object (all processors must participate)
  RTOpPack::SPMD_apply_op(
    locallyReplicated ? NULL : spmdSpc.getComm().get(), // comm
    pri_op, // op
    col_rng.size(), // num_cols
    multi_vecs.size(), // multi_vecs.size()
    multi_vecs.size() ? &sub_multi_vecs[0] : NULL, // sub_multi_vecs
    targ_multi_vecs.size(), // targ_multi_vecs.size()
    targ_multi_vecs.size() ? &targ_sub_multi_vecs[0] : NULL, // targ_sub_multi_vecs
    reduct_objs.size() ? &reduct_objs_ptr[0] : 0 // reduct_objs
    );

  // Free and commit the local data
  for(int k = 0; k < multi_vecs.size(); ++k ) {
    sub_multi_vecs[k].setGlobalOffset(local_rng.lbound());
    multi_vecs[k]->releaseDetachedView( &sub_multi_vecs[k] );
  }
  for(int k = 0; k < targ_multi_vecs.size(); ++k ) {
    targ_sub_multi_vecs[k].setGlobalOffset(local_rng.lbound());
    targ_multi_vecs[k]->commitDetachedView( &targ_sub_multi_vecs[k] );
  }

  // Flag that we are leaving applyOp()
  in_applyOp_ = false;

}


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::acquireDetachedMultiVectorViewImpl(
  const Range1D &rowRng_in,
  const Range1D &colRng_in,
  RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
  ) const
{
  using Teuchos::outArg;
  const Range1D rowRng = validateRowRange(rowRng_in);
  const Range1D colRng = validateColRange(colRng_in);
  if( rowRng.lbound() < localOffset_ || localOffset_+localSubDim_-1 < rowRng.ubound() ) {
    // rng consists of off-processor elements so use the default implementation!
    MultiVectorDefaultBase<Scalar>::acquireDetachedMultiVectorViewImpl(
      rowRng_in,colRng_in,sub_mv
      );
    return;
  }
  ArrayRCP<const Scalar> localValues;
  Ordinal leadingDim = 0;
  this->getLocalData(outArg(localValues), outArg(leadingDim));
  sub_mv->initialize(
    rowRng.lbound(), // globalOffset
    rowRng.size(), // subDim
    colRng.lbound(), // colOffset
    colRng.size(), // numSubCols
    localValues
    +(rowRng.lbound()-localOffset_)
    +colRng.lbound()*leadingDim, // values
    leadingDim // leadingDim
    );
}


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::releaseDetachedMultiVectorViewImpl(
  RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
  ) const
{
  if(
    sub_mv->globalOffset() < localOffset_ 
    ||
    localOffset_+localSubDim_ < sub_mv->globalOffset()+sub_mv->subDim()
    )
  {
    // Let the default implementation handle it!
    MultiVectorDefaultBase<Scalar>::releaseDetachedMultiVectorViewImpl(sub_mv);
    return;
  }
  sub_mv->uninitialize();
}


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::acquireNonconstDetachedMultiVectorViewImpl(
  const Range1D &rowRng_in,
  const Range1D &colRng_in,
  RTOpPack::SubMultiVectorView<Scalar> *sub_mv
  )
{
  using Teuchos::outArg;
  const Range1D rowRng = validateRowRange(rowRng_in);
  const Range1D colRng = validateColRange(colRng_in);
  if(
    rowRng.lbound() < localOffset_
    ||
    localOffset_+localSubDim_-1 < rowRng.ubound()
    )
  {
    // rng consists of off-processor elements so use the default implementation!
    MultiVectorDefaultBase<Scalar>::acquireNonconstDetachedMultiVectorViewImpl(
      rowRng_in, colRng_in, sub_mv
      );
    return;
  }
  ArrayRCP<Scalar> localValues;
  Ordinal leadingDim = 0;
  this->getNonconstLocalData(outArg(localValues), outArg(leadingDim));
  sub_mv->initialize(
    rowRng.lbound() // globalOffset
    ,rowRng.size() // subDim
    ,colRng.lbound() // colOffset
    ,colRng.size() // numSubCols
    ,localValues
    +(rowRng.lbound()-localOffset_)
    +colRng.lbound()*leadingDim // values
    ,leadingDim // leadingDim
    );
}


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::commitNonconstDetachedMultiVectorViewImpl(
  RTOpPack::SubMultiVectorView<Scalar>* sub_mv
  )
{
  if(
    sub_mv->globalOffset() < localOffset_
    ||
    localOffset_+localSubDim_ < sub_mv->globalOffset()+sub_mv->subDim()
    )
  {
    // Let the default implementation handle it!
    MultiVectorDefaultBase<Scalar>::commitNonconstDetachedMultiVectorViewImpl(sub_mv);
    return;
  }
  sub_mv->uninitialize();
}


// Protected functions overridden from MultiVectorAdapterBase


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::euclideanApply(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  Teuchos::Time timerTotal("dummy",true);
  Teuchos::Time timer("dummy");
#endif

  //
  // This function performs one of two operations.
  //
  // The first operation (M_trans == NOTRANS) is:
  //
  // Y = beta * Y + alpha * M * X
  //
  // where Y and M have compatible (distributed?) range vector
  // spaces and X is a locally replicated serial multi-vector. This
  // operation does not require any global communication.
  //
  // The second operation (M_trans == TRANS) is:
  //
  // Y = beta * Y + alpha * M' * X
  //
  // where M and X have compatible (distributed?) range vector spaces
  // and Y is a locally replicated serial multi-vector. This operation
  // requires a local reduction.
  //

  //
  // Get spaces and validate compatibility
  //

  // Get the SpmdVectorSpace
  const SpmdVectorSpaceBase<Scalar> &spmdSpc = *this->spmdSpace();

  // Get the Spmd communicator
  const RCP<const Teuchos::Comm<Ordinal> > comm = spmdSpc.getComm();

#ifdef TEUCHOS_DEBUG
  const VectorSpaceBase<Scalar>
    &Y_range = *Y->range(),
    &X_range = *X.range();
//	std::cout << "SpmdMultiVectorBase<Scalar>::apply(...): comm = " << comm << std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(
    ( globalDim_ > localSubDim_ ) && is_null(comm), std::logic_error
    ,"SpmdMultiVectorBase<Scalar>::apply(...MultiVectorBase<Scalar>...): Error!"
    );
  // ToDo: Write a good general validation function that I can call that will replace
  // all of these TEUCHOS_TEST_FOR_EXCEPTION(...) uses

  TEUCHOS_TEST_FOR_EXCEPTION(
    real_trans(M_trans)==NOTRANS && !spmdSpc.isCompatible(Y_range), Exceptions::IncompatibleVectorSpaces
    ,"SpmdMultiVectorBase<Scalar>::apply(...MultiVectorBase<Scalar>...): Error!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    real_trans(M_trans)==TRANS && !spmdSpc.isCompatible(X_range), Exceptions::IncompatibleVectorSpaces
    ,"SpmdMultiVectorBase<Scalar>::apply(...MultiVectorBase<Scalar>...): Error!"
    );
#endif

  //
  // Get explicit (local) views of Y, M and X
  //

#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.start();
#endif
 
  DetachedMultiVectorView<Scalar>
    Y_local(
      *Y,
      real_trans(M_trans)==NOTRANS ? Range1D(localOffset_,localOffset_+localSubDim_-1) : Range1D(),
      Range1D()
      );
  ConstDetachedMultiVectorView<Scalar>
    M_local(
      *this,
      Range1D(localOffset_,localOffset_+localSubDim_-1),
      Range1D()
      );
  ConstDetachedMultiVectorView<Scalar>
    X_local(
      X
      ,real_trans(M_trans)==NOTRANS ? Range1D() : Range1D(localOffset_,localOffset_+localSubDim_-1)
      ,Range1D()
      );
#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.stop();
  std::cout << "\nSpmdMultiVectorBase<Scalar>::apply(...): Time for getting view = " << timer.totalElapsedTime() << " seconds\n";
#endif
#ifdef TEUCHOS_DEBUG		
  TEUCHOS_TEST_FOR_EXCEPTION(
    real_trans(M_trans)==NOTRANS && ( M_local.numSubCols() != X_local.subDim() || X_local.numSubCols() != Y_local.numSubCols() )
    , Exceptions::IncompatibleVectorSpaces
    ,"SpmdMultiVectorBase<Scalar>::apply(...MultiVectorBase<Scalar>...): Error!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    real_trans(M_trans)==TRANS && ( M_local.subDim() != X_local.subDim() || X_local.numSubCols() != Y_local.numSubCols() )
    , Exceptions::IncompatibleVectorSpaces
    ,"SpmdMultiVectorBase<Scalar>::apply(...MultiVectorBase<Scalar>...): Error!"
    );
#endif

  //
  // If nonlocal (i.e. M_trans==TRANS) then create temporary storage
  // for:
  //
  // Y_local_tmp = alpha * M(local) * X(local) : on nonroot processes
  //
  // or
  //
  // Y_local_tmp = beta*Y_local + alpha * M(local) * X(local) : on root process (localOffset_==0)
  // 
  // and set
  //
  // localBeta = ( localOffset_ == 0 ? beta : 0.0 )
  //
  // Above, we choose localBeta such that we will only perform
  // Y_local = beta * Y_local + ... on one process (the root
  // process where localOffset_==0x). Then, when we add up Y_local
  // on all of the processors and we will get the correct result.
  //
  // If strictly local (i.e. M_trans == NOTRANS) then set:
  //
  // Y_local_tmp = Y_local
  // localBeta = beta
  //

#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.start();
#endif
 
  const bool isNonLocalAdjoint = 
    (
      real_trans(M_trans) == TRANS
      &&
      (globalDim_ > localSubDim_  || (nonnull(comm) && comm->getSize() > 1))
      );

  Workspace<Scalar> Y_local_tmp_store(wss, Y_local.subDim()*Y_local.numSubCols(), false);
  RTOpPack::SubMultiVectorView<Scalar> Y_local_tmp;
  Scalar localBeta;
  if (isNonLocalAdjoint) {
    // Nonlocal
    Y_local_tmp.initialize(
      0, Y_local.subDim(),
      0, Y_local.numSubCols(),
      Teuchos::arcpFromArrayView(Y_local_tmp_store()),
      Y_local.subDim() // leadingDim == subDim (columns are adjacent)
      );
    if( localOffset_ == 0 ) {
      // Root process: Must copy Y_local into Y_local_tmp
      for( int j = 0; j < Y_local.numSubCols(); ++j ) {
        Scalar *Y_local_j = Y_local.values() + Y_local.leadingDim()*j;
        std::copy( Y_local_j, Y_local_j + Y_local.subDim(), Y_local_tmp.values() + Y_local_tmp.leadingDim()*j );
      }
      localBeta = beta;
    }
    else {
      // Not the root process
      localBeta = 0.0;
    }
  }
  else {
    // Local
    Y_local_tmp = Y_local.smv(); // Shallow copy only!
    localBeta = beta;
  }

#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.stop();
  std::cout << "\nSpmdMultiVectorBase<Scalar>::apply(...): Time for setting up Y_local_tmp and localBeta = " << timer.totalElapsedTime() << " seconds\n";
#endif
 
  //
  // Perform the local multiplication:
  //
  // Y(local) = localBeta * Y(local) + alpha * op(M(local)) * X(local)
  //
  // or in BLAS lingo:
  //
  // C = beta * C + alpha * op(A) * op(B)
  //

#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.start();
#endif
  Teuchos::ETransp t_transp;
  if(ST::isComplex) {
    switch(M_trans) {
      case NOTRANS: t_transp = Teuchos::NO_TRANS; break;
      case TRANS: t_transp = Teuchos::TRANS; break;
      case CONJTRANS: t_transp = Teuchos::CONJ_TRANS; break;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  else {
    switch(real_trans(M_trans)) {
      case NOTRANS: t_transp = Teuchos::NO_TRANS; break;
      case TRANS: t_transp = Teuchos::TRANS; break;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  if (M_local.numSubCols() > 0) {
    // AGS: Added std::max on ld? below, following what is done in
    // Epetra_MultiVector Multiply use of GEMM. Allows for 0 length.
    blas_.GEMM(
      t_transp // TRANSA
      ,Teuchos::NO_TRANS // TRANSB
      ,Y_local.subDim() // M
      ,Y_local.numSubCols() // N
      ,real_trans(M_trans)==NOTRANS ? M_local.numSubCols() : M_local.subDim() // K
      ,alpha // ALPHA
      ,const_cast<Scalar*>(M_local.values()) // A
      ,std::max((int) M_local.leadingDim(),1) // LDA
      ,const_cast<Scalar*>(X_local.values()) // B
      ,std::max((int) X_local.leadingDim(),1) // LDB
      ,localBeta // BETA
      ,Y_local_tmp.values().get() // C
      ,std::max((int) Y_local_tmp.leadingDim(),1) // LDC
      );
  }
  else {
    std::fill( Y_local_tmp.values().begin(), Y_local_tmp.values().end(),
      ST::zero() );
  }
#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.stop();
  std::cout
    << "\nSpmdMultiVectorBase<Scalar>::apply(...): Time for GEMM = "
    << timer.totalElapsedTime() << " seconds\n";
#endif

  if (nonnull(comm)) {
 
    //
    // Perform the global reduction of Y_local_tmp back into Y_local
    //
 
    if (isNonLocalAdjoint) {
      // Contiguous buffer for final reduction
      Workspace<Scalar> Y_local_final_buff(wss,Y_local.subDim()*Y_local.numSubCols(),false);
      // Perform the reduction
      Teuchos::reduceAll<Ordinal,Scalar>(
        *comm,Teuchos::REDUCE_SUM,Y_local_final_buff.size(),Y_local_tmp.values().get(),
        &Y_local_final_buff[0]
        );
      // Load Y_local_final_buff back into Y_local
      const Scalar *Y_local_final_buff_ptr = &Y_local_final_buff[0];
      for( int j = 0; j < Y_local.numSubCols(); ++j ) {
        Scalar *Y_local_ptr = Y_local.values() + Y_local.leadingDim()*j;
        for( int i = 0; i < Y_local.subDim(); ++i ) {
          (*Y_local_ptr++) = (*Y_local_final_buff_ptr++);
        }
      }
    }
  }
  else {

    // When you get here the view Y_local will be committed back to Y
    // in the destructor to Y_local

  }

#ifdef THYRA_SPMD_MULTI_VECTOR_BASE_PRINT_TIMES
  timer.stop();
  std::cout 
    << "\nSpmdMultiVectorBase<Scalar>::apply(...): Total time = "
    << timerTotal.totalElapsedTime() << " seconds\n";
#endif

}


// Protected functions for subclasses to call


template<class Scalar>
void SpmdMultiVectorBase<Scalar>::updateSpmdSpace()
{
  if(globalDim_ == 0) {
    const SpmdVectorSpaceBase<Scalar> *l_spmdSpace = this->spmdSpace().get();
    if(l_spmdSpace) {
      globalDim_ = l_spmdSpace->dim();
      localOffset_ = l_spmdSpace->localOffset();
      localSubDim_ = l_spmdSpace->localSubDim();
      numCols_ = this->domain()->dim();
    }
    else {
      globalDim_ = 0;
      localOffset_ = -1;
      localSubDim_ = 0;
      numCols_ = 0;
    }
  }
}


template<class Scalar>
Range1D SpmdMultiVectorBase<Scalar>::validateRowRange( const Range1D &rowRng_in ) const
{
  const Range1D rowRng = Teuchos::full_range(rowRng_in,0,globalDim_-1);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= rowRng.lbound() && rowRng.ubound() < globalDim_ ), std::invalid_argument
    ,"SpmdMultiVectorBase<Scalar>::validateRowRange(rowRng): Error, the range rowRng = ["
    <<rowRng.lbound()<<","<<rowRng.ubound()<<"] is not "
    "in the range [0,"<<(globalDim_-1)<<"]!"
    );
#endif
  return rowRng;
}


template<class Scalar>
Range1D SpmdMultiVectorBase<Scalar>::validateColRange( const Range1D &colRng_in ) const
{
  const Range1D colRng = Teuchos::full_range(colRng_in,0,numCols_-1);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= colRng.lbound() && colRng.ubound() < numCols_), std::invalid_argument
    ,"SpmdMultiVectorBase<Scalar>::validateColRange(colRng): Error, the range colRng = ["
    <<colRng.lbound()<<","<<colRng.ubound()<<"] is not "
    "in the range [0,"<<(numCols_-1)<<"]!"
    );
#endif
  return colRng;
}


} // end namespace Thyra


#endif // THYRA_SPMD_MULTI_VECTOR_BASE_DEF_HPP
