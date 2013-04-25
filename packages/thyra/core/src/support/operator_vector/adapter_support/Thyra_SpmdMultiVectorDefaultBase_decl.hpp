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

#ifndef THYRA_SPMD_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_SPMD_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_MultiVectorAdapterBase_decl.hpp"
#include "Teuchos_BLAS.hpp"


namespace Thyra {


/** \brief Base node implementation class for SPMD multi-vectors.
 *
 * By inheriting from this base class, multi-vector implementations allow
 * their multi-vector objects to be seamlessly combined with other SPMD
 * multi-vector objects (of different concrete types) in <tt>applyOp()</tt>
 * and <tt>apply()</tt>.  A big part of this protocol is that every
 * multi-vector object can expose an <tt>SpmdVectorSpaceBase</tt> object
 * through the virtual function <tt>spmdSpace()</tt>.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must override only five functions:
 * <tt>spmdSpace()</tt>, <tt>getLocalData(const Scalar**,Ordinal*)</tt>,
 * <tt>freeLocalData(const Scalar**,Ordinal*)</tt>,
 * <tt>getLocalData(Scalar**,Ordinal*)</tt>, and
 * <tt>commitLocalData(Scalar**,Ordinal*)</tt>.  Note that overriding the
 * <tt>spmdSpace()</tt> function requires implementing or using a
 * pre-implemented concrete <tt>SpmdVectorSpaceBase</tt> subclass.
 *
 * If the <tt>acquireDetachedView()</tt> functions are ever called with index
 * ranges outside of those of the local process, then the default
 * implementations in <tt>MultiVectorBase</tt> of all of the functions
 * (<tt>const</tt>) <tt>MultiVectorBase::acquireDetachedView()</tt>,
 * <tt>MultiVectorBase::releaseDetachedView()</tt>, (non-<tt>const</tt>)
 * <tt>MultiVectorBase::acquireDetachedView()</tt> and
 * <tt>MultiVectorBase::commitDetachedView()</tt> are called in instead.
 * Alternatively, a subclass could provide more specialized implementations of
 * these functions (for more efficient gather/scatter operations) if desired
 * but this should not be needed for most use cases.
 *
 * It is interesting to note that in the above use case that the
 * explicit subvector access functions call on its default
 * implementation defined in <tt>MultiVectorBase</tt> (which calls on
 * <tt>applyOp()</tt>) and the operator will be properly applied since
 * the version of <tt>applyOp()</tt> implemented in this class will
 * only request local vector data and hence there will only be two
 * levels of recursion for any call to an explicit subvector access
 * function.  This is a truly elegant result.
 *
 * Note that a multi-vector subclass derived from this base class must only be
 * directly used in SPMD mode for this to work properly.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_support_grp
 */
template<class Scalar>
class SpmdMultiVectorDefaultBase
  : virtual public SpmdMultiVectorBase<Scalar>,
    virtual public MultiVectorAdapterBase<Scalar>
{
public:

  /** @name  Constructors / initializers / accessors */
  //@{

  /** \brief . */
  SpmdMultiVectorDefaultBase();

  //@}

  /** @name Overridden public functions from MultiVectorAdapterBase */
  //@{

  /** \brief Returns <tt>spmdSpace</tt>. */
  RCP< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;

  //@}

protected:

  /** @name Protected funtions overridden from SpmdMultiVectorBase. */
  //@{

  /** \brief . */
  RTOpPack::SubMultiVectorView<Scalar> getNonconstLocalSubMultiVectorImpl();

  /** \brief . */
  RTOpPack::ConstSubMultiVectorView<Scalar> getLocalSubMultiVectorImpl() const;

  //@}

  /** @name Protected functions overridden from MultiVectorBase */
  //@{
  /** \brief . */
  void mvMultiReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
    const Ordinal primary_global_offset
    ) const;
  /** \brief . */
  void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng
    ,RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
    ) const;
  /** \brief . */
  void releaseDetachedMultiVectorViewImpl(
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const;
  /** \brief . */
  void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar> *sub_mv
    );
  /** \brief . */
  void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    );
  //@}

  /** @name Protected functions overridden from MultiVectorAdapterBase */
  //@{

  /** \brief Uses <tt>GEMM()</tt> and <tt>Teuchos::reduceAll()</tt> to
   * implement.
   *
   * ToDo: Finish documentation!
   */
  void euclideanApply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

  /** @name Protected functions for subclasses to call. */
  //@{

  /** \brief Subclasses should call whenever the structure of the
   * VectorSpaceBase changes.
   *
   * <b>WARNING!</b> This function can be overridden by subclasses but this
   * particular function implementation must be called back from within any
   * override (i.e. call
   * <tt>SpmdMultiVectorDefaultBase<Scalar>::updateSpmdSpace();</tt>).
   */
  virtual void updateSpmdSpace();

  /** \brief Validate and resize the row range.
   *
   * This function throws an exception if the input range is invalid
   */
  Range1D validateRowRange( const Range1D& rowRng ) const;

  /** \brief Validate and resize the column range.
   *
   * This function throws an exception if the input range is invalid
   */
  Range1D validateColRange( const Range1D& rowCol ) const;

  //@}
  
private:
  
  // ///////////////////////////////////////
  // Private data members
  
  mutable bool in_applyOp_;

  mutable Teuchos::BLAS<int,Scalar> blas_;

  // cached
  Ordinal  globalDim_;
  Ordinal  localOffset_;
  Ordinal  localSubDim_;
  Ordinal  numCols_;
  
}; // end class SpmdMultiVectorDefaultBase


} // end namespace Thyra


#endif // THYRA_SPMD_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
