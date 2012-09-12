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

#ifndef THYRA_SPMD_MULTI_VECTOR_BASE_DECL_HPP
#define THYRA_SPMD_MULTI_VECTOR_BASE_DECL_HPP

#include "Thyra_MultiVectorAdapterBase_decl.hpp"
#include "Teuchos_BLAS.hpp"


namespace Thyra {


/** \brief . */
template<class Scalar> class SpmdVectorSpaceBase;


/** \brief Base class for SPMD multi-vectors.
 *
 * By inheriting from this base class, multi-vector implementations allow
 * their multi-vector objects to be seamlessly combined with other SPMD
 * multi-vector objects (of different concrete types) in <tt>applyOp()</tt>
 * and <tt>apply()</tt>.  A big part of this protocol is that every
 * multi-vector object can expose an <tt>SpmdVectorSpaceBase</tt> object
 * through the virtual function <tt>spmdSpace()</tt>.
 *
 * This base class contains an implementation of <tt>applyOp()</tt> that
 * relies on implementations of the <tt>const</tt> functions
 * <tt>acquireDetachedView()</tt> and <tt>releaseDetachedView()</tt>, and the
 * non-<tt>const</tt> functions <tt>acquireDetachedView()</tt> and
 * <tt>commitDetachedView()</tt> (which all have default implementations in
 * this subclass).  In essence, this implementation will only call the
 * <tt>acquireDetachedView()</tt> functions using a range of (global) indexes
 * for elements that exist in the local process.  As long as the number of
 * local elements in each process is fairly large, the virtual function call
 * overhead will be minimal and this will result in a near optimal
 * implementation.
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
class SpmdMultiVectorBase
  : virtual public MultiVectorAdapterBase<Scalar>
{
public:

  /** @name  Constructors / initializers / accessors */
  //@{

  /** \brief . */
  SpmdMultiVectorBase();

  //@}

  /** @name Pure virtual functions to be overridden by subclasses */
  //@{

  /** \brief Returns the SPMD vector space object for the range of
   * <tt>*this</tt> multi-vector.
   */
  virtual RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const = 0;

  /** \brief Returns a non-<tt>const</tt> pointer to a Fortran-style view of
   * the local multi-vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to
   * the first element in the first column of the local multi-vector stored as
   * a column-major dense Fortran-style matrix.
   *
   * \param leadingDim [out] On output <tt>*leadingDim</tt> gives the leading
   * dimension of the Fortran-style local multi-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>leadingDim!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*leadingDim!=0</tt>
   * </ul>
   *
   * The function <tT>commitLocalData()</tt> must be called to
   * commit changes to the data.
   */
  void getNonconstLocalData(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    )
    {
      getNonconstLocalDataImpl(localValues, leadingDim);
    }

  /** \brief Returns a <tt>const</tt> pointer to a Fortran-style view of the
   * local multi-vector data.
   *
   * \param localValues [out] On output <tt>*localValues</tt> will point to
   * the first element in the first column of the local multi-vector stored as
   * a column-major dense Fortran-style matrix.
   *
   * \param leadingDim [out] On output <tt>*leadingDim</tt> gives the leading
   * dimension of the Fortran-style local multi-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>leadingDim!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*leadingDim!=0</tt>
   * </ul>
   */
  void getLocalData(
    const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const
    {
      getLocalDataImpl(localValues, leadingDim);
    }

  //@}

  /** @name Overridden public functions from MultiVectorAdapterBase */
  //@{

  /** \brief Returns <tt>spmdSpace</tt>. */
  RCP< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;

  //@}

#ifndef THYRA_HIDE_DEPRECATED_CODE
  /** \name Deprecated. */
  //@{

  /** \brief Deprecated (use ArrayRCP version of getLocalData(...). */
  THYRA_DEPRECATED void getLocalData( Scalar **localValues_out, Ordinal *leadingDim_out )
    {
      using Teuchos::outArg;
      ArrayRCP<Scalar> localValues;
      this->getNonconstLocalData(outArg(localValues), outArg(*leadingDim_out));
      *localValues_out = localValues.getRawPtr();
    }
  
  /** \brief Deprecated. */
  THYRA_DEPRECATED void commitLocalData( Scalar *localValues )
    {
      // Do nothing!
    }

  /** \brief Deprecated. */
  THYRA_DEPRECATED void getLocalData(
    const Scalar **localValues_out, Ordinal *leadingDim_out
    ) const
    {
      using Teuchos::outArg;
      ArrayRCP<const Scalar> localValues;
      this->getLocalData(outArg(localValues), outArg(*leadingDim_out));
      *localValues_out = localValues.getRawPtr();

    }

  /** \brief Deprecated. */
  THYRA_DEPRECATED void freeLocalData( const Scalar *localValues ) const
    {
      // Do nothing!
    }

  //@}

#endif // THYRA_HIDE_DEPRECATED_CODE
protected:

  /** @name Virtual functions to be overridden by sublcasses. */
  //@{

  /** \brief Virtual implementation for getNonconstLocalData(). */
  virtual void getNonconstLocalDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) = 0;

  /** \brief Virtual implementation for getLocalData(). */
  virtual void getLocalDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const = 0;

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
   * <tt>SpmdMultiVectorBase<Scalar>::updateSpmdSpace();</tt>).
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
  
}; // end class SpmdMultiVectorBase

} // end namespace Thyra

#endif // THYRA_SPMD_MULTI_VECTOR_BASE_DECL_HPP
