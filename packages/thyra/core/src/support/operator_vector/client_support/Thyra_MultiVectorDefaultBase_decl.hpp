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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpDefaultBase_decl.hpp"


namespace Thyra {


/** \brief Node subclass that uses a default <tt>MultiVectorBase</tt>
 * implementation</tt> to provide default implementations for as many other
 * functions in <tt>MultiVectorBase</tt> interface the as is reasonable.
 *
 * <b>Notes to subclass developers</b>
 *
 * Only three function overrides are required in order to create a concrete
 * <tt>MultiVectorBase</tt> subclass: <tt>range()</tt>, <tt>domain()</tt> and
 * the non-const version of <tt>col()</tt>.  All of the other functions have
 * default implementations.  However, a good implementation will provide
 * optimized overrides of at least the functions <tt>apply()</tt> and
 * <tt>applyTranspose()</tt>.  The non-const versions of <tt>subView()</tt>
 * should be overridden if subviews are important.  The default implementation
 * will not achieve near-optimal performance in many cases.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class MultiVectorDefaultBase
  : virtual public MultiVectorBase<Scalar>
  , virtual protected LinearOpDefaultBase<Scalar>
{
public:

  /** \name Overridden public member functions from MultiVectorBase */
  //@{

  /** \brief .
   *
   * This implementation uses the vector space to create a new multi-vector
   * object and then uses a transformation operator to assign the vector
   * elements.  A subclass should only override this function if it can do
   * something more sophisticated (i.e. lazy evaluation) but in general, this
   * is not needed.
   */
  virtual RCP<MultiVectorBase<Scalar> > clone_mv() const;

  //@}

protected:

  /** \name Overridden protected member functions from MultiVectorBase */
  //@{

  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  contigSubViewImpl( const Range1D& colRng ) const;

  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl( const Range1D& colRng );

  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  nonContigSubViewImpl( const ArrayView<const int> &cols ) const;

  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstNonContigSubViewImpl( const ArrayView<const int> &cols );

  /** \brief .
   *
   * This implementation calls <tt>VectorBase::applyOp()</tt> on each column
   * <tt>this->col(j)</tt> for <tt>j = 0 ... this->range()->dim()-1</tt>.
   */
  virtual void mvMultiReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
    const Ordinal primary_global_offset
    ) const;

  /** \brief .
   *
   * This implementation calls <tt>applyOp()</tt> where an array of reduction
   * objects is taken.
   */
  virtual void mvSingleReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const RTOpPack::RTOpT<Scalar> &secondary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal primary_global_offset
    ) const;

  /** \brief .
   *
   * This implementation is based on the vector operation
   * <tt>VectorBase::acquireDetachedView()</tt> called on the non-changeable vector
   * objects returned from <tt>col()</tt>.  Note that the footprint of the
   * reduction object (both internal and external state) will be
   * O(<tt>rowRng.size()*colRng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, however, this is a terrible implementation and
   * must be overridden if <tt>rowRng.size()</tt> is large at all.  Although,
   * this function should not even be used in cases where the multi-vector is
   * very large.  If a subclass does override this function, it must also
   * override <tt>releaseDetachedView()</tt> which has an implementation
   * which is a companion to this function's implementation.
   */
  virtual void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
    ) const;

  /** \brief .
   *
   * This implementation is a companion to the implementation for
   * <tt>acquireDetachedView()</tt>.  If <tt>acquireDetachedView()</tt> is
   * overridden by a subclass then this function must be overridden also!
   */
  virtual void releaseDetachedMultiVectorViewImpl(
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const;

  /** \brief .
   *
   * This implementation is based on the vector operation
   * <tt>VectorBase::acquireDetachedView()</tt> called on the changeable vector
   * objects returned from <tt>col()</tt>.  Note that the footprint of the
   * reduction object (both internal and external state) will be
   * O(<tt>rowRng.size()*colRng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, however, this is a terrible implementation and
   * must be overridden if <tt>rowRng.size()</tt> is large at all.  Although,
   * this function should not even be used in case where the multi-vector is
   * very large.  If a subclass does override this function, it must also
   * override <tt>commitDetachedView()</tt> which has an implementation
   * which is a companion to this function's implementation.
   */
  virtual void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar> *sub_mv
    );

  /** \brief .
   *
   * This implementation is a companion to the default implementation for
   * <tt>acquireDetachedView()</tt>.  If <tt>acquireDetachedView()</tt> is
   * overridden by a subclass then this function must be overridden also!
   */
  virtual void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    );

  //@}

};


} // namespace Thyra


#define THYRA_ASSERT_MV_COLS(FUNCNAME, cols) \
  { \
    const int numCols = cols.size(); \
    const Thyra::Ordinal dimDomain = this->domain()->dim(); \
    const std::string msgErr = this->description()+"::"+FUNCNAME; \
    TEST_FOR_EXCEPTION( !( 1 <= numCols && numCols <= dimDomain ), \
      std::invalid_argument, msgErr<<"Error!"); \
    for (int k = 0; k < numCols; ++k) { \
      const int col_k = cols[k]; \
      TEST_FOR_EXCEPTION( \
        !( 0<= col_k && col_k < dimDomain ), std::out_of_range, \
        msgErr<<": col["<<k<<"] = " << col_k \
        << " is not in the range [0,"<<(dimDomain-1)<<"]!" \
        ); \
    } \
  } \
  typedef int THYRA_ASSERT_MV_COLS_t 

#ifdef THYRA_DEBUG
#  define THYRA_DEBUG_ASSERT_MV_COLS(FUNCNAME, cols) \
  THYRA_ASSERT_MV_COLS(FUNCNAME, cols)
#else
#  define THYRA_DEBUG_ASSERT_MV_COLS(FUNCNAME, cols)
#endif

#endif // THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
