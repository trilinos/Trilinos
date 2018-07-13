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

#ifndef THYRA_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorDefaultBase_decl.hpp"


namespace Thyra {


/** \brief Convenient node subclass for concrete <tt>VectorBase</tt>
 * subclasses that relies on a default <tt>MultiVectorBase</tt>
 * implementation.
 *
 * <b>Notes for subclass developers</b>
 *
 * In order to create a concrete subclass of this interface, only two
 * operations must be overridden: <tt>space()</tt> and <tt>applyOp()</tt>.
 * Overriding the <tt>space()</tt> operation requires defining a concrete
 * <tt>VectorSpaceBase</tt> class (which has only three pure virtual
 * operations if using <tt>VectorSpaceDefaultBase</tt>).
 *
 * Note that all of the inherited <tt>LinearOpBase</tt> and
 * <tt>MultiVectorBase</tt> functions are overridden in this subclass and are
 * given perfectly good implementations.  Therefore, a concrete subclass of
 * <tt>%VectorDefaultBase</tt> should not have to re-override any of these
 * functions.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class VectorDefaultBase
  : virtual public VectorBase<Scalar>,
    virtual protected MultiVectorDefaultBase<Scalar>
{
public:

  /** @name Public functions overridden from Teuchos::Describable */
  //@{

  /** \brief Default description that gives the label, type, and dimenstion . */
  virtual std::string description() const;

  /** \brief Generates a default outputting for all vectors.
   *
   * Calls on the <tt>this->describe(void)</tt> function for the name of the
   * class (and possibly its instance name) and then if
   * <tt>verbLevel>=VERB_HIGH</tt>, then the vector elements themselves are
   * printed as well.  The format of the output is is shown below:
   
   \verbatim

   type = 'this->description()', size = n
     0:x1
     1:x2
     .
     .
     .
     n-1:xn
   \endverbatim
   *
   * The <tt>index:value</tt> lines are offset an additional
   * <tt>Teuchos::OSTab</tt> amount.  A newline is printed after the last
   * <tt>n-1:xn</tt> entry.
   */
  virtual void describe(
    Teuchos::FancyOStream                &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Returns <tt>this->space()</tt>. */
  virtual RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief Returns a <tt>DefaultSerialVectorSpace</tt> object with dimension
   * 1.
   */
  virtual RCP< const VectorSpaceBase<Scalar> > domain() const;

  //@}

  /** @name Overridden from MultiVectorBase */
  //@{

  /** \brief Returns <tt>this->clone_v()</tt>. */
  RCP<MultiVectorBase<Scalar> > clone_mv() const;

  //@}

  /** \name Overridden from VectorBase */
  //@{

  /** \brief Simply creates a new vector and copies the contents from
   * <tt>*this</tt>.
   */
  RCP<VectorBase<Scalar> > clone_v() const;

  //@}

protected:

  /** @name Overridden protected functions from MultiVectorBase */
  //@{

  /** \brief Default implementation of assign(vector) using RTOps. */
  virtual void assignVecImpl(const VectorBase<Scalar>& x);

  /** \brief Default implementation of randomize using RTOps. */
  virtual void randomizeImpl(Scalar l, Scalar u);

  /** \brief Default implementation of abs using RTOps. */
  virtual void absImpl(const VectorBase<Scalar>& x);

  /** \brief Default implementation of reciprocal using RTOps. */
  virtual void reciprocalImpl(const VectorBase<Scalar>& x);

  /** \brief Default implementation of ele_wise_scale using RTOps. */
  virtual void eleWiseScaleImpl(const VectorBase<Scalar>& x);

  /** \brief Default implementation of update using RTOps. */
  virtual void updateVecImpl(
    Scalar alpha,
    const VectorBase<Scalar>& x);

  /** \brief Default implementation of linear_combination using RTOps. */
  virtual void linearCombinationVecImpl(
    const ArrayView<const Scalar>& alpha,
    const ArrayView<const Ptr<const VectorBase<Scalar> > >& x,
    const Scalar& beta
    );

  /** \brief Default implementation of dot using RTOps. */
  virtual Scalar dotImpl(const VectorBase<Scalar>& x) const;

  /** \brief Default implementation of norm_1 using RTOps. */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  norm1Impl() const;

  /** \brief Default implementation of norm_2 using RTOps. */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  norm2Impl() const;

  /** \brief Default implementation of norm_2 (weighted) using RTOps. */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  norm2WeightedImpl(const VectorBase<Scalar>& x) const;

  /** \brief Default implementation of norm_inf using RTOps. */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  normInfImpl() const;

  /** \brief Returns <tt>Teuchos::rcp(this,false)</tt>. */
  virtual RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j);
  /** \brief Returns <tt>Teuchos::rcp(this,false)</tt>. */
  virtual RCP<const MultiVectorBase<Scalar> >
  contigSubViewImpl( const Range1D& col_rng ) const;
  /** \brief Returns <tt>Teuchos::rcp(this,false)</tt>. */
  virtual RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl( const Range1D& col_rng );
  /** \brief Returns <tt>Teuchos::rcp(this,false)</tt>. */
  virtual RCP<const MultiVectorBase<Scalar> >
  nonContigSubViewImpl( const ArrayView<const int> &cols ) const;
  /** \brief Returns <tt>Teuchos::rcp(this,false)</tt>. */
  virtual RCP<MultiVectorBase<Scalar> >
  nonconstNonContigSubViewImpl( const ArrayView<const int> &cols );
  /** \brief Implemented in terms of <tt>this->acquireDetachedView()</tt>. */
  virtual void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
    ) const;
  /** \brief Implemented in terms of <tt>this->releaseDetachedView()</tt>. */
  virtual void releaseDetachedMultiVectorViewImpl(
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const;
  /** \brief Implemented in terms of <tt>this->acquireDetachedView()</tt>. */
  virtual void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar> *sub_mv
    );
  /** \brief Implemented in terms of <tt>this->commitDetachedView()</tt>. */
  virtual void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    );

  //@}

  /** @name Overridden protected functions from VectorBase */
  //@{

  /** \brief .
   *
   * This implementation is based on a vector reduction operator class (see
   * <tt>RTOpPack::ROpGetSubVector</tt>) and calls <tt>applyOp()</tt>.  Note
   * that the footprint of the reduction object (both internal and external
   * state) will be O(<tt>rng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, however, this is a terrible implementation and
   * must be overridden if <tt>rng.size()</tt> is large at all.  Although,
   * this function should not even be used in case where the vector is very
   * large.  If a subclass does override this function, it must also override
   * <tt>releaseDetachedView()</tt> which has a implementation which is a companion
   * to this function's implementation.
   */
  virtual void acquireDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief .
   *
   * This implementation is a companion to the implementation for the
   * non-<tt>const</tt> version of <tt>acquireDetachedView()</tt>.  If
   * <tt>acquireDetachedView()</tt> is overridden by a subclass then this function
   * must be overridden also!
   */
  virtual void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief .
   *
   * This implementation is based on a vector reduction operator class (see
   * <tt>RTOpPack::ROpGetSubVector</tt>) and calls <tt>applyOp()</tt>.  Note
   * that the footprint of the reduction object (both internal and external
   * state) will be O(<tt>rng.size()</tt>).  For serial applications this is
   * fairly reasonable and will not be a major performance penalty.  For
   * parallel applications, this will be a terrible thing to do and must be
   * overridden if <tt>rng.size()</tt> is large at all.  If a subclass does
   * override this function, it must also override <tt>commitDetachedView()</tt>
   * which has a implementation which is a companion to this function's
   * implementation.
   */
  virtual void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief .
   *
   * This function has an implementation which is a companion to the
   * implementation for <tt>acquireDetachedView()</tt>.  If <tt>acquireDetachedView()</tt>
   * is overridden by a subclass then this function must be overridden also!
   */
  virtual void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief .
   *
   * This implementation uses a transformation operator class (see
   * <tt>RTOpPack::TOpSetSubVector</tt>) and calls <tt>applyOp()</tt>.  Be
   * forewarned however, that the operator objects state data (both internal
   * and external) will be order O(<tt>sub_vec.subNz()</tt>).  For serial
   * applications, this is entirely adequate.  For parallel applications this
   * may be bad!
   */
  virtual void setSubVectorImpl( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec );

  //@}

  /** @name Overridden protected functions from LinearOpBase */
  //@{

  /** \brief For complex <tt>Scalar</tt> types returns <tt>true</tt> for
   * <tt>NOTRANS</tt> and <tt>CONJTRANS</tt> and for real types returns true
   * for all values of <tt>M_trans</tt>.
   */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief. Applies vector or its adjoint (transpose) as a linear
   * operator.
   */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private data members

  mutable RCP<const VectorSpaceBase<Scalar> > domain_;
  // Above only initialized if *this is used as a MultiVectorBase

  // /////////////////////////////////////
  // Private member functions

  void validateColRng( const Range1D &rowRng ) const;
  void validateColIndexes( const ArrayView<const int> &cols ) const;

};


} // end namespace Thyra


#endif  // THYRA_VECTOR_DEFAULT_BASE_DECL_HPP
