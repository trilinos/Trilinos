// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_Spmd_MULTI_VECTOR_STD_DECL_HPP
#define THYRA_Spmd_MULTI_VECTOR_STD_DECL_HPP

#include "Thyra_SpmdMultiVectorDefaultBase_decl.hpp"

namespace Thyra {

/** \brief Efficient concrete implementation subclass for SPMD multi-vectors.
 *
 * This subclass provides a very efficient and very general concrete
 * implementation of a <tt>Thyra::MultiVectorBase</tt> object.
 *
 * Objects of this type generally should not be constructed directly by a
 * client but instead by using the concrete vector space subclass
 * <tt>Thyra::DefaultSpmdVectorSpace</tt> and using the function
 * <tt>Thyra::createMembers()</tt>.
 *
 * The storage type can be anything since a <tt>Teuchos::ArrayRCP</tt> is used
 * to pass in the local values pointer into the constructor and
 * <tt>initialize()</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template <class Scalar>
class DefaultSpmdMultiVector : virtual public SpmdMultiVectorDefaultBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct to uninitialized
  DefaultSpmdMultiVector();

  /// Calls <tt>initialize()</tt>
  DefaultSpmdMultiVector(
      const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace,
      const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace);

  /// Calls <tt>initialize()</tt>
  DefaultSpmdMultiVector(
      const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace,
      const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
      const ArrayRCP<Scalar> &localValues,
      const Ordinal leadingDim = -1);

  /** \brief Initialize only with vector spaces where storage is allocated
   * internally..
   *
   * \param spmdRangeSpace [in] Smart pointer to <tt>SpmdVectorSpaceBase</tt>
   * object that defines the data distribution for <tt>spmdSpace()</tt> and
   * <tt>range()</tt>.
   *
   * \param domainSpace [in] Smart pointer to <tt>VectorSpaceBase</tt> object
   * that defines <tt>domain()</tt> space.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>spmdRangeSpace.get()!=NULL</tt>
   * <li><tt>domainSpace.get()!=NULL</tt>
   * </ul>
   */
  void initialize(
      const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace,
      const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace);

  /** \brief Initialize using externally allocated storage.
   *
   * \param spmdRangeSpace [in] Smart pointer to <tt>SpmdVectorSpaceBase</tt>
   * object that defines the data distribution for <tt>spmdSpace()</tt> and
   * <tt>range()</tt>.
   *
   * \param domainSpace [in] Smart pointer to <tt>VectorSpaceBase</tt> object
   * that defines <tt>domain()</tt> space.
   *
   * \param localValues [in] Smart pointer to beginning of Fortran-style
   * column-major array that defines the local localValues in the
   * multi-vector.  This array must be at least of dimension
   * <tt>spmdRangeSpace->localSubDim()*domainSpace->dim()</tt> and
   * <tt>(&*localValues)[ i + j*leadingDim ]</tt> gives the local value of the
   * zero-based entry <tt>(i,j)</tt> where
   * <tt>i=0...spmdSpace()->localSubDim()-1</tt> and
   * <tt>j=0...domainSpace->dim()-1</tt>.
   *
   * \param leadingDim [in] The leading dimension of the multi-vector.  If -1,
   * then is taken to be spmdRangeSpace->localSubDim().
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>spmdRangeSpace.get()!=NULL</tt>
   * <li><tt>domainSpace.get()!=NULL</tt>
   * <li><tt>localValues.get()!=NULL</tt>
   * <li><tt>leadingDim >= spmdRangeSpace->localSubDim()</tt>
   * </ul>
   */
  void initialize(
      const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace,
      const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
      const ArrayRCP<Scalar> &localValues,
      const Ordinal leadingDim = -1);

  /** \brief Set to an uninitialized state.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->spmdSpace().get() == NULL</tt>.
   */
  void uninitialize(
      RCP<const SpmdVectorSpaceBase<Scalar> > *spmdRangeSpace    = NULL,
      RCP<const ScalarProdVectorSpaceBase<Scalar> > *domainSpace = NULL,
      ArrayRCP<Scalar> *localValues                              = NULL,
      Ordinal *leadingDim                                        = NULL);

  /** \brief . */
  RCP<const ScalarProdVectorSpaceBase<Scalar> >
  domainScalarProdVecSpc() const;

  //@}

 protected:
  /** @name Overridden protected functions from MultiVectorBase */
  //@{
  /** \brief . */
  RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j);
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  contigSubViewImpl(const Range1D &colRng) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl(const Range1D &colRng);
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  nonContigSubViewImpl(const ArrayView<const int> &cols) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstNonContigSubViewImpl(const ArrayView<const int> &cols);
  //@}

  /** @name Overridden protected functions from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpaceImpl() const;
  /** \brief . */
  void getNonconstLocalMultiVectorDataImpl(
      const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim);
  /** \brief . */
  void getLocalMultiVectorDataImpl(
      const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim) const;
  //@}

 private:
  // ///////////////////////////////////////
  // Private data members

  RCP<const SpmdVectorSpaceBase<Scalar> > spmdRangeSpace_;
  RCP<const ScalarProdVectorSpaceBase<Scalar> > domainSpace_;
  ArrayRCP<Scalar> localValues_;
  Ordinal leadingDim_;

  // ///////////////////////////////////////
  // Private member functions

  ArrayRCP<Scalar> createContiguousCopy(const ArrayView<const int> &cols) const;

 public:
#ifdef THYRA_DEBUG
  // Unit testing sensing varaible
  static int numSkipCopyBack;
#endif

};  // end class DefaultSpmdMultiVector

template <class Scalar>
RCP<DefaultSpmdMultiVector<Scalar> >
defaultSpmdMultiVector(
    const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdRangeSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const ArrayRCP<Scalar> &localValues,
    const Ordinal leadingDim = -1) {
  return Teuchos::rcp(
      new DefaultSpmdMultiVector<Scalar>(
          spmdRangeSpace, domainSpace, localValues, leadingDim));
}

}  // end namespace Thyra

#endif  // THYRA_Spmd_MULTI_VECTOR_STD_DECL_HPP
