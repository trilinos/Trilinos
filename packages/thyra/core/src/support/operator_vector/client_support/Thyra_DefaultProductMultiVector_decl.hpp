// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DECL_HPP
#define THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DECL_HPP

#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief . */
template <class Scalar>
class DefaultProductVectorSpace;

/** \brief Concrete implementation of a product multi-vector.
 *
 * Note that clients should almost never be creating objects of this
 * type explicitly and should instead use <tt>DefaultProductVectorSpace</tt>
 * as a factory.
 *
 * ToDo: Finish documentation!
 *
 * The default constructor is made private to avoid accidental default
 * construction.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template <class Scalar>
class DefaultProductMultiVector
  : virtual public ProductMultiVectorBase<Scalar>,
    virtual protected MultiVectorDefaultBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  DefaultProductMultiVector();

  /** \brief . */
  void initialize(
      const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
      const int numMembers);

  /** \brief . */
  void initialize(
      const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
      const ArrayView<const RCP<MultiVectorBase<Scalar> > > &multiVecs);

  /** \brief . */
  void initialize(
      const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
      const ArrayView<const RCP<const MultiVectorBase<Scalar> > > &multiVecs);

  /** \brief Uninitialize.
   *
   * ToDo: Finish documentation.
   */
  void uninitialize();

  //@}

  /** @name Overridden public functions from Teuchos::Describable */
  //@{

  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
      Teuchos::FancyOStream &out,
      const Teuchos::EVerbosityLevel verbLevel) const;

  //@}

  /** \name Overridden public functions from ProductMultiVectorBase */
  //@{

  /** \brief . */
  RCP<const ProductVectorSpaceBase<Scalar> >
  productSpace() const;
  /** \brief . */
  bool blockIsConst(const int k) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  getNonconstMultiVectorBlock(const int k);
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  getMultiVectorBlock(const int k) const;

  //@}

  /** \name Overriden public functions from MultiVectorBase */
  //@{
  /** \brief . */
  RCP<MultiVectorBase<Scalar> > clone_mv() const;

  //@}

  /** \name Overriden from LinearOpBase */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> >
  range() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> >
  domain() const;

  //@}

 protected:
  /** \name Overriden protected functions from MultiVectorBase */
  //@{

  /** \brief . */
  virtual void assignImpl(Scalar alpha);
  /** \brief . */
  virtual void assignMultiVecImpl(const MultiVectorBase<Scalar> &mv);
  /** \brief . */
  virtual void scaleImpl(Scalar alpha);
  /** \brief . */
  virtual void updateImpl(
      Scalar alpha,
      const MultiVectorBase<Scalar> &mv);
  /** \brief . */
  virtual void linearCombinationImpl(
      const ArrayView<const Scalar> &alpha,
      const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &mv,
      const Scalar &beta);
  /** \brief . */
  virtual void dotsImpl(
      const MultiVectorBase<Scalar> &mv,
      const ArrayView<Scalar> &prods) const;
  /** \brief . */
  virtual void norms1Impl(
      const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms) const;
  /** \brief . */
  virtual void norms2Impl(
      const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms) const;
  /** \brief . */
  virtual void normsInfImpl(
      const ArrayView<typename ScalarTraits<Scalar>::magnitudeType> &norms) const;
  /** \brief . */
  RCP<const VectorBase<Scalar> > colImpl(Ordinal j) const;
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
  /** \brief . */
  void mvMultiReductApplyOpImpl(
      const RTOpPack::RTOpT<Scalar> &primary_op,
      const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
      const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
      const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
      const Ordinal primary_global_offset) const;
  /** \brief . */
  void acquireDetachedMultiVectorViewImpl(
      const Range1D &rowRng,
      const Range1D &colRng,
      RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv) const;
  /** \brief . */
  void releaseDetachedMultiVectorViewImpl(
      RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv) const;
  /** \brief . */
  void acquireNonconstDetachedMultiVectorViewImpl(
      const Range1D &rowRng,
      const Range1D &colRng,
      RTOpPack::SubMultiVectorView<Scalar> *sub_mv);
  /** \brief . */
  void commitNonconstDetachedMultiVectorViewImpl(
      RTOpPack::SubMultiVectorView<Scalar> *sub_mv);

  //@}

  /** \name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
      const EOpTransp M_trans,
      const MultiVectorBase<Scalar> &X,
      const Ptr<MultiVectorBase<Scalar> > &Y,
      const Scalar alpha,
      const Scalar beta) const;

  //@}

 public:
 private:
  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<MultiVectorBase<Scalar> > CNMVC;

  // //////////////////////////////
  // Private data members

  RCP<const DefaultProductVectorSpace<Scalar> > productSpace_;
  Teuchos::Array<CNMVC> multiVecs_;
  // cache
  int numBlocks_;

  // //////////////////////////////
  // Private member functions

  template <class MultiVectorType>
  void initializeImpl(
      const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
      const ArrayView<const RCP<MultiVectorType> > &multiVecs);

  void assertInitialized() const;

  void validateColIndex(const int j) const;
};

/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template <class Scalar>
RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector();

/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template <class Scalar>
RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const int numMembers);

/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template <class Scalar>
RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const ArrayView<const RCP<MultiVectorBase<Scalar> > > &multiVecs);

/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template <class Scalar>
RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const ArrayView<const RCP<const MultiVectorBase<Scalar> > > &multiVecs);

/** \brief Dynamic cast to a const product multi-vector or create a new
 * product multi-vector with one component if the input multi-vector is not a
 * product vector.
 *
 * \relates DefaultProductMultiVector
 */
template <class Scalar>
RCP<const ProductMultiVectorBase<Scalar> >
castOrCreateSingleBlockProductMultiVector(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const RCP<const MultiVectorBase<Scalar> > &mv);

/** \brief Dynamic cast to a const product multi-vector or create a new
 * product multi-vector with one component if the input multi-vector is not a
 * product vector.
 *
 * \relates DefaultProductMultiVector
 */
template <class Scalar>
RCP<ProductMultiVectorBase<Scalar> >
nonconstCastOrCreateSingleBlockProductMultiVector(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const RCP<MultiVectorBase<Scalar> > &mv);

// /////////////////////////
// Inline members

#ifndef TEUCHOS_DEBUG

template <class Scalar>
inline void DefaultProductMultiVector<Scalar>::assertInitialized() const {}

template <class Scalar>
inline void DefaultProductMultiVector<Scalar>::validateColIndex(const int /* j */) const {}

#endif  // TEUCHOS_DEBUG

}  // namespace Thyra

#endif  // THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DECL_HPP
