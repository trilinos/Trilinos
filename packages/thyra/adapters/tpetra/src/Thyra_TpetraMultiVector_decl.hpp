// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_TPETRA_MULTIVECTOR_DECL_HPP
#define THYRA_TPETRA_MULTIVECTOR_DECL_HPP

#include "Thyra_SpmdMultiVectorDefaultBase.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Concrete implementation of Thyra::MultiVector in terms of
 * Tpetra::MultiVector.
 *
 * \todo Finish documentation!
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraMultiVector
  : virtual public SpmdMultiVectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct to uninitialized
  TpetraMultiVector();

  /** \brief Initialize.
   */
  void initialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
    );

  /** \brief Initialize.
   */
  void constInitialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
    );

  /** \brief Extract the underlying non-const Tpetra::MultiVector object.*/
  RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMultiVector();

  /** \brief Extract the underlying const Tpetra::MultiVector object.*/
  RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector() const;

  //@}

  /** @name Overridden public functions form MultiVectorAdapterBase */
  //@{
  /** \brief . */
  RCP< const ScalarProdVectorSpaceBase<Scalar> >
  domainScalarProdVecSpc() const;
  //@}

protected:

  /** @name Overridden protected functions from MultiVectorBase */
  //@{
  /** \brief . */
  virtual void assignImpl(Scalar alpha);

  /** \brief . */
  virtual void assignMultiVecImpl(const MultiVectorBase<Scalar>& mv);

  /** \brief . */
  virtual void scaleImpl(Scalar alpha);

  /** \brief . */
  virtual void updateImpl(
    Scalar alpha,
    const MultiVectorBase<Scalar>& mv
    );

  /** \brief . */
  virtual void linearCombinationImpl(
    const ArrayView<const Scalar>& alpha,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
    const Scalar& beta
    );

  /** \brief . */
  virtual void dotsImpl(
    const MultiVectorBase<Scalar>& mv,
    const ArrayView<Scalar>& prods
    ) const;

  /** \brief . */
  virtual void norms1Impl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const;

  /** \brief . */
  virtual void norms2Impl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const;

  /** \brief . */
  virtual void normsInfImpl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const;

  /** \brief . */
  RCP<const VectorBase<Scalar> > colImpl(Ordinal j) const;
  /** \brief . */
  RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j);

  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  contigSubViewImpl(const Range1D& colRng) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl(const Range1D& colRng);
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  nonContigSubViewImpl(const ArrayView<const int>& cols_in) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstNonContigSubViewImpl(const ArrayView<const int>& cols_in);

  /** \brief . */
  virtual void mvMultiReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
    const Ordinal primary_global_offset
    ) const;

  /** \brief . */
  void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const;

  /** \brief . */
  void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    );

  /** \brief . */
  void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    );

//  /** \brief . */
//  RCP<const MultiVectorBase<Scalar> >
//  nonContigSubViewImpl( const ArrayView<const int> &cols ) const;
//  /** \brief . */
//  RCP<MultiVectorBase<Scalar> >
//  nonconstNonContigSubViewImpl( const ArrayView<const int> &cols );
  //@}

  /** @name Overridden protected functions from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpaceImpl() const;
  /** \brief . */
  void getNonconstLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    );
  /** \brief . */
  void getLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const;

  //@}

  /** @name Overridden protected functions from MultiVectorAdapterBase */
  //@{
  /** \brief . */
  virtual void euclideanApply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members

  RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraVectorSpace_;
  RCP<const ScalarProdVectorSpaceBase<Scalar> > domainSpace_;
  Teuchos::ConstNonconstObjectContainer<Tpetra::MultiVector<Scalar, LocalOrdinal,GlobalOrdinal,Node> >
  tpetraMultiVector_;

  // ////////////////////////////////////
  // Private member functions

  template<class TpetraMultiVector_t>
  void initializeImpl(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const RCP<TpetraMultiVector_t> &tpetraMultiVector
    );

  // Non-throwing Tpetra MultiVector extraction methods.
  // Return null if casting failed.
  RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMultiVector(const RCP<MultiVectorBase<Scalar> >& mv) const;

  RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const;

};


/** \brief Nonmember constructor for non-const TpetraMultiVector.
 *
 * \relates TpetraMultiVector.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraMultiVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
  )
{
  RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmv =
    Teuchos::rcp(new TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  tmv->initialize(tpetraVectorSpace, domainSpace, tpetraMultiVector);
  return tmv;
}


/** \brief Nonmember constructor for const TpetraMultiVector.
 *
 * \relates TpetraMultiVector.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
constTpetraMultiVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
  )
{
  RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmv =
    Teuchos::rcp(new TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  tmv->constInitialize(tpetraVectorSpace, domainSpace, tpetraMultiVector);
  return tmv;
}


} // end namespace Thyra


#endif // THYRA_TPETRA_MULTIVECTOR_DECL_HPP
