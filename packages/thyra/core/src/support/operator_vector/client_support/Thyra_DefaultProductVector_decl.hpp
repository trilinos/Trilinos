// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_PRODUCT_VECTOR_DECL_HPP
#define THYRA_DEFAULT_PRODUCT_VECTOR_DECL_HPP

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_as.hpp"


namespace Thyra {


/** \brief . */
template <class Scalar> class DefaultProductVectorSpace;


/** \brief Concrete implementation of a product vector.
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
template<class Scalar>
class DefaultProductVector
  : virtual public ProductVectorBase<Scalar>,
    virtual protected VectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  DefaultProductVector();

  /** \brief Constructs to initialized (calls <tt>initialize()</tt>). */
  DefaultProductVector(
    const RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const ArrayView<const RCP<VectorBase<Scalar> > > &vecs
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const ArrayView<const RCP<const VectorBase<Scalar> > > &vecs
    );

  /** \brief Uninitialize.
   *
   * ToDo: Finish documentation.
   */
  void uninitialize();

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
 
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** @name Extensions to ProductVectorBase suitable for physically-blocked vectors */
  //@{

  /** \brief . */
  void setBlock(int i, const RCP<const VectorBase<Scalar> >& b); 
  /** \brief . */
  void setNonconstBlock(int i, const RCP<VectorBase<Scalar> >& b); 
  //@}

  /** @name Overridden from ProductVectorBase */
  //@{

  /** \brief . */
  RCP<VectorBase<Scalar> > getNonconstVectorBlock(const int k); 
  /** \brief . */
  RCP<const VectorBase<Scalar> > getVectorBlock(const int k) const;

  //@}

  /** @name Overridden public functions from ProductMultiVectorBase */
  //@{

  /** \brief . */
  RCP<const ProductVectorSpaceBase<Scalar> > productSpace() const;
  /** \brief . */
  bool blockIsConst(const int k) const; 
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  getNonconstMultiVectorBlock(const int k);
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  getMultiVectorBlock(const int k) const;

  //@}

  /** @name Overridden public functions from VectorBase */
  //@{

  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > space() const;

  //@}

protected:

  /** @name Overridden protected functions from VectorBase */
  //@{

  /** \brief . */
  //virtual void randomizeImpl(Scalar l, Scalar u);
  /** \brief . */
  virtual void absImpl(const VectorBase<Scalar>& x);
  /** \brief . */
  virtual void reciprocalImpl(const VectorBase<Scalar>& x);
  /** \brief . */
  virtual void eleWiseScaleImpl(const VectorBase<Scalar>& x);
  /** \brief . */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  norm2WeightedImpl(const VectorBase<Scalar>& x) const;
  /** \brief . */
  void applyOpImpl(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const;
  /** \brief . */
  void acquireDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief . */
  void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief . */
  void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief . */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief . */
  void setSubVectorImpl(
    const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
    );

  //@}

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

  //#}

public:

private:

  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<VectorBase<Scalar> > CNVC;

  // //////////////////////////////
  // Private data members

  RCP<const DefaultProductVectorSpace<Scalar> > productSpace_;
  Array<CNVC> vecs_;
  // cache
  int numBlocks_;

};


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
inline
RCP<DefaultProductVector<Scalar> >
defaultProductVector()
{
  return Teuchos::rcp(new DefaultProductVector<Scalar>);
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
inline
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace
  )
{
  return Teuchos::rcp(
    new DefaultProductVector<Scalar>(productSpace)
    );
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const ArrayView<const RCP<VectorBase<Scalar> > > &vecs
  )
{
  RCP<DefaultProductVector<Scalar> > pv = defaultProductVector<Scalar>();
  pv->initialize(productSpace, vecs);
  return pv;
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const ArrayView<const RCP<const VectorBase<Scalar> > > &vecs
  )
{
  RCP<DefaultProductVector<Scalar> > pv = defaultProductVector<Scalar>();
  pv->initialize(productSpace, vecs);
  return pv;
}


/** \brief Return a casted non-const ProductVectorBase object or create a new
 * DefaultProductVector object with one component.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<ProductVectorBase<Scalar> >
castOrCreateNonconstProductVectorBase(const RCP<VectorBase<Scalar> > v);


/** \brief Return a casted const ProductVectorBase object or create a new
 * DefaultProductVector object with one component.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<const ProductVectorBase<Scalar> >
castOrCreateProductVectorBase(const RCP<const VectorBase<Scalar> > v);


} // namespace Thyra


#endif // THYRA_DEFAULT_PRODUCT_VECTOR_DECL_HPP
