// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_DECL_HPP
#define THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_DECL_HPP

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_VectorDefaultBase.hpp"

namespace Thyra {

/** \brief . */
template <class Scalar> class DefaultClusteredSpmdProductVectorSpace;

/** \brief Concrete implementation of a clustered Spmd-based product vector.
 *
 * ToDo: Finish documentation!
 *
 * The default constructor is made private to avoid accidental default
 * construction.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class Scalar>
class DefaultClusteredSpmdProductVector
  : virtual public ProductVectorBase<Scalar>
  , virtual protected VectorDefaultBase<Scalar>
{
public:

#ifndef _MSC_VER
  /** \brief . */
  using ProductVectorBase<Scalar>::applyOp;
#endif
  
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized (see postconditions for <tt>uninitialize()</tt>). */
  DefaultClusteredSpmdProductVector();

  /** \brief Constructs to initialized (calls <tt>initialize()</tt>). */
  DefaultClusteredSpmdProductVector(
    const Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  &productSpace
    ,const Teuchos::RCP<VectorBase<Scalar> >                                   vecs[]
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  &productSpace
    ,const Teuchos::RCP<VectorBase<Scalar> >                                   vecs[]
    );

  /** \brief Uninitialize.
   *
   * ToDo: Finish documentation.
   */
  void uninitialize(
    Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >  *productSpace = NULL
    ,Teuchos::RCP<VectorBase<Scalar> >                                   vecs[]        = NULL
    );

  //@}

  /** @name Overridden from ProductVectorBase */
  //@{

  /** \brief . */
  Teuchos::RCP<VectorBase<Scalar> > getNonconstVectorBlock(const int k); 
  /** \brief . */
  Teuchos::RCP<const VectorBase<Scalar> > getVectorBlock(const int k) const;

  //@}

  /** @name Overridden from ProductMultiVectorBase */
  //@{

  /** \brief . */
  Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > productSpace() const;
  /** \brief . */
  bool blockIsConst(const int k) const; 
  /** \brief . */
  Teuchos::RCP<MultiVectorBase<Scalar> > getNonconstMultiVectorBlock(const int k); 
  /** \brief . */
  Teuchos::RCP<const MultiVectorBase<Scalar> > getMultiVectorBlock(const int k) const;

  //@}

  /** @name Overridden from VectorBase */
  //@{

  /** \brief . */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > space() const;

  //@}

protected:

  /** @name Overridden protected members from VectorBase */
  //@{

  /** \brief . */
  void applyOpImpl(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const;

  //@}

private:

  // //////////////////////////////
  // Private data members

  Teuchos::RCP<const DefaultClusteredSpmdProductVectorSpace<Scalar> >   productSpace_;
  std::vector<Teuchos::RCP<VectorBase<Scalar> > >                       vecs_;

};

} // namespace Thyra

#endif // THYRA_DEFAULT_CLUSTERED_SPMD_PRODUCT_VECTOR_DECL_HPP
