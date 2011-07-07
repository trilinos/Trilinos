// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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
