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

#ifndef THYRA_DEFAULT_CLUSTERED_MPI_PRODUCT_VECTOR_DECL_HPP
#define THYRA_DEFAULT_CLUSTERED_MPI_PRODUCT_VECTOR_DECL_HPP

#include "Thyra_ProductVectorBase.hpp"

namespace Thyra {

/** \brief . */
template <class Scalar> class DefaultClusteredMPIProductVectorSpace;

/** \brief . */
template <class Scalar> class MPIVectorBase;

/** \brief Concrete implementation of a clustered MPI-based product vector.
 *
 * ToDo: Finish documentation!
 *
 * The default constructor is made private to avoid accidental default
 * construction.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_support_grp
 */
template<class Scalar>
class DefaultClusteredMPIProductVector
  : virtual public ProductVectorBase<Scalar>
  , virtual protected VectorDefaultBase<Scalar>
{
public:

  /** \brief . */
  using ProductVectorBase<Scalar>::applyOp;

  /** @name Constructors/initializers/accessors */
  //@{

  /// Constructs to initialized (calls <tt>initialize()</tt>).
  DefaultClusteredMPIProductVector(
    const Teuchos::RefCountPtr<const DefaultClusteredMPIProductVectorSpace<Scalar> >  &productSpace
    ,const Teuchos::RefCountPtr<VectorBase<Scalar> >                                  vecs[]
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const Teuchos::RefCountPtr<const DefaultClusteredMPIProductVectorSpace<Scalar> >  &productSpace
    ,const Teuchos::RefCountPtr<VectorBase<Scalar> >                                  vecs[]
    );

  /** \brief Uninitialize.
   *
   * ToDo: Finish documentation.
   */
  void uninitialize(
    Teuchos::RefCountPtr<const DefaultClusteredMPIProductVectorSpace<Scalar> >  *productSpace = NULL
    ,Teuchos::RefCountPtr<VectorBase<Scalar> >                                  vecs[]        = NULL
    );

  //@}

  /** @name Overridden from ProductVectorBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > productSpace() const;
  /** \brief . */
  Teuchos::RefCountPtr<VectorBase<Scalar> > getBlock(const int k); 
  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > getBlock(const int k) const;

  //@}

  /** @name Overridden from VectorBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > space() const;
  /** \brief . */
  void applyOp(
    const RTOpPack::RTOpT<Scalar>    &op
    ,const int                       num_vecs
    ,const VectorBase<Scalar>*const  vecs[]
    ,const int                       num_targ_vecs
    ,VectorBase<Scalar>*const        targ_vecs[]
    ,RTOpPack::ReductTarget          *reduct_obj
    ,const Index                     first_ele
    ,const Index                     sub_dim
    ,const Index                     global_offset
    ) const;

  //@}

private:

  // //////////////////////////////
  // Private data members

  Teuchos::RefCountPtr<const DefaultClusteredMPIProductVectorSpace<Scalar> >       productSpace_;
  std::vector<Teuchos::RefCountPtr<VectorBase<Scalar> > >                          vecs_;

protected:

  // //////////////////////////////
  // Protected member functions
  // Added to allow TSFExtended DefaultClusteredMPIProductVector to derive from this.
  DefaultClusteredMPIProductVector();

};

} // namespace Thyra

#endif // THYRA_DEFAULT_CLUSTERED_MPI_PRODUCT_VECTOR_DECL_HPP
