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

#ifndef THYRA_EUCLIDEAN_LINEAR_OP_DECL_HPP
#define THYRA_EUCLIDEAN_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpBaseDecl.hpp"

namespace Thyra {

/** \brief Base interface for Euclidean linear operators.
 *
 * Most generic subclass implementations should derive from this base
 * interface since it allows an application-specific definition of the
 * scalar product.  Note that almost every concrete implementation of
 * <tt>MultiVectorBase</tt> should derive from this interface.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class EuclideanLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** \brief . */
  using LinearOpBase<Scalar>::apply;

  /** @name Pure virtual functions to override in subclasses */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const = 0;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const = 0;

  /** \brief Apply the linear operator to a vector with respect to a
   * Euclidean vector space where the scalar product is the dot
   * product.
   */
  virtual void euclideanApply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief Apply the linear operator to a multi-vector with respect
   * to a Euclidean vector space where the scalar product is the dot
   * product.
   *
   * The default implementation calls the single-vector version
   * <tt>this->euclideanApply()</tt>.  A subclass should only override
   * this version if it can do something more intelligent with
   * multi-vectors that simply applying the operator one column at a
   * time.
   */
  virtual void euclideanApply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

  /** @name Overridden functions from OpBase */
  //@{
  /// Returns <tt>this->rangeScalarProdVecSpc()</tt>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > range() const;
  /// Returns <tt>this->domainScalarProdVecSpc()</tt>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > domain() const;
  //@}

  /** @name Overridden functions from LinearOpBase */
  //@{
  /** \brief Apply the linear operator to a vector using an
   * application-specific definition of the scalar product.
   *
   * ToDo: Finish Documentation!
   */
  void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha
    ,const Scalar                beta
    ) const;
  /** \brief Apply the linear operator to a multi-vector using an
   * application-specific definition of the scalar product.
   *
   * ToDo: Finish Documentation!
   */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  //@}

protected:
  
  void euclidean_apply_impl(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

}; // end class EuclideanLinearOpBase<Scalar>

} // namespace Thyra

#endif // THYRA_EUCLIDEAN_LINEAR_OP_DECL_HPP
