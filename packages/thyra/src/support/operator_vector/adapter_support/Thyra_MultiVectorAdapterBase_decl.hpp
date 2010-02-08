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

#ifndef THYRA_MULTI_VECTOR_ADAPTER_BASE_DECL_HPP
#define THYRA_MULTI_VECTOR_ADAPTER_BASE_DECL_HPP

#include "Thyra_MultiVectorDefaultBase.hpp"


namespace Thyra {


/** \brief Forward decl. */
template<class Scalar> class ScalarProdVectorSpaceBase;


/** \brief Node subclass for MultiVectorBase subclasses that allows the
 * insertion of an application defined scalar product.
 *
 * Most concrete MultiVector adapter subclasses should derive from this base
 * subclass in order to allow for the incorporate of application-defined
 * scalar products.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class MultiVectorAdapterBase : virtual public MultiVectorDefaultBase<Scalar>
{
public:

  /** @name Pure virtual functions to override in subclasses */
  //@{

  /** \brief . */
  virtual RCP<const ScalarProdVectorSpaceBase<Scalar> >
  rangeScalarProdVecSpc() const = 0;

  /** \brief . */
  virtual RCP<const ScalarProdVectorSpaceBase<Scalar> >
  domainScalarProdVecSpc() const = 0;

  /** \brief Apply the linear operator to a multi-vector with respect
   * to a Euclidean vector space where the scalar product is the dot
   * product.
   *
   * Preconditions:<ul>
   * <li><tt>this->applySupports(conj)==true</tt>
   * </ul>
   */
  virtual void euclideanApply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar>   &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const = 0;

  //@}

  /** @name Overridden functions from LinearOp */
  //@{
  /// Returns <tt>this->rangeScalarProdVecSpc()</tt>
  RCP<const VectorSpaceBase<Scalar> > range() const;
  /// Returns <tt>this->domainScalarProdVecSpc()</tt>
  RCP<const VectorSpaceBase<Scalar> > domain() const;
  //@}

protected:

  /** @name Overridden protected functions from LinearOpBase */
  //@{
  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief .  */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;
  //@}

};


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_ADAPTER_BASE_DECL_HPP
