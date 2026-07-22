// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
