// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP
#define THYRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP


#include "Thyra_ScalarProdBase_decl.hpp"


namespace Thyra {

/** \brief Concrete implementation of a scalar product for a Euclidean vector
 * space (i.e. using the dot product).
 *
 * Because this subclass is implemented using an RTOp, it will work with any
 * <tt>VectorBase</tt> or <tt>MultiVectorBase</tt> implementation no matter
 * what.
 *
 * \ingroup Thyra_Op_Vec_basic_adapter_support_grp
 */
template<class Scalar>
class EuclideanScalarProd : public ScalarProdBase<Scalar> {
protected:
  
  /** @name Overridden from ScalarProdBase */
  //@{

  /** \brief Returns <tt>true</tt>. */
  virtual bool isEuclideanImpl() const;

  /** \brief Simply calls <tt>dots(X,Y,scalar_prods)</tt>. */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds
    ) const;

  //@}

};


} // end namespace Thyra


#endif  // THYRA_EUCLIDEAN_SCALAR_PROD_DECL_HPP
