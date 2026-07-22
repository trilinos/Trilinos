// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP
#define THYRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP

#include "Thyra_EuclideanScalarProd_decl.hpp"
#include "Thyra_ScalarProdBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"


namespace Thyra {


template<class Scalar>
bool EuclideanScalarProd<Scalar>::isEuclideanImpl() const
{
  return true;
}


template<class Scalar>
void EuclideanScalarProd<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out
  ) const
{
  dots(X, Y, scalarProds_out);
}


} // end namespace Thyra


#endif  // THYRA_EUCLIDEAN_SCALAR_PROD_DEF_HPP
