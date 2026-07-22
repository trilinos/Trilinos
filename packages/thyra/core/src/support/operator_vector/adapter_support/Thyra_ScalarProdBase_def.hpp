// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SCALAR_PROD_BASE_DEF_HPP
#define THYRA_SCALAR_PROD_BASE_DEF_HPP

#include "Thyra_ScalarProdBase_decl.hpp"
#include "Thyra_VectorBase.hpp"


namespace Thyra {


// Protected virtual functions


template<class Scalar>
Scalar ScalarProdBase<Scalar>::scalarProdImpl(
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& y
  ) const
{
  Tuple<Scalar,1> scalarProds_out;
  this->scalarProds(x, y, scalarProds_out());
  return scalarProds_out[0];
}


} // end namespace Thyra


#endif  // THYRA_SCALAR_PROD_BASE_DEF_HPP
