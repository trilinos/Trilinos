// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SCALAR_PROD_VECTOR_SPACE_BASE_DEF_HPP
#define THYRA_SCALAR_PROD_VECTOR_SPACE_BASE_DEF_HPP

#include "Thyra_ScalarProdVectorSpaceBase_decl.hpp"
#include "Thyra_VectorSpaceDefaultBase.hpp"
#include "Thyra_EuclideanScalarProd.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// Constructors / initializers


template<class Scalar>
ScalarProdVectorSpaceBase<Scalar>::ScalarProdVectorSpaceBase()
  :scalarProd_(Teuchos::rcp(new EuclideanScalarProd<Scalar>()))
{}

  
template<class Scalar>
ScalarProdVectorSpaceBase<Scalar>::ScalarProdVectorSpaceBase(
  const Teuchos::RCP<const ScalarProdBase<Scalar> > &scalarProd_in
  )
  :scalarProd_(scalarProd_in.assert_not_null())
{}


template<class Scalar>
void ScalarProdVectorSpaceBase<Scalar>::setScalarProd(
  const Teuchos::RCP<const ScalarProdBase<Scalar> > &scalarProd_in
  )
{
  scalarProd_ = scalarProd_in.assert_not_null();
}


template<class Scalar>
Teuchos::RCP<const ScalarProdBase<Scalar> >
ScalarProdVectorSpaceBase<Scalar>::getScalarProd() const
{
  return scalarProd_;
}


// Overridden from VectorSpaceBase


template<class Scalar>
bool ScalarProdVectorSpaceBase<Scalar>::isEuclidean() const
{
  return scalarProd_->isEuclidean();
}


template<class Scalar>
Scalar ScalarProdVectorSpaceBase<Scalar>::scalarProd( 
  const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_VEC_SPACES("ScalarProdVectorSpaceBase<Scalar>::scalarProd(...)",
    *x.space(), *this);
  THYRA_ASSERT_VEC_SPACES("ScalarProdVectorSpaceBase<Scalar>::scalarProd(...)",
    *y.space(), *this);
#endif
  return scalarProd_->scalarProd(x,y);
}


template<class Scalar>
void ScalarProdVectorSpaceBase<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_VEC_SPACES("ScalarProdVectorSpaceBase<Scalar>::scalarProds(...)",
    *X.range(), *this);
  THYRA_ASSERT_VEC_SPACES("ScalarProdVectorSpaceBase<Scalar>::scalarProds(...)",
    *Y.range(), *this);
  THYRA_ASSERT_VEC_SPACES("ScalarProdVectorSpaceBase<Scalar>::scalarProds(...)",
    *X.domain(), *Y.domain());
#endif
  scalarProd_->scalarProds(X, Y, scalarProds_out);
}


} // end namespace Thyra


#endif  // THYRA_SCALAR_PROD_VECTOR_SPACE_BASE_DEF_HPP
