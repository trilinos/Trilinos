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
