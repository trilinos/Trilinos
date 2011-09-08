// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
