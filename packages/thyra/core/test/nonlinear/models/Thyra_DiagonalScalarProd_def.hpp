/*
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
*/

#ifndef THYRA_DIAGONAL_SCALAR_PROD_DEF_HPP
#define THYRA_DIAGONAL_SCALAR_PROD_DEF_HPP


#include "Thyra_DiagonalScalarProd_decl.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


// Consturctors/Initializers/Accessors


template<class Scalar>
DiagonalScalarProd<Scalar>::DiagonalScalarProd()
{}


template<class Scalar>
void DiagonalScalarProd<Scalar>::initialize(
  const RCP<const VectorBase<Scalar> > &s_diag )
{
  s_diag_ = s_diag.assert_not_null();
}


// Overridden from ScalarProdBase


template<class Scalar>
bool DiagonalScalarProd<Scalar>::isEuclideanImpl() const
{
  return false;
}


template<class Scalar>
void DiagonalScalarProd<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out ) const
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Ordinal Ordinal;

  const Ordinal m = X.domain()->dim();

#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_VEC_SPACES( "DiagonalScalarProd<Scalar>::scalarProds(X,Y,sclarProds)",
    *s_diag_->space(), *Y.range() );
  THYRA_ASSERT_VEC_SPACES( "DiagonalScalarProd<Scalar>::scalarProds(X,Y,sclarProds)",
    *X.range(), *Y.range() );
  THYRA_ASSERT_VEC_SPACES( "DiagonalScalarProd<Scalar>::scalarProds(X,Y,sclarProds)",
    *X.domain(), *Y.domain() );
  TEUCHOS_ASSERT_EQUALITY( as<Ordinal>(scalarProds_out.size()), m );
#endif

  const ConstDetachedSpmdVectorView<Scalar> s_diag(s_diag_);

  const RCP<const Teuchos::Comm<Ordinal> > comm = s_diag.spmdSpace()->getComm();
  
  for (Ordinal j = 0; j < m; ++j) {

    const ConstDetachedSpmdVectorView<Scalar> x(X.col(j));
    const ConstDetachedSpmdVectorView<Scalar> y(Y.col(j));
    
    Scalar scalarProd_j = ST::zero();

    for (Ordinal i = 0; i < x.subDim(); ++i) {
      scalarProd_j += ST::conjugate(x[i]) * s_diag[i] * y[i];
    }

    if (!is_null(comm)) {
      Scalar g_scalarProd_j = 0.0;
      Teuchos::reduceAll<Ordinal,Scalar>(
        *comm, Teuchos::REDUCE_SUM,
        scalarProd_j,
        Teuchos::outArg(g_scalarProd_j)
        );
      scalarProds_out[j] = g_scalarProd_j;
    }
    else {
      scalarProds_out[j] = scalarProd_j;
    }

  }

}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DiagonalScalarProd<Scalar>::getLinearOpImpl() const
{
  TEST_FOR_EXCEPT_MSG(true, "ToDo: Implement when needed!")
  return Teuchos::null;
}


} // end namespace Thyra


#endif  // THYRA_DIAGONAL_SCALAR_PROD_DEF_HPP
