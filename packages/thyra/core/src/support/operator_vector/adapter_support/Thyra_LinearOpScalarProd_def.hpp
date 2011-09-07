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

#ifndef THYRA_LINEAR_OP_SCALAR_PROD_DEF_HPP
#define THYRA_LINEAR_OP_SCALAR_PROD_DEF_HPP

#include "Thyra_LinearOpScalarProd_decl.hpp"
#include "Thyra_ScalarProdBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace Thyra {


// Constructors, initializers, accessors


template<class Scalar>
LinearOpScalarProd<Scalar>::LinearOpScalarProd()
{}


template<class Scalar>
LinearOpScalarProd<Scalar>::LinearOpScalarProd(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &op_in )
{
  this->initialize(op_in);
}


template<class Scalar>
void LinearOpScalarProd<Scalar>::initialize(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &op_in
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(op_in));
#endif
  op_ = op_in;
}


template<class Scalar>
void LinearOpScalarProd<Scalar>::uninitialize(
  const Ptr<RCP<const LinearOpBase<Scalar> > > &op_out
  )
{
  if (!is_null(op_out)) *op_out = op_;
  op_ = Teuchos::null;
}


// Overridden from ScalarProdBase


template<class Scalar>
bool LinearOpScalarProd<Scalar>::isEuclideanImpl() const
{
  return false;
}


template<class Scalar>
void LinearOpScalarProd<Scalar>::scalarProdsImpl(
  const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
  const ArrayView<Scalar> &scalarProds_out
  ) const
{
  Teuchos::RCP<MultiVectorBase<Scalar> >
    T = createMembers(Y.range() ,Y.domain()->dim());
  Thyra::apply(*op_, NOTRANS,Y, T.ptr());
  dots(X, *T, scalarProds_out);
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
LinearOpScalarProd<Scalar>::getLinearOpImpl() const
{
  return op_;
}


} // end namespace Thyra


#endif  // THYRA_LINEAR_OP_SCALAR_PROD_DEF_HPP
