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

#ifndef THYRA_DIAGONAL_LINEAR_OP_DEF_HPP
#define THYRA_DIAGONAL_LINEAR_OP_DEF_HPP


#include "Thyra_DefaultDiagonalLinearOp_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_AssertOp.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp()
{}


template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp(
  const RCP<const VectorSpaceBase<Scalar> > &space
  )
{
  initialize(space);
}


template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp(
  const RCP<VectorBase<Scalar> > &diag
  )
{
  initialize(diag);
}


template<class Scalar>
DefaultDiagonalLinearOp<Scalar>::DefaultDiagonalLinearOp(
  const RCP<const VectorBase<Scalar> > &diag
  )
{
  initialize(diag);
}


template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::initialize(
  const RCP<const VectorSpaceBase<Scalar> > &space
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(space.get()==NULL);
#endif
  initialize(createMember(space)); // Note that the space is guaranteed to be remembered here!
}


template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::initialize(
  const RCP<VectorBase<Scalar> > &diag
  )
{
  diag_.initialize(diag);
}


template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::initialize(
  const RCP<const VectorBase<Scalar> > &diag
  )
{
  diag_.initialize(diag);
}


template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::uninitialize()
{
  diag_.uninitialize();
}


// Overridden from DiagonalLinearOpBase


template<class Scalar>
bool DefaultDiagonalLinearOp<Scalar>::isDiagConst() const
{
  return diag_.isConst();
}


template<class Scalar>
RCP<VectorBase<Scalar> > 
DefaultDiagonalLinearOp<Scalar>::getNonconstDiag()
{
  return diag_.getNonconstObj();
}


template<class Scalar>
RCP<const VectorBase<Scalar> > 
DefaultDiagonalLinearOp<Scalar>::getDiag() const
{
  return diag_.getConstObj();
}


// Overridden from LinearOpBase


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultDiagonalLinearOp<Scalar>::range() const
{
  return diag_.getConstObj()->space();
}


template<class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DefaultDiagonalLinearOp<Scalar>::domain() const
{
  return diag_.getConstObj()->space();
}


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
DefaultDiagonalLinearOp<Scalar>::clone() const
{
  return Teuchos::rcp(new DefaultDiagonalLinearOp<Scalar>(diag_.getConstObj()->clone_v()));
}


// protected


// Protected functions overridden from LinearOpBase


template<class Scalar>
bool DefaultDiagonalLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  return true;
}


template<class Scalar>
void DefaultDiagonalLinearOp<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultDiagonalLinearOp<Scalar>::apply(...)",*this, M_trans, X, &*Y
    );
#endif // TEUCHOS_DEBUG

  // Y = beta * Y

  if( beta != ST::one() ) scale<Scalar>(beta, Y);

  // Y += alpha *op(M) * X

  const Ordinal m = X.domain()->dim();

  for (Ordinal col_j = 0; col_j < m; ++col_j) {
    const RCP<const VectorBase<Scalar> > x = X.col(col_j);
    const RCP<VectorBase<Scalar> > y = Y->col(col_j);
    if (ST::isComplex) {
      if ( M_trans==NOTRANS || M_trans==TRANS ) {
        ele_wise_prod( alpha, *diag_.getConstObj(), *x, y.ptr() );
      }
      else {
        ele_wise_conj_prod( alpha, *diag_.getConstObj(), *x, y.ptr() );
      }
    }
    else {
      ele_wise_prod( alpha, *diag_.getConstObj(), *x, y.ptr() );
    }
  }

}


}	// end namespace Thyra


#endif	// THYRA_DIAGONAL_LINEAR_OP_DEF_HPP
