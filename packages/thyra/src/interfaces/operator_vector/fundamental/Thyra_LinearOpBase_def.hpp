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

#ifndef THYRA_LINEAR_OP_BASE_DEF_HPP
#define THYRA_LINEAR_OP_BASE_DEF_HPP

#include "Thyra_LinearOpBase_decl.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace Thyra {


// Public interface functions


template<class Scalar>
RCP<const LinearOpBase<Scalar> > 
LinearOpBase<Scalar>::clone() const
{
  return Teuchos::null;
}


// Deprecated


template<class Scalar>
bool LinearOpBase<Scalar>::applySupports(
  const EConj conj
  ) const
{
  return Thyra::opSupported(*this, applyConjToTrans(conj));
}

template<class Scalar>
void LinearOpBase<Scalar>::apply(
  const EConj conj,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply(*this, applyConjToTrans(conj), X, Teuchos::ptr(Y), alpha, beta);
}


template<class Scalar>
bool LinearOpBase<Scalar>::applyTransposeSupports(
  const EConj conj
  ) const
{
  return Thyra::opSupported(*this, applyTransposeConjToTrans(conj));
}


template<class Scalar>
void LinearOpBase<Scalar>::applyTranspose(
  const EConj conj,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply(*this, applyTransposeConjToTrans(conj), X, Teuchos::ptr(Y), alpha, beta);
}


}	// end namespace Thyra


// ToDo: You can move this back to the decl file after you have refactored
// apply(...) to not use raw pointers.  Otherwise the Y.ptr() call needs to have
// the definition of MultiVectorBase.


template<class Scalar>
void Thyra::apply(
  const LinearOpBase<Scalar> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  )
{
  M.apply(M_trans, X, Y, alpha, beta);
}


//
// Explicit instantiation macro
//

#define THYRA_LINEAR_OP_BASE_INSTANT(SCALAR) \
  \
  template class LinearOpBase<SCALAR >; \
  \
  template void apply(  \
    const LinearOpBase<SCALAR > &M,  \
    const EOpTransp M_trans,  \
    const MultiVectorBase<SCALAR > &X,  \
    const Ptr<MultiVectorBase<SCALAR > > &Y,  \
    const SCALAR alpha,  \
    const SCALAR beta  \
    );


#endif // THYRA_LINEAR_OP_BASE_DEF_HPP
