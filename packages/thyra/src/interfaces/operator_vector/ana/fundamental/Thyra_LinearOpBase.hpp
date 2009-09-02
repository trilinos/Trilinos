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

#ifndef THYRA_LINEAR_OP_BASE_HPP
#define THYRA_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBaseDecl.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace Thyra {


// Virtual functions with default implementations


template<class RangeScalar, class DomainScalar>
bool LinearOpBase<RangeScalar,DomainScalar>::applySupports(
  const EConj conj
  ) const
{
  return (
    Teuchos::ScalarTraits<RangeScalar>::isComplex ? conj==NONCONJ_ELE : true
    );
}


template<class RangeScalar, class DomainScalar>
bool LinearOpBase<RangeScalar,DomainScalar>::applyTransposeSupports(
  const EConj conj
  ) const
{
  return false;
}


template<class RangeScalar, class DomainScalar>
void LinearOpBase<RangeScalar,DomainScalar>::applyTranspose(
  const EConj conj,
  const MultiVectorBase<RangeScalar> &X,
  MultiVectorBase<DomainScalar> *Y,
  const DomainScalar alpha,
  const DomainScalar beta
  ) const
{
  const std::string
    &typeName = TypeNameTraits<LinearOpBase<RangeScalar,DomainScalar> >::name();
  TEST_FOR_EXCEPTION(
    true, std::logic_error,
    typeName << "::applyTranspose(...): " <<
    "Error, the concrete subclass described as { " << this->description() << " } "
    " with this->applyTransposeSupports("<<toString(conj)<<")="<<this->applyTransposeSupports(conj)
    << " did not override this function and does not support transposes."
    );
}


template<class RangeScalar, class DomainScalar>
RCP<const LinearOpBase<RangeScalar,DomainScalar> > 
LinearOpBase<RangeScalar,DomainScalar>::clone() const
{
  return Teuchos::null;
}


}	// end namespace Thyra


// ToDo: You can move this back to the decl file after you have refactored
// apply(...) to not use raw pointers.  Otherwise the &*Y call needs to have
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
  if(real_trans(M_trans)==NOTRANS) {
    M.apply(transToConj(M_trans), X, &*Y, alpha, beta);
  }
  else {
    M.applyTranspose(transToConj(M_trans), X, &*Y, alpha, beta);
  }
}



#endif // THYRA_LINEAR_OP_BASE_HPP
