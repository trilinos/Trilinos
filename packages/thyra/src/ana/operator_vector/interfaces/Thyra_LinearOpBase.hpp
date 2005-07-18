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
#include "Thyra_AssertOp.hpp"
#include "Thyra_describeLinearOp.hpp"

namespace Thyra {

// Virtual functions with default implementations

template<class RangeScalar, class DomainScalar>
bool LinearOpBase<RangeScalar,DomainScalar>::applySupports( const EConj conj ) const
{
  return ( Teuchos::ScalarTraits<Scalar>::isComplex ? conj==NONCONJ_ELE : true );
}

template<class RangeScalar, class DomainScalar>
bool LinearOpBase<RangeScalar,DomainScalar>::applyTransposeSupports( const EConj conj ) const
{
  return false;
}

template<class RangeScalar, class DomainScalar>
void LinearOpBase<RangeScalar,DomainScalar>::applyTranspose(
  const EConj                            conj
  ,const MultiVectorBase<RangeScalar>    &X
  ,MultiVectorBase<DomainScalar>         *Y
  ,const DomainScalar                    alpha
  ,const DomainScalar                    beta
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"LinearOpBase<"<<Teuchos::ScalarTraits<RangeScalar>::name()<<","<<Teuchos::ScalarTraits<DomainScalar>::name()<<">::applyTranspose(...): "
    "Error, the concrete subclass described as { " << this->description() << " } "
    " with this->applyTransposeSupports("<<toString(conj)<<")="<<this->applyTransposeSupports(conj)
    << " did not override this function and does not support transposes."
    );
}

template<class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> > 
LinearOpBase<RangeScalar,DomainScalar>::clone() const
{
  return Teuchos::null;
}

// Overridden from Teuchos::Describable

template<class RangeScalar, class DomainScalar>
std::ostream& LinearOpBase<RangeScalar,DomainScalar>::describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const
{
  return describeLinearOp(*this,out,verbLevel,leadingIndent,indentSpacer);
}

}	// end namespace Thyra

#endif // THYRA_LINEAR_OP_BASE_HPP
