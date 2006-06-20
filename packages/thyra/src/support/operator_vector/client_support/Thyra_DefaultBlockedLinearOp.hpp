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

#ifndef THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP
#define THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP

#include "Thyra_DefaultBlockedLinearOpDecl.hpp"

namespace Thyra {

// Constructors

template<class Scalar>
DefaultBlockedLinearOp<Scalar>::DefaultBlockedLinearOp()
{}

// Overridden from PhysicallyBlockedLinearOpBase

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill()
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::beginBlockFill(
  const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >  &productRange
  ,const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > &productDomain
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockFillIsActive() const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::acceptsBlock(
  const int i, const int j
  ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::setNonconstBlock(
  const int i, const int j
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &block
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::setBlock(
  const int i, const int j
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &block
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::endBlockFill()
{
  TEST_FOR_EXCEPT(true);
}

// Overridden from BlockedLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::productRange() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::productDomain() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockExists(
  const int i, const int j
  ) const
{
  TEST_FOR_EXCEPT(true);
} 

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::blockIsConst(
  const int i, const int j
  ) const
{
  TEST_FOR_EXCEPT(true);
} 

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::getNonconstBlock(const int i, const int j)
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
} 

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::getBlock(const int i, const int j) const
{
  TEST_FOR_EXCEPT(true);
} 

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::range() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::domain() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultBlockedLinearOp<Scalar>::clone() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultBlockedLinearOp<Scalar>::description() const
{
  TEST_FOR_EXCEPT(true);
  return "";
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  TEST_FOR_EXCEPT(true);
}

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool DefaultBlockedLinearOp<Scalar>::opSupported(
  ETransp M_trans
  ) const
{
  TEST_FOR_EXCEPT(true);
  return false;
}

template<class Scalar>
void DefaultBlockedLinearOp<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  TEST_FOR_EXCEPT(true);
}

} // namespace Thyra

#endif	// THYRA_DEFAULT_BLOCKED_LINEAR_OP_HPP
