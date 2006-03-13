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

// Define to use DefaultColumnwiseMultiVector instead of SerialMultiVector
//#define THYRA_VECTOR_SPACE_USE_MULTI_VECTOR_COLS

#ifndef THYRA_SERIAL_VECTOR_SPACE_STD_HPP
#define THYRA_SERIAL_VECTOR_SPACE_STD_HPP

#include "Thyra_DefaultSerialVectorSpaceDecl.hpp"
#include "Thyra_SerialVectorSpaceBase.hpp"
#include "Thyra_DefaultSerialVector.hpp"
#ifndef THYRA_VECTOR_SPACE_USE_MULTI_VECTOR_COLS
#include "Thyra_DefaultSerialMultiVector.hpp"
#endif
#include "Teuchos_TestForException.hpp"

namespace Thyra {

template<class Scalar>
DefaultSerialVectorSpace<Scalar>::DefaultSerialVectorSpace( int dim )
{
  initialize(dim);
}

template<class Scalar>
void DefaultSerialVectorSpace<Scalar>::initialize( int dim )
{
  dim_ = dim;
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultSerialVectorSpace<Scalar>::description() const
{
  return (std::string("DefaultSerialVectorSpace<") + Teuchos::ScalarTraits<Scalar>::name() + std::string(">"));
}

// Overridden from VectorSpaceBase

template<class Scalar>
Index DefaultSerialVectorSpace<Scalar>::dim() const
{
  return dim_;
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
DefaultSerialVectorSpace<Scalar>::createMember() const
{
  return Teuchos::rcp(new DefaultSerialVector<Scalar>(Teuchos::rcp(this,false)));
}

template<class Scalar>
Teuchos::RefCountPtr< MultiVectorBase<Scalar> >
DefaultSerialVectorSpace<Scalar>::createMembers(int numMembers) const
{
#ifndef THYRA_VECTOR_SPACE_USE_MULTI_VECTOR_COLS
  return Teuchos::rcp(
    new DefaultSerialMultiVector<Scalar>(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(this->smallVecSpcFcty()->createVecSpc(numMembers),true)
      )
    );
#else
  return VectorSpaceBase<Scalar>::createMembers(numMembers);
#endif
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultSerialVectorSpace<Scalar>::clone() const
{
  return Teuchos::rcp(new DefaultSerialVectorSpace<Scalar>(*this));
}

} // end namespace Thyra

#endif // THYRA_SERIAL_VECTOR_SPACE_STD_HPP
