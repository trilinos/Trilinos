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

#ifndef THYRA_VECTOR_SERIAL_STD_HPP
#define THYRA_VECTOR_SERIAL_STD_HPP

#include "Thyra_DefaultSerialVectorDecl.hpp"
#include "Thyra_SerialVectorBase.hpp"
#include "Thyra_DefaultSerialVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Thyra {

template<class Scalar>
DefaultSerialVector<Scalar>::DefaultSerialVector(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc
  )
{
  this->initialize(vecSpc);
}

template<class Scalar>
DefaultSerialVector<Scalar>::DefaultSerialVector(
  const Index dim
  )
{
  this->initialize(dim);
}

template<class Scalar>
DefaultSerialVector<Scalar>::DefaultSerialVector(
  const Teuchos::RefCountPtr<Scalar>                          &v
  ,const Index                                                vs
  ,const Index                                                dim
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc
  )
{
  this->initialize(v,vs,dim,vecSpc);
}

template<class Scalar>
void DefaultSerialVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc
  )
{
  const Index dim = vecSpc->dim();
  this->initialize(
    Teuchos::rcp( new Scalar[dim], Teuchos::DeallocArrayDelete<Scalar>(), true )
    ,1
    ,dim
    ,vecSpc
    );
}

template<class Scalar>
void DefaultSerialVector<Scalar>::initialize(
  const Index dim
  )
{
  this->initialize(
    Teuchos::rcp( new Scalar[dim], Teuchos::DeallocArrayDelete<Scalar>(), true )
    ,1
    ,dim
    );
}

template<class Scalar>
void DefaultSerialVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<Scalar>                          &v
  ,const Index                                                vs
  ,const Index                                                dim
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vecSpc
  )
{
  if(vecSpc.get()) {
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION( vecSpc.get()!=NULL && dim != vecSpc->dim(), std::invalid_argument, "DefaultSerialVector<Scalar>::initialize(...): Error!" );
#endif
    space_serial_ = vecSpc;
  }
  else {
    space_serial_ = Teuchos::rcp(new DefaultSerialVectorSpace<Scalar>(dim));
  }
  v_       = v;
  vs_      = vs;
  dim_     = dim;
}

// Overridden from SerialVectorBase

template<class Scalar>
void DefaultSerialVector<Scalar>::getData( Scalar** values, Index* stride )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(values==NULL || stride==NULL);
#endif
  *values = v_.get();
  *stride = vs_;
}

template<class Scalar>
void DefaultSerialVector<Scalar>::commitData( Scalar** values )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( values==NULL || *values==NULL );
#endif
  // There is nothing to commit, the client had direct pointer access to internal data!
  *values = NULL;
}

template<class Scalar>
void DefaultSerialVector<Scalar>::getData( const Scalar** values, Index* stride ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(values==NULL || stride==NULL);
#endif
  *values = v_.get();
  *stride = vs_;
}

template<class Scalar>
void DefaultSerialVector<Scalar>::freeData( const Scalar** values ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( values==NULL || *values==NULL );
#endif
  // There is nothing to free!
  *values = NULL;
}

// Overridden from VectorBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultSerialVector<Scalar>::space() const
{
  return space_serial_;
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultSerialVector<Scalar>::description() const
{
  return (std::string("Thyra::DefaultSerialVector<") + Teuchos::ScalarTraits<Scalar>::name() + std::string(">"));
}

} // end namespace Thyra

#endif // THYRA_VECTOR_SERIAL_STD_HPP
