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

#ifndef THYRA_SERIAL_MULTI_VECTOR_STD_HPP
#define THYRA_SERIAL_MULTI_VECTOR_STD_HPP

// Define to make some verbose output
//#define THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT

#include "Thyra_SerialMultiVectorStdDecl.hpp"
#include "Thyra_SerialMultiVectorBase.hpp"
#include "Thyra_SerialVectorStd.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
SerialMultiVectorStd<Scalar>::SerialMultiVectorStd()
  :leadingDim_(0)
{}

template<class Scalar>
SerialMultiVectorStd<Scalar>::SerialMultiVectorStd(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  )
{
  initialize(range,domain);
}

template<class Scalar>
SerialMultiVectorStd<Scalar>::SerialMultiVectorStd(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  ,const Teuchos::RefCountPtr<Scalar>                                    &values
  ,const Index                                                           leadingDim
  )
{
  initialize(range,domain,values,leadingDim);
}

template<class Scalar>
void SerialMultiVectorStd<Scalar>::initialize(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
#endif
  const Index numRows = range->dim(), numCols = domain->dim();
  initialize(
    range, domain
    ,Teuchos::rcp(new Scalar[numRows*numCols],Teuchos::DeallocArrayDelete<Scalar>(),true)
    ,numRows
    );
}

template<class Scalar>
void SerialMultiVectorStd<Scalar>::initialize(
  const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
  ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
  ,const Teuchos::RefCountPtr<Scalar>                                    &values
  ,const Index                                                           leadingDim
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
  TEST_FOR_EXCEPT(values.get()==NULL);
  TEST_FOR_EXCEPT(leadingDim < range->dim());
#endif
  range_       = range;
  domain_      = domain;
  values_      = values;
  leadingDim_  = leadingDim;
  numRows_     = range->dim();
  numCols_     = domain->dim();
  this->updateSpace();
}

template<class Scalar>
void SerialMultiVectorStd<Scalar>::uninitialize(
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >        *range
  ,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >       *domain
  ,Teuchos::RefCountPtr<Scalar>                                         *values
  ,Index                                                                *leadingDim
  )
{
  if(range)         *range         = range_;
  if(domain)        *domain        = domain_;
  if(values)        *values        = values_;
  if(leadingDim)    *leadingDim    = leadingDim_;

  range_       = Teuchos::null;
  domain_      = Teuchos::null;
  values_      = Teuchos::null;
  leadingDim_  = 0;

  this->updateSpace();
}

// Overridden from Teuchos::Describable

template<class Scalar>
std::string SerialMultiVectorStd<Scalar>::description() const
{
  return (std::string("SerialMultiVectorStd<") + Teuchos::ScalarTraits<Scalar>::name() + std::string(">"));
}

// Overridden from EuclideanLinearOpBase

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
SerialMultiVectorStd<Scalar>::rangeScalarProdVecSpc() const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::range() const called!\n";
#endif
  return range_;
}

template<class Scalar>
Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> >
SerialMultiVectorStd<Scalar>::domainScalarProdVecSpc() const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::domain() const called!\n";
#endif
  return domain_;
}

// Overridden from MultiVectorBase

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
SerialMultiVectorStd<Scalar>::col(Index j)
{
  using Teuchos::rcp;
#ifdef _DEBUG
  TEST_FOR_EXCEPT( j < 1 || numCols_ < j );
#endif
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::col() called!\n";
#endif
  return rcp(new SerialVectorStd<Scalar>(rcp((&*values_)+(j-1)*leadingDim_,false),1,numRows_,range_));
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
SerialMultiVectorStd<Scalar>::subView( const Range1D& col_rng_in )
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::subView() called!\n";
#endif
  const Range1D colRng = this->validateColRange(col_rng_in);
  return Teuchos::rcp(
    new SerialMultiVectorStd<Scalar>(
      range_
      ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(range_->smallVecSpcFcty()->createVecSpc(colRng.size()),true)
      ,Teuchos::rcp( (&*values_) + (colRng.lbound()-1)*leadingDim_, false )
      ,leadingDim_
      )
    );
}

// Overridden from MPIMultiVectorBase

template<class Scalar>
void SerialMultiVectorStd<Scalar>::getData( const Scalar **values, Index *leadingDim ) const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::getData() const called!\n";
#endif
  *values     = &*values_;
  *leadingDim = leadingDim_;
}

template<class Scalar>
void SerialMultiVectorStd<Scalar>::freeData( const Scalar *values ) const
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::freeData() called!\n";
#endif
}

template<class Scalar>
void SerialMultiVectorStd<Scalar>::getData( Scalar **values, Index *leadingDim )
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::getData() called!\n";
#endif
  *values     = &*values_;
  *leadingDim = leadingDim_;
}

template<class Scalar>
void SerialMultiVectorStd<Scalar>::commitData( Scalar *values )
{
#ifdef THYRA_SERIAL_MULTI_VECTOR_STD_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nSerialMultiVectorStd<Scalar>::commitData() called!\n";
#endif
}

} // end namespace Thyra

#endif // THYRA_SERIAL_MULTI_VECTOR_STD_HPP
