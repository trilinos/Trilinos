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

#ifndef THYRA_DEFAULT_SPMD_VECTOR_DEF_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_DEF_HPP


#include "Thyra_DefaultSpmdVector_decl.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template<class Scalar>
DefaultSpmdVector<Scalar>::DefaultSpmdVector()
  :stride_(0)
{}


template<class Scalar>
DefaultSpmdVector<Scalar>::DefaultSpmdVector(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdSpace_in,
  const ArrayRCP<Scalar> &localValues,
  const Ordinal stride
  )
{
  initialize(spmdSpace_in, localValues, stride);
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::initialize(
  const RCP<const SpmdVectorSpaceBase<Scalar> > &spmdSpace_in
  ,const ArrayRCP<Scalar> &localValues
  ,const Ordinal stride
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(spmdSpace_in));
  TEST_FOR_EXCEPT(spmdSpace_in->localSubDim() > 0 && localValues.get()==NULL);
  TEST_FOR_EXCEPT(stride==0);
#endif
  spmdSpace_ = spmdSpace_in;
  localValues_ = localValues;
  stride_ = stride;
  this->updateSpmdSpace();
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::uninitialize(
  RCP<const SpmdVectorSpaceBase<Scalar> > *spmdSpace_in
  ,ArrayRCP<Scalar> *localValues
  ,Ordinal *stride
  )
{
  if(spmdSpace_in) *spmdSpace_in = spmdSpace_;
  if(localValues) *localValues = localValues_;
  if(stride) *stride = stride_;

  spmdSpace_ = Teuchos::null;
  localValues_ = Teuchos::null;
  stride_ = 0;

  this->updateSpmdSpace();
}


// Overridden from SpmdVectorBase


template<class Scalar>
RCP<const SpmdVectorSpaceBase<Scalar> >
DefaultSpmdVector<Scalar>::spmdSpace() const
{
  return spmdSpace_;
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::getNonconstLocalDataImpl(
  const Ptr<ArrayRCP<Scalar> > &localValues )
{
  *localValues = localValues_;
}


template<class Scalar>
void DefaultSpmdVector<Scalar>::getLocalDataImpl(
  const Ptr<ArrayRCP<const Scalar> > &localValues ) const
{
  *localValues = localValues_;
}


} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_VECTOR_DEF_HPP
