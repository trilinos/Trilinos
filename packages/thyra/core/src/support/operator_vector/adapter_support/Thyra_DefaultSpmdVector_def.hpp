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
