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

#ifndef THYRA_DEFUALT_LINEAR_OP_SOURCE_HPP
#define THYRA_DEFUALT_LINEAR_OP_SOURCE_HPP


#include "Thyra_DefaultLinearOpSource_decl.hpp"
#include "Thyra_LinearOpBase.hpp"


namespace Thyra {


// Constructors/initializers/accessors


template <class Scalar>
DefaultLinearOpSource<Scalar>::DefaultLinearOpSource()
{}


template <class Scalar>
DefaultLinearOpSource<Scalar>::DefaultLinearOpSource(
  const Teuchos::RCP<LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
DefaultLinearOpSource<Scalar>::DefaultLinearOpSource(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
void DefaultLinearOpSource<Scalar>::initialize(
  const Teuchos::RCP<LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
void DefaultLinearOpSource<Scalar>::initialize(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &op
  )
{
  op_.initialize(op);
}


template <class Scalar>
void DefaultLinearOpSource<Scalar>::uninitialize()
{
  op_.uninitialize();
}


// Overridden from LinearOpSourceBase


template <class Scalar>
bool DefaultLinearOpSource<Scalar>::isOpConst() const
{
  return op_.isConst();
}


template <class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultLinearOpSource<Scalar>::getNonconstOp()
{
  return op_.getNonconstObj();
}


template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultLinearOpSource<Scalar>::getOp() const
{
  return op_.getConstObj();
}


} // namespace Thyra


#endif // THYRA_DEFUALT_LINEAR_OP_SOURCE_HPP
