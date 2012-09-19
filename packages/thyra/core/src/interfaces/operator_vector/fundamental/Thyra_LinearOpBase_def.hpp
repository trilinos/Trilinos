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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
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
