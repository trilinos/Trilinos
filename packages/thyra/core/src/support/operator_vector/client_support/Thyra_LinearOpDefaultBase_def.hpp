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

#ifndef THYRA_LINEAR_OP_DEFAULT_BASE_HPP
#define THYRA_LINEAR_OP_DEFAULT_BASE_HPP

#include "Thyra_LinearOpDefaultBase_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_describeLinearOp.hpp"

namespace Thyra {


// Overridden from Teuchos::Describable


template<class Scalar>
std::string LinearOpDefaultBase<Scalar>::description() const
{
  std::ostringstream oss;
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >
    l_range = this->range();
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >
    l_domain = this->domain();
  oss << Teuchos::Describable::description();
  if(!l_range.get()) {
    oss << "{range=NULL,domain=NULL}"; 
  }
  else {
    const Ordinal dimDomain = l_domain->dim(), dimRange = l_range->dim();
    oss
      << "{rangeDim=" << dimRange
      << ",domainDim=" << dimDomain << "}";
  }
  return oss.str();
}


template<class Scalar>
void LinearOpDefaultBase<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  describeLinearOp(*this, out, verbLevel);
}


}	// end namespace Thyra

#endif // THYRA_LINEAR_OP_DEFAULT_BASE_HPP
