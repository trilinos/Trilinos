// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"

// Define this to throw exceptions when any Teuchos::ScalarTraits function
// encounters a NaN or an Inf.
//#define TEUCHOS_SCALAR_TRAITS_THROW_NAN_INF_ERR


namespace {

// These functions exist to trick the compiler into not returning a warning
// message for 0.0/0.0 or refusing to compile the code.  If a compiler gets
// too smart, we can put these definitions into a different *.cpp file such
// that most compilers would not be able to know at compile-time if a NaN or
// an Inf was being created.

float returnFloatZero() { return 0.0; }

double returnDoubleZero() { return 0.0; }

} // namespace


void Teuchos::throwScalarTraitsNanInfError( const std::string &errMsg )
{
#ifdef TEUCHOS_SCALAR_TRAITS_THROW_NAN_INF_ERR
  TEST_FOR_EXCEPTION( true, std::runtime_error, errMsg );
#endif
}

#ifdef HAVE_TEUCHOS_GNU_MP
gmp_randclass Teuchos::gmp_rng ( gmp_randinit_default );
#endif

#ifdef HAVE_TEUCHOS_QD
bool Teuchos::operator&&(const dd_real &a, const dd_real &b) {
  return !a.is_zero() && !b.is_zero();
}
bool Teuchos::operator&&(const qd_real &a, const qd_real &b) {
  return !a.is_zero() && !b.is_zero();
}
#endif

#ifndef __sun
// This is an intentional computation of NaN.
const float  Teuchos::flt_nan = +returnFloatZero()/returnFloatZero();
const double Teuchos::dbl_nan = +returnDoubleZero()/returnDoubleZero();
#endif
