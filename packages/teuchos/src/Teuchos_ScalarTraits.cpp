// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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
