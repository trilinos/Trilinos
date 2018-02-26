/*
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
*/

#include "Teuchos_UnitTestHarness.hpp"

#include "Teuchos_PrintDouble.hpp"

#include <cstdlib>
#include <sstream>

namespace {

static void test_print_double(const char* s, bool& success, std::ostream& out) {
  auto v = std::atof(s);
  std::stringstream ss;
  Teuchos::print_double(ss, v);
  auto s2 = ss.str();
  TEST_EQUALITY_CONST( s2, s );
}

TEUCHOS_UNIT_TEST( PrintDouble, Basic )
{
  test_print_double("1.", success, out);
  test_print_double("-1.", success, out);
  test_print_double("0.", success, out);
  test_print_double("4.", success, out);
  test_print_double("55.", success, out);
  test_print_double("3.14159", success, out);
  test_print_double("5404319552844595.", success, out);
  test_print_double("1e300", success, out);
  test_print_double("1e2", success, out);
  test_print_double("10.", success, out);
  test_print_double(".003", success, out);
  /* in order to ensure this gets converted back to the right
     value even if people change their IEEE rounding mode,
     this long representation is needed */
  test_print_double("9.999999999999999e22", success, out);
  test_print_double("inf", success, out);
  test_print_double("-inf", success, out);
  test_print_double("nan", success, out);
}

} // namespace
