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

#include "Trilinos_LinearSolverSetupFailure.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <type_traits>

namespace { // (anonymous)

  void throws_LSSF () {
    throw Trilinos::LinearSolverSetupFailure ("Not an std::runtime_error");
  }

  void does_not_throw () {
  }

  void throws_runtime_error () {
    throw std::runtime_error ("Not a LinearSolverSetupFailure");
  }

  void catches_LSSF_0 (bool& success, Teuchos::FancyOStream& out) {
    out << "Test that LinearSolverSetupFailure has a what() method that "
      "returns something convertible to std::string, and test that the "
      "resulting message is as expected" << std::endl;
    Teuchos::OSTab tab1 (out);
    bool threw = false;
    try {
      throws_LSSF ();
    }
    catch (Trilinos::LinearSolverSetupFailure& e) {
      threw = true;
      const std::string msg (e.what ());
      TEST_ASSERT( msg == "Not an std::runtime_error" );
    }
    TEST_ASSERT( threw );
  }

  void catches_LSSF_1 (bool& success, Teuchos::FancyOStream& out) {
    out << "Test that LinearSolverSetupFailure is a subclass of "
      "std::runtime_error" << std::endl;
    Teuchos::OSTab tab1 (out);

    bool first_threw = false;
    try {
      throws_LSSF ();
    }
    catch (std::runtime_error& e) {
      first_threw = true;
      const std::string msg (e.what ());
      TEST_ASSERT(msg == "Not an std::runtime_error");
    }
    TEST_ASSERT( first_threw );
    static_assert (std::has_virtual_destructor<Trilinos::LinearSolverSetupFailure>::value,
      "LinearSolverSetupFailure does not have a virtual destructor.  "
      "This is bad for memory safety, because LinearSolverSetupFailure "
      "inherits from a class with virtual methods.");

    bool second_threw = false;
    try {
      throws_runtime_error ();
    }
    catch (std::runtime_error& e) {
      second_threw = true;
      const std::string msg (e.what ());
      TEST_ASSERT(msg == "Not a LinearSolverSetupFailure");
    }
    TEST_ASSERT( second_threw );
  }

  TEUCHOS_UNIT_TEST( LinearSolverSetupFailure, Test0 )
  {
    using std::endl;

    out << "LinearSolverSetupFailure Test0" << endl;
    Teuchos::OSTab tab1 (out);
    TEST_NOTHROW( does_not_throw () ); // sanity-check TEST_NOTHROW macro itself
    TEST_NOTHROW( catches_LSSF_0 (success, out) );
    TEST_NOTHROW( catches_LSSF_1 (success, out) );
  }
} // (anonymous)
