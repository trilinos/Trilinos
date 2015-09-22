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


#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_dyn_cast.hpp"


namespace {


TEUCHOS_UNIT_TEST( UnitTest, verbose )
{
  // This test checks to see that the 'verbose' bool is returned correctly.
  // This test uses knowlege of the internals of
  // Teuchos::UnitTestRepository::runUnitTests(...) to determine if this is
  // set correctly according to the --details input option.  This test is
  // *very* closely tied to the interneral implemetation of
  // Teuchos::UnitTestRepository.
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::FancyOStream;
  bool verbose = Teuchos::UnitTestRepository::verboseUnitTests();
  const Teuchos::RCP<std::ostream> wrappedOut = out.getOStream();
  if (verbose) {
    TEST_THROW(dyn_cast<std::ostringstream>(*wrappedOut), Teuchos::m_bad_cast);
  }
  else {
    TEST_NOTHROW(dyn_cast<std::ostringstream>(*wrappedOut));
  }
}


} // namespace
