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
#include "Teuchos_TimeMonitor.hpp"


namespace {


int someFunction()
{
  TEUCHOS_FUNC_TIME_MONITOR("someFunction1");
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("someFunction2", diff);
  return 5;
}


TEUCHOS_UNIT_TEST( TimeMonitor, someFunction_timed )
{
  const int rtn = someFunction();
  TEST_EQUALITY(rtn, 5);
}

TEUCHOS_UNIT_TEST( TimeMonitor, formatting)
{
    using namespace Teuchos;

    // Zero out previous timers
    TimeMonitor::zeroOutTimers();

    std::ostringstream out1, out2;
    const std::string filter = "";

    // Check std::fixed formatting
    out1 << std::fixed;
    TimeMonitor::report(out1, filter);
    bool test1 = (out1.str().find("someFunction1    0.0000 (0)") != std::string::npos);
    TEST_EQUALITY(test1, true);

    // Check std::scientific formatting
    out2 << std::scientific;
    TimeMonitor::report(out2, filter);
    bool test2 = (out2.str().find("someFunction1    0.0000e+00 (0)") != std::string::npos);
    TEST_EQUALITY(test2, true);
}


} // namespace
