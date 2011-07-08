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

#include "Teuchos_TimeMonitor.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


void func_time_monitor1()
{
  TEUCHOS_FUNC_TIME_MONITOR("FUNC_TIME_MONITOR1");
}


void func_time_monitor2()
{
  TEUCHOS_FUNC_TIME_MONITOR("FUNC_TIME_MONITOR2");
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("FUNC_TIME_MONITOR2_inner", inner);
  }
}


} // namespace


namespace Teuchos {


TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR  )
{
  func_time_monitor1();
  std::ostringstream oss;
  TimeMonitor::summarize(oss);
  out << oss.str() << "\n";
  const size_t substr_i = oss.str().find("FUNC_TIME_MONITOR1");
  TEST_INEQUALITY(substr_i, std::string::npos);
}


TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR_tested  )
{
  func_time_monitor2();
  std::ostringstream oss;
  TimeMonitor::summarize(oss);
  out << oss.str() << "\n";
  const size_t substr_i = oss.str().find("FUNC_TIME_MONITOR2");
  TEST_INEQUALITY(substr_i, std::string::npos);
  const size_t substr_inner_i = oss.str().find("FUNC_TIME_MONITOR2_inner");
  TEST_INEQUALITY(substr_inner_i, std::string::npos);
}


} // namespace Teuchos
