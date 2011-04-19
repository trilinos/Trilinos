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
