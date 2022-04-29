// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "gtest/gtest.h"
#include "stk_util/environment/RuntimeMessage.hpp"  // for MessageCode
#include "stk_util/environment/RuntimeWarning.hpp"  // for RuntimeWarningAdHoc
#include "stk_util/util/ReportHandler.hpp"          // for set_report_handler
#include <iostream>                                 // for operator<<, endl, basic_ostream, basi...
#include <string>                                   // for char_traits, string


static std::ostringstream s_os;

void my_report_handler(const char* message, int type) {
    s_os << type<<": "<<message << "; ";
}

namespace {

TEST(UnitTestRuntimeWarning, Throttle)
{
  stk::set_report_handler(my_report_handler);
  s_os.str(std::string());
  stk::reset_warning_count();

  size_t throttleLimit = 3;
  stk::MessageCode id(throttleLimit);
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XY" << std::endl;

  std::string warningstring = s_os.str();
  std::string expected("0: runtime-warning XX\n; 0: runtime-warning XX\n; 0: runtime-warning XX\n; -2147483648: Maximum count for this unknown (previous message) has been exceeded and will no longer be displayed; ");
  EXPECT_EQ(expected, warningstring);

  EXPECT_EQ(4u, stk::get_warning_count());
  EXPECT_EQ(throttleLimit, stk::get_warning_printed_count());
  EXPECT_EQ(throttleLimit, stk::get_warning_printed_count(id));
}

TEST(UnitTestRuntimeWarning, ThrottlePrinted)
{
  stk::set_report_handler(my_report_handler);
  s_os.str(std::string());
  stk::reset_warning_count();

  size_t throttle_limit = 3;
  stk::MessageCode id(throttle_limit);
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;

  stk::RuntimeWarning() << "runtime-warning default not throttled" << std::endl;
  stk::RuntimeWarning() << "runtime-warning default not throttled" << std::endl;
  stk::RuntimeWarning() << "runtime-warning default not throttled" << std::endl;
  stk::RuntimeWarning() << "runtime-warning default not throttled" << std::endl;

  EXPECT_EQ(10u, stk::get_warning_count());

  EXPECT_EQ(throttle_limit, stk::get_warning_printed_count(id));
  EXPECT_EQ(throttle_limit + 4u, stk::get_warning_printed_count());

  stk::MessageCode defaultMessageCode = stk::MessageCode::s_defaultMessageCode;
  EXPECT_EQ(4u, stk::get_warning_printed_count(defaultMessageCode));
}

TEST(UnitTestRuntimeWarningP0, WarningPrintedCountParallel)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  stk::set_report_handler(my_report_handler);
  s_os.str(std::string());
  stk::reset_warning_count();

  stk::RuntimeWarningP0() << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningP0() << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningP0() << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningP0() << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningP0() << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningP0() << "runtime-warning XX" << std::endl;

  unsigned expectedWarningCount = (stk::parallel_machine_rank(MPI_COMM_WORLD)==0) ? 6 : 0;

  EXPECT_EQ(expectedWarningCount, stk::get_warning_count());
  EXPECT_EQ(expectedWarningCount, stk::get_warning_printed_count());
}

}

