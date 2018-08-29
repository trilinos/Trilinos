// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <iostream>                     // for operator<<, ostringstream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <gtest/gtest.h>
#include <string>                       // for operator==, basic_string, etc


static std::ostringstream s_os;

void my_report_handler(const char* message, int type) {
    s_os << type<<": "<<message << "; ";
}

namespace {

TEST(UnitTestRuntimeWarning, Throttle)
{
  stk::set_report_handler(my_report_handler);

  int throttle_limit = 3;
  stk::MessageCode id(throttle_limit);
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XX" << std::endl;
  stk::RuntimeWarningAdHoc(id) << "runtime-warning XY" << std::endl;

  std::string warningstring = s_os.str();
  std::string expected("0: runtime-warning XX\n; 0: runtime-warning XX\n; 0: runtime-warning XX\n; -2147483648: Maximum count for this unknown (previous message) has been exceeded and will no longer be displayed; ");
  EXPECT_EQ(expected, warningstring);
}

}

