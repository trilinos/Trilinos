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

//-BEGIN
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_util/environment/EnvData.hpp>  // for EnvData
#include <stk_util/environment/FileUtils.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <string>                       // for string, allocator, etc
#include <utility>                      // for make_pair
#include "boost/program_options/variables_map.hpp"  // for variable_value, etc

namespace
{
  TEST(StkUtilHowTo, useFilenameSubstitutionWithNoCommandLineOptions)
  {
    const std::string default_base_filename = "stdin";
    const std::string numProcsString = "1";
    const std::string expected_filename = default_base_filename + "-" + numProcsString + ".e";

    std::string file_name = "%B-%P.e";
    stk::util::filename_substitution(file_name);
    EXPECT_EQ(expected_filename, file_name);
  }

  void setFilenameInCommandLineOptions(const std::string &filename)
  {
    boost::program_options::variables_map &command_line_options = stk::get_variables_map();
    command_line_options.insert(std::make_pair("input-deck", boost::program_options::variable_value(filename, false)));
    stk::EnvData::instance().m_inputFile = filename;
  }
  TEST(StkUtilHowTo, useFilenameSubstitutionWithFileComingFromCommandLineOptions)
  {
    const std::string base_filename = "myfile";
    const std::string full_filename = "/path/to/" + base_filename + ".g";
    setFilenameInCommandLineOptions(full_filename);

    const std::string numProcsString = "1";
    const std::string expected_filename = base_filename + "-" + numProcsString + ".e";

    std::string file_name = "%B-%P.e";
    stk::util::filename_substitution(file_name);

    EXPECT_EQ(expected_filename, file_name);
  }
}
//-END
