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

//BEGINFilenameSubstitution
#include "gtest/gtest.h"
#include "stk_util/environment/EnvData.hpp"         // for EnvData
#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/FileUtils.hpp"       // for filename_substitution
#include "stk_util/environment/ParsedOptions.hpp"   // for ParsedOptions
#include "stk_util/environment/ProgramOptions.hpp"  // for get_parsed_options
#include <string>                                   // for allocator, operator+, string, char_tr...

namespace
{
TEST(StkUtilHowTo, useFilenameSubstitutionWithNoCommandLineOptions)
{
  const std::string default_base_filename = "stdin";
  const int numProcs = stk::parallel_machine_size(sierra::Env::parallel_comm());
  const std::string numProcsString = std::to_string(numProcs);
  const std::string expected_filename = default_base_filename + "-" + numProcsString + ".e";

  std::string file_name = "%B-%P.e";
  stk::util::filename_substitution(file_name);
  EXPECT_EQ(expected_filename, file_name);
}

void setFilenameInCommandLineOptions(const std::string &filename)
{
  stk::get_parsed_options().insert("input-deck", filename);
  stk::EnvData::instance().m_inputFile = filename;
}

TEST(StkUtilHowTo, useFilenameSubstitutionWithFileComingFromCommandLineOptions)
{
  const std::string base_filename = "myfile";
  const std::string full_filename = "/path/to/" + base_filename + ".g";
  setFilenameInCommandLineOptions(full_filename);

  const int numProcs = stk::parallel_machine_size(sierra::Env::parallel_comm());
  const std::string numProcsString = std::to_string(numProcs);
  const std::string expected_filename = base_filename + "-" + numProcsString + ".e";

  std::string file_name = "%B-%P.e";
  stk::util::filename_substitution(file_name);

  EXPECT_EQ(expected_filename, file_name);
}
}
//ENDFilenameSubstitution
