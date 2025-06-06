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

#include <gtest/gtest.h>
#include "stk_unit_test_utils/getOption.h"            // for get_command_lin...
#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"            // for ThrowRequireMsg
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_transfer_util/TransferMainOptions.hpp"

namespace {

TEST(TransferMainOptions, basic)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"exe", "--from-mesh", "foo.exo", "--to-mesh", "bar.exo"});
  stk::transfer_util::TransferMainOptions options(comm, args.argc(), args.argv());

  EXPECT_EQ("foo.exo", options.get_fromMesh_filename());
  EXPECT_EQ("bar.exo", options.get_toMesh_filename());
}

TEST(TransferMainOptions, help)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"exe", "--help"});
  stk::transfer_util::TransferMainOptions options(comm, args.argc(), args.argv());

  EXPECT_EQ(stk::CommandLineParser::ParseHelpOnly, options.get_parse_state());
  std::ostringstream os;
  options.print_usage_help(os);
  std::string usageHelp = os.str();
  bool found = usageHelp.find("--from-mesh") != std::string::npos;
  EXPECT_TRUE(found);
  found = usageHelp.find("--to-mesh") != std::string::npos;
  EXPECT_TRUE(found);
}

}
