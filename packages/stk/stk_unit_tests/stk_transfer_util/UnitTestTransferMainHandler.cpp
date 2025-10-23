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

#include "gtest/gtest.h"
#include "stk_transfer_util/TransferMainHandler.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/GeneratedMeshToFile.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include <stk_util/parallel/OutputStreams.hpp>

namespace {

void build_serial_mesh(std::string filename) {
  stk::unit_test_util::generated_mesh_to_file_in_serial("1x1x4", filename);
}

TEST(TransferMainHandler, basic)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string fromMesh = "source.exo";
  std::string toMesh = "target.exo";

  build_serial_mesh(fromMesh);
  build_serial_mesh(toMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--from-mesh", fromMesh, "--to-mesh", toMesh});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerBasic(comm, args.argc(), args.argv());
  handlerBasic.run();

  EXPECT_EQ(handlerBasic.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  unlink(fromMesh.c_str());
  unlink(toMesh.c_str());
}

TEST(TransferMainHandler, all)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string fromMesh = "source.exo";
  std::string toMesh = "target.exo";

  build_serial_mesh(fromMesh);
  build_serial_mesh(toMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--from-mesh", fromMesh, "--to-mesh", toMesh,
                                 "--field-list", "field_1", "--extrapolate-option", "PROJECT" });
  
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerAll(comm, args.argc(), args.argv());
  handlerAll.run();

  EXPECT_EQ(handlerAll.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  unlink(fromMesh.c_str());
  unlink(toMesh.c_str());
}

TEST(TransferMainHandler, noArgs)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer"});

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerNoArgs(comm, args.argc(), args.argv());
  handlerNoArgs.run();
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerNoArgs.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);
}

TEST(TransferMainHandler, missingRequiredFromMesh)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string toMesh = "target.exo";
  build_serial_mesh(toMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--to-mesh", toMesh});

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerNoFrom(comm, args.argc(), args.argv());
  handlerNoFrom.run();
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerNoFrom.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);

  unlink(toMesh.c_str());
}

TEST(TransferMainHandler, missingRequiredToMesh)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string fromMesh = "source.exo";
  build_serial_mesh(fromMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--from-mesh", fromMesh});

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerNoTo(comm, args.argc(), args.argv());
  handlerNoTo.run();
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerNoTo.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);

  unlink(fromMesh.c_str());
}

TEST(TransferMainHandler, invalidExtrapolateOption)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string fromMesh = "source.exo";
  std::string toMesh = "target.exo";

  build_serial_mesh(fromMesh);
  build_serial_mesh(toMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--from-mesh", fromMesh, "--to-mesh", toMesh,
                                 "--field-list", "field_1", "--extrapolate_option", "FAKE" });

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerAll(comm, args.argc(), args.argv());
  handlerAll.run(); 
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerAll.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);

  unlink(fromMesh.c_str());
  unlink(toMesh.c_str());
}

TEST(TransferMainHandler, sameFromAndToMesh)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string fromMesh = "source.exo";
  build_serial_mesh(fromMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--from-mesh", fromMesh, "--to-mesh", fromMesh});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerSameFromTo(comm, args.argc(), args.argv());
  handlerSameFromTo.run();

  EXPECT_EQ(handlerSameFromTo.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  unlink(fromMesh.c_str());
}

TEST(TransferMainHandler, help)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--help"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerHelp(comm, args.argc(), args.argv());
  handlerHelp.run();

  EXPECT_EQ(handlerHelp.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ONLY);
}

TEST(TransferMainHandler, version)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--version"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerVersion(comm, args.argc(), args.argv());
  handlerVersion.run();

  EXPECT_EQ(handlerVersion.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ONLY);
}

}
