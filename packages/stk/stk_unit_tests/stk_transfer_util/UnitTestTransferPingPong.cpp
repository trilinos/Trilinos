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
#include "stk_transfer_util/TransferMainBroker.hpp"
#include "stk_search_util/MasterElementProviderIntrepid2.hpp"
#include "stk_unit_test_utils/FieldEvaluator.hpp"
#include "stk_unit_test_utils/meshCreationHelpers.hpp"
#include "stk_unit_test_utils/TransferPingPong.hpp"
#include "stk_unit_test_utils/stk_mesh_fixtures/Hex20Fixture.hpp"
#include "stk_io/WriteMesh.hpp"

namespace {

void write_hex20_mesh_file(MPI_Comm comm, size_t nx, size_t ny, size_t nz,
                           const std::string& exoFileName)
{
  stk::mesh::fixtures::Hex20Fixture fixture(comm, nx, ny, nz);
  fixture.generate_mesh();
  stk::unit_test_util::scale_to_unit_bbox(fixture.m_bulk_data);
  stk::io::write_mesh(exoFileName, fixture.m_bulk_data);
}

TEST(TestTransferPingPong, intrepid2_linear_checkError)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  constexpr unsigned spatialDim = 3;
  stk::unit_test_util::LinearFieldEvaluator fieldEvalLinear(spatialDim);
  const int numSteps = 10;

  const double maxErr =
    stk::unit_test_util::run_ping_pong_transfer(comm, numSteps, fieldEvalLinear,
            std::make_shared<stk::search::MasterElementProviderIntrepid2>());

  EXPECT_GT(1.e-9, maxErr);
  std::cout<<"maxErr = "<<maxErr<<std::endl;
}

TEST(TestTransferPingPong, intrepid2_exponential_default_hex8_checkError)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  constexpr unsigned spatialDim = 3;
  stk::unit_test_util::ExponentialFieldEvaluator fieldEvalExp(spatialDim);
  const int numSteps = 10;

  const double maxErr =
    stk::unit_test_util::run_ping_pong_transfer(comm, numSteps, fieldEvalExp,
            std::make_shared<stk::search::MasterElementProviderIntrepid2>());

  EXPECT_GT(1.9, maxErr);
  std::cout<<"maxErr = "<<maxErr<<std::endl;
}

TEST(TestTransferPingPong, intrepid2_exponential_hex20_checkError)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMeshName("hex20_3x3x3.exo");
  std::string recvMeshName("hex20_5x5x5.exo");
  write_hex20_mesh_file(comm, 3, 3, 3, sendMeshName);
  write_hex20_mesh_file(comm, 5, 5, 5, recvMeshName);

  constexpr unsigned spatialDim = 3;
  stk::unit_test_util::ExponentialFieldEvaluator fieldEvalExp(spatialDim);
  const int numSteps = 10;

  const double maxErr =
    stk::unit_test_util::run_ping_pong_transfer(comm, numSteps, fieldEvalExp,
            std::make_shared<stk::search::MasterElementProviderIntrepid2>(),
            sendMeshName, recvMeshName);

  EXPECT_GT(0.05, maxErr);
  std::cout<<"maxErr = "<<maxErr<<std::endl;
}

}

