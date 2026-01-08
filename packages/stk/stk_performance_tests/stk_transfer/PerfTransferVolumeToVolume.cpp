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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_mesh/base/BulkData.hpp>                      // for BulkData
#include "stk_mesh/base/Field.hpp"                         // for Field
#include "stk_search_util/MasterElementProvider.hpp"
#include "stk_transfer_util/spmd/GeometricTransfer.hpp"
#include "stk_transfer_util/spmd/SimpleTransfer.hpp"
#include <stk_unit_test_utils/getOption.h>
#include "stk_unit_test_utils/BuildMesh.hpp"               // for build_mesh
#include "stk_unit_test_utils/FieldEvaluator.hpp"
#include "stk_unit_test_utils/UnitTestSearchUtils.hpp"
#include "stk_unit_test_utils/UnitTestTransferUtils.hpp"

#include <gtest/gtest.h>
#include "mpi.h"            // for MPI_COMM_WORLD

#include <algorithm>                                       // for max, copy
#include <cstddef>                                         // for size_t
#include <cstdint>                                         // for int64_t
#include <iostream>                                        // for operator<<
#include <memory>                                          // for shared_ptr
#include <stdexcept>                                       // for runtime_error
#include <string>                                          // for string
#include <utility>                                         // for pair
#include <vector>                                          // for vector, swap

#include <stk_unit_test_utils/timer.hpp>

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace {

namespace impl {
void setup_mesh_with_field(stk::mesh::BulkData& bulk,
                           const std::string& meshSpec,
                           const std::string& fieldName,
                           stk::mesh::EntityRank fieldRank,
                           unsigned numCopies = 1)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  unsigned spatialDim = meta.spatial_dimension();
  std::vector<double> initialValue(spatialDim, 0.0);

  stk::mesh::Field<double>& field = meta.declare_field<double>(fieldRank, fieldName);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), spatialDim, numCopies, initialValue.data());
  stk::io::fill_mesh(meshSpec, bulk);
}

void run_perf_master_element_to_gauss_point_transfer(const std::string& sendMeshSpec, const std::string& recvMeshSpec)
{
  unsigned spatialDim(3u);
  stk::unit_test_util::LinearFieldEvaluator eval(spatialDim);

  std::string sendFieldName{"field1"};
  std::string recvFieldName{"field1"};

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::search::ObjectOutsideDomainPolicy extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::transfer_util::MasterElementProvider>();

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  unsigned numSendCopies = 1;
  unsigned numRecvCopies = masterElemProvider->num_integration_points(hex8Topo);
  EXPECT_EQ(8u, numRecvCopies);

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 20;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_mesh(spatialDim, comm);
    std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(spatialDim, comm);

    stk::mesh::MetaData& sendMeta = sendBulk->mesh_meta_data();
    stk::mesh::MetaData& recvMeta = recvBulk->mesh_meta_data();

    setup_mesh_with_field(*sendBulk, sendMeshSpec, sendFieldName, sendRank, numSendCopies);
    setup_mesh_with_field(*recvBulk, recvMeshSpec, recvFieldName, recvRank, numRecvCopies);

    stk::mesh::FieldBase* sendField = sendMeta.get_field(sendRank, sendFieldName);
    stk::mesh::FieldBase* recvField = recvMeta.get_field(recvRank, recvFieldName);

    ASSERT_TRUE(nullptr != sendField);
    ASSERT_TRUE(nullptr != recvField);

    stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

    stk::transfer::spmd::SimpleTransfer transfer("TransferTest");
    transfer.add_send_field(*sendField);
    transfer.add_recv_field(*recvField);

    std::string transferPartName{"block_1"};
    transfer.add_send_part_name(transferPartName);
    transfer.add_recv_part_name(transferPartName);

    transfer.setup_master_element_transfer(*sendBulk, *recvBulk,
                                           stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                           masterElemProvider, extrapolateOption);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      transfer.initialize();
      transfer.apply();
    }
    batchTimer.stop_batch_timer();

  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

void run_perf_patch_recovery_to_gauss_point_transfer(const std::string& sendMeshSpec, const std::string& recvMeshSpec)
{
  unsigned spatialDim(3u);
  stk::unit_test_util::LinearFieldEvaluator eval(spatialDim);

  std::string sendFieldName{"field1"};
  std::string recvFieldName{"field1"};

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::search::ObjectOutsideDomainPolicy extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::transfer_util::MasterElementProvider>();

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  unsigned numSendCopies = 1;
  unsigned numRecvCopies = masterElemProvider->num_integration_points(hex8Topo);
  EXPECT_EQ(8u, numRecvCopies);

  MPI_Comm comm = MPI_COMM_WORLD;
  const unsigned NUM_RUNS = 5;
  const unsigned NUM_ITERS = 20;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();

  for (unsigned run = 0; run < NUM_RUNS; ++run) {

    std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_mesh(spatialDim, comm);
    std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(spatialDim, comm);

    stk::mesh::MetaData& sendMeta = sendBulk->mesh_meta_data();
    stk::mesh::MetaData& recvMeta = recvBulk->mesh_meta_data();

    setup_mesh_with_field(*sendBulk, sendMeshSpec, sendFieldName, sendRank, numSendCopies);
    setup_mesh_with_field(*recvBulk, recvMeshSpec, recvFieldName, recvRank, numRecvCopies);

    stk::mesh::FieldBase* sendField = sendMeta.get_field(sendRank, sendFieldName);
    stk::mesh::FieldBase* recvField = recvMeta.get_field(recvRank, recvFieldName);

    ASSERT_TRUE(nullptr != sendField);
    ASSERT_TRUE(nullptr != recvField);

    stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

    stk::transfer::spmd::SimpleTransfer transfer("TransferTest");
    transfer.add_send_field(*sendField);
    transfer.add_recv_field(*recvField);

    std::string transferPartName{"block_1"};
    transfer.add_send_part_name(transferPartName);
    transfer.add_recv_part_name(transferPartName);

    transfer.setup_patch_recovery_transfer(*sendBulk, *recvBulk,
                                           stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                           masterElemProvider, extrapolateOption);

    batchTimer.start_batch_timer();
    for (unsigned i = 0; i < NUM_ITERS; ++i) {
      transfer.initialize();
      transfer.apply();
    }
    batchTimer.stop_batch_timer();

  }
  batchTimer.print_batch_timing(NUM_ITERS);
}
}

TEST(StkTransfer_VolumeToVolume, masterElementToGaussPoint)
{
  std::string meshSpec = stk::unit_test_util::get_option("-i", "generated:40x80x20");
  impl::run_perf_master_element_to_gauss_point_transfer(meshSpec, meshSpec);
}

TEST(StkTransfer_VolumeToVolume, patchRecoveryToGaussPoint)
{
  std::string meshSpec = stk::unit_test_util::get_option("-i", "generated:20x20x10");
  impl::run_perf_patch_recovery_to_gauss_point_transfer(meshSpec, meshSpec);
}
}

