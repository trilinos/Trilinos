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
#include "stk_io/FillMesh.hpp"
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

#include <iostream>                                        // for operator<<
#include <memory>                                          // for shared_ptr
#include <stdexcept>                                       // for runtime_error
#include <string>                                          // for string
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace {
//BEGIN_supporting_functions
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

void test_expected_value_at_location(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field,
                                     unsigned fieldDimension, unsigned numCopies, stk::mesh::Entity& entity,
                                     const stk::unit_test_util::FieldEvaluator& eval, const double* location,
                                     double tolerance = 1.0e-6)
{
  const unsigned nDim = bulk.mesh_meta_data().spatial_dimension();
  const double* data = static_cast<const double *>(stk::mesh::field_data(*field, entity));
  for(unsigned i = 0; i < numCopies; ++i) {
    double x = location[nDim * i + 0];
    double y = location[nDim * i + 1];
    double z = (nDim == 2 ? 0.0 : location[nDim * i + 2]);

    double expectedValue = eval(entity, x, y, z);

    for(unsigned j = 0; j < fieldDimension; ++j) {
      EXPECT_NEAR(expectedValue, data[fieldDimension * i + j], tolerance)
          << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entity)
          << " at location (" << x << "," << y << "," << z << ")";
    }
  }
}

void test_node_transfer(const stk::mesh::BulkData& recvBulk,
                        const std::string& recvPartName, const std::string& recvFieldName,
                        stk::unit_test_util::FieldEvaluator& eval, stk::mesh::EntityRank recvRank)
{
  const stk::mesh::MetaData& recvMeta = recvBulk.mesh_meta_data();
  const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();

  const stk::mesh::FieldBase* recvField = recvMeta.get_field(stk::topology::NODE_RANK, recvFieldName);
  ASSERT_TRUE(nullptr != recvField);

  const stk::mesh::Part* part = recvMeta.get_part(recvPartName);
  ASSERT_TRUE(nullptr != part);

  stk::mesh::Selector selector(*part);
  stk::mesh::EntityVector entities;

  stk::mesh::get_selected_entities(selector, recvBulk.buckets(recvRank), entities);
  unsigned numCopies = 1;

  for(stk::mesh::Entity entity : entities) {
    const double* location = static_cast<const double *>(stk::mesh::field_data(*recvCoordinateField, entity));
    test_expected_value_at_location(recvBulk, recvField, recvMeta.spatial_dimension(), numCopies, entity, eval, location);
  }
}
//END_supporting_functions

//BEGIN_simple_transfer
TEST(StkTransferHowTo, useSimpleTransfer)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nProc = stk::parallel_machine_size(comm);

  if ((nProc != 1) && (nProc != 2)) {
    GTEST_SKIP();
  }

  unsigned spatialDim(3u);
  stk::unit_test_util::LinearFieldEvaluator eval(spatialDim);

  const std::string sendFieldName{"field1"};
  const std::string recvFieldName{"field1"};

  const std::string sendMeshSpec{"generated:2x2x2|bbox:0,0,0,1,1,1"};
  const std::string recvMeshSpec{"generated:4x4x4|bbox:0,0,0,1,1,1"};

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  auto extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;
  auto masterElemProvider = std::make_shared<stk::unit_test_util::MasterElementProvider>();

  unsigned numSendCopies = 1;
  unsigned numRecvCopies = 1;

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

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", comm);
  transfer.add_send_field(*sendField);
  transfer.add_recv_field(*recvField);

  std::string transferPartName{"block_1"};
  transfer.add_send_part_name(transferPartName);
  transfer.add_recv_part_name(transferPartName);

  transfer.setup_master_element_transfer(*sendBulk, *recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         masterElemProvider, extrapolateOption);

  transfer.initialize();
  transfer.apply();

  test_node_transfer(*recvBulk, transferPartName, recvFieldName, eval, recvRank);
}
//END_simple_transfer

}
