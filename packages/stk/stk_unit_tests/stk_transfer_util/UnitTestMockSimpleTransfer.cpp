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
#include <stk_io/FillMesh.hpp>                             // for fill_mesh
#include "stk_io/WriteMesh.hpp"                            // for write_mesh
#include <stk_mesh/base/BulkData.hpp>                      // for BulkData
#include "stk_mesh/base/EntityKey.hpp"                     // for EntityKey
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/base/Field.hpp"                         // for Field
#include "stk_mesh/base/FieldBLAS.hpp"                     // for field_copy
#include "stk_search_util/MasterElementProvider.hpp"
#include "stk_transfer_util/spmd/GeometricTransfer.hpp"
#include "stk_transfer_util/spmd/SimpleTransfer.hpp"
#include <stk_unit_test_utils/getOption.h>
#include "stk_unit_test_utils/BuildMesh.hpp"               // for build_mesh
#include "stk_unit_test_utils/FieldEvaluator.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"                // for get_full_t...
#include "stk_unit_test_utils/UnitTestSearchUtils.hpp"
#include "stk_unit_test_utils/UnitTestTransferUtils.hpp"

#include <gtest/gtest.h>
#include "mpi.h"            // for MPI_COMM_WORLD

#include <algorithm>                                       // for max, copy
#include <cmath>                                           // for sqrt, pow
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

double transform_scale_by_half(double value)    { return value * 0.5; }
double transform_scale_by_two(double value)     { return value * 2.0; }
double transform_scale_by_four(double value)    { return value * 4.0; }
double transform_scale_by_eight(double value)   { return value * 8.0; }

}

// TODO:
// Optimize field size retrieval

class MockSimpleTransferTest : public ::testing::Test
{
 public:

  MockSimpleTransferTest(stk::ParallelMachine comm)
  : m_comm(comm)
  { }

  MockSimpleTransferTest() = default;

  template <typename T, stk::mesh::FieldAccessTag FieldAccess, typename Alg>
  inline
  void field_data_test(const stk::mesh::FieldBase& fieldBase, const stk::mesh::Entity entity, Alg&& alg)
  {
    if (fieldBase.host_data_layout() == stk::mesh::Layout::Right) {
      auto fieldData = fieldBase.data<T, FieldAccess, stk::ngp::HostSpace, stk::mesh::Layout::Right>();
      alg(fieldData, entity);
    }
    else if (fieldBase.host_data_layout() == stk::mesh::Layout::Left) {
      auto fieldData = fieldBase.data<T, FieldAccess, stk::ngp::HostSpace, stk::mesh::Layout::Left>();
      alg(fieldData, entity);
    }
    else {
      STK_ThrowErrorMsg("Unsupported host Field data layout: " << fieldBase.host_data_layout());
    }
  }

  void test_expected_indexed_value_at_location(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field, unsigned recvIndex,
                                               stk::mesh::Entity& entity, const stk::unit_test_util::FieldEvaluator& eval,
                                               const double* loc, double tolerance = 1.0e-6)
  {
    // Coordinate data at one point e.g centroid or node
    const unsigned nDim = bulk.mesh_meta_data().spatial_dimension();

    auto alg = [&](auto& fieldData, stk::mesh::Entity entityArg) {
                 auto data = fieldData.entity_values(entityArg);

                 double x = loc[0];
                 double y = loc[1];
                 double z = (nDim != 3 ? 0.0 : loc[2]);
                 double expectedValue = eval(entity, x, y, z);

                 for (stk::mesh::CopyIdx copy : data.copies()) {
                   for(stk::mesh::ComponentIdx comp : data.components()) {
                     double value = data(copy, comp);

                     unsigned index = comp;
                     if((recvIndex == 0) || index == (recvIndex - 1)) {
                       EXPECT_NEAR(expectedValue, value, tolerance)
                           << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                           << " at location (" << x << "," << y << "," << z << ")";
                     }
                     else {
                       EXPECT_EQ(0, value)
                           << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                           << " at location (" << x << "," << y << "," << z << ")";
                     }
                   }
                 }
      };

    field_data_test<double, stk::mesh::ReadOnly>(*field, entity, alg);
  }

  void test_expected_indexed_value_at_locations(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field, unsigned recvIndex,
                                                stk::mesh::Entity& entity, const stk::unit_test_util::FieldEvaluator& eval,
                                                const double* loc, const unsigned locCompStride, const unsigned locCopyStride,
                                                double tolerance = 1.0e-6)
  {
    // Coordinate data per copy e.g gauss points
    const unsigned nDim = bulk.mesh_meta_data().spatial_dimension();

    auto alg = [&](auto& fieldData, stk::mesh::Entity entityArg) {
                 auto data = fieldData.entity_values(entityArg);

                 for (stk::mesh::CopyIdx copy : data.copies()) {
                   unsigned copyIndex = copy;
                   double x = loc[copyIndex*locCopyStride + 0*locCompStride];
                   double y = loc[copyIndex*locCopyStride + 1*locCompStride];
                   double z = (nDim != 3 ? 0.0 : loc[copyIndex*locCopyStride + 2*locCompStride]);
                   double expectedValue = eval(entityArg, x, y, z);

                   for(stk::mesh::ComponentIdx comp : data.components()) {
                     double value = data(copy, comp);

                     unsigned index = comp;
                     if((recvIndex == 0) || index == (recvIndex - 1)) {
                       EXPECT_NEAR(expectedValue, value, tolerance)
                           << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                           << " at location (" << x << "," << y << "," << z << ")";
                     }
                     else {
                       EXPECT_EQ(0, value)
                           << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                           << " at location (" << x << "," << y << "," << z << ")";
                     }
                   }
                 }
      };

    field_data_test<double, stk::mesh::ReadOnly>(*field, entity, alg);
  }

  void test_expected_value_at_location(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field,
                                       stk::mesh::Entity& entity, const stk::unit_test_util::FieldEvaluator& eval,
                                       unsigned numComponentsUsed,
                                       const double* loc, double tolerance = 1.0e-6)
  {
    // Coordinate data at one point e.g centroid or node
    const unsigned nDim = bulk.mesh_meta_data().spatial_dimension();

    auto alg = [&](auto& fieldData, stk::mesh::Entity entityArg) {
                 auto data = fieldData.entity_values(entityArg);

                 double x = loc[0];
                 double y = loc[1];
                 double z = (nDim != 3 ? 0.0 : loc[2]);
                 double expectedValue = eval(entityArg, x, y, z);

                 for (stk::mesh::CopyIdx copy : data.copies()) {
                   for (stk::mesh::ComponentIdx comp=0_comp; comp < static_cast<stk::mesh::ComponentIdx>(numComponentsUsed); ++comp) {
                     double value = data(copy, comp);

                     EXPECT_NEAR(expectedValue, value, tolerance)
                         << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                         << " at loc (" << x << "," << y << "," << z << ")";
                   }
                 }
    };

    field_data_test<double, stk::mesh::ReadOnly>(*field, entity, alg);
  }

  void test_expected_value_at_locations(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field,
                                        stk::mesh::Entity& entity, const stk::unit_test_util::FieldEvaluator& eval,
                                        unsigned numComponentsUsed,
                                        const double* loc, const unsigned locCompStride, const unsigned locCopyStride,
                                        double tolerance = 1.0e-6)
  {
    // Coordinate data per copy e.g gauss points
    const unsigned nDim = bulk.mesh_meta_data().spatial_dimension();

    auto alg = [&](auto& fieldData, stk::mesh::Entity entityArg) {
                 auto data = fieldData.entity_values(entityArg);

                 for (stk::mesh::CopyIdx copy : data.copies()) {
                   unsigned copyIndex = copy;
                   double x = loc[copyIndex*locCopyStride + 0*locCompStride];
                   double y = loc[copyIndex*locCopyStride + 1*locCompStride];
                   double z = (nDim != 3 ? 0.0 : loc[copyIndex*locCopyStride + 2*locCompStride]);
                   double expectedValue = eval(entityArg, x, y, z);

                   for (stk::mesh::ComponentIdx comp=0_comp; comp < static_cast<stk::mesh::ComponentIdx>(numComponentsUsed); ++comp) {
                     double value = data(copy, comp);

                     EXPECT_NEAR(expectedValue, value, tolerance)
                         << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                         << " at loc (" << x << "," << y << "," << z << ")";
                   }
                 }
    };

    field_data_test<double, stk::mesh::ReadOnly>(*field, entity, alg);
  }

  void setup_mesh_with_field(stk::mesh::BulkData& bulk,
                             const std::string& meshSpec,
                             const std::string& fieldName,
                             stk::mesh::EntityRank fieldRank,
                             unsigned numCopies = 1)
  {
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    unsigned spatialDim = meta.spatial_dimension();
    std::vector<double> initialValue(spatialDim*numCopies, 0.0);

    stk::mesh::Field<double>& field = meta.declare_field<double>(fieldRank, fieldName);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), spatialDim, numCopies, initialValue.data());
    stk::io::fill_mesh(meshSpec, bulk);
  }

  void internal_setup_transfer(stk::transfer::spmd::SimpleTransfer& transfer,
                               stk::mesh::EntityRank sendRank,
                               stk::mesh::EntityRank recvRank)
  {
    m_sendBulk = stk::unit_test_util::build_mesh(m_sendSpatialDim, m_comm);
    m_recvBulk = stk::unit_test_util::build_mesh(m_recvSpatialDim, m_comm);

    stk::mesh::MetaData& sendMeta = m_sendBulk->mesh_meta_data();
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();

    setup_mesh_with_field(*m_sendBulk, m_sendMeshSpec, m_sendFieldName, sendRank, m_numSendCopies);
    setup_mesh_with_field(*m_recvBulk, m_recvMeshSpec, m_recvFieldName, recvRank, m_numRecvCopies);

    m_sendField = sendMeta.get_field(sendRank, m_sendFieldName);
    m_recvField = recvMeta.get_field(recvRank, m_recvFieldName);

    ASSERT_TRUE(nullptr != m_sendField);
    ASSERT_TRUE(nullptr != m_recvField);

    transfer.add_send_field(*m_sendField);
    transfer.add_recv_field(*m_recvField);

    transfer.add_send_part_name(m_sendPartName);
    transfer.add_recv_part_name(m_recvPartName);
  }

  void setup_copies(unsigned numCopies,
                    stk::search::SearchTopology* topo = nullptr,
                    unsigned expectedNumRecvCopies = 0)
  {
    m_masterElemProvider = std::make_shared<stk::unit_test_util::MasterElementProvider>();

    m_numSendCopies = numCopies;
    m_numRecvCopies = (nullptr != topo) ? m_masterElemProvider->num_integration_points(*topo)*numCopies : numCopies;

    if(nullptr != topo) {
      EXPECT_EQ(expectedNumRecvCopies, m_numRecvCopies);
    }
  }

  void setup_transfer(stk::transfer::spmd::SimpleTransfer& transfer,
                      stk::mesh::EntityRank sendRank,
                      stk::mesh::EntityRank recvRank,
                      unsigned numCopies,
                      stk::search::SearchTopology& topo,
                      unsigned expectedNumRecvCopies)
  {
    setup_copies(numCopies, &topo, expectedNumRecvCopies);

    internal_setup_transfer(transfer, sendRank, recvRank);
  }

  void setup_transfer(stk::transfer::spmd::SimpleTransfer& transfer,
                      stk::mesh::EntityRank sendRank,
                      stk::mesh::EntityRank recvRank,
                      unsigned numCopies = 1)
  {
    setup_copies(numCopies);

    internal_setup_transfer(transfer, sendRank, recvRank);
  }

  void test_indexed_node_transfer(stk::unit_test_util::FieldEvaluator& eval,
                                  stk::mesh::EntityRank recvRank, unsigned recvIndex)
  {
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();
    const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
    const stk::mesh::Part* part = recvMeta.get_part(m_recvPartName);
    ASSERT_TRUE(nullptr != part);
    stk::mesh::Selector selector(*part);
    stk::mesh::EntityVector entities;

    stk::mesh::get_selected_entities(selector, m_recvBulk->buckets(recvRank), entities);

    std::vector<double> location;
    for(stk::mesh::Entity entity : entities) {
      stk::search::determine_centroid(m_sendSpatialDim, entity, *recvCoordinateField, location);
      test_expected_indexed_value_at_location(*m_recvBulk, m_recvField, recvIndex, entity, eval, location.data());
    }
  }

  void test_node_transfer(stk::unit_test_util::FieldEvaluator& eval,
                          stk::mesh::EntityRank recvRank,
                          double tolerance = 1.0e-6)
  {
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();
    const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
    const stk::mesh::Part* part = recvMeta.get_part(m_recvPartName);
    ASSERT_TRUE(nullptr != part);
    stk::mesh::Selector selector(*part);
    stk::mesh::EntityVector entities;

    stk::mesh::get_selected_entities(selector, m_recvBulk->buckets(recvRank), entities);

    std::vector<double> location;
    for(stk::mesh::Entity entity : entities) {
      stk::search::determine_centroid(m_recvSpatialDim, entity, *recvCoordinateField, location);
      test_expected_value_at_location(*m_recvBulk, m_recvField, entity, eval, std::min(m_sendSpatialDim, m_recvSpatialDim), location.data(), tolerance);
    }
  }

  void test_indexed_centroid_transfer(stk::unit_test_util::FieldEvaluator& eval,
                                      stk::mesh::EntityRank recvRank, unsigned recvIndex)
  {
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();
    const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
    const stk::mesh::Part* part = recvMeta.get_part(m_recvPartName);
    ASSERT_TRUE(nullptr != part);
    stk::mesh::Selector selector(*part);
    stk::mesh::EntityVector entities;

    stk::mesh::get_selected_entities(selector, m_recvBulk->buckets(recvRank), entities);

    std::vector<double> location;
    for(stk::mesh::Entity entity : entities) {
      stk::search::determine_centroid(m_sendSpatialDim, entity, *recvCoordinateField, location);
      test_expected_indexed_value_at_location(*m_recvBulk, m_recvField, recvIndex, entity, eval, location.data());
    }
  }

  void test_centroid_transfer(stk::unit_test_util::FieldEvaluator& eval, stk::mesh::EntityRank recvRank, [[maybe_unused]] unsigned numCopies = 1)
  {
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();
    const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
    const stk::mesh::Part* part = recvMeta.get_part(m_recvPartName);
    ASSERT_TRUE(nullptr != part);
    stk::mesh::Selector selector(*part);
    stk::mesh::EntityVector entities;

    stk::mesh::get_selected_entities(selector, m_recvBulk->buckets(recvRank), entities);

    std::vector<double> location;
    for(stk::mesh::Entity entity : entities) {
      stk::search::determine_centroid(m_recvSpatialDim, entity, *recvCoordinateField, location);
      test_expected_value_at_location(*m_recvBulk, m_recvField, entity, eval, std::min(m_sendSpatialDim, m_recvSpatialDim), location.data());
    }
  }

  void test_indexed_gauss_point_transfer(stk::unit_test_util::FieldEvaluator& eval,
                                         stk::mesh::EntityRank recvRank, unsigned recvIndex)
  {
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();
    const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
    const stk::mesh::Part* part = recvMeta.get_part(m_recvPartName);
    ASSERT_TRUE(nullptr != part);
    stk::mesh::Selector selector(*part);
    stk::mesh::EntityVector entities;

    stk::mesh::get_selected_entities(selector, m_recvBulk->buckets(recvRank), entities);

    unsigned coordCompStride = 1;
    unsigned coordCopyStride = 3;
    std::vector<double> gpCoordinates;
    for(stk::mesh::Entity entity : entities) {
      stk::search::determine_gauss_points(*m_recvBulk, entity, *m_masterElemProvider, *recvCoordinateField, gpCoordinates);
      test_expected_indexed_value_at_locations(*m_recvBulk, m_recvField, recvIndex, entity, eval,
                                               gpCoordinates.data(), coordCompStride, coordCopyStride);
    }
  }

  void test_gauss_point_transfer(stk::unit_test_util::FieldEvaluator& eval, stk::mesh::EntityRank recvRank)
  {
    stk::mesh::MetaData& recvMeta = m_recvBulk->mesh_meta_data();
    const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
    const stk::mesh::Part* part = recvMeta.get_part(m_recvPartName);
    ASSERT_TRUE(nullptr != part);
    stk::mesh::Selector selector(*part);
    stk::mesh::EntityVector entities;

    stk::mesh::get_selected_entities(selector, m_recvBulk->buckets(recvRank), entities);

    unsigned coordCompStride = 1;
    unsigned coordCopyStride = 3;
    std::vector<double> gpCoordinates;
    for(stk::mesh::Entity entity : entities) {
      stk::search::determine_gauss_points(*m_recvBulk, entity, *m_masterElemProvider, *recvCoordinateField, gpCoordinates);
      test_expected_value_at_locations(*m_recvBulk, m_recvField, entity, eval, std::min(m_sendSpatialDim, m_recvSpatialDim),
                                       gpCoordinates.data(), coordCompStride, coordCopyStride);
    }
  }

  void initialize_volume_transfer_variables(unsigned nx = 1, unsigned ny = 1, unsigned nz = 1)
  {
    // Create surfaces and faces
    std::ostringstream oss;

    oss << "generated:";
    oss << nx << "x";
    oss << ny << "x";
    oss << nz;

    initialize_transfer_variables(oss.str(), oss.str(), "block_1", "block_1");
  }

  void initialize_surface_transfer_variables(unsigned nx = 1, unsigned ny = 1, unsigned nz = 1)
  {
    // Create surfaces and faces
    std::ostringstream oss;

    oss << "generated:";
    oss << nx << "x";
    oss << ny << "x";
    oss << nz << "|sideset:xXyYzZ";

    initialize_transfer_variables(oss.str(), oss.str(), "surface_1", "surface_1");
  }

  void initialize_transfer_variables(const std::string& sendMeshSpec, const std::string& recvMeshSpec,
                                     const std::string& sendPartName, const std::string& recvPartName,
                                     int sendSpatialDimension=3,    int recvSpatialDimension=3)
  {
    m_sendMeshSpec = sendMeshSpec;
    m_recvMeshSpec = recvMeshSpec;
    m_sendPartName = sendPartName;
    m_recvPartName = recvPartName;

    m_sendSpatialDim = sendSpatialDimension;
    m_recvSpatialDim = recvSpatialDimension;
  }

 protected:
  stk::ParallelMachine m_comm{MPI_COMM_WORLD};

  unsigned m_sendSpatialDim{3u};
  unsigned m_recvSpatialDim{3u};

  std::string m_sendFieldName{"field1"};
  std::string m_recvFieldName{"field1"};

  // Initialization for volume transfer
  std::string m_sendPartName{"block_1"};
  std::string m_recvPartName{"block_1"};

  // Initialization for volume transfer
  std::string m_sendMeshSpec{"generated:1x1x1"};
  std::string m_recvMeshSpec{"generated:1x1x1"};

  std::shared_ptr<stk::mesh::BulkData> m_sendBulk;
  std::shared_ptr<stk::mesh::BulkData> m_recvBulk;

  stk::mesh::FieldBase* m_sendField{nullptr};
  stk::mesh::FieldBase* m_recvField{nullptr};

  unsigned m_numSendCopies{1};
  unsigned m_numRecvCopies{1};

  stk::search::ObjectOutsideDomainPolicy m_extrapolateOption{stk::search::ObjectOutsideDomainPolicy::ABORT};
  std::shared_ptr<stk::search::MasterElementProviderInterface> m_masterElemProvider;
};

TEST_F(MockSimpleTransferTest, faceSendMesh_masterElementInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  initialize_surface_transfer_variables();

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_masterElementInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_masterElementInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  initialize_surface_transfer_variables();

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::FACE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::FACE_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_masterElementInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_masterElementInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  initialize_surface_transfer_variables();

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::FACE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  stk::search::SearchTopology quad4Topo(stk::topology::QUAD_4);
  setup_transfer(transfer, sendRank, recvRank, 1u, quad4Topo, 4u);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_masterElementInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_patchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_surface_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::FACE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_patchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_patchRecoveryInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_surface_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::FACE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::FACE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::FACE_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_patchRecoveryInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_patchRecoveryInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_surface_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::FACE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::FACE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  stk::search::SearchTopology quad4Topo(stk::topology::QUAD_4);
  setup_transfer(transfer, sendRank, recvRank, 1u, quad4Topo, 4u);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_patchRecoveryInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, nodeSendMesh_copyToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::NODE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_copyToFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  initialize_surface_transfer_variables();

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::FACE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::FACE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::FACE_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_copyToElement)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}


TEST_F(MockSimpleTransferTest, elementSendMesh_copyToElement_multiFieldCopy)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  unsigned numCopies = 3;
  setup_transfer(transfer, sendRank, recvRank, numCopies);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank, numCopies);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_copyToElement_GaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  unsigned numCopies = 8;
  setup_transfer(transfer, sendRank, recvRank, numCopies);
  stk::unit_test_util::set_entity_field_gauss_point(*m_sendBulk, *m_sendField, eval, m_masterElemProvider);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(eval, recvRank);
}

using MockSimpleIndexedTransferTest = MockSimpleTransferTest;

class SimpleIndexedTransfer : public stk::transfer::spmd::SimpleTransfer {
 public:
  SimpleIndexedTransfer(const std::string& transferName, unsigned sendIndex, unsigned recvIndex)
  : SimpleTransfer(transferName)
  , m_sendIndex(sendIndex)
  , m_recvIndex(recvIndex)
  { }

  ~SimpleIndexedTransfer() = default;

  using stk::transfer::spmd::SimpleTransfer::add_send_field;
  using stk::transfer::spmd::SimpleTransfer::add_recv_field;

  void add_send_field(stk::mesh::FieldBase & field) override
  {
    STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
    add_field(m_sendFieldSpecs, field, m_sendIndex);
  }
  void add_recv_field(stk::mesh::FieldBase & field) override
  {
    STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
    add_field(m_recvFieldSpecs, field, m_recvIndex);
  }

 protected:
  unsigned m_sendIndex{0};
  unsigned m_recvIndex{0};

  void add_field(std::vector<stk::transfer::FieldSpec> & fieldSpecs, stk::mesh::FieldBase & field, unsigned fieldIndex)
  {
    stk::transfer::FieldSpec fieldSpec(field.name(), field.state(), fieldIndex);
    fieldSpecs.push_back(fieldSpec);

    // all fields need same entity rank
    const stk::mesh::MetaData& meta = field.mesh_meta_data();

    stk::mesh::FieldBase* firstField = meta.get_field(field.entity_rank(), fieldSpecs.front().name);
    STK_ThrowRequire(nullptr != firstField);
  }
};

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_masterElementInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_node_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_masterElementInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_centroid_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_masterElementInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_gauss_point_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_patchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK,stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_node_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_patchRecoveryInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_centroid_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_patchRecoveryInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_gauss_point_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, nodeSendMesh_copyToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::NODE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_node_transfer(eval, recvRank, recvIndex);
}

TEST_F(MockSimpleIndexedTransferTest, elementSendMesh_copyToElement)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // One based indexing
  unsigned sendIndex = 1;
  unsigned recvIndex = 3;

  SimpleIndexedTransfer transfer("TransferTest", sendIndex, recvIndex);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_indexed_centroid_transfer(eval, recvRank, recvIndex);
}


using MockSimpleBoundedTransferTest = MockSimpleTransferTest;

class SimpleBoundedTransfer : public stk::transfer::spmd::SimpleTransfer {
 public:
  SimpleBoundedTransfer(const std::string& transferName, double lowerBound, double upperBound)
  : SimpleTransfer(transferName)
  , m_lowerBound(lowerBound)
  , m_upperBound(upperBound)
  { }

  ~SimpleBoundedTransfer() = default;

  using stk::transfer::spmd::SimpleTransfer::add_send_field;
  using stk::transfer::spmd::SimpleTransfer::add_recv_field;

  void add_send_field(stk::mesh::FieldBase & field) override
  {
    STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
    add_field(m_sendFieldSpecs, field);
  }
  void add_recv_field(stk::mesh::FieldBase & field) override
  {
    STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
    add_field(m_recvFieldSpecs, field);
  }

 protected:
  double m_lowerBound{0};
  double m_upperBound{0};

  void add_field(std::vector<stk::transfer::FieldSpec> & fieldSpecs, stk::mesh::FieldBase & field)
  {
    stk::transfer::FieldSpec fieldSpec(field.name(), field.state());
    fieldSpec.lowerBound = m_lowerBound;
    fieldSpec.upperBound = m_upperBound;
    fieldSpecs.push_back(fieldSpec);

    // all fields need same entity rank
    const stk::mesh::MetaData& meta = field.mesh_meta_data();

    stk::mesh::FieldBase* firstField = meta.get_field(field.entity_rank(), fieldSpecs.front().name);
    STK_ThrowRequire(nullptr != firstField);
  }
};

TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_masterElementInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // Eval has a range of [1,3] on unit cube
  double lowerBound = 1.5;
  double upperBound = 2.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_node_transfer(boundedEval, recvRank);
}


TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_masterElementInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Eval has a range of [1,3] on unit cube
  double lowerBound = 1.5;
  double upperBound = 2.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_centroid_transfer(boundedEval, recvRank);
}

TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_masterElementInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Eval has a range of [1,3] on unit cube
  double lowerBound = 1.5;
  double upperBound = 2.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_gauss_point_transfer(boundedEval, recvRank);
}

TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_patchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // Eval has a range of [1,5] on 2x2x2 cube
  double lowerBound = 2.5;
  double upperBound = 3.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_node_transfer(boundedEval, recvRank);
}

TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_patchRecoveryInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Eval has a range of [1,5] on 2x2x2 cube
  double lowerBound = 2.5;
  double upperBound = 3.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_centroid_transfer(boundedEval, recvRank);
}

TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_patchRecoveryInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Eval has a range of [1,5] on 2x2x2 cube
  double lowerBound = 2.5;
  double upperBound = 3.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_gauss_point_transfer(boundedEval, recvRank);
}

TEST_F(MockSimpleBoundedTransferTest, nodeSendMesh_copyToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // Eval has a range of [1,3] on unit cube
  double lowerBound = 1.5;
  double upperBound = 2.5;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::NODE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_node_transfer(boundedEval, recvRank);
}

TEST_F(MockSimpleBoundedTransferTest, elementSendMesh_copyToElementCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Eval has a range of [1,3] on unit cube
  // Centroid value is 2.0 so we set a lower bound
  double lowerBound = 2.5;
  double upperBound = 3.0;

  SimpleBoundedTransfer transfer("TransferTest", lowerBound, upperBound);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  stk::unit_test_util::BoundedLinearFieldEvaluator boundedEval(m_sendSpatialDim, lowerBound, upperBound);
  test_centroid_transfer(boundedEval, recvRank);
}


using MockSimpleTransformTransferTest = MockSimpleTransferTest;

class SimpleTransformTransfer : public stk::transfer::spmd::SimpleTransfer {
 public:
  SimpleTransformTransfer(const std::string& transferName,
                              stk::transfer::FieldTransform preTransform,
                              stk::transfer::FieldTransform postTransform)
  : SimpleTransfer(transferName)
  , m_preTransform(preTransform)
  , m_postTransform(postTransform)
  { }

  ~SimpleTransformTransfer() = default;

  using stk::transfer::spmd::SimpleTransfer::add_send_field;
  using stk::transfer::spmd::SimpleTransfer::add_recv_field;

  void add_send_field(stk::mesh::FieldBase & field) override
  {
    STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
    add_field(m_sendFieldSpecs, field);
  }
  void add_recv_field(stk::mesh::FieldBase & field) override
  {
    STK_ThrowRequireMsg(!m_committed, "Transfer named: " << m_transferName << " has already been committed");
    add_field(m_recvFieldSpecs, field);
  }

 protected:
  stk::transfer::FieldTransform m_preTransform;
  stk::transfer::FieldTransform m_postTransform;

  void add_field(std::vector<stk::transfer::FieldSpec> & fieldSpecs, stk::mesh::FieldBase & field)
  {
    stk::transfer::FieldSpec fieldSpec(field.name(), field.state());
    fieldSpec.preTransform = m_preTransform;
    fieldSpec.postTransform = m_postTransform;
    fieldSpecs.push_back(fieldSpec);

    // all fields need same entity rank
    const stk::mesh::MetaData& meta = field.mesh_meta_data();

    stk::mesh::FieldBase* firstField = meta.get_field(field.entity_rank(), fieldSpecs.front().name);
    STK_ThrowRequire(nullptr != firstField);
  }
};

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_masterElementInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 4.0, 4.0, 4.0, 4.0); // x4

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // Pre-transform by 0.5, then post-transform by 8 => final scale by 4
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_eight);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(goldEval, recvRank);
}


TEST_F(MockSimpleTransformTransferTest, elementSendMesh_masterElementInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 8.0, 8.0, 8.0, 8.0); // x8

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Pre-transform by 2, then post-transform by 4 => final scale by 8
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_two, impl::transform_scale_by_four);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_masterElementInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 4.0, 4.0, 4.0, 4.0); // x4

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Pre-transform by 0.5, then post-transform by 8 => final scale by 4
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_eight);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_patchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 4.0, 4.0, 4.0, 4.0); // x4

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // Pre-transform by 0.5, then post-transform by 8 => final scale by 4
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_eight);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_patchRecoveryInterpolationToCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 4.0, 4.0, 4.0, 4.0); // x4

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Pre-transform by 0.5, then post-transform by 8 => final scale by 4
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_eight);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_patchRecoveryInterpolationToGaussPoint)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 2.0, 2.0, 2.0, 2.0); // x2

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Pre-transform by 4, then post-transform by 0.5 => final scale by 2
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_four, impl::transform_scale_by_half);

  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
  setup_transfer(transfer, sendRank, recvRank, 1u, hex8Topo, 8u);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_gauss_point_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, nodeSendMesh_copyToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 4.0, 4.0, 4.0, 4.0); // x4

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  // Pre-transform by 0.5, then post-transform by 8 => final scale by 4
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_eight);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::NODE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_copyToElement)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 2.0, 2.0, 2.0, 2.0); // x2

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Pre-transform by 0.5, then post-transform by 4 => final scale by 2
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_four);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransformTransferTest, elementSendMesh_copyToElement_multiFieldCopy)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 2.0, 2.0, 2.0, 2.0); // x2

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  // Pre-transform by 0.5, then post-transform by 4 => final scale by 2
  SimpleTransformTransfer transfer("TransferTest", impl::transform_scale_by_half, impl::transform_scale_by_four);

  unsigned numCopies = 3;
  setup_transfer(transfer, sendRank, recvRank, numCopies);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.setup_copy_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                       stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                       m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(goldEval, recvRank, numCopies);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_linearPatchRecoveryInterpolationToNode_mls)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Linear patch recovery in 3D requires a minimum mesh of 2x2x2
  initialize_volume_transfer_variables(2,2,2);

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.set_patch_recovery_evaluation_type(stk::transfer::PatchRecoveryEvaluationType::LINEAR_MOVING_LEAST_SQUARES);
  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_quadraticPatchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Quadratic patch recovery in 3D requires a minimum send mesh of 3x3x3
  // Create a 1x1x1 recv mesh centered in this send mesh so that interpolation
  // is based on the central hex as the seed element
  m_sendMeshSpec = "generated:3x3x3";
  m_recvMeshSpec = "generated:1x1x1|bbox:1.05,1.05,1.05,1.95,1.95,1.95";

  std::vector<double> xCoeffs{ 1.0, 1.0, 1.0 };
  std::vector<double> yCoeffs{ 1.0, 1.0, 1.0 };
  std::vector<double> zCoeffs{ 1.0, 1.0, 1.0 };
  stk::unit_test_util::QuadraticFieldEvaluatorWithCoefficients eval(xCoeffs, yCoeffs, zCoeffs);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.set_patch_recovery_evaluation_type(stk::transfer::PatchRecoveryEvaluationType::QUADRATIC_LEAST_SQUARES);
  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_quadraticPatchRecoveryInterpolationToNode_mls)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Quadratic patch recovery in 3D requires a minimum send mesh of 3x3x3
  // Create a 1x1x1 recv mesh centered in this send mesh so that interpolation
  // is based on the central hex as the seed element
  m_sendMeshSpec = "generated:3x3x3";
  m_recvMeshSpec = "generated:1x1x1|bbox:1.05,1.05,1.05,1.95,1.95,1.95";

  std::vector<double> xCoeffs{ 1.0, 1.0, 1.0 };
  std::vector<double> yCoeffs{ 1.0, 1.0, 1.0 };
  std::vector<double> zCoeffs{ 1.0, 1.0, 1.0 };
  stk::unit_test_util::QuadraticFieldEvaluatorWithCoefficients eval(xCoeffs, yCoeffs, zCoeffs);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.set_patch_recovery_evaluation_type(stk::transfer::PatchRecoveryEvaluationType::QUADRATIC_MOVING_LEAST_SQUARES);
  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_cubicPatchRecoveryInterpolationToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Cubic patch recovery in 3D requires a minimum send mesh of 4x4x4
  // Create a 1x1x1 recv mesh centered in this send mesh so that interpolation
  // is based on the central hex as the seed element
  m_sendMeshSpec = "generated:5x5x5";
  m_recvMeshSpec = "generated:1x1x1|bbox:2.05,2.05,2.05,2.95,2.95,2.95";

  std::vector<double> xCoeffs{ 1.0, 1.0, 1.0, 1.0 };
  std::vector<double> yCoeffs{ 1.0, 1.0, 1.0, 1.0 };
  std::vector<double> zCoeffs{ 1.0, 1.0, 1.0, 1.0 };
  stk::unit_test_util::CubicFieldEvaluatorWithCoefficients eval(xCoeffs, yCoeffs, zCoeffs);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.set_patch_recovery_evaluation_type(stk::transfer::PatchRecoveryEvaluationType::CUBIC_LEAST_SQUARES);
  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  double tolerance = 1.0e-4;
  test_node_transfer(eval, recvRank, tolerance);
}

TEST_F(MockSimpleTransferTest, elementSendMesh_cubicPatchRecoveryInterpolationToNode_mls)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  // Cubic patch recovery in 3D requires a minimum send mesh of 4x4x4
  // Create a 1x1x1 recv mesh centered in this send mesh so that interpolation
  // is based on the central hex as the seed element
  m_sendMeshSpec = "generated:5x5x5";
  m_recvMeshSpec = "generated:1x1x1|bbox:2.05,2.05,2.05,2.95,2.95,2.95";

  std::vector<double> xCoeffs{ 1.0, 1.0, 1.0, 1.0 };
  std::vector<double> yCoeffs{ 1.0, 1.0, 1.0, 1.0 };
  std::vector<double> zCoeffs{ 1.0, 1.0, 1.0, 1.0 };
  stk::unit_test_util::CubicFieldEvaluatorWithCoefficients eval(xCoeffs, yCoeffs, zCoeffs);

  stk::mesh::EntityRank sendRank = stk::topology::ELEM_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_entity_field(*m_sendBulk, *m_sendField, eval);

  transfer.set_patch_recovery_evaluation_type(stk::transfer::PatchRecoveryEvaluationType::CUBIC_MOVING_LEAST_SQUARES);
  transfer.setup_patch_recovery_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEMENT_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  double tolerance = 1.0e-4;
  test_node_transfer(eval, recvRank, tolerance);
}

TEST_F(MockSimpleTransferTest, nodeSendMesh_sumToNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator     eval(m_sendSpatialDim, 1.0, 1.0, 1.0, 1.0); // x1
  stk::unit_test_util::LinearFieldEvaluator goldEval(m_sendSpatialDim, 2.0, 2.0, 2.0, 2.0); // x2

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  stk::transfer::spmd::SimpleTransfer transfer("TransferTest");

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  // Initialize recv mesh values to send mesh values so that sum results in a 2x
  stk::unit_test_util::set_node_field(*m_recvBulk, *m_recvField, eval);

  transfer.setup_sum_nearest_transfer(*m_sendBulk, *m_recvBulk,
                                      stk::topology::NODE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                      m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(goldEval, recvRank);
}

TEST_F(MockSimpleTransferTest, ThreeDElementToTwoDNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(2);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  std::string sendMesh = "generated:2x2x2|bbox:0,0,0,2,2,2";
  std::string recvMesh = "textmesh:"
                         "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                         "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                         "0,4,QUAD_4_2D,5,6,9,8,block_1\n"
                         "|coordinates:"
                         "0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2|dimension:2";

  initialize_transfer_variables(sendMesh, recvMesh, "block_1", "block_1", 3, 2);
  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  auto coordTransform = std::make_shared<stk::search::CoordTransformAddZ>("z=0");
  transfer.add_recv_mesh_coord_transform(coordTransform);
  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, faceSendMesh_ThreeDElementToTwoDNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(2);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  std::string sendMesh = "generated:2x2x2|sideset:Z|bbox:0,0,0,2,2,2";
  std::string recvMesh = "textmesh:"
                         "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                         "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                         "0,4,QUAD_4_2D,5,6,9,8,block_1\n"
                         "|coordinates:"
                         "0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2|dimension:2";

  initialize_transfer_variables(sendMesh, recvMesh, "surface_1", "block_1", 3, 2);
  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  auto coordTransform = std::make_shared<stk::search::CoordTransformAddZ>("z=2");
  transfer.add_recv_mesh_coord_transform(coordTransform);
  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, ThreeDElementToTwoDCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(2);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::ELEM_RANK;

  std::string sendMesh = "generated:2x2x2|sideset:z|bbox:0,0,0,2,2,2";
  std::string recvMesh = "textmesh:"
                         "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                         "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                         "0,4,QUAD_4_2D,5,6,9,8,block_1\n"
                         "|coordinates:"
                         "0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2|dimension:2";

  initialize_transfer_variables(sendMesh, recvMesh, "surface_1", "block_1", 3, 2);
  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  auto coordTransform = std::make_shared<stk::search::CoordTransformAddZ>("z=0");
  transfer.add_recv_mesh_coord_transform(coordTransform);
  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::FACE_RANK, stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, TwoDElementTo3DNode)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvRank = stk::topology::NODE_RANK;

  std::string sendMesh = "textmesh:"
                         "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                         "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                         "0,4,QUAD_4_2D,5,6,9,8,block_1\n"
                         "|coordinates:"
                         "0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2|dimension:2";
  std::string recvMesh = "generated:2x2x2|sideset:z";

  initialize_transfer_variables(sendMesh, recvMesh, "block_1", "surface_1", 2, 3);
  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendRank, recvRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  auto coordTransform = std::make_shared<stk::search::CoordTransformRemoveZ>();
  transfer.add_recv_mesh_coord_transform(coordTransform);
  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         stk::topology::ELEM_RANK, stk::transfer::spmd::RecvMeshType::NODE,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_node_transfer(eval, recvRank);
}

TEST_F(MockSimpleTransferTest, TwoDElementToThreeDCentroid)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  stk::unit_test_util::LinearFieldEvaluator eval(m_sendSpatialDim);

  stk::mesh::EntityRank sendSearchEntityRank = stk::topology::ELEM_RANK;
  stk::transfer::spmd::RecvMeshType recvSearchEntityType = stk::transfer::spmd::RecvMeshType::FACE_CENTROID;
  stk::mesh::EntityRank sendFieldRank = stk::topology::NODE_RANK;
  stk::mesh::EntityRank recvFieldRank = stk::topology::FACE_RANK;

  std::string sendMesh = "textmesh:"
                         "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                         "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                         "0,4,QUAD_4_2D,5,6,9,8,block_1\n"
                         "|coordinates:"
                         "0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2|dimension:2";
  std::string recvMesh = "generated:2x2x2|sideset:z";


  initialize_transfer_variables(sendMesh, recvMesh, "block_1", "surface_1", 2, 3);
  stk::transfer::spmd::SimpleTransfer transfer("TransferTest", m_comm);

  setup_transfer(transfer, sendFieldRank, recvFieldRank);
  stk::unit_test_util::set_node_field(*m_sendBulk, *m_sendField, eval);

  auto coordTransform = std::make_shared<stk::search::CoordTransformRemoveZ>();
  transfer.add_recv_mesh_coord_transform(coordTransform);
  transfer.setup_master_element_transfer(*m_sendBulk, *m_recvBulk,
                                         sendSearchEntityRank, recvSearchEntityType,
                                         m_masterElemProvider, m_extrapolateOption);
  transfer.initialize();
  transfer.apply();

  test_centroid_transfer(eval, recvFieldRank);
}

}

