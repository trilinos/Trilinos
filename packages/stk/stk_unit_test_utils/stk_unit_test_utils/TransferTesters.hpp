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

#ifndef stk_unit_test_utils_transfertesters_hpp
#define stk_unit_test_utils_transfertesters_hpp

#include <stk_util/stk_config.h>
#include "stk_util/environment/Env.hpp"
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_search_util/MasterElementProvider.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/FieldEvaluator.hpp>
#ifdef STK_HAVE_INTREPID2
#include <Intrepid2_CellTools.hpp>
#endif
#include <Ioss_Field.h>
#include <string>
#include <memory>
#include <gtest/gtest.h>

namespace stk
{
namespace unit_test_util
{

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

inline void test_expected_value_at_location(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field,
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
               double expectedValue = eval(entityArg, x, y, z);

               for (stk::mesh::CopyIdx copy : data.copies()) {
                 for(stk::mesh::ComponentIdx comp : data.components()) {
                   double value = data(copy, comp);

                   EXPECT_NEAR(expectedValue, value, tolerance)
                       << " for processor " << bulk.parallel_rank() << " with entity " << bulk.entity_key(entity)
                       << " at loc (" << x << "," << y << "," << z << ")";
                 }
               }
  };

  field_data_test<double, stk::mesh::ReadOnly>(*field, entity, alg);
}

 inline void test_expected_value_at_locations(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* field,
                                        stk::mesh::Entity& entity, const stk::unit_test_util::FieldEvaluator& eval,
                                        const double* loc, const unsigned locCompStride, const unsigned locCopyStride,
                                        double tolerance = 1.0e-6)
{
  // Coordinate data per copy e.g gauss points
  const unsigned nDim = bulk.mesh_meta_data().spatial_dimension();

  auto alg = [&](auto& fieldData, stk::mesh::Entity entityArg) {
               auto data = fieldData.entity_values(entity);

               for (stk::mesh::CopyIdx copy : data.copies()) {
                 unsigned copyIndex = copy;
                 double x = loc[copyIndex*locCopyStride + 0*locCompStride];
                 double y = loc[copyIndex*locCopyStride + 1*locCompStride];
                 double z = (nDim != 3 ? 0.0 : loc[copyIndex*locCopyStride + 2*locCompStride]);

                 double expectedValue = eval(entityArg, x, y, z);

                 for(stk::mesh::ComponentIdx comp : data.components()) {
                   double value = data(copy, comp);

                   EXPECT_NEAR(expectedValue, value, tolerance)
                       << " for processor " << sierra::Env::parallel_rank() << " with entity " << bulk.entity_key(entityArg)
                       << " at loc (" << x << "," << y << "," << z << ")";
                 }
               }
  };

  field_data_test<double, stk::mesh::ReadOnly>(*field, entity, alg);
}

inline void test_node_transfer(stk::unit_test_util::FieldEvaluator& eval,
                        std::shared_ptr<stk::mesh::BulkData> recvBulk,
                        stk::mesh::EntityRank recvRank,
                        double tolerance = 1.0e-6)
{
  stk::mesh::MetaData& recvMeta = recvBulk->mesh_meta_data();
  const unsigned recvDim = recvMeta.spatial_dimension();

  auto recvField = stk::mesh::get_field_by_name("field_1", recvMeta);
  const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
  const stk::mesh::Part* part = recvMeta.get_part("block_1");
  ASSERT_TRUE(nullptr != part);
  stk::mesh::Selector selector(*part);
  stk::mesh::EntityVector entities;

  stk::mesh::get_selected_entities(selector, recvBulk->buckets(recvRank), entities);

  std::vector<double> location;
  for(stk::mesh::Entity entity : entities) {
    stk::search::determine_centroid(recvDim, entity, *recvCoordinateField, location);
    test_expected_value_at_location(*recvBulk, recvField, entity, eval, location.data(), tolerance);
  }
}

inline void test_centroid_transfer(stk::unit_test_util::FieldEvaluator& eval,
                            std::shared_ptr<stk::mesh::BulkData> recvBulk,
                            stk::mesh::EntityRank recvRank)
{
  stk::mesh::MetaData& recvMeta = recvBulk->mesh_meta_data();
  const unsigned recvDim = recvMeta.spatial_dimension();

  auto recvField = stk::mesh::get_field_by_name("field_1", recvMeta);
  const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
  const stk::mesh::Part* part = recvMeta.get_part("block_1");
  ASSERT_TRUE(nullptr != part);
  stk::mesh::Selector selector(*part);
  stk::mesh::EntityVector entities;

  stk::mesh::get_selected_entities(selector, recvBulk->buckets(recvRank), entities);

  std::vector<double> location;
  for(stk::mesh::Entity entity : entities) {
    stk::search::determine_centroid(recvDim, entity, *recvCoordinateField, location);
    test_expected_value_at_location(*recvBulk, recvField, entity, eval, location.data());
  }
}

inline void test_gauss_point_transfer(stk::unit_test_util::FieldEvaluator& eval,
                               std::shared_ptr<stk::mesh::BulkData> recvBulk,
                               stk::mesh::EntityRank recvRank,
                               std::string recvFieldName,
                               std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
{
  stk::mesh::MetaData& recvMeta = recvBulk->mesh_meta_data();
  auto recvField = stk::mesh::get_field_by_name(recvFieldName, recvMeta);
  const stk::mesh::FieldBase* recvCoordinateField = recvMeta.coordinate_field();
  const stk::mesh::Part* part = recvMeta.get_part("block_1");
  ASSERT_TRUE(nullptr != part);
  stk::mesh::Selector selector(*part);
  stk::mesh::EntityVector entities;

  stk::mesh::get_selected_entities(selector, recvBulk->buckets(recvRank), entities);

  unsigned coordCompStride = 1;
  unsigned coordCopyStride = 3;
  std::vector<double> gpCoordinates;
  for(stk::mesh::Entity entity : entities) {
    stk::search::determine_gauss_points(*recvBulk, entity, *masterElemProvider, *recvCoordinateField, gpCoordinates);
    test_expected_value_at_locations(*recvBulk, recvField, entity, eval,
                                     gpCoordinates.data(), coordCompStride, coordCopyStride);
  }
}

inline
void test_find_parametric_coordinates(const std::string& meshDesc,
                                      const stk::mesh::EntityRank rank,
                                      const stk::mesh::EntityId id,
                                      std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
                                      const std::vector<double>& point,
                                      const double expectedParametricDistance,
                                      const double parametricTolerance = 1.0e-6)
{
  MPI_Comm comm = stk::parallel_machine_world();

  unsigned spatialDim = point.size();
  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(spatialDim).create();
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  stk::unit_test_util::setup_text_mesh(mesh, meshDesc);

  stk::mesh::Entity entity = mesh.get_entity(rank, id);
  ASSERT_TRUE(mesh.is_valid(entity));

  std::vector<double> elemNodeCoords;
  stk::search::SearchField coordField(meta.coordinate_field());
  stk::search::spmd::EntityKeyPair entityKeyPair(entity, mesh.entity_key(entity));
  unsigned numNodes = 0;

  masterElemProvider->nodal_field_data(entityKeyPair, coordField,
                                       spatialDim, numNodes, elemNodeCoords);

  stk::topology stkTopo = mesh.bucket(entity).topology();

  std::vector<double> parCoords = {99.9, 99.9, 99.9};
  stk::search::SearchTopology srchTopo(stkTopo);
  double parametricDistance = 999.9;
  masterElemProvider->find_parametric_coordinates(srchTopo, spatialDim,
                                             elemNodeCoords, point, parCoords,
                                             parametricDistance);

  EXPECT_NEAR(expectedParametricDistance, parametricDistance, parametricTolerance);
}

}
}

#endif

