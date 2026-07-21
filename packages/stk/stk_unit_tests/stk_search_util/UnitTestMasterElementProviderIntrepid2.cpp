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
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include "stk_mesh/base/Entity.hpp"
#include <stk_mesh/base/BulkData.hpp>                      // for BulkData
#include "stk_mesh/base/EntityKey.hpp"                     // for EntityKey
#include "stk_mesh/base/Field.hpp"                         // for Field
#include "stk_mesh/base/FieldBLAS.hpp"                     // for field_copy
#include "stk_search_util/MasterElementProviderIntrepid2.hpp"
#include "stk_search_util/Intrepid2_HasParametricDistance.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"               // for build_mesh
#include "stk_unit_test_utils/TextMesh.hpp"                // for get_full_t...
#include "stk_unit_test_utils/UnitTestSearchUtils.hpp"

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

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace {

TEST(MasterElementProviderIntrepid2Test, Hex_8)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::HEX_8, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));
  unsigned numIntegrationPts = masterElemProvider->num_integration_points(meTopo);
  EXPECT_EQ(8u, numIntegrationPts);

  std::vector<double> integrationPoints;
  masterElemProvider->integration_points(meTopo, integrationPoints);

  for(unsigned i=0; i<numIntegrationPts; i++) {
    unsigned offset = 3*i;
    EXPECT_TRUE(integrationPoints[offset+0] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+0] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+1] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+1] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+2] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+2] <  1.0);
  }

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  EXPECT_EQ(numFieldComponents, meta.spatial_dimension());
  EXPECT_EQ(numNodes, topo.num_nodes());

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, elementCentroid,
                                                  paramCoords, paramDistance);

  EXPECT_NEAR(0.0, paramCoords[0], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[1], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[2], 1.0e-6);
  EXPECT_NEAR(0.0, paramDistance, 1.0e-6);

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> result(numFieldComponents, 0.);
  masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result);
  EXPECT_NEAR(5.0, result[0], 1.0e-6);
  EXPECT_NEAR(5.0, result[1], 1.0e-6);
  EXPECT_NEAR(5.0, result[2], 1.0e-6);

  std::vector<double> coordCenter;
  masterElemProvider->coordinate_center(meTopo, coordCenter);
  EXPECT_NEAR(0.0, coordCenter[0], 1.0e-6);
  EXPECT_NEAR(0.0, coordCenter[1], 1.0e-6);
  EXPECT_NEAR(0.0, coordCenter[2], 1.0e-6);

}

TEST(MasterElementProviderIntrepid2Test, Tet_4)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::string meshDesc = "0,1,TET_4,1,2,3,4"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1";
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::TET_4, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));
  unsigned numIntegrationPts = masterElemProvider->num_integration_points(meTopo);
  EXPECT_EQ(4u, numIntegrationPts);

  std::vector<double> integrationPoints;
  masterElemProvider->integration_points(meTopo, integrationPoints);

  for(unsigned i=0; i<numIntegrationPts; i++) {
    unsigned offset = 3*i;

    EXPECT_TRUE(integrationPoints[offset+0] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+0] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+1] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+1] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+2] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+2] <  1.0);
  }

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  EXPECT_EQ(numFieldComponents, meta.spatial_dimension());
  EXPECT_EQ(numNodes, topo.num_nodes());

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, elementCentroid,
                                                  paramCoords, paramDistance);

  EXPECT_NEAR(0.25, paramCoords[0], 1.0e-6);
  EXPECT_NEAR(0.25, paramCoords[1], 1.0e-6);
  EXPECT_NEAR(0.25, paramCoords[2], 1.0e-6);
  EXPECT_NEAR(0.0, paramDistance, 1.0e-6);

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> result(numFieldComponents, 0.);
  masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result);
  EXPECT_NEAR(5.0, result[0], 1.0e-6);
  EXPECT_NEAR(5.0, result[1], 1.0e-6);
  EXPECT_NEAR(5.0, result[2], 1.0e-6);

  std::vector<double> coordCenter;
  masterElemProvider->coordinate_center(meTopo, coordCenter);
  EXPECT_NEAR(0.25, coordCenter[0], 1.0e-6);
  EXPECT_NEAR(0.25, coordCenter[1], 1.0e-6);
  EXPECT_NEAR(0.25, coordCenter[2], 1.0e-6);

}

TEST(MasterElementProviderIntrepid2Test, Pyramid_5)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0.5,0.5,1";
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::PYRAMID_5, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));
  unsigned numIntegrationPts = masterElemProvider->num_integration_points(meTopo);
  EXPECT_EQ(8u, numIntegrationPts);

  std::vector<double> integrationPoints;
  masterElemProvider->integration_points(meTopo, integrationPoints);

  for(unsigned i=0; i<numIntegrationPts; i++) {
    unsigned offset = 3*i;

    EXPECT_TRUE(integrationPoints[offset+0] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+0] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+1] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+1] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+2] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+2] <  1.0);
  }

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  EXPECT_EQ(numFieldComponents, meta.spatial_dimension());
  EXPECT_EQ(numNodes, topo.num_nodes());

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, elementCentroid,
                                                  paramCoords, paramDistance);

  EXPECT_NEAR(0.0, paramCoords[0], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[1], 1.0e-6);
  EXPECT_NEAR(0.2, paramCoords[2], 1.0e-6);
  EXPECT_NEAR(0.4, paramDistance, 1.0e-6);

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> result(numFieldComponents, 0.);
  masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result);
  EXPECT_NEAR(5.0, result[0], 1.0e-6);
  EXPECT_NEAR(5.0, result[1], 1.0e-6);
  EXPECT_NEAR(5.0, result[2], 1.0e-6);

  std::vector<double> coordCenter;
  masterElemProvider->coordinate_center(meTopo, coordCenter);
  EXPECT_NEAR(0.0, coordCenter[0], 1.0e-6);
  EXPECT_NEAR(0.0, coordCenter[1], 1.0e-6);
  EXPECT_NEAR(0.25, coordCenter[2], 1.0e-6);

}

TEST(MasterElementProviderIntrepid2Test, Wedge_6)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::string meshDesc = "0,1,WEDGE_6,1,2,3,4,5,6"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,1,1, 1,1,1";

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::WEDGE_6, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));
  unsigned numIntegrationPts = masterElemProvider->num_integration_points(meTopo);
  EXPECT_EQ(6u, numIntegrationPts);

  std::vector<double> integrationPoints;
  masterElemProvider->integration_points(meTopo, integrationPoints);

  for(unsigned i=0; i<numIntegrationPts; i++) {
    unsigned offset = 3*i;

    EXPECT_TRUE(integrationPoints[offset+0] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+0] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+1] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+1] <  1.0);

    EXPECT_TRUE(integrationPoints[offset+2] > -1.0);
    EXPECT_TRUE(integrationPoints[offset+2] <  1.0);
  }

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  EXPECT_EQ(numFieldComponents, meta.spatial_dimension());
  EXPECT_EQ(numNodes, topo.num_nodes());

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, elementCentroid,
                                                  paramCoords, paramDistance);

  EXPECT_NEAR(1.0 / 3.0, paramCoords[0], 1.0e-6);
  EXPECT_NEAR(1.0 / 3.0, paramCoords[1], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[2], 1.0e-6);
  EXPECT_NEAR(0.0, paramDistance, 1.0e-6);

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> result(numFieldComponents, 0.);
  masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result);
  EXPECT_NEAR(5.0, result[0], 1.0e-6);
  EXPECT_NEAR(5.0, result[1], 1.0e-6);
  EXPECT_NEAR(5.0, result[2], 1.0e-6);

  std::vector<double> coordCenter;
  masterElemProvider->coordinate_center(meTopo, coordCenter);
  EXPECT_NEAR(1.0 / 3.0, coordCenter[0], 1.0e-6);
  EXPECT_NEAR(1.0 / 3.0, coordCenter[1], 1.0e-6);
  EXPECT_NEAR(0.0, coordCenter[2], 1.0e-6);

}

TEST(MasterElementProviderIntrepid2Test, UnsupportedTopology)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,SHELL_QUAD_4,1,2,3,4"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0";

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::SHELL_QUAD_4, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));
  EXPECT_NO_THROW(masterElemProvider->num_integration_points(meTopo));

  std::vector<double> integrationPoints;
  EXPECT_NO_THROW(masterElemProvider->integration_points(meTopo, integrationPoints));

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  if constexpr (stk::search::intrepid2_has_param_dist) {
    EXPECT_NO_THROW(masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords,
                                                       elementCentroid, paramCoords, paramDistance));
  }
  else {
    EXPECT_ANY_THROW(masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords,
                                                       elementCentroid, paramCoords, paramDistance));
  }

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> result(numFieldComponents, 0.);
  if constexpr (stk::search::intrepid2_has_param_dist) {
    EXPECT_NO_THROW(masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result));
  }
  else {
    EXPECT_ANY_THROW(masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result));
  }

  std::vector<double> coordCenter;
  EXPECT_NO_THROW(masterElemProvider->coordinate_center(meTopo, coordCenter));
}

TEST(MasterElementProviderIntrepid2Test, FindParametricCoords_Throws)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
    stk::topology topo = bulk->bucket(element).topology();
  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));

  std::vector<double> inputCoords;
  stk::search::determine_centroid(3, element, *coords, inputCoords);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  std::vector<double> elementCoordsTooSmall;
  std::vector<double> inputCoordsTooSmall;
  EXPECT_ANY_THROW(masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, inputCoordsTooSmall,
                                                  paramCoords, paramDistance));
  EXPECT_ANY_THROW(masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoordsTooSmall, inputCoords,
                                                  paramCoords, paramDistance));
  EXPECT_NO_THROW(masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, inputCoords,
                                                  paramCoords, paramDistance));
}

TEST(MasterElementProviderIntrepid2Test, FindParametricCoords_LargerInputs)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
    stk::topology topo = bulk->bucket(element).topology();
  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));

  std::vector<double> inputCoords;
  stk::search::determine_centroid(3, element, *coords, inputCoords);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  elementCoords.push_back(1234.56);
  inputCoords.push_back(6543.21);
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, inputCoords,
                                                  paramCoords, paramDistance);

  EXPECT_NEAR(0.0, paramCoords[0], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[1], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[2], 1.0e-6);
  EXPECT_NEAR(0.0, paramDistance, 1.0e-6);
}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_Throws)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
    stk::topology topo = bulk->bucket(element).topology();
  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));

  std::vector<double> inputCoords;
  stk::search::determine_centroid(3, element, *coords, inputCoords);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, inputCoords,
                                                  paramCoords, paramDistance);

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> fieldDataTooSmall;
  std::vector<double> paramCoordsTooSmall;
  std::vector<double> result(numFieldComponents, 0.);

  EXPECT_ANY_THROW(masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldDataTooSmall, result));
  EXPECT_ANY_THROW(masterElemProvider->evaluate_field(meTopo, paramCoordsTooSmall, numFieldComponents, fieldData, result));
  EXPECT_NO_THROW(masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result));

}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_LargerInputs)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
    stk::topology topo = bulk->bucket(element).topology();
  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));

  std::vector<double> inputCoords;
  stk::search::determine_centroid(3, element, *coords, inputCoords);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, inputCoords,
                                                  paramCoords, paramDistance);

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  fieldData.push_back(1234.56);
  paramCoords.push_back(6543.21);
  std::vector<double> result(numFieldComponents, 0.);

  masterElemProvider->evaluate_field(meTopo, paramCoords, numFieldComponents, fieldData, result);
  EXPECT_NEAR(5.0, result[0], 1.0e-6);
  EXPECT_NEAR(5.0, result[1], 1.0e-6);
  EXPECT_NEAR(5.0, result[2], 1.0e-6);
}

TEST(MasterElementProviderIntrepid2Test, Hex_8_TwoEvalPoints)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) {
    GTEST_SKIP();
  }

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::HEX_8, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  unsigned numFieldComponents;
  unsigned numNodes;
  std::vector<double> elementCoords;

  masterElemProvider->nodal_field_data(elementKey, coords, numFieldComponents, numNodes, elementCoords);

  double paramDistance{0.0};
  std::vector<double> paramCoords;
  masterElemProvider->find_parametric_coordinates(meTopo, numFieldComponents, elementCoords, elementCentroid,
                                                  paramCoords, paramDistance);


  const unsigned numEvalPoints = 2;
  std::vector<double> paramCoordsMultPoints = {paramCoords[0], paramCoords[1], paramCoords[2],
                                               paramCoords[0], paramCoords[1], paramCoords[2]};

  std::vector<double> fieldData(numFieldComponents*numNodes, 5.0);
  std::vector<double> result(numFieldComponents*numEvalPoints, 0.);
  masterElemProvider->evaluate_field(meTopo, numEvalPoints, paramCoordsMultPoints, numFieldComponents, fieldData, result);
  EXPECT_NEAR(5.0, result[0], 1.0e-6);
  EXPECT_NEAR(5.0, result[1], 1.0e-6);
  EXPECT_NEAR(5.0, result[2], 1.0e-6);
  EXPECT_NEAR(5.0, result[3], 1.0e-6);
  EXPECT_NEAR(5.0, result[4], 1.0e-6);
  EXPECT_NEAR(5.0, result[5], 1.0e-6);

}

void initialize_nodal_field_using_nodeIds(const stk::mesh::BulkData& mesh,
                                          const stk::mesh::FieldBase& field)
{
  auto fieldData = field.template data<double,stk::mesh::ReadWrite>();
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  stk::mesh::for_each_entity_run(mesh, stk::topology::NODE_RANK, meta.universal_part(),
  [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity node) {
    auto fieldDataValues = fieldData.entity_values(node);
    for(stk::mesh::ScalarIdx scalarIdx : fieldDataValues.scalars()) {
      fieldDataValues(scalarIdx) = bulk.identifier(node);
    }
  });
}

void test_evaluate_field_nodal(const std::string& meshDesc,
                         const std::vector<double>& point,
                         double expectedNodalScalarFieldValue,
                         const std::vector<double>& expectedNodalVectorFieldValues)
{
  MPI_Comm comm = stk::parallel_machine_world();
  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  auto& nodalScalarField = meta.declare_field<double>(stk::topology::NODE_RANK, "nodalScalarField");
  auto& nodalVectorField = meta.declare_field<double>(stk::topology::NODE_RANK, "nodalVectorField");
  stk::mesh::put_field_on_mesh(nodalScalarField, meta.universal_part(), 1, nullptr);
  const size_t numComponents = expectedNodalVectorFieldValues.size();
  stk::mesh::put_field_on_mesh(nodalVectorField, meta.universal_part(), numComponents, nullptr);

  stk::unit_test_util::setup_text_mesh(mesh, meshDesc);

  initialize_nodal_field_using_nodeIds(mesh, nodalScalarField);
  initialize_nodal_field_using_nodeIds(mesh, nodalVectorField);

  //hard-coded! assume element 1
  stk::mesh::Entity elem1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  ASSERT_TRUE(mesh.is_valid(elem1));

  unsigned spatialDim = meta.spatial_dimension();

  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  std::vector<double> elemNodeCoords;
  stk::search::SearchField coordField(meta.coordinate_field());
  stk::search::spmd::EntityKeyPair entityKeyPair(elem1, mesh.entity_key(elem1));
  unsigned numNodes = 0;

  masterElemProvider->nodal_field_data(entityKeyPair, coordField,
                                       spatialDim, numNodes, elemNodeCoords);

  stk::topology stkTopo = mesh.bucket(elem1).topology();

  std::vector<double> parCoords = {99.9, 99.9, 99.9};
  stk::search::SearchTopology srchTopo(stkTopo);
  double parametricDistance = 999.9;

  masterElemProvider->find_parametric_coordinates(srchTopo, spatialDim,
                                             elemNodeCoords, point, parCoords,
                                             parametricDistance);
  //scalar field
  stk::search::SearchField scalarField(&nodalScalarField);
  unsigned numFieldComponents = 0;
  std::vector<double> scalarFieldData;

  masterElemProvider->nodal_field_data(entityKeyPair, scalarField,
                                       numFieldComponents, numNodes, scalarFieldData);
  std::vector<double> resultScalarFieldData = {0.0};

  masterElemProvider->evaluate_field(srchTopo, parCoords,
                                     1, scalarFieldData, resultScalarFieldData);

  EXPECT_DOUBLE_EQ(expectedNodalScalarFieldValue, resultScalarFieldData[0]);

  //vector field
  stk::search::SearchField vectorField(&nodalVectorField);
  std::vector<double> vectorFieldData;

  masterElemProvider->nodal_field_data(entityKeyPair, vectorField,
                                       numFieldComponents, numNodes, vectorFieldData);
  std::vector<double> resultVectorFieldData(expectedNodalVectorFieldValues.size());

  masterElemProvider->evaluate_field(srchTopo, parCoords,
                                     numComponents, vectorFieldData,
                                     resultVectorFieldData);

  for(size_t i=0; i<expectedNodalVectorFieldValues.size(); ++i) {
    EXPECT_DOUBLE_EQ(expectedNodalVectorFieldValues[i], resultVectorFieldData[i]);
  }
}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_hex8)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";

  std::vector<double> point = {0.99, 0.99, 0.99};

  double expectedNodalScalarFieldValue = 6.9598;
  std::vector<double> expectedNodalVectorFieldValues = {6.9598, 6.9598, 6.9598};

  test_evaluate_field_nodal(meshDesc, point,
                            expectedNodalScalarFieldValue,
                            expectedNodalVectorFieldValues);
}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_tet4)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TET_4,1,2,3,4"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1";

  std::vector<double> pointInside = {0.2, 0.2, 0.2};
  std::vector<double> pointOnFace = {0.0, 0.5, 0.0};

  double expectedScalarValueInside = 2.2;
  double expectedScalarValueOnFace = 2.0;
  std::vector<double> expectedVectorValuesInside = {2.2, 2.2, 2.2};
  std::vector<double> expectedVectorValuesOnFace = {2.0, 2.0, 2.0};

  test_evaluate_field_nodal(meshDesc, pointInside, expectedScalarValueInside,
                            expectedVectorValuesInside);
  test_evaluate_field_nodal(meshDesc, pointOnFace, expectedScalarValueOnFace,
                            expectedVectorValuesOnFace);
}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_pyramid5)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0.5,0.5,1";

  std::vector<double> pointInside = {0.2, 0.2, 0.2};
  std::vector<double> pointOnTip = {0.5, 0.5, 1.0};

  double expectedScalarValueInside = 2.175;
  double expectedScalarValueOnTip = 5.0;
  std::vector<double> expectedVectorValuesInside = {2.175, 2.175, 2.175};
  std::vector<double> expectedVectorValuesOnTip = {5.0, 5.0, 5.0};

  test_evaluate_field_nodal(meshDesc, pointInside, expectedScalarValueInside,
                            expectedVectorValuesInside);
  test_evaluate_field_nodal(meshDesc, pointOnTip, expectedScalarValueOnTip,
                            expectedVectorValuesOnTip);
}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_wedge6)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,WEDGE_6,1,2,3,4,5,6"
        "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,0,1, 1,0,1, 1,1,1";

  std::vector<double> pointInside = {0.7, 0.5, 0.5};
  std::vector<double> pointOnFace = {0.5, 0.5, 0.5};

  double expectedScalarValueInside = 3.7;
  double expectedScalarValueOnFace = 3.5;
  std::vector<double> expectedVectorValuesInside = {3.7, 3.7, 3.7};
  std::vector<double> expectedVectorValuesOnFace = {3.5, 3.5, 3.5};

  test_evaluate_field_nodal(meshDesc, pointInside, expectedScalarValueInside,
                            expectedVectorValuesInside);
  test_evaluate_field_nodal(meshDesc, pointOnFace, expectedScalarValueOnFace,
                            expectedVectorValuesOnFace);
}

void test_evaluate_field_nodal_two_points(const std::string& meshDesc,
                         const std::vector<double>& point1,
                         const std::vector<double>& point2,
                         double expectedNodalScalarFieldValue1,
                         double expectedNodalScalarFieldValue2)
{
  MPI_Comm comm = stk::parallel_machine_world();
  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  auto& nodalScalarField = meta.declare_field<double>(stk::topology::NODE_RANK, "nodalScalarField");
  stk::mesh::put_field_on_mesh(nodalScalarField, meta.universal_part(), 1, nullptr);

  stk::unit_test_util::setup_text_mesh(mesh, meshDesc);

  initialize_nodal_field_using_nodeIds(mesh, nodalScalarField);

  //hard-coded! assume element 1
  stk::mesh::Entity elem1 = mesh.get_entity(stk::topology::ELEM_RANK, 1);
  ASSERT_TRUE(mesh.is_valid(elem1));

  unsigned spatialDim = meta.spatial_dimension();

  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::search::MasterElementProviderIntrepid2>();

  std::vector<double> elemNodeCoords;
  stk::search::SearchField coordField(meta.coordinate_field());
  stk::search::spmd::EntityKeyPair entityKeyPair(elem1, mesh.entity_key(elem1));
  unsigned numNodes = 0;

  masterElemProvider->nodal_field_data(entityKeyPair, coordField,
                                       spatialDim, numNodes, elemNodeCoords);

  stk::topology stkTopo = mesh.bucket(elem1).topology();

  std::vector<double> parCoords1 = {99.9, 99.9, 99.9};
  stk::search::SearchTopology srchTopo(stkTopo);
  double parametricDistance1 = 999.9;

  masterElemProvider->find_parametric_coordinates(srchTopo, spatialDim,
                                             elemNodeCoords, point1, parCoords1,
                                             parametricDistance1);

  std::vector<double> parCoords2 = {99.9, 99.9, 99.9};
  double parametricDistance2 = 999.9;

  masterElemProvider->find_parametric_coordinates(srchTopo, spatialDim,
                                             elemNodeCoords, point2, parCoords2,
                                             parametricDistance2);
  //scalar field
  stk::search::SearchField scalarField(&nodalScalarField);
  unsigned numFieldComponents = 0;
  std::vector<double> scalarFieldData;

  masterElemProvider->nodal_field_data(entityKeyPair, scalarField,
                                       numFieldComponents, numNodes, scalarFieldData);

  std::vector<double> resultScalarFieldData = {0.0};
  std::vector<double> parCoordsCombined = {parCoords1[0], parCoords1[1], parCoords1[2],
                                           parCoords2[0], parCoords2[1], parCoords2[2]};

  masterElemProvider->evaluate_field(srchTopo, 2, parCoordsCombined,
                                     1, scalarFieldData, resultScalarFieldData);

  EXPECT_DOUBLE_EQ(expectedNodalScalarFieldValue1, resultScalarFieldData[0]);
  EXPECT_DOUBLE_EQ(expectedNodalScalarFieldValue2, resultScalarFieldData[1]);

}

TEST(MasterElementProviderIntrepid2Test, EvaluateField_tet4_twoEvalPoints)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshDesc = "0,1,TET_4,1,2,3,4"
        "|coordinates:   0,0,0, 1,0,0, 0,1,0, 0,0,1";

  std::vector<double> pointInside = {0.2, 0.2, 0.2};
  std::vector<double> pointOnFace = {0.0, 0.5, 0.0};

  double expectedScalarValueInside = 2.2;
  double expectedScalarValueOnFace = 2.0;
  std::vector<double> expectedVectorValuesInside = {2.2, 2.2, 2.2};
  std::vector<double> expectedVectorValuesOnFace = {2.0, 2.0, 2.0};

  test_evaluate_field_nodal_two_points(meshDesc, pointInside, pointOnFace,
                                       expectedScalarValueInside, expectedScalarValueOnFace);
}

}
