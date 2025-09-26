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
#include "stk_mesh/base/Field.hpp"                         // for Field
#include "stk_mesh/base/FieldBLAS.hpp"                     // for field_copy
#include "stk_search/FilterCoarseSearch.hpp"               // for ObjectOuts...
#include "stk_search_util/spmd/GeometricSearch.hpp"
#include "stk_search_util/spmd/GeometricSearchDispatch.hpp"
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

TEST(MockSearchTest, expandBoundingBox)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    double expansionFactor = 1.5;
    double expansionSum = 0.5;

    double radius = 1.0;
    stk::search::Point<double> center(0, 0, 0);
    stk::search::Sphere<double> sphere(center, radius);

    stk::search::scale_by(sphere, expansionFactor, expansionSum);
    double expectedRadius = 2.0;
    EXPECT_EQ(expectedRadius, sphere.radius());

    stk::search::Point<double> minCorner(-1, -1, -1);
    stk::search::Point<double> maxCorner(1, 1, 1);
    stk::search::Box<double> box(minCorner, maxCorner);

    stk::search::Point<double> expectedMinCorner(-2, -2, -2);
    stk::search::Point<double> expectedMaxCorner(2, 2, 2);
    stk::search::Box<double> expectedBox(expectedMinCorner, expectedMaxCorner);
    stk::search::scale_by(box, expansionFactor, expansionSum);
    EXPECT_EQ(expectedBox, box);
  }
}

TEST(MockSearchTest, computeCentroid)
{
  double radius = 1.0;
  stk::search::Point<double> center(0.5, 1.0, 2.0);
  stk::search::Sphere<double> sphere(center, radius);

  stk::search::Point<double> centroid = stk::search::spmd::impl::compute_centroid(sphere);
  EXPECT_EQ(center, centroid);

  stk::search::Point<double> minCorner(-1, -1, -1);
  stk::search::Point<double> maxCorner(1, 1, 1);
  stk::search::Box<double> box(minCorner, maxCorner);

  stk::search::Point<double> expectedCentroid(0, 0, 0);
  centroid = stk::search::spmd::impl::compute_centroid(box);
  EXPECT_EQ(expectedCentroid, centroid);
}

TEST(MockSearchTest, comparisonLessThan)
{
  double d1 = 8.5;
  double d2 = 9;
  double epsilon = 1.0;

  EXPECT_FALSE(stk::search::less_than(d1, d2, epsilon));
  EXPECT_FALSE(stk::search::less_than(d2, d1, epsilon));

  epsilon = 0.1;
  EXPECT_TRUE(stk::search::less_than(d1, d2, epsilon));
  EXPECT_FALSE(stk::search::less_than(d2, d1, epsilon));
}

TEST(MockSearchTest, createSearchDispatch)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  using STKMockSearch   = stk::search::spmd::GeometricSearch<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;
  using STKMockDispatch = stk::search::spmd::GeometricSearchDispatch<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;

  double parametricTolerance = 0.001;
  double geometricTolerance = 0.1;
  double expansionSum = geometricTolerance;
  double expansionFactor = 0.0;
  double slantHeight = 1.0;
  const stk::search::SearchMethod search_method = stk::search::KDTREE;

  std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_slanted_single_hex_bulk(slantHeight);
  auto sendMesh = stk::unit_test_util::construct_hex_send_mesh(*sendBulk, parametricTolerance);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_single_point_bulk(0.5, 0.5, 0.5);
  auto recvMesh = stk::unit_test_util::construct_node_recv_mesh(*recvBulk,
                                                                parametricTolerance, geometricTolerance);

  sendMesh->set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy::ABORT);

  auto searchObj = std::make_shared<STKMockSearch>(sendMesh, recvMesh, "searchObj",
                                                   MPI_COMM_WORLD, expansionFactor,
                                                   expansionSum, search_method);

  auto searchDispatch = std::make_shared<STKMockDispatch>(searchObj);

  searchDispatch->set_print_search_warnings(false);
  searchDispatch->set_do_initial_search_expansion(false);
  searchDispatch->set_output_stream(std::cerr);
  searchDispatch->set_expansion_factor(0.1);

  EXPECT_NO_THROW(searchObj->initialize());
  EXPECT_NO_THROW(searchObj->coarse_search());

  EXPECT_TRUE(searchObj->get_unpaired_recv_entities().empty());
  EXPECT_EQ(1u, searchObj->get_range_to_domain().size());
}

TEST(MockSearchTest, unpairedEntities)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  using STKMockSearch = stk::search::spmd::GeometricSearch<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;

  double parametricTolerance = 0.001;
  double geometricTolerance = 0.1;
  double expansionSum = geometricTolerance;
  double expansionFactor = 0.0;
  double slantHeight = 1.0;
  const stk::search::SearchMethod search_method = stk::search::KDTREE;

  std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_slanted_single_hex_bulk(slantHeight);
  auto sendMesh = stk::unit_test_util::construct_hex_send_mesh(*sendBulk, parametricTolerance);

  std::shared_ptr<stk::mesh::BulkData> recvBulkInside = stk::unit_test_util::build_single_point_bulk(0.5, 0.5, 0.5);
  auto destMeshInside = stk::unit_test_util::construct_node_recv_mesh(*recvBulkInside,
                                                                      parametricTolerance, geometricTolerance);

  std::shared_ptr<stk::mesh::BulkData> recvBulkInTolerance = stk::unit_test_util::build_single_point_bulk(0.5, 0.5, 1.09);
  auto destMeshInTolerance = stk::unit_test_util::construct_node_recv_mesh(*recvBulkInTolerance,
                                                                           parametricTolerance, geometricTolerance);

  std::shared_ptr<stk::mesh::BulkData> recvBulkOutside = stk::unit_test_util::build_single_point_bulk(0.5, 0.5, 1.5);
  auto destMeshOutside = stk::unit_test_util::construct_node_recv_mesh(*recvBulkOutside,
                                                                       parametricTolerance, geometricTolerance);

  sendMesh->set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy::EXTRAPOLATE);

  STKMockSearch searchObjInDomain(sendMesh, destMeshInside, "searchObjInDomain", MPI_COMM_WORLD, expansionFactor,
                                  expansionSum, search_method);

  searchObjInDomain.set_print_search_warnings(false);
  searchObjInDomain.coarse_search_for_objects_inside_domain();
  EXPECT_TRUE(searchObjInDomain.get_unpaired_recv_entities().empty());
  EXPECT_EQ(1u, searchObjInDomain.get_range_to_domain().size());

  STKMockSearch searchObjInTolerance(sendMesh, destMeshInTolerance, "searchObjInTolerance", MPI_COMM_WORLD,
                                     expansionFactor, expansionSum, search_method);

  searchObjInTolerance.set_print_search_warnings(false);
  searchObjInTolerance.coarse_search_for_objects_inside_domain();
  EXPECT_TRUE(searchObjInTolerance.get_unpaired_recv_entities().empty());
  EXPECT_EQ(1u, searchObjInTolerance.get_range_to_domain().size());

  STKMockSearch searchObjOutside(sendMesh, destMeshOutside, "searchObjOutside", MPI_COMM_WORLD, expansionFactor,
                                 expansionSum, search_method);

  searchObjOutside.set_print_search_warnings(false);
  searchObjOutside.coarse_search_for_objects_inside_domain();
  EXPECT_EQ(1u, searchObjOutside.get_unpaired_recv_entities().size());
  EXPECT_EQ(0u, searchObjOutside.get_range_to_domain().size());

  searchObjOutside.coarse_search_for_objects_outside_domain();
  EXPECT_EQ(0u, searchObjOutside.get_unpaired_recv_entities().size());
  EXPECT_EQ(1u, searchObjOutside.get_range_to_domain().size());
}

TEST(MockSearchTest, filterToNearestWithAbort)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  using STKMockSearch = stk::search::spmd::GeometricSearch<stk::search::spmd::ElementSendMesh, stk::search::spmd::NodeRecvMesh>;

  double x = 1.5;
  double y = 0.5;
  double z = 0.5;
  double expansionFactor = 0.0;
  double parametricTolerance = 0.001;
  double geometricTolerance = 0.1;

  std::shared_ptr<stk::mesh::BulkData> sendBulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *sendBulk);

  auto sendMesh = stk::unit_test_util::construct_hex_send_mesh(*sendBulk, parametricTolerance);
  sendMesh->set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy::ABORT);
  sendMesh->set_name("sendMesh");

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_single_point_bulk(x, y, z);
  auto destMesh = stk::unit_test_util::construct_node_recv_mesh(*recvBulk, parametricTolerance, geometricTolerance);
  destMesh->set_name("recvMesh");

  STKMockSearch search(sendMesh, destMesh, "mocktest", MPI_COMM_WORLD, expansionFactor, 0.0);
  search.set_print_search_warnings(false);
  search.initialize();

  search.coarse_search();
  EXPECT_THROW(search.local_search(), std::runtime_error);
}

TEST(MockSearchTest, centroidPointEvaluator)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator =
      std::make_shared<stk::search::CentroidEvaluator>(*bulk, coords);

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::HEX_8, topo);

  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  EXPECT_NEAR(0.5, elementCentroid[0], 1.0e-6);
  EXPECT_NEAR(0.5, elementCentroid[1], 1.0e-6);
  EXPECT_NEAR(0.5, elementCentroid[2], 1.0e-6);

  std::vector<double> evalPoint;
  size_t numPts = pointEvaluator->num_points(elementKey, topo);
  EXPECT_EQ(1u, numPts);

  pointEvaluator->coordinates(elementKey, 0, evalPoint);

  EXPECT_NEAR(elementCentroid[0], evalPoint[0], 1.0e-6);
  EXPECT_NEAR(elementCentroid[1], evalPoint[1], 1.0e-6);
  EXPECT_NEAR(elementCentroid[2], evalPoint[2], 1.0e-6);
}

TEST(MockSearchTest, gaussPointEvaluator)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>(0);
  std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator =
      std::make_shared<stk::search::MasterElementGaussPointEvaluator>(*bulk, coords, masterElemProvider);

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));

  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));
  stk::topology topo = bulk->bucket(element).topology();
  EXPECT_EQ(stk::topology::HEX_8, topo);

  stk::search::SearchTopology meTopo(topo, elementKey, bulk->bucket_ptr(element));
  unsigned numGaussPts = masterElemProvider->num_integration_points(meTopo);
  EXPECT_EQ(8u, numGaussPts);

  std::vector<double> gaussPoints;
  masterElemProvider->integration_points(meTopo, gaussPoints);

  for(unsigned i=0; i<numGaussPts; i++) {
    // Bounded by [-1,1] x [-1,1] x [-1,1]
    unsigned offset = 3*i;

    EXPECT_TRUE(gaussPoints[offset+0] > -1.0);
    EXPECT_TRUE(gaussPoints[offset+0] <  1.0);

    EXPECT_TRUE(gaussPoints[offset+1] > -1.0);
    EXPECT_TRUE(gaussPoints[offset+1] <  1.0);

    EXPECT_TRUE(gaussPoints[offset+2] > -1.0);
    EXPECT_TRUE(gaussPoints[offset+2] <  1.0);
  }

  std::vector<double> evalPoint;

  size_t numPts = pointEvaluator->num_points(elementKey, topo);
  EXPECT_EQ(numGaussPts, numPts);

  for(size_t i=0; i<numPts; i++) {
    pointEvaluator->coordinates(elementKey, i, evalPoint);

    // Point is bounded by [0,1] x [0,1] x [0,1]
    EXPECT_TRUE(evalPoint[0] > 0.0);
    EXPECT_TRUE(evalPoint[0] < 1.0);

    EXPECT_TRUE(evalPoint[1] > 0.0);
    EXPECT_TRUE(evalPoint[1] < 1.0);

    EXPECT_TRUE(evalPoint[2] > 0.0);
    EXPECT_TRUE(evalPoint[2] < 1.0);
  }
}

TEST(MockSearchTest, findParametricCoordinates)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    GTEST_SKIP();
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, MPI_COMM_WORLD);;
  stk::io::fill_mesh("generated:1x1x1", *bulk);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  double parametricTolerance = 0.01;

  const stk::mesh::FieldBase* coords = meta.coordinate_field();
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider =
      std::make_shared<stk::unit_test_util::MasterElementProvider>(0);
  std::shared_ptr<stk::search::MasterElementParametricCoordsFinder> paramCoordsFinder =
      std::make_shared<stk::search::MasterElementParametricCoordsFinder>(*bulk, coords, masterElemProvider, parametricTolerance);

  stk::mesh::Entity element = bulk->get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulk->is_valid(element));
  stk::search::spmd::EntityKeyPair elementKey(stk::search::spmd::make_entity_key_pair(bulk,element));

  bool isWithinParametricTolerance{false};
  double paramDistance{0.0};
  std::vector<double> paramCoords;
  std::vector<double> elementCentroid;
  stk::search::determine_centroid(3, element, *coords, elementCentroid);

  // Interval [-1,1] x [-1,1] x [-1,1] -> centroid: (xi,eta,zeta) = (0,0,0)
  paramCoordsFinder->find_parametric_coords(elementKey, elementCentroid,
                                            paramCoords, paramDistance, isWithinParametricTolerance);

  EXPECT_TRUE(isWithinParametricTolerance);
  EXPECT_NEAR(0.0, paramCoords[0], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[1], 1.0e-6);
  EXPECT_NEAR(0.0, paramCoords[2], 1.0e-6);
}

class LinearPeriodicSearch : public ::testing::Test
{
 public:
  using STKNodeToNodeSearch = stk::search::spmd::GeometricSearch<stk::search::spmd::NodeSendMesh, stk::search::spmd::NodeRecvMesh>;
  using NodeToNodeRelationVec = STKNodeToNodeSearch::EntityProcRelationVec;
  using NodeToNodeRelation = STKNodeToNodeSearch::EntityProcRelation;

  LinearPeriodicSearch(stk::ParallelMachine comm, unsigned elemsPerDim,
                       double parametricTolerance = 1.0e-12,
                       double coarseSearchTolerance = 1.0e-12)
  : m_comm(comm)
  , m_numElemPerAxis(elemsPerDim)
  , m_parametricTolerance(parametricTolerance)
  , m_coarseSearchTolerance(coarseSearchTolerance)
  { }

  LinearPeriodicSearch() = default;

  void displace_mesh(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase const* coord, const std::vector<double>& disp)
  {
    const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    STK_ThrowRequireMsg(disp.size() >= meta.spatial_dimension(), "Displacement vector is of insufficient length");

    stk::mesh::EntityVector entities;
    stk::mesh::Selector selector(meta.universal_part());
    stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK), entities);

    for(stk::mesh::Entity entity : entities) {
      double* coordData = static_cast<double*>(stk::mesh::field_data(*coord, entity));
      for(unsigned i = 0; i < meta.spatial_dimension(); ++i) {
        coordData[i] += disp[i];
      }
    }
  }

  std::string get_3d_hex8_generated_mesh_spec()
  {
    unsigned nProc = stk::parallel_machine_size(m_comm);
    std::ostringstream oss;
    oss << "generated:"
        << m_numElemPerAxis << "x" << m_numElemPerAxis << "x" << m_numElemPerAxis * nProc
        << "|bbox:0,0,0,1,1," << nProc
        << "|sideset:xX";
    std::string meshSpec(oss.str());
    return meshSpec;
  }

  std::shared_ptr<stk::mesh::BulkData> create_send_bulk()
  {
    const std::string meshSpec(get_3d_hex8_generated_mesh_spec());

    std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, m_comm);
    stk::io::fill_mesh(meshSpec, *bulk);

    return bulk;
  }

  std::shared_ptr<stk::mesh::BulkData> create_recv_bulk(const std::string& periodicCoordsName)
  {
    const std::string meshSpec(get_3d_hex8_generated_mesh_spec());

    std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, m_comm);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Field<double> &field = meta.declare_field<double>(stk::topology::NODE_RANK, periodicCoordsName, 1);

    const std::vector<double> init(meta.spatial_dimension(), 0);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), meta.spatial_dimension(), init.data());

    stk::io::fill_mesh(meshSpec, *bulk);

    return bulk;
  }

  std::shared_ptr<stk::search::spmd::NodeSendMesh> create_send_mesh()
  {
    const stk::mesh::MetaData& srcMeta = m_sendBulk->mesh_meta_data();
    const stk::mesh::FieldBase* srcCoords = srcMeta.coordinate_field();

    stk::mesh::Part* srcPart = srcMeta.get_part("surface_1");
    EXPECT_TRUE(srcPart != nullptr);
    stk::mesh::PartVector srcParts{ srcPart };
    std::shared_ptr<stk::search::spmd::NodeSendMesh> srcMesh(
        new stk::search::spmd::NodeSendMesh(m_sendBulk.get(), srcCoords, srcParts, m_comm, m_parametricTolerance));
    srcMesh->set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy::IGNORE);

    return srcMesh;
  }

  std::shared_ptr<stk::search::spmd::NodeRecvMesh>
  create_recv_mesh(const std::string& periodicCoordsName)
  {
    const stk::mesh::MetaData& destMeta = m_recvBulk->mesh_meta_data();
    stk::mesh::Part* destPart = destMeta.get_part("surface_2");
    EXPECT_TRUE(destPart != nullptr);
    stk::mesh::PartVector destParts{ destPart };
    const stk::mesh::FieldBase* destCoords = destMeta.coordinate_field();

    stk::mesh::FieldBase* transformedCoords = destMeta.get_field(stk::topology::NODE_RANK, periodicCoordsName);
    stk::mesh::field_copy(*destCoords, *transformedCoords);

    // Shift mesh so that recv.surface_2 aligns with send.surface_1
    displace_mesh(*m_recvBulk, transformedCoords, std::vector<double>{ -1.0, 0.0, 0.0 });

    std::shared_ptr<stk::search::spmd::NodeRecvMesh>
    destMesh(new stk::search::spmd::NodeRecvMesh(m_recvBulk.get(), transformedCoords, destParts, m_comm,
                                                 m_parametricTolerance, m_coarseSearchTolerance));

    return destMesh;
  }

  void create_search(const std::string& periodicCoordsName)
  {
    m_sendBulk = create_send_bulk();
    m_sendMesh = create_send_mesh();

    m_recvBulk = create_recv_bulk(periodicCoordsName);
    m_recvMesh = create_recv_mesh(periodicCoordsName);

    m_search = std::make_shared<STKNodeToNodeSearch>(m_sendMesh, m_recvMesh, "searchTest", m_comm, 0.0, 1e-6);
  }

 protected:
  stk::ParallelMachine m_comm{MPI_COMM_WORLD};
  unsigned m_numElemPerAxis{1};
  double m_parametricTolerance = 1.0e-6;
  double m_coarseSearchTolerance = 1.0e-6;

  std::shared_ptr<stk::mesh::BulkData> m_sendBulk;
  std::shared_ptr<stk::mesh::BulkData> m_recvBulk;

  std::shared_ptr<stk::search::spmd::NodeSendMesh> m_sendMesh;
  std::shared_ptr<stk::search::spmd::NodeRecvMesh> m_recvMesh;

  std::shared_ptr<STKNodeToNodeSearch> m_search;
};

TEST_F(LinearPeriodicSearch, periodicSearch)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    return;
  }

  const std::string periodicCoordsName("transformed_coords");
  m_numElemPerAxis = 5;

  create_search(periodicCoordsName);

  m_search->set_closest_bounding_box_using_nearest_node(false);
  m_search->initialize();
  m_search->coarse_search();
  EXPECT_NO_THROW(m_search->local_search());

  EXPECT_EQ( (m_numElemPerAxis+1)*(m_numElemPerAxis+1), m_search->get_range_to_domain().size());
  EXPECT_EQ( 0u, m_search->get_unpaired_recv_entities().size() );

  const NodeToNodeRelationVec& relation = m_search->get_range_to_domain();

  for(const NodeToNodeRelation& thisRelation : relation) {
    EXPECT_EQ(thisRelation.first.id().rank(), stk::topology::NODE_RANK);
    EXPECT_EQ(thisRelation.second.id().rank(), stk::topology::NODE_RANK);
    unsigned offset = (thisRelation.first.id().id() - thisRelation.second.id().id());
    EXPECT_EQ(m_numElemPerAxis, offset);
  }
}

}
