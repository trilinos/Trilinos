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
#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/timer.hpp>

namespace perf_test {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using HostSpace = Kokkos::DefaultHostExecutionSpace;

using ElemIdentProc = stk::search::IdentProc<unsigned,int>;
using NodeIdentProc = stk::search::IdentProc<stk::mesh::EntityId,int>;
using SphereIdentProc = stk::search::BoxIdentProc<stk::search::Sphere<double>,ElemIdentProc>;
using PointIdentProc = stk::search::BoxIdentProc<stk::search::Point<double>,NodeIdentProc>;
using Intersection = stk::search::IdentProcIntersection<ElemIdentProc,NodeIdentProc>;

using DomainViewType = Kokkos::View<SphereIdentProc*,ExecSpace>;
using RangeViewType = Kokkos::View<PointIdentProc*,ExecSpace>;
using ResultViewType = Kokkos::View<Intersection*,ExecSpace>;

using FastMeshIndicesViewType = Kokkos::View<stk::mesh::FastMeshIndex*,ExecSpace>;

constexpr unsigned maxNumNeighbors = 16; //we're only expecting 8 per element

DomainViewType create_elem_spheres(const stk::mesh::BulkData& mesh, double radius)
{
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  const unsigned numLocalElems = stk::mesh::count_entities(mesh, stk::topology::ELEM_RANK, meta.locally_owned_part());
  DomainViewType elemSpheres("elemSpheres", numLocalElems);

  const stk::mesh::FieldBase* coordField = meta.coordinate_field();
  const stk::mesh::NgpField<double>& ngpCoords = stk::mesh::get_updated_ngp_field<double>(*coordField);
  const stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);
  const int myRank = mesh.parallel_rank();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part(),
  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elemIndex) {
      stk::mesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemIndex);
      stk::mesh::Entity elemEntity = ngpMesh.get_entity(stk::topology::ELEM_RANK, elemIndex);
      unsigned elemLocalId = ngpMesh.local_id(elemEntity);
      stk::search::Point<double> center(0,0,0);
      for(unsigned j=0; j<nodes.size(); ++j) {
        stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[j]);
        stk::mesh::EntityFieldData<double> coords = ngpCoords(nodeIndex);
        center[0] += coords[0];
        center[1] += coords[1];
        center[2] += coords[2];
      }

      center[0] /= nodes.size();
      center[1] /= nodes.size();
      center[2] /= nodes.size();

      elemSpheres(elemLocalId) = SphereIdentProc{stk::search::Sphere<double>(center, radius), ElemIdentProc(elemEntity.local_offset(), myRank)};
    }
  );
  
  return elemSpheres;
}

RangeViewType create_node_points(const stk::mesh::BulkData& mesh)
{
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  const unsigned numLocalNodes = stk::mesh::count_entities(mesh, stk::topology::NODE_RANK, meta.locally_owned_part());
  const int myRank = mesh.parallel_rank();
  RangeViewType nodePoints("nodePoints", numLocalNodes);

  const stk::mesh::FieldBase* coordField = meta.coordinate_field();
  const stk::mesh::NgpField<double>& ngpCoords = stk::mesh::get_updated_ngp_field<double>(*coordField);
  const stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(mesh);

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, mesh.mesh_meta_data().locally_owned_part(),
  KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& nodeIndex) {
      stk::mesh::EntityFieldData<double> coords = ngpCoords(nodeIndex);
      stk::mesh::Entity node = ngpMesh.get_entity(stk::topology::NODE_RANK, nodeIndex);
      unsigned nodeLocalId = ngpMesh.local_id(node);
      nodePoints(nodeLocalId) = PointIdentProc{stk::search::Point<double>(coords[0], coords[1], coords[2]), NodeIdentProc(ngpMesh.identifier(node), myRank)};
    }
  );
  
  return nodePoints;
}

void ghost_node_neighbors_to_elements(stk::mesh::BulkData& mesh, const ResultViewType& searchResults)
{
  auto hostSpace = HostSpace{};
  auto hostSearchResults = Kokkos::create_mirror_view_and_copy(hostSpace, searchResults);

  mesh.modification_begin();
  stk::mesh::Ghosting& neighborGhosting = mesh.create_ghosting("neighbors");
  std::vector<stk::mesh::EntityProc> nodesToGhost;

  const int myRank = mesh.parallel_rank();

  for(size_t i=0; i<hostSearchResults.size(); ++i) {
    auto result = hostSearchResults(i);
    if (result.domainIdentProc.proc() != myRank && result.rangeIdentProc.proc() == myRank) {
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, result.rangeIdentProc.id());
      nodesToGhost.emplace_back(node, result.domainIdentProc.proc());
    }
  }

  mesh.change_ghosting(neighborGhosting, nodesToGhost);
  mesh.modification_end();
}

void unpack_search_results_into_field(stk::mesh::BulkData& mesh,
                                      stk::mesh::Field<unsigned>& neighborField,
                                      const ResultViewType& searchResults)
{
  auto hostSpace = HostSpace{};
  auto hostSearchResults = Kokkos::create_mirror_view_and_copy(hostSpace, searchResults);

  const int myRank = mesh.parallel_rank();

  for(size_t i=0; i<hostSearchResults.size(); ++i) {
    auto result = hostSearchResults(i);
    if (result.domainIdentProc.proc() == myRank) {
      stk::mesh::Entity elem(result.domainIdentProc.id());
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, result.rangeIdentProc.id());
      ASSERT_TRUE(mesh.is_valid(node));
      unsigned* neighborData = stk::mesh::field_data(neighborField, elem);
      unsigned& numNeighbors = neighborData[0];
      ASSERT_TRUE(numNeighbors < maxNumNeighbors);
      neighborData[1+numNeighbors] = node.local_offset();
      ++numNeighbors;
    }
  }
}

void verify_8_neighbors_per_element(const stk::mesh::BulkData& mesh,
                                    const stk::mesh::Field<unsigned>& neighborField)
{
  stk::mesh::for_each_entity_run(mesh, stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part(),
  [&](const stk::mesh::BulkData& /*bulk*/, stk::mesh::Entity elem) {
    const unsigned* neighborData = stk::mesh::field_data(neighborField, elem);
    EXPECT_EQ(8u, neighborData[0]);
  });
}

TEST(PerfNgpSearch, elemNodeNeighbors)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) > 4) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD)
                                     .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                     .set_spatial_dimension(3)
                                     .set_maintain_local_ids(true)
                                     .create();
  stk::mesh::MetaData& meta = meshPtr->mesh_meta_data();

  std::string meshSpec("generated:80x80x80|bbox:-20,-20,-20,20,20,20");
  const double radius = 0.5;
  //elems are 0.5 cubes, so radius 0.5 from the center should find 8 nodes per hex

  stk::mesh::Field<unsigned>& neighborField = meta.declare_field<unsigned>(stk::topology::ELEM_RANK, "nodeNeighbors");
  stk::mesh::put_field_on_mesh(neighborField, meta.universal_part(), maxNumNeighbors+1, nullptr);

  stk::io::fill_mesh(meshSpec, *meshPtr);

  std::vector<size_t> globalCounts;
  stk::mesh::comm_mesh_counts(*meshPtr, globalCounts);
  if (stk::parallel_machine_rank(comm) == 0) {
    std::cout<<"using mesh " << meshSpec << " with " << globalCounts[stk::topology::ELEM_RANK] << " elements" << std::endl;
  }

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  const unsigned NUM_RUNS = 50;

  const unsigned numLocalElems = stk::mesh::count_entities(*meshPtr, stk::topology::ELEM_RANK, meta.locally_owned_part());
  const unsigned numLocalOwnedNodes = stk::mesh::count_entities(*meshPtr, stk::topology::NODE_RANK, meta.locally_owned_part());
  stk::mesh::Selector sharedAndOwned = meta.globally_shared_part() & meta.locally_owned_part();
  const unsigned numSharedAndOwnedNodes = stk::mesh::count_entities(*meshPtr, stk::topology::NODE_RANK, sharedAndOwned);
  const bool enforceSearchResultSymmetry = true;
  ExecSpace execSpace = ExecSpace{};
  
  ResultViewType searchResults;
  for(unsigned i=0; i<NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();
    
    DomainViewType elemSpheres = create_elem_spheres(*meshPtr, radius);
    RangeViewType nodePoints = create_node_points(*meshPtr);
  
    EXPECT_EQ(elemSpheres.size(), numLocalElems);
    EXPECT_EQ(nodePoints.size(), numLocalOwnedNodes);
  
    stk::search::SearchMethod searchMethod = stk::search::MORTON_LBVH;
  
    stk::search::coarse_search(elemSpheres, nodePoints, searchMethod, comm, searchResults, execSpace, enforceSearchResultSymmetry);
  
    batchTimer.stop_batch_timer();
  }

  constexpr unsigned numNodesPerElement = 8;
  unsigned expectedNumResults = numLocalElems * numNodesPerElement;
  if (enforceSearchResultSymmetry) {
    EXPECT_GE(searchResults.size(), expectedNumResults+numSharedAndOwnedNodes);
  }
  else {
    EXPECT_EQ(searchResults.size(), expectedNumResults);
  }

  ghost_node_neighbors_to_elements(*meshPtr, searchResults);

  unpack_search_results_into_field(*meshPtr, neighborField, searchResults);

  verify_8_neighbors_per_element(*meshPtr, neighborField);

  batchTimer.print_batch_timing(NUM_RUNS);
}

}  // namespace perf_test

