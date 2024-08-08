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
#include "searchMockMesh.hpp"

namespace doc_test {

//BEGINfilter_coarse_search
TEST(StkSearchHowTo, useFilterCoarseSearch)
{
  using Relation = std::pair<SinglePointMesh::EntityProc, Hex8SourceMesh::EntityProc>;
  using RelationVec = std::vector<Relation>;

  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  // Build 8 element cube
  const std::string meshSpec("generated:2x2x2");
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::io::fill_mesh(meshSpec, *mesh);

  // Point in element 1
  double x = 0.5, y = 0.5, z = 0.5;
  double geometricTolerance = 0.1;
  double parametricTolerance = 0.001;
  stk::mesh::EntityKey expectedSendKey(stk::topology::ELEM_RANK, 1u);

  // Create recv mesh
  auto recvMesh = std::make_shared<SinglePointMesh>(communicator, x, y, z, parametricTolerance, geometricTolerance);

  // Create send mesh
  stk::mesh::Part* part = meta.get_part("block_1");
  STK_ThrowRequireMsg(nullptr != part, "Error: block_1 does not exist");
  stk::mesh::PartVector parts{part};
  auto sendMesh = std::make_shared<Hex8SourceMesh>(*mesh, parts, mesh->parallel(), parametricTolerance);

  RelationVec relationVec;

  // Get single recv point
  SinglePointMesh::EntityKey expectedRecvKey(1);
  SinglePointMesh::EntityProc rangeEntry(expectedRecvKey, 0);

  // Load all elements as coarse search candidates
  stk::mesh::BucketVector const& buckets = mesh->get_buckets(stk::topology::ELEM_RANK, meta.universal_part());
  for(auto&& ib : buckets) {
    stk::mesh::Bucket& b = *ib;

    for(auto elem : b) {
      stk::mesh::EntityKey domainKey = mesh->entity_key(elem);
      Hex8SourceMesh::EntityProc domainEntry(domainKey, 0);

      relationVec.emplace_back(rangeEntry, domainEntry);
    }
  }

  EXPECT_EQ(8u, relationVec.size());

  bool useNearestNodeForClosestBoundingBox{false};
  bool useCentroidForGeometricProximity{false};
  bool verbose{false};
  auto extrapolateOption = stk::search::ObjectOutsideDomainPolicy::ABORT;

  stk::search::FilterCoarseSearchOptions options(std::cout, extrapolateOption,
                                                 useNearestNodeForClosestBoundingBox,
                                                 useCentroidForGeometricProximity, verbose);
  stk::search::FilterCoarseSearchResultVector<SinglePointMesh> searchResults;
  stk::search::filter_coarse_search("filter", relationVec, *sendMesh, *recvMesh, options, searchResults);

  EXPECT_EQ(1u, relationVec.size());

  auto relation = relationVec[0];
  const SinglePointMesh::EntityKey recvEntityKey = relation.first.id();
  const Hex8SourceMesh::EntityKey sendEntityKey = relation.second.id();

  EXPECT_EQ(expectedRecvKey, recvEntityKey);
  EXPECT_EQ(expectedSendKey, sendEntityKey);
}
//ENDfilter_coarse_search

}
