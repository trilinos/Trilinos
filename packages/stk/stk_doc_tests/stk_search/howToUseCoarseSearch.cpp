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
#include "stk_search/CoarseSearch.hpp"

namespace doc_test
{
template <typename SENDMESH, typename RECVMESH>
struct CoarseSearchTrait {
  using SendMesh = SENDMESH;
  using RecvMesh = RECVMESH;
  using SendEntityProc = typename SENDMESH::EntityProc;
  using RecvEntityProc = typename RECVMESH::EntityProc;
  using SendBoundingBox = typename SENDMESH::BoundingBox;
  using RecvBoundingBox = typename RECVMESH::BoundingBox;

  using EntityProcRelation = std::pair<RecvEntityProc, SendEntityProc>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;
};

template <typename BoundingBoxType>
struct BoundingBoxCompare {
  bool operator()(const BoundingBoxType& a, const BoundingBoxType& b) const { return a.second.id() < b.second.id(); }
};

template <typename ForwardIterator, typename Compare>
bool local_is_sorted(ForwardIterator first, ForwardIterator last, Compare comparison)
{
  if (first == last) return true;
  ForwardIterator next = first;
  while (++next != last) {
    if (comparison(*next, *first)) return false;
    ++first;
  }
  return true;
}

template <typename T, typename U>
inline void inflate_bounding_box(stk::search::Sphere<T>& s, U const& mult_fact, U const& add_fact)
{
  s.set_radius(s.radius() * mult_fact + add_fact);
}

template <typename T, typename U>
inline void inflate_bounding_box(stk::search::Box<T>& b, U const& mult_fact, U const& add_fact)
{
  const double zero = 0.0;

  STK_ThrowRequire(mult_fact >= zero);
  STK_ThrowRequire(add_fact >= zero);

  stk::search::Point<T>& min_corner = b.min_corner();
  stk::search::Point<T>& max_corner = b.max_corner();

  T diag = 0.0;
  for (int i = 0; i < 3; ++i) {
    const T dd = max_corner[i] - min_corner[i];
    diag += dd * dd;
  }
  diag = std::sqrt(diag);

  const T d = mult_fact * diag + add_fact;

  for (int i = 0; i < 3; ++i) {
    min_corner[i] -= d;
    max_corner[i] += d;
  }
}

//BEGINcoarse_search
template <typename CoarseSearchType>
void do_coarse_search(typename CoarseSearchType::SendMesh& sendMesh,
    typename CoarseSearchType::RecvMesh& recvMesh,
    const double expansionFactor,
    const double expansionSum,
    typename CoarseSearchType::EntityProcRelationVec& coarseSearchResult)
{
  using SendBoundingBox = typename CoarseSearchType::SendBoundingBox;
  using RecvBoundingBox = typename CoarseSearchType::RecvBoundingBox;

  std::vector<SendBoundingBox> domain_vector;
  std::vector<RecvBoundingBox> range_vector;

  sendMesh.bounding_boxes(domain_vector);
  recvMesh.bounding_boxes(range_vector);

  if (!local_is_sorted(domain_vector.begin(), domain_vector.end(), BoundingBoxCompare<SendBoundingBox>()))
    std::sort(domain_vector.begin(), domain_vector.end(), BoundingBoxCompare<SendBoundingBox>());

  if (!local_is_sorted(range_vector.begin(), range_vector.end(), BoundingBoxCompare<RecvBoundingBox>()))
    std::sort(range_vector.begin(), range_vector.end(), BoundingBoxCompare<RecvBoundingBox>());

  for (SendBoundingBox& i : domain_vector) {
    inflate_bounding_box(i.first, expansionFactor, expansionSum);
  }

  stk::search::coarse_search(range_vector, domain_vector, stk::search::KDTREE, sendMesh.comm(), coarseSearchResult);

  std::sort(coarseSearchResult.begin(), coarseSearchResult.end());
}
//ENDcoarse_search

TEST(StkSearchHowTo, useCoarseSearch)
{
//BEGINuse_coarse_search
  using CoarseSearchType = CoarseSearchTrait<Hex8SourceMesh, SinglePointMesh>;
  using Relation = typename CoarseSearchType::EntityProcRelation;
  using RelationVec = typename CoarseSearchType::EntityProcRelationVec;

  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) {
    GTEST_SKIP();
  }

  // Build 8 element cube
  const std::string meshSpec("generated:2x2x2");
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::io::fill_mesh(meshSpec, *mesh);

  // Point in element 1
  double x = 0.5, y = 1, z = 1;
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

  RelationVec coarseSearchResult;

  // Get single recv point
  SinglePointMesh::EntityKey expectedRecvKey(1);
  SinglePointMesh::EntityProc rangeEntry(expectedRecvKey, 0);

  double expansionFactor = 0.01;
  double expansionSum = 0.005;

  do_coarse_search<CoarseSearchType>(*sendMesh, *recvMesh, expansionFactor, expansionSum, coarseSearchResult);

  EXPECT_EQ(4u, coarseSearchResult.size());
//ENDuse_coarse_search

  stk::mesh::EntityIdVector expectedSendIds{1, 3, 5, 7};

  for (Relation& result : coarseSearchResult) {
    const Hex8SourceMesh::EntityKey sendEntityKey = result.second.id();
    EXPECT_EQ(stk::topology::ELEM_RANK, sendEntityKey.rank());

    auto found = std::binary_search(expectedSendIds.begin(), expectedSendIds.end(), sendEntityKey.id());
    EXPECT_TRUE(found);

    const int recvEntityKey = result.first.id();
    EXPECT_EQ(1, recvEntityKey);
  }
}

}  // namespace doc_test
