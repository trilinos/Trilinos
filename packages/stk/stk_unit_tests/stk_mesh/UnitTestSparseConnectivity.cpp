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
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/SparseConnectivity.hpp>

namespace {

const stk::mesh::EntityRank nodeRank = stk::topology::NODE_RANK;
const stk::mesh::EntityRank edgeRank = stk::topology::EDGE_RANK;
const stk::mesh::EntityRank faceRank = stk::topology::FACE_RANK;
const stk::mesh::EntityRank elemRank = stk::topology::ELEM_RANK;
const stk::mesh::EntityRank constraintRank = stk::topology::CONSTRAINT_RANK;

TEST(SparseConnectivity, hasPermutation)
{
  EXPECT_FALSE(stk::mesh::SparseConnectivity::has_permutation(nodeRank, elemRank));
  EXPECT_TRUE (stk::mesh::SparseConnectivity::has_permutation(edgeRank, elemRank));
  EXPECT_TRUE (stk::mesh::SparseConnectivity::has_permutation(faceRank, elemRank));
  EXPECT_TRUE (stk::mesh::SparseConnectivity::has_permutation(edgeRank, faceRank));
  EXPECT_FALSE(stk::mesh::SparseConnectivity::has_permutation(nodeRank, edgeRank));
  EXPECT_FALSE(stk::mesh::SparseConnectivity::has_permutation(nodeRank, faceRank));
  EXPECT_FALSE(stk::mesh::SparseConnectivity::has_permutation(nodeRank, constraintRank));
  EXPECT_FALSE(stk::mesh::SparseConnectivity::has_permutation(faceRank, constraintRank));
  EXPECT_FALSE(stk::mesh::SparseConnectivity::has_permutation(elemRank, constraintRank));
}

TEST(SparseConnectivity, addAndRemoveConnectivity)
{
  stk::mesh::SparseConnectivity sparseConnectivity;

  unsigned nodeEntityOffset = 1;
  unsigned elemEntityOffset = 9;
  sparseConnectivity.update_size_of_entity_index_space(10);

  stk::mesh::Entity node(nodeEntityOffset);
  stk::mesh::Entity elem(elemEntityOffset);
  stk::mesh::ConnectivityOrdinal ordinal = 0;
  stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;

  EXPECT_TRUE(sparseConnectivity.add_connectivity(nodeRank, node, elemRank, elem, ordinal, perm));
  bool addAlreadyExistingShouldBeFalse = sparseConnectivity.add_connectivity(nodeRank, node, elemRank, elem, ordinal, perm);
  EXPECT_FALSE(addAlreadyExistingShouldBeFalse);

  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node, elemRank));

  stk::mesh::Entity elem2(elemEntityOffset+1);
  sparseConnectivity.add_connectivity(nodeRank, node, elemRank, elem2, ordinal+1, perm);

  EXPECT_EQ(2u, sparseConnectivity.num_connectivity(node, elemRank));

  sparseConnectivity.remove_connectivity(nodeRank, node, elemRank, elem, ordinal);

  stk::mesh::PairIterEntity conn = sparseConnectivity.get_connectivity(node, elemRank);
  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node, elemRank));
  EXPECT_EQ(1u, conn.size());

  const stk::mesh::Entity* elems = sparseConnectivity.begin_connectivity(node, elemRank);
  EXPECT_EQ(elem2, elems[0]);
  EXPECT_EQ(elems, conn.begin());

  sparseConnectivity.remove_connectivity(nodeRank, node, elemRank, elem2, ordinal+1);

  EXPECT_EQ(0u, sparseConnectivity.num_connectivity(node, elemRank));
  conn = sparseConnectivity.get_connectivity(node, elemRank);
  EXPECT_TRUE(conn.empty());
}

TEST(SparseConnectivity, replaceOrAddAndRemoveConnectivity)
{
  stk::mesh::SparseConnectivity sparseConnectivity;

  unsigned nodeEntityOffset = 1;
  unsigned elemEntityOffset = 9;
  sparseConnectivity.update_size_of_entity_index_space(10);

  stk::mesh::Entity node(nodeEntityOffset);
  stk::mesh::Entity elem(elemEntityOffset);
  stk::mesh::Entity elem2(elemEntityOffset+1);
  stk::mesh::ConnectivityOrdinal ordinal = 0;

  EXPECT_TRUE(sparseConnectivity.replace_or_add_connectivity(nodeRank, node, elemRank, elem2, ordinal));
  bool replaceAlreadyExistingShouldBeFalse = sparseConnectivity.replace_or_add_connectivity(nodeRank, node, elemRank, elem2, ordinal);
  EXPECT_FALSE(replaceAlreadyExistingShouldBeFalse);

  bool replaceWithDifferentElemShouldBeTrue = sparseConnectivity.replace_or_add_connectivity(nodeRank, node, elemRank, elem, ordinal);
  EXPECT_TRUE(replaceWithDifferentElemShouldBeTrue);

  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node, elemRank));

  sparseConnectivity.replace_or_add_connectivity(nodeRank, node, elemRank, elem2, ordinal+1);

  EXPECT_EQ(2u, sparseConnectivity.num_connectivity(node, elemRank));

  sparseConnectivity.remove_connectivity(nodeRank, node, elemRank, elem, ordinal);

  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node, elemRank));

  const stk::mesh::Entity* elems = sparseConnectivity.begin_connectivity(node, elemRank);
  EXPECT_EQ(elem2, elems[0]);

  sparseConnectivity.remove_connectivity(nodeRank, node, elemRank, elem2, ordinal+1);

  EXPECT_EQ(0u, sparseConnectivity.num_connectivity(node, elemRank));
}

TEST(SparseConnectivity, beginEnd)
{
  stk::mesh::SparseConnectivity sparseConnectivity;

  unsigned nodeEntityOffset = 1;
  unsigned elemEntityOffset = 9;
  sparseConnectivity.update_size_of_entity_index_space(10);
  stk::mesh::Entity node(nodeEntityOffset);
  stk::mesh::Entity elem(elemEntityOffset);
  stk::mesh::ConnectivityOrdinal ordinal = 0;
  stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;

  sparseConnectivity.add_connectivity(nodeRank, node, elemRank, elem, ordinal, perm);

  stk::mesh::Entity elem2(elemEntityOffset+1);
  sparseConnectivity.add_connectivity(nodeRank, node, elemRank, elem2, ordinal+1, perm);

  const stk::mesh::Entity* beginElems = sparseConnectivity.begin_connectivity(node, elemRank);
  const stk::mesh::Entity* endElems = sparseConnectivity.end_connectivity(node, elemRank);

  EXPECT_EQ(2u, std::distance(beginElems, endElems));

  EXPECT_EQ(elem, beginElems[0]);
  EXPECT_EQ(elem2, beginElems[1]);

  const stk::mesh::ConnectivityOrdinal* beginElemOrds = sparseConnectivity.begin_ordinals(node, elemRank);
  const stk::mesh::ConnectivityOrdinal* endElemOrds = sparseConnectivity.end_ordinals(node, elemRank);

  EXPECT_EQ(2u, std::distance(beginElemOrds, endElemOrds));

  EXPECT_EQ(ordinal, beginElemOrds[0]);
  EXPECT_EQ(ordinal+1, beginElemOrds[1]);

  EXPECT_EQ(nullptr, sparseConnectivity.begin_permutations(node, elemRank));
  EXPECT_EQ(nullptr, sparseConnectivity.end_permutations(node, elemRank));
}

TEST(SparseConnectivity, permutations)
{
  stk::mesh::SparseConnectivity sparseConnectivity;

  unsigned faceEntityOffset = 1;
  unsigned elemEntityOffset = 9;
  sparseConnectivity.update_size_of_entity_index_space(10);
  stk::mesh::Entity face(faceEntityOffset);
  stk::mesh::Entity elem(elemEntityOffset);
  stk::mesh::ConnectivityOrdinal ordinal = 0;
  stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(1);
  stk::mesh::Permutation perm2 = static_cast<stk::mesh::Permutation>(-1);

  sparseConnectivity.add_connectivity(faceRank, face, elemRank, elem, ordinal, perm);

  stk::mesh::Entity elem2(elemEntityOffset+1);
  sparseConnectivity.add_connectivity(faceRank, face, elemRank, elem2, ordinal+1, perm2);

  const stk::mesh::Permutation* beginElemPerms = sparseConnectivity.begin_permutations(face, elemRank);
  const stk::mesh::Permutation* endElemPerms = sparseConnectivity.end_permutations(face, elemRank);

  EXPECT_EQ(2u, std::distance(beginElemPerms, endElemPerms));

  EXPECT_EQ(perm, beginElemPerms[0]);
  EXPECT_EQ(perm2, beginElemPerms[1]);
}

TEST(SparseConnectivity, addConnMultiple)
{
  stk::mesh::SparseConnectivity sparseConnectivity;

  unsigned nodeEntityOffset = 1;
  unsigned node2EntityOffset = 2;
  unsigned elemEntityOffset = 3;
  unsigned faceEntityOffset = 10;
  sparseConnectivity.update_size_of_entity_index_space(4);

  stk::mesh::Entity node(nodeEntityOffset);
  stk::mesh::Entity node2(node2EntityOffset);
  stk::mesh::Entity face(faceEntityOffset);
  stk::mesh::Entity elem(elemEntityOffset);
  stk::mesh::ConnectivityOrdinal nodeElemOrd = 5;
  stk::mesh::ConnectivityOrdinal node2ElemOrd = 6;
  stk::mesh::ConnectivityOrdinal nodeFaceOrd = 2;
  stk::mesh::ConnectivityOrdinal elemFaceOrd = 3;
  stk::mesh::Permutation invalidPerm = stk::mesh::Permutation::INVALID_PERMUTATION;
  stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);

  sparseConnectivity.add_connectivity(nodeRank, node, elemRank, elem, nodeElemOrd, invalidPerm);
  const stk::mesh::Connectivity& nodeConn = sparseConnectivity.get_connectivity(node);
  (void)nodeConn;

  sparseConnectivity.add_connectivity(nodeRank, node, faceRank, face, nodeFaceOrd, invalidPerm);
  sparseConnectivity.add_connectivity(elemRank, elem, faceRank, face, elemFaceOrd, perm);

  sparseConnectivity.update_size_of_entity_index_space(11);

  sparseConnectivity.add_connectivity(nodeRank, node2, elemRank, elem, node2ElemOrd, invalidPerm);

  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node, elemRank));
  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node, faceRank));
  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(elem, faceRank));
  EXPECT_EQ(1u, sparseConnectivity.num_connectivity(node2, elemRank));

  const stk::mesh::Permutation* perms = sparseConnectivity.begin_permutations(elem, faceRank);
  const stk::mesh::Permutation* endPerms = sparseConnectivity.end_permutations(elem, faceRank);
  EXPECT_TRUE(perms != nullptr);
  EXPECT_TRUE(endPerms != nullptr);
  EXPECT_EQ(1u, std::distance(perms, endPerms));
  EXPECT_EQ(perm, perms[0]);

  EXPECT_TRUE(sparseConnectivity.remove_connectivity(nodeRank, node, elemRank, elem, nodeElemOrd));
  EXPECT_TRUE(sparseConnectivity.remove_connectivity(nodeRank, node2, elemRank, elem, node2ElemOrd));

  EXPECT_EQ(0u, sparseConnectivity.num_connectivity(node, elemRank));
  EXPECT_EQ(0u, sparseConnectivity.num_connectivity(node2, elemRank));
}

}

