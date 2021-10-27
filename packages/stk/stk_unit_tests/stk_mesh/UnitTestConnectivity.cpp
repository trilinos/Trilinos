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
#include <stk_mesh/base/Connectivity.hpp>

namespace {

const stk::mesh::EntityRank nodeRank = stk::topology::NODE_RANK;
const stk::mesh::EntityRank faceRank = stk::topology::FACE_RANK;
const stk::mesh::EntityRank elemRank = stk::topology::ELEM_RANK;
const stk::mesh::EntityRank constraintRank = stk::topology::CONSTRAINT_RANK;

class ConnectivityTester : public stk::mesh::Connectivity {
public:
  void expect_null_ptrs()
  {
    EXPECT_EQ(nullptr, m_connectivity);
    EXPECT_EQ(nullptr, m_offsetsAndOrdinals);
  }

  void expect_zeros_and_empty()
  {
    const unsigned numOffsets = stk::mesh::Connectivity::MAX_NUM_RANKS+1;
    for(unsigned r=0; r<numOffsets; ++r) {
      EXPECT_EQ(0u, m_offsetsAndOrdinals[r]);
    }
  }

  void expect_num_per_rank(const std::vector<unsigned>& numPerRank)
  {
    EXPECT_TRUE(stk::mesh::Connectivity::MAX_NUM_RANKS == numPerRank.size());
    for(unsigned r=0; r<numPerRank.size(); ++r) {
      EXPECT_EQ(numPerRank[r], num_connectivity(static_cast<stk::mesh::EntityRank>(r)));
    }
  }

  void expect_connectivity(const std::vector<stk::mesh::Entity>& conn)
  {
    EXPECT_EQ(conn.size(), m_offsetsAndOrdinals[stk::mesh::Connectivity::MAX_NUM_RANKS]);
    for(unsigned i=0; i<conn.size(); ++i) {
      EXPECT_EQ(m_connectivity[i], conn[i]);
    }
  }

  void expect_ordinals(const std::vector<stk::mesh::ConnectivityOrdinal>& ords)
  {
    EXPECT_EQ(ords.size(), m_offsetsAndOrdinals[stk::mesh::Connectivity::MAX_NUM_RANKS]);
    for(unsigned i=0; i<ords.size(); ++i) {
      EXPECT_EQ(m_offsetsAndOrdinals[stk::mesh::Connectivity::m_firstOrdinalIndex + i], ords[i]);
    }
  }
};

class TestConnectivity : public ::testing::Test
{
public:
  TestConnectivity()
  : node1(1), node2(2), node3(3),
    face1(11), face2(22), face3(33),
    elem1(68), elem2(67), elem3(66),
    constraint1(5555555), highest1(9999),
    ordinal1(0), ordinal2(1), ordinal3(2)
  {
  }

  void expect_num_per_rank(stk::mesh::Connectivity& connObject, const std::vector<unsigned>& numPerRank)
  {
    ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connObject);
    connTester.expect_num_per_rank(numPerRank);
  }

  void expect_num_per_rank(const std::vector<unsigned>& numPerRank)
  {
    expect_num_per_rank(connectivity, numPerRank);
  }

  void expect_connectivity(stk::mesh::Connectivity& connObject, const std::vector<stk::mesh::Entity>& conn)
  {
    ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connObject);
    connTester.expect_connectivity(conn);
  }

  void expect_connectivity(const std::vector<stk::mesh::Entity>& conn)
  {
    expect_connectivity(connectivity, conn);
  }

  void expect_ordinals(const std::vector<stk::mesh::ConnectivityOrdinal>& ords)
  {
    ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connectivity);
    connTester.expect_ordinals(ords);
  }

  void add_2nodes_2elems()
  {
    const bool allowDuplicateOrdinals = true;

    EXPECT_EQ(0, connectivity.add_connectivity(nodeRank, node1, ordinal1, allowDuplicateOrdinals));
    EXPECT_EQ(1, connectivity.add_connectivity(nodeRank, node2, ordinal2, allowDuplicateOrdinals));
    EXPECT_EQ(0, connectivity.add_connectivity(elemRank, elem1, ordinal1, allowDuplicateOrdinals));
    EXPECT_EQ(1, connectivity.add_connectivity(elemRank, elem2, ordinal2, allowDuplicateOrdinals));
  }

protected:
  stk::mesh::Connectivity connectivity;
  stk::mesh::Entity node1, node2, node3;
  stk::mesh::Entity face1, face2, face3;
  stk::mesh::Entity elem1, elem2, elem3;
  stk::mesh::Entity constraint1;
  stk::mesh::Entity highest1;
  stk::mesh::ConnectivityOrdinal ordinal1, ordinal2, ordinal3;
};

TEST_F(TestConnectivity, defaultConstruct_zerosAndEmpty)
{
  ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connectivity);
  connTester.expect_zeros_and_empty();
}

TEST_F(TestConnectivity, addConnectivity)
{
  add_2nodes_2elems();
  expect_num_per_rank({2, 0, 0, 2, 0, 0});
  expect_connectivity({node1, node2, elem1, elem2});
  expect_ordinals({ordinal1, ordinal2, ordinal1, ordinal2});
}
  
TEST_F(TestConnectivity, addConnectivity_then_clear)
{
  add_2nodes_2elems();
  expect_num_per_rank({2, 0, 0, 2, 0, 0});
  connectivity.clear();
  expect_num_per_rank({0, 0, 0, 0, 0, 0});
}
  
TEST_F(TestConnectivity, assignment)
{
  add_2nodes_2elems();
  expect_num_per_rank({2, 0, 0, 2, 0, 0});

  stk::mesh::Connectivity connCopy = connectivity;
  expect_num_per_rank(connCopy, {2, 0, 0, 2, 0, 0});  
  expect_num_per_rank(connectivity, {2, 0, 0, 2, 0, 0});  

  connCopy = std::move(connectivity);
  expect_num_per_rank(connCopy, {2, 0, 0, 2, 0, 0});  
  ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connectivity);
  connTester.expect_null_ptrs();
}

TEST_F(TestConnectivity, replaceOrAddConnectivity_not_pre_existing)
{
  EXPECT_EQ(0, connectivity.replace_or_add_connectivity(nodeRank, node1, ordinal1));
  EXPECT_EQ(1, connectivity.replace_or_add_connectivity(nodeRank, node2, ordinal2));
  EXPECT_EQ(0, connectivity.replace_or_add_connectivity(elemRank, elem1, ordinal1));
  EXPECT_EQ(1, connectivity.replace_or_add_connectivity(elemRank, elem2, ordinal2));
  expect_num_per_rank({2, 0, 0, 2, 0, 0});
  expect_connectivity({node1, node2, elem1, elem2});
  expect_ordinals({ordinal1, ordinal2, ordinal1, ordinal2});
}
  
TEST_F(TestConnectivity, replaceOrAddConnectivity_pre_existing)
{
  const bool allowDuplicateOrdinals = true;
  EXPECT_EQ(0, connectivity.add_connectivity(nodeRank, node1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(1, connectivity.add_connectivity(nodeRank, node2, ordinal2, allowDuplicateOrdinals));
  EXPECT_EQ(0, connectivity.add_connectivity(elemRank, elem1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(1, connectivity.add_connectivity(elemRank, elem2, ordinal2, allowDuplicateOrdinals));
  EXPECT_EQ(-1, connectivity.replace_or_add_connectivity(nodeRank, node1, ordinal1));
  EXPECT_EQ(-1, connectivity.replace_or_add_connectivity(nodeRank, node2, ordinal2));
  EXPECT_EQ(-1, connectivity.replace_or_add_connectivity(elemRank, elem1, ordinal1));
  EXPECT_EQ(-1, connectivity.replace_or_add_connectivity(elemRank, elem2, ordinal2));
  expect_num_per_rank({2, 0, 0, 2, 0, 0});
  expect_connectivity({node1, node2, elem1, elem2});
  expect_ordinals({ordinal1, ordinal2, ordinal1, ordinal2});
}
  
TEST_F(TestConnectivity, addAndRemoveConnectivity)
{
  const bool allowDuplicateOrdinals = true;

  EXPECT_EQ(0, connectivity.add_connectivity(nodeRank, node1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(1, connectivity.add_connectivity(nodeRank, node2, ordinal2, allowDuplicateOrdinals));
  EXPECT_EQ(0, connectivity.add_connectivity(elemRank, elem1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(1, connectivity.add_connectivity(elemRank, elem2, ordinal2, allowDuplicateOrdinals));
  expect_num_per_rank({2, 0, 0, 2, 0, 0});
  expect_connectivity({node1, node2, elem1, elem2});
  expect_ordinals({ordinal1, ordinal2, ordinal1, ordinal2});
  
  EXPECT_EQ(0, connectivity.remove_connectivity(elemRank, elem1, ordinal1));
  expect_num_per_rank({2, 0, 0, 1, 0, 0});
  expect_connectivity({node1, node2, elem2});
  expect_ordinals({ordinal1, ordinal2, ordinal2});

  stk::mesh::PairIterEntity conn = connectivity.get_connectivity(elemRank);
  const stk::mesh::Entity* elems = connectivity.begin_connectivity(elemRank);
  const stk::mesh::ConnectivityOrdinal* ords = connectivity.begin_ordinals(elemRank);
  EXPECT_EQ(elem2, elems[0]);
  EXPECT_EQ(1u, conn.size());
  EXPECT_EQ(conn.begin(), elems);
  EXPECT_EQ(conn.end(), elems+1);
  EXPECT_EQ(ordinal2, ords[0]);

  EXPECT_EQ(0, connectivity.remove_connectivity(elemRank, elem2, ordinal2));
  expect_num_per_rank({2, 0, 0, 0, 0, 0});
  expect_connectivity({node1, node2});
  expect_ordinals({ordinal1, ordinal2});
  conn = connectivity.get_connectivity(elemRank);
  EXPECT_EQ(connectivity.begin_connectivity(elemRank), connectivity.end_connectivity(elemRank));
  EXPECT_TRUE(conn.empty());
  EXPECT_EQ(conn.begin(), conn.end());
  EXPECT_EQ(connectivity.begin_ordinals(elemRank), connectivity.end_ordinals(elemRank));
}

TEST_F(TestConnectivity, addAndRemoveMultipleRanks)
{
  const bool allowDuplicateOrdinals = true;

  EXPECT_EQ(0, connectivity.add_connectivity(nodeRank, node1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(0, connectivity.add_connectivity(faceRank, face1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(0, connectivity.add_connectivity(elemRank, elem1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(1, connectivity.add_connectivity(elemRank, elem2, ordinal2, allowDuplicateOrdinals));
  EXPECT_EQ(0, connectivity.add_connectivity(constraintRank, constraint1, ordinal1, allowDuplicateOrdinals));
  expect_num_per_rank({1, 0, 1, 2, 1, 0});
  
  EXPECT_EQ(0, connectivity.remove_connectivity(elemRank, elem1, ordinal1));
  expect_num_per_rank({1, 0, 1, 1, 1, 0});

  const stk::mesh::Entity* elems = connectivity.begin_connectivity(elemRank);
  const stk::mesh::ConnectivityOrdinal* ords = connectivity.begin_ordinals(elemRank);
  EXPECT_EQ(elem2, elems[0]);
  EXPECT_EQ(ordinal2, ords[0]);

  EXPECT_EQ(0, connectivity.remove_connectivity(elemRank, elem2, ordinal2));
  expect_num_per_rank({1, 0, 1, 0, 1, 0});

  EXPECT_EQ(connectivity.begin_connectivity(elemRank), connectivity.end_connectivity(elemRank));
  EXPECT_EQ(connectivity.begin_ordinals(elemRank), connectivity.end_ordinals(elemRank));
  EXPECT_NE(connectivity.begin_connectivity(faceRank), connectivity.end_connectivity(faceRank));
  EXPECT_NE(connectivity.begin_ordinals(faceRank), connectivity.end_ordinals(faceRank));
}

TEST_F(TestConnectivity, addRemoveAddAddRemoveAdd)
{
  const bool allowDuplicateOrdinals = true;

  EXPECT_EQ(0, connectivity.add_connectivity(elemRank, elem1, ordinal1, allowDuplicateOrdinals));
  expect_num_per_rank({0, 0, 0, 1, 0, 0});

  EXPECT_EQ(0, connectivity.remove_connectivity(elemRank, elem1, ordinal1));
  expect_num_per_rank({0, 0, 0, 0, 0, 0});

  EXPECT_EQ(0, connectivity.add_connectivity(elemRank, elem1, ordinal1, allowDuplicateOrdinals));
  expect_num_per_rank({0, 0, 0, 1, 0, 0});

  EXPECT_EQ(1, connectivity.add_connectivity(elemRank, elem2, ordinal2, allowDuplicateOrdinals));
  expect_num_per_rank({0, 0, 0, 2, 0, 0});

  EXPECT_EQ(1, connectivity.remove_connectivity(elemRank, elem2, ordinal2));
  expect_num_per_rank({0, 0, 0, 1, 0, 0});

  EXPECT_EQ(1, connectivity.add_connectivity(elemRank, elem2, ordinal2, allowDuplicateOrdinals));
  expect_num_per_rank({0, 0, 0, 2, 0, 0});

  const stk::mesh::Entity* elems = connectivity.begin_connectivity(elemRank);
  const stk::mesh::ConnectivityOrdinal* ords = connectivity.begin_ordinals(elemRank);
  EXPECT_EQ(elem1, elems[0]);
  EXPECT_EQ(elem2, elems[1]);
  EXPECT_EQ(ordinal1, ords[0]);
  EXPECT_EQ(ordinal2, ords[1]);
}

TEST_F(TestConnectivity, addAndRemoveMultipleConnectivity)
{
  const bool allowDuplicateOrdinals = true;

  EXPECT_EQ(0, connectivity.add_connectivity(faceRank, face1, ordinal1, allowDuplicateOrdinals));
  EXPECT_EQ(1, connectivity.add_connectivity(faceRank, face2, ordinal2, allowDuplicateOrdinals));
  EXPECT_EQ(2, connectivity.add_connectivity(faceRank, face3, ordinal3, allowDuplicateOrdinals));
  expect_num_per_rank({0, 0, 3, 0, 0, 0});

  {
    const stk::mesh::Entity* faces = connectivity.begin_connectivity(faceRank);
    const stk::mesh::ConnectivityOrdinal* ordinals = connectivity.begin_ordinals(faceRank);

    EXPECT_EQ(face1, faces[0]);
    EXPECT_EQ(face2, faces[1]);
    EXPECT_EQ(face3, faces[2]);
    EXPECT_EQ(ordinal1, ordinals[0]);
    EXPECT_EQ(ordinal2, ordinals[1]);
    EXPECT_EQ(ordinal3, ordinals[2]);
  }

  EXPECT_EQ(1, connectivity.remove_connectivity(faceRank, face2, ordinal2));
  expect_num_per_rank({0, 0, 2, 0, 0, 0});

  int indexOfRemovedEntity = connectivity.remove_connectivity(faceRank, face2, ordinal2);
  bool duplicateRemoveCausedNoChange = (indexOfRemovedEntity == -1);
  EXPECT_TRUE(duplicateRemoveCausedNoChange);

  expect_num_per_rank({0, 0, 2, 0, 0, 0});
  expect_connectivity({face1, face3});
  expect_ordinals({ordinal1, ordinal3});

  EXPECT_EQ(1, connectivity.remove_connectivity(faceRank, face3, ordinal3));
  expect_num_per_rank({0, 0, 1, 0, 0, 0});
  expect_connectivity({face1});
  expect_ordinals({ordinal1});

  EXPECT_EQ(0, connectivity.remove_connectivity(faceRank, face1, ordinal1));
  expect_num_per_rank({0, 0, 0, 0, 0, 0});

  EXPECT_EQ(connectivity.begin_connectivity(faceRank), connectivity.end_connectivity(faceRank));
  EXPECT_EQ(connectivity.begin_ordinals(faceRank), connectivity.end_ordinals(faceRank));
}

TEST_F(TestConnectivity, moveInVector)
{
  const bool allowDuplicateOrdinals = true;
  std::vector<stk::mesh::Connectivity> connectivities(1);
  connectivities[0].add_connectivity(nodeRank, node1, ordinal1, allowDuplicateOrdinals);
  {
    ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connectivities[0]);
    connTester.expect_num_per_rank({1, 0, 0, 0, 0, 0});
    connTester.expect_connectivity({node1});
    connTester.expect_ordinals({ordinal1});
  }

  connectivities.resize(2);
  {
    ConnectivityTester& connTester = static_cast<ConnectivityTester&>(connectivities[0]);
    connTester.expect_num_per_rank({1, 0, 0, 0, 0, 0});
    connTester.expect_connectivity({node1});
    connTester.expect_ordinals({ordinal1});
  }
}

class MyClass {
public:
  stk::mesh::Entity* entityPtr;
  stk::mesh::ConnectivityOrdinal* ordPtr;
};

TEST_F(TestConnectivity, simpleStructSize)
{
  EXPECT_EQ(16u, sizeof(MyClass));
  EXPECT_EQ(sizeof(MyClass), sizeof(stk::mesh::Connectivity));
}

}

