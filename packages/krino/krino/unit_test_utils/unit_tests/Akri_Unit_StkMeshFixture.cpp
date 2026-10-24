// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <vector>

namespace krino {

class FourRightTrisFixture : public StkMeshTriFixture
{
protected:
  // Owner of each element, on every processor.
  std::vector<int> get_element_owners() const
  {
    const size_t numElements = meshSpec.allElementConn.size();
    std::vector<int> localOwners(numElements, -1);
    for (size_t i=0; i<numElements; ++i)
      if (is_owned_element(i))
        localOwners[i] = parallel_rank();
    std::vector<int> owners(numElements);
    stk::all_reduce_max(mComm, localOwners.data(), owners.data(), numElements);
    return owners;
  }

  FourRightTris meshSpec;
};

class TwoRightTrisFixture : public StkMeshTriFixture
{
protected:
  TwoRightTris meshSpec;
};

static std::vector<int> round_robin_owners(const size_t numElements, const int numProcs)
{
  std::vector<int> owners;
  for (size_t i=0; i<numElements; ++i)
    owners.push_back(i % numProcs);
  return owners;
}

TEST_F(FourRightTrisFixture, buildMeshFromSpec_defaultOwnersAreRoundRobin)
{
  build_mesh_from_spec(meshSpec);
  if (!is_valid_proc_size_for_test()) return;

  EXPECT_EQ(meshSpec.allElementConn.size(), get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  EXPECT_EQ(meshSpec.nodeLocs.size(), get_global_num_entities(mMesh, stk::topology::NODE_RANK));
  EXPECT_EQ(round_robin_owners(meshSpec.allElementConn.size(), parallel_size()), get_element_owners());
}

TEST_F(FourRightTrisFixture, isOwnedElementAndNode_eachIsOwnedOnExactlyOneProcAndLocalWhereOwned)
{
  build_mesh_from_spec(meshSpec);
  if (!is_valid_proc_size_for_test()) return;

  for (size_t i=0; i<meshSpec.allElementConn.size(); ++i)
  {
    if (is_owned_element(i))
    {
      EXPECT_TRUE(is_element_present_on_proc(i));
    }
    EXPECT_EQ(1, stk::get_global_sum(mComm, is_owned_element(i) ? 1 : 0)) << "element " << i;
  }
  for (size_t i=0; i<meshSpec.nodeLocs.size(); ++i)
  {
    if (is_owned_node(i))
    {
      EXPECT_TRUE(is_node_present_on_proc(i));
    }
    EXPECT_EQ(1, stk::get_global_sum(mComm, is_owned_node(i) ? 1 : 0)) << "node " << i;
  }
}

TEST_F(FourRightTrisFixture, buildMeshFromSpec_givenOwnersForThisProcCountAreUsed)
{
  build_mesh_from_spec(meshSpec, {{2, {0,0,1,1}}, {4, {3,2,1,0}}});
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<int> goldOwners =
      (parallel_size() == 2) ? std::vector<int>{0,0,1,1} :
      (parallel_size() == 4) ? std::vector<int>{3,2,1,0} :
      round_robin_owners(meshSpec.allElementConn.size(), parallel_size());
  EXPECT_EQ(goldOwners, get_element_owners());
}

TEST_F(FourRightTrisFixture, buildMeshFromSpec_ownersForOtherProcCountsDoNotChangeTheDefault)
{
  const int otherProcCount = parallel_size() + 1;
  build_mesh_from_spec(meshSpec, {{otherProcCount, std::vector<int>(meshSpec.allElementConn.size(), 0)}});
  if (!is_valid_proc_size_for_test()) return;

  EXPECT_EQ(round_robin_owners(meshSpec.allElementConn.size(), parallel_size()), get_element_owners());
}

TEST_F(FourRightTrisFixture, buildMeshFromSpec_invalidOwnersThrow)
{
  const int np = parallel_size();
  EXPECT_ANY_THROW(build_mesh_from_spec(meshSpec, {{np, {0,0,0}}}));     // wrong number of owners
  EXPECT_ANY_THROW(build_mesh_from_spec(meshSpec, {{np, {0,0,0,np}}}));  // owner is not a processor
  EXPECT_ANY_THROW(build_mesh_from_spec(meshSpec, {{np, {0,0,0,-1}}}));
}

TEST_F(TwoRightTrisFixture, buildMeshFromSpec_validProcSizesDefaultToUpToTheNumberOfElements)
{
  build_mesh_from_spec(meshSpec); // does not throw for more processors than elements
  EXPECT_EQ(parallel_size() <= 2, is_valid_proc_size_for_test());
  if (parallel_size() > 2)
  {
    EXPECT_EQ(0u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(TwoRightTrisFixture, buildMeshFromSpec_ownersGivenForMoreProcsThanElementsAreUsed)
{
  build_mesh_from_spec(meshSpec, {{4, {3,1}}});
  if (parallel_size() == 4)
  {
    EXPECT_EQ(2u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
  }
}

TEST_F(FourRightTrisFixture, buildMeshFromSpecInBlocks_elementsAreInTheirBlocks)
{
  build_mesh_from_spec_in_blocks(meshSpec, {1,2,2,1});
  if (!is_valid_proc_size_for_test()) return;

  EXPECT_EQ(round_robin_owners(meshSpec.allElementConn.size(), parallel_size()), get_element_owners());
  const stk::mesh::Part * block1 = mMesh.mesh_meta_data().get_part("block_1");
  const stk::mesh::Part * block2 = mMesh.mesh_meta_data().get_part("block_2");
  ASSERT_TRUE(block1 != nullptr && block2 != nullptr);
  EXPECT_EQ(2u, get_global_num_entities(mMesh, *block1));
  EXPECT_EQ(2u, get_global_num_entities(mMesh, *block2));
  for (unsigned i : {0u, 3u})
    if (mMesh.is_valid(get_assigned_element_for_index(i)))
    {
      EXPECT_TRUE(mMesh.bucket(get_assigned_element_for_index(i)).member(*block1));
    }
  for (unsigned i : {1u, 2u})
    if (mMesh.is_valid(get_assigned_element_for_index(i)))
    {
      EXPECT_TRUE(mMesh.bucket(get_assigned_element_for_index(i)).member(*block2));
    }
}

TEST_F(FourRightTrisFixture, buildMeshFromSpecInBlocks_wrongNumberOfBlockIdsThrows)
{
  EXPECT_ANY_THROW(build_mesh_from_spec_in_blocks(meshSpec, {1,2}));
}

TEST_F(FourRightTrisFixture, buildMeshFromSpec_excludedProcSizeBuildsNothing)
{
  set_valid_proc_sizes_for_test({2});
  build_mesh_from_spec(meshSpec);
  EXPECT_EQ(parallel_size() == 2, is_valid_proc_size_for_test());
  EXPECT_EQ((parallel_size() == 2) ? meshSpec.allElementConn.size() : 0u, get_global_num_entities(mMesh, stk::topology::ELEMENT_RANK));
}

TEST_F(TwoRightTrisFixture, setValidProcSizes_overridesTheDefault)
{
  set_valid_proc_sizes_for_test({1});
  build_mesh_from_spec(meshSpec);
  EXPECT_EQ(parallel_size() == 1, is_valid_proc_size_for_test());
}

TEST_F(TwoRightTrisFixture, validProcSizes_neitherSetNorBuiltFromSpec_throws)
{
  EXPECT_ANY_THROW(is_valid_proc_size_for_test());
}

} // namespace krino
