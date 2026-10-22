#include <Akri_AllReduce.hpp>
#include <Akri_MeshHelpers.hpp>
#include <gtest/gtest.h>

#include <Akri_MeshSpecs.hpp>
#include <Akri_SideAttachedElements.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino {

class SideAttachedTriElements : public StkMeshTriFixture
{
public:
  SideAttachedTriElements() {}

  void build_quad_split_4tri_with_sideset(const std::vector<unsigned> &elementBlockIDs = {1,1,1,1})
  {
    mBuilder.create_sideset_part(1);
    QuadSplit4Tri meshSpec;
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, mBuilder.get_processor_distribution_for_num_elements(meshSpec.allElementConn.size()));
  }

  void build_quad_split_4tri_and_quad_split_2tri(const std::vector<unsigned> &elementBlockIDs = {1,1,1,1,1,1})
  {
    QuadSplit4TriAndQuadSplit2Tri meshSpec;
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, mBuilder.get_processor_distribution_for_num_elements(meshSpec.allElementConn.size()));
  }

  stk::mesh::Entity get_tri_side_with_node_indices(const unsigned node0, const unsigned node1) const { return mBuilder.get_side_with_nodes({get_assigned_node_for_index(node0), get_assigned_node_for_index(node1)}); }
  void add_side_with_node_indices_to_sideset(const unsigned node0, const unsigned node1, const unsigned sidesetId)
  {
    stk::mesh::Entity side;
    if (mMesh.is_valid(get_assigned_node_for_index(node0)) && mMesh.is_valid(get_assigned_node_for_index(node1)))
      side = get_tri_side_with_node_indices(node0, node1);
    if (mMesh.is_valid(side) && mMesh.bucket(side).owned())
      mBuilder.add_sides_to_sidesets({side}, {{sidesetId}});
    else
      mBuilder.add_sides_to_sidesets({}, {});
  }

  void test_num_unattached(const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const size_t goldNumUnattached)
  {
    const std::vector<stk::mesh::Entity> unattachedElems = get_selected_owned_side_unattached_elements(get_mesh(), elementSelector, sideSelector);
    size_t numUnattached = unattachedElems.size();
    all_reduce_sum(get_mesh().parallel(), numUnattached);
    EXPECT_EQ(goldNumUnattached, numUnattached);
  }

  void test_num_elements_not_in_largest_group_of_selected_side_attached_elements(const stk::mesh::Selector & elementSelector, const size_t goldNumNotInLargestGroup)
  {
    const std::vector<stk::mesh::Entity> elemsNotInLargestGroup = find_owned_elements_that_are_not_in_the_largest_group_of_selected_side_attached_elements(get_mesh(), elementSelector);
    size_t numNotInLargestGroup = elemsNotInLargestGroup.size();
    all_reduce_sum(get_mesh().parallel(), numNotInLargestGroup);
    EXPECT_EQ(goldNumNotInLargestGroup, numNotInLargestGroup);
  }



protected:


};

TEST_F(SideAttachedTriElements, emptySideSelector_allElementsAreUnAttached)
{
  build_quad_split_4tri_with_sideset();

  stk::mesh::Selector elemSelector = *mBuilder.get_block_parts()[0];
  stk::mesh::Selector emptySideSelector;

  test_num_unattached(elemSelector, emptySideSelector, 4);
}

TEST_F(SideAttachedTriElements, sideSelectorWithSideAttachedToAllElements_noElementsAreUnAttached)
{
  build_quad_split_4tri_with_sideset();
  const unsigned sidesetId = 1;
  add_side_with_node_indices_to_sideset(0,1, sidesetId);

  stk::mesh::Selector elemSelector = *mBuilder.get_block_parts()[0];
  stk::mesh::Selector sideSelector = mBuilder.get_sideset_part(sidesetId);
  test_num_unattached(elemSelector, sideSelector, 0);
}

TEST_F(SideAttachedTriElements, sideSelectorWithSideAttachedToOneOfTwoElementsInBlockThatAreEdgeConnected_oneElementsIsUnAttached)
{
  build_quad_split_4tri_with_sideset({1,2,1,2});
  const unsigned sidesetId = 1;
  add_side_with_node_indices_to_sideset(0,1, sidesetId);

  stk::mesh::Selector elemSelector = *mBuilder.get_block_parts()[0];
  stk::mesh::Selector sideSelector = mBuilder.get_sideset_part(sidesetId);

  test_num_unattached(elemSelector, sideSelector, 1);
}

TEST_F(SideAttachedTriElements, meshWithClusterOf2AndClusterOf4tris_numNotInLargesGroup_equal2)
{
  build_quad_split_4tri_and_quad_split_2tri();

  stk::mesh::Selector elemSelector = *mBuilder.get_block_parts()[0];

  test_num_elements_not_in_largest_group_of_selected_side_attached_elements(elemSelector, 2);
}


}
