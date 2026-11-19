// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_SharpFeature.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

template <typename MESHSPEC> class SharpFeatureFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
protected:
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;

  void find_sharp_features()
  {
    mySharpFeatureInfo.find_sharp_features(mMesh,
        mMesh.mesh_meta_data().coordinate_field(),
        mMesh.mesh_meta_data().universal_part(),
        AuxMetaData::get(mMesh.mesh_meta_data()).exposed_boundary_part(),
        myCosFeatureAngle);
  }

  void build_mesh_and_find_sharp_features()
  {
    this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    find_sharp_features();
  }

  void test_is_node_pinned(const stk::mesh::Entity node, const bool goldIsNodePinned)
  {
    const SharpFeatureConstraint * constraint = mySharpFeatureInfo.get_constraint(node);
    if (goldIsNodePinned)
    {
      EXPECT_TRUE(constraint != nullptr && constraint->is_pinned());
    }
    else
    {
      EXPECT_TRUE(constraint == nullptr || !constraint->is_pinned());
    }
  }

  bool is_node_in_assigned_nodes_for_indices(const stk::mesh::BulkData & /*mesh*/,
      const stk::mesh::Entity node,
      const std::vector<unsigned> & nodeIndices)
  {
    for (auto nodeIndex : nodeIndices)
      if (this->get_assigned_node_for_index(nodeIndex) == node) return true;
    return false;
  }

  void test_are_nodes_pinned(const std::vector<unsigned> & goldPinnedNodeIndices)
  {
    std::vector<stk::mesh::Entity> ownedNodes;
    stk::mesh::get_selected_entities(mMesh.mesh_meta_data().locally_owned_part(),
        mMesh.buckets(stk::topology::NODE_RANK),
        ownedNodes);

    for (auto && node : ownedNodes)
    {
      const bool goldIsNodePinned =
          is_node_in_assigned_nodes_for_indices(mMesh, node, goldPinnedNodeIndices);
      test_is_node_pinned(node, goldIsNodePinned);
    }
  }

  void test_is_node_constrained_on_edge(const unsigned sharpEdgeNodeIndex,
      const std::array<unsigned, 2> & goldSharpEdgeNodeNbrIndices)
  {
    stk::mesh::Entity sharpEdgeNode = this->get_assigned_node_for_index(sharpEdgeNodeIndex);
    if (mMesh.is_valid(sharpEdgeNode) &&
        mMesh.parallel_owner_rank(sharpEdgeNode) == mMesh.parallel_rank())
    {
      const SharpFeatureConstraint * constraint = mySharpFeatureInfo.get_constraint(sharpEdgeNode);
      ASSERT_TRUE(constraint != nullptr && constraint->is_constrained_on_edge());
      const std::array<stk::mesh::Entity, 2> sharpEdgeNodes = constraint->get_sharp_edge_nodes();
      for (auto goldSharpEdgeNodeNbrIndex : goldSharpEdgeNodeNbrIndices)
      {
        stk::mesh::Entity goldSharpEdgeNodeNbr =
            this->get_assigned_node_for_index(goldSharpEdgeNodeNbrIndex);
        EXPECT_TRUE(
            sharpEdgeNodes[0] == goldSharpEdgeNodeNbr || sharpEdgeNodes[1] == goldSharpEdgeNodeNbr);
      }
    }
  }

  void test_are_all_nodes_pinned()
  {
    std::vector<stk::mesh::Entity> ownedNodes;
    stk::mesh::get_selected_entities(mMesh.mesh_meta_data().locally_owned_part(),
        mMesh.buckets(stk::topology::NODE_RANK),
        ownedNodes);

    for (auto && node : ownedNodes)
      test_is_node_pinned(node, true);
  }

  MESHSPEC meshSpec;
  double myCosFeatureAngle{std::cos(M_PI / 180. * 135.0)};
  SharpFeatureInfo mySharpFeatureInfo;
};

typedef SharpFeatureFixture<RegularTet> SharpFeatureRegularTetFixture;

TEST_F(SharpFeatureRegularTetFixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<RightTet> SharpFeatureRightTetFixture;

TEST_F(SharpFeatureRightTetFixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<RegularTri> SharpFeatureRegularTriFixture;

TEST_F(SharpFeatureRegularTriFixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<Tri306090> SharpFeatureTri306090Fixture;

TEST_F(SharpFeatureTri306090Fixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<TwoTri306090> SharpFeatureTwoTri306090Fixture;

TEST_F(
    SharpFeatureTwoTri306090Fixture, meshWithCornerNodesAndUnconstrainedNode_constraintsAreCorrect)
{
  if (stk::parallel_machine_size(mComm) <= 2)
  {
    this->build_mesh_from_spec(meshSpec);

    find_sharp_features();

    const std::vector<unsigned> goldPinnedNodeIndices{1, 2, 3};
    test_are_nodes_pinned(goldPinnedNodeIndices);
  }
}

typedef SharpFeatureFixture<FourRightTets> SharpFeatureFourRightTetsFixture;

TEST_F(
    SharpFeatureFourRightTetsFixture, meshWithCornerNodesAndUnconstrainedNode_constraintsAreCorrect)
{
  if (stk::parallel_machine_size(mComm) <= 4)
  {
    this->build_mesh_from_spec(meshSpec, {{2, {0, 0, 1, 1}}, {3, {0, 1, 2, 2}}});

    find_sharp_features();

    const std::vector<unsigned> goldPinnedNodeIndices{1, 2, 3, 4, 5};
    test_are_nodes_pinned(goldPinnedNodeIndices);
  }
}

typedef SharpFeatureFixture<TwoRightTets> SharpFeatureTwoRightTetsFixture;

TEST_F(SharpFeatureTwoRightTetsFixture, meshWithCornerNodesAndEdgeNode_constraintsAreCorrect)
{
  if (stk::parallel_machine_size(mComm) <= 2)
  {
    this->build_mesh_from_spec(meshSpec);

    find_sharp_features();

    const std::vector<unsigned> goldPinnedNodeIndices{1, 2, 3, 4};
    test_are_nodes_pinned(goldPinnedNodeIndices);

    test_is_node_constrained_on_edge(0, {{1, 3}});
  }
}

} // namespace krino
