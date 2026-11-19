// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <string>
#include <vector>

namespace krino
{

class NodalBoundingBox : public StkMeshTriFixture
{
public:
  NodalBoundingBox()
  {
    set_valid_proc_sizes_for_test({1, 2, 3, 4});
    build_mesh_from_spec_in_blocks(meshSpec, {1, 2, 2, 1});
  }

protected:
  stk::mesh::Part & block_part(const unsigned blockId) const
  {
    return *mMesh.mesh_meta_data().get_part("block_" + std::to_string(blockId));
  }
  FieldRef coords() const { return mMesh.mesh_meta_data().coordinate_field(); }

  // The box of the spec locations of the given nodes that are on this proc.
  BoundingBox expected_local_bbox(const std::vector<unsigned> & nodeIndices) const
  {
    BoundingBox bbox;
    for (unsigned n : nodeIndices)
      if (is_node_present_on_proc(n))
        bbox.accommodate(stk::math::Vector3d(meshSpec.nodeLocs[n][0], meshSpec.nodeLocs[n][1], 0.));
    return bbox;
  }

  BoundingBox global_bbox(const BoundingBox & localBBox) const
  {
    BoundingBox::VecType min, max;
    stk::all_reduce_min(mMesh.parallel(), localBBox.get_min().data(), min.data(), 3);
    stk::all_reduce_max(mMesh.parallel(), localBBox.get_max().data(), max.data(), 3);
    return BoundingBox(min, max);
  }

  static void expect_bbox_2d(const BoundingBox & bbox,
      const std::array<double, 2> & min,
      const std::array<double, 2> & max)
  {
    for (unsigned d = 0; d < 2; ++d)
    {
      EXPECT_EQ(min[d], bbox.get_min()[d]);
      EXPECT_EQ(max[d], bbox.get_max()[d]);
    }
  }

  static void expect_same_bbox_2d(const BoundingBox & gold, const BoundingBox & bbox)
  {
    for (unsigned d = 0; d < 2; ++d)
    {
      EXPECT_EQ(gold.get_min()[d], bbox.get_min()[d]);
      EXPECT_EQ(gold.get_max()[d], bbox.get_max()[d]);
    }
  }

  FourRightTris meshSpec;
};

TEST_F(NodalBoundingBox, selector_boundsSelectedNodesOnThisProc)
{
  if (!is_valid_proc_size_for_test()) return;
  const stk::mesh::MetaData & meta = mMesh.mesh_meta_data();

  // Nodes 1-4 are in block 2, nodes 0, 1, 3, 5 in block 1.
  const BoundingBox block2 = compute_nodal_bbox(mMesh, block_part(2), coords());
  expect_same_bbox_2d(expected_local_bbox({1, 2, 3, 4}), block2);
  expect_bbox_2d(global_bbox(block2), {{0., 0.}}, {{1., 2.}});

  const BoundingBox ownedBlock1 =
      compute_nodal_bbox(mMesh, block_part(1) & meta.locally_owned_part(), coords());
  expect_bbox_2d(global_bbox(ownedBlock1), {{0., -1.}}, {{1., 1.}});

  const BoundingBox all = compute_nodal_bbox(mMesh, meta.universal_part(), coords());
  expect_same_bbox_2d(expected_local_bbox({0, 1, 2, 3, 4, 5}), all);
  expect_bbox_2d(global_bbox(all), {{0., -1.}}, {{1., 2.}});
}

TEST_F(NodalBoundingBox, nodeList_boundsListedNodesIncludingSharedAndGhosted)
{
  if (!is_valid_proc_size_for_test()) return;

  const std::vector<unsigned> nodeIndices{2, 4, 5};
  std::vector<stk::mesh::Entity> nodes;
  for (unsigned n : nodeIndices)
    if (is_node_present_on_proc(n)) nodes.push_back(get_assigned_node_for_index(n));

  const BoundingBox bbox = compute_nodal_bbox(mMesh, coords(), nodes);
  expect_same_bbox_2d(expected_local_bbox(nodeIndices), bbox);
  expect_bbox_2d(global_bbox(bbox), {{0., -1.}}, {{1., 2.}});

  EXPECT_FALSE(compute_nodal_bbox(mMesh, coords(), std::vector<stk::mesh::Entity>{}).valid());
}

} // namespace krino
