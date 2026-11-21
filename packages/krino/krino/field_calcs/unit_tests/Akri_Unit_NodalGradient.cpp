// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_NodalGradient.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <functional>
#include <vector>

namespace krino
{

// FourRightTris with node 4 moved to (0,3) and node 5 to (0,-2), so elements 2 and 3 have area 1
// and elements 0 and 1 area 1/2: a volume-weighted average differs from a plain one.
class NodalGradientFixture : public StkMeshTriFixture
{
protected:
  explicit NodalGradientFixture(const stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    set_valid_proc_sizes_for_test({1, 2, 3, 4});
    mMesh.set_automatic_aura_option(auraOption);
    mScalar = mBuilder.get_aux_meta().register_field("scalar",
        FieldType::REAL,
        stk::topology::NODE_RANK,
        1,
        1,
        mMesh.mesh_meta_data().universal_part());
    register_nodal_gradient_for_scalar_field(mMesh.mesh_meta_data(), mScalar);
    build_mesh_from_spec_in_blocks(meshSpec, {1, 2, 2, 1});

    for (const auto & loc : meshSpec.nodeLocs)
      mNodeLocs.emplace_back(loc[0], loc[1], 0.);
    mNodeLocs[4] = stk::math::Vector3d(0., 3., 0.);
    mNodeLocs[5] = stk::math::Vector3d(0., -2., 0.);
    if (is_valid_proc_size_for_test())
      for (unsigned n : {4u, 5u})
        if (has_node(n)) set_node_coordinates(get_assigned_node_for_index(n), mNodeLocs[n]);
  }

  void set_scalar(const std::function<double(const stk::math::Vector3d &)> & fn)
  {
    auto scalarData = mScalar.data<double, stk::mesh::ReadWrite>();
    for (auto * bucketPtr :
        mMesh.get_buckets(stk::topology::NODE_RANK, mMesh.mesh_meta_data().universal_part()))
      for (auto node : *bucketPtr)
        scalarData.entity_values(node)() = fn(get_node_coordinates(node));
  }

  void update()
  {
    update_nodal_gradient(mMesh, mMesh.mesh_meta_data().coordinate_field(), mScalar);
  }

  bool has_node(const unsigned nodeIndex) const { return is_node_present_on_proc(nodeIndex); }
  stk::math::Vector3d gradient_at(const unsigned nodeIndex) const
  {
    const FieldRef grad = get_nodal_gradient_for_scalar_field(mMesh.mesh_meta_data(), mScalar);
    return get_vector_field(mMesh, grad, get_assigned_node_for_index(nodeIndex), 2);
  }

  double element_area(const unsigned elemIndex) const
  {
    const auto & conn = meshSpec.allElementConn[elemIndex];
    const stk::math::Vector3d e1 = mNodeLocs[conn[1]] - mNodeLocs[conn[0]];
    const stk::math::Vector3d e2 = mNodeLocs[conn[2]] - mNodeLocs[conn[0]];
    return 0.5 * std::abs(e1[0] * e2[1] - e1[1] * e2[0]);
  }

  // The gradient of the linear interpolant of fn over element elemIndex.
  stk::math::Vector3d element_gradient(
      const unsigned elemIndex, const std::function<double(const stk::math::Vector3d &)> & fn) const
  {
    const auto & conn = meshSpec.allElementConn[elemIndex];
    const stk::math::Vector3d & x0 = mNodeLocs[conn[0]];
    const stk::math::Vector3d e1 = mNodeLocs[conn[1]] - x0, e2 = mNodeLocs[conn[2]] - x0;
    const double d1 = fn(mNodeLocs[conn[1]]) - fn(x0), d2 = fn(mNodeLocs[conn[2]]) - fn(x0);
    const double det = e1[0] * e2[1] - e1[1] * e2[0];
    return stk::math::Vector3d(
        (d1 * e2[1] - d2 * e1[1]) / det, (e1[0] * d2 - e2[0] * d1) / det, 0.);
  }

  // Sum of A_e g_e / sum of A_e over the given elements.
  stk::math::Vector3d area_weighted_gradient(const std::vector<unsigned> & elemIndices,
      const std::function<double(const stk::math::Vector3d &)> & fn) const
  {
    stk::math::Vector3d sum = stk::math::Vector3d::ZERO;
    double area = 0.;
    for (unsigned e : elemIndices)
    {
      sum += element_area(e) * element_gradient(e, fn);
      area += element_area(e);
    }
    return sum / area;
  }

  FourRightTris meshSpec;
  std::vector<stk::math::Vector3d> mNodeLocs;
  FieldRef mScalar;
};

// The averaging divides by a nodal volume summed over all of a node's elements, so in parallel it
// needs the aura.
class NodalGradient : public NodalGradientFixture
{
protected:
  NodalGradient() : NodalGradientFixture(stk::mesh::BulkData::AUTO_AURA) {}
};

class NodalGradientNoAura : public NodalGradientFixture
{
protected:
  NodalGradientNoAura() : NodalGradientFixture(stk::mesh::BulkData::NO_AUTO_AURA) {}
};

static double piecewise_linear_test_field(const stk::math::Vector3d & x)
{
  return x[0] * x[0] + 2. * x[0] * x[1] - 3. * x[1] * x[1];
}

TEST_F(NodalGradient, linearField_exactAtEveryNode)
{
  if (!is_valid_proc_size_for_test()) return;
  const stk::math::Vector3d gold(0.3, -0.7, 0.);
  set_scalar([&](const stk::math::Vector3d & x) { return 0.2 + Dot(gold, x); });
  update();

  for (unsigned n = 0; n < mNodeLocs.size(); ++n)
    if (has_node(n)) expect_near_absolute(gold, gradient_at(n), 1.e-14);
}

TEST_F(NodalGradient, piecewiseLinearField_areaWeightedAverageOfAttachedElements)
{
  if (!is_valid_proc_size_for_test()) return;
  set_scalar(piecewise_linear_test_field);
  update();

  const std::vector<std::vector<unsigned>> nodeElements{
      {0, 3}, {0, 1, 3}, {1, 2}, {0, 1, 2}, {2}, {3}};
  for (unsigned n = 0; n < nodeElements.size(); ++n)
    if (has_node(n))
      expect_near_absolute(area_weighted_gradient(nodeElements[n], piecewise_linear_test_field),
          gradient_at(n),
          1.e-13);
}

TEST_F(NodalGradient, inactiveElement_doesNotContribute)
{
  if (!is_valid_proc_size_for_test()) return;
  std::vector<stk::mesh::Entity> elemsToDeactivate;
  if (is_owned_element(1)) elemsToDeactivate.push_back(get_assigned_element_for_index(1));
  mMesh.batch_change_entity_parts(elemsToDeactivate,
      std::vector<stk::mesh::PartVector>(elemsToDeactivate.size()),
      std::vector<stk::mesh::PartVector>(
          elemsToDeactivate.size(), {&mBuilder.get_aux_meta().active_part()}));

  set_scalar(piecewise_linear_test_field);
  update();

  const std::vector<std::vector<unsigned>> nodeElements{{0, 3}, {0, 3}, {2}, {0, 2}, {2}, {3}};
  for (unsigned n = 0; n < nodeElements.size(); ++n)
    if (has_node(n))
      expect_near_absolute(area_weighted_gradient(nodeElements[n], piecewise_linear_test_field),
          gradient_at(n),
          1.e-13);
}

TEST_F(NodalGradientNoAura, throwsInParallel)
{
  if (!is_valid_proc_size_for_test()) return;
  set_scalar(piecewise_linear_test_field);
  if (mMesh.parallel_size() == 1)
  {
    EXPECT_NO_THROW(update());
  }
  else
  {
    EXPECT_ANY_THROW(update());
  }
}

} // namespace krino
