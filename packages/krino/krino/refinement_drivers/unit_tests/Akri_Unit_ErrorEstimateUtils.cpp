// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_ErrorEstimateUtils.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <cmath>
#include <functional>
#include <vector>

namespace krino {

// Errors are accumulated on owned elements from the normals of their side neighbors, which may be ghosts.  So the
// tests set normals on every element present on the processor (as the estimators' callers do) and check each
// element's error on its owner.
template <typename MESHSPEC>
class ErrorEstimateFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
protected:
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;

  ErrorEstimateFixture()
  {
    this->build_mesh_from_spec(meshSpec);
  }

  stk::mesh::Entity elem(const unsigned index) const { return this->get_assigned_element_for_index(index); }
  stk::mesh::Selector all() const { return mMesh.mesh_meta_data().universal_part(); }

  std::vector<stk::math::Vector3d> no_normals() const { return create_vector_indexable_by_entity_offset(mMesh, stk::topology::ELEMENT_RANK, stk::math::Vector3d::ZERO); }
  std::vector<double> zero_errors() const { return create_vector_indexable_by_entity_offset(mMesh, stk::topology::ELEMENT_RANK, 0.); }

  void set_normal(std::vector<stk::math::Vector3d> & normals, const unsigned index, const stk::math::Vector3d & normal) const
  {
    if (this->is_element_present_on_proc(index))
      normals[elem(index).local_offset()] = normal;
  }

  // Value of the owned element with this index, on every processor.
  double reduce_from_owner(const unsigned index, const std::function<double(stk::mesh::Entity)> & valueOfElem) const
  {
    return stk::get_global_sum(mComm, this->is_owned_element(index) ? valueOfElem(elem(index)) : 0.);
  }

  double error(const std::vector<double> & errors, const unsigned index) const
  {
    return reduce_from_owner(index, [&errors](const stk::mesh::Entity e) { return errors[e.local_offset()]; });
  }

  double rms_edge_length(const unsigned index) const
  {
    return reduce_from_owner(index, [this](const stk::mesh::Entity e)
    {
      std::vector<stk::math::Vector3d> coords;
      fill_element_node_coordinates(mMesh, e, FieldRef(this->get_coordinates_field()), coords);
      return compute_simplex_RMS_edge_lengths(coords);
    });
  }

  void accumulate_angle(const std::vector<stk::math::Vector3d> & normals, std::vector<double> & errors) const
  {
    accumulate_facet_angle_errors_from_interface_normals(mMesh, all(), normals, errors);
  }

  void accumulate_chordal(const std::vector<stk::math::Vector3d> & normals, std::vector<double> & errors) const
  {
    accumulate_facet_chordal_errors_from_interface_normals(mMesh, FieldRef(this->get_coordinates_field()), all(), normals, errors);
  }

  MESHSPEC meshSpec;
};

// FourRightTris: the side neighbors form a chain, element 3 - element 0 - element 1 - element 2.  With round-robin
// owners on 2 or more processors, every neighbor pair is on different processors.
class FourRightTrisErrorEstimate : public ErrorEstimateFixture<FourRightTris>
{
protected:
  void set_normal(std::vector<stk::math::Vector3d> & normals, const unsigned index, const double angle) const
  {
    ErrorEstimateFixture<FourRightTris>::set_normal(normals, index, stk::math::Vector3d(std::cos(angle), std::sin(angle), 0.));
  }
};

static constexpr double tol = 1.e-12;

TEST_F(FourRightTrisErrorEstimate, facetAngle_flatInterfaceHasNoError)
{
  if (!is_valid_proc_size_for_test()) return;

  std::vector<stk::math::Vector3d> normals = no_normals();
  for (unsigned i=0; i<4; ++i)
    set_normal(normals, i, 0.3);
  std::vector<double> errors = zero_errors();
  accumulate_angle(normals, errors);
  for (unsigned i=0; i<4; ++i)
    EXPECT_NEAR(0., error(errors, i), tol);
}

TEST_F(FourRightTrisErrorEstimate, facetAngle_isMeanAngleToNeighborsWithNormals)
{
  if (!is_valid_proc_size_for_test()) return;

  // Normals on elements 3, 0, 1 (element 2 has none): element 0 has two neighbors with normals, elements 3 and 1 one each.
  const double theta1 = 0.2;
  const double theta3 = 0.5;
  std::vector<stk::math::Vector3d> normals = no_normals();
  set_normal(normals, 0, 0.);
  set_normal(normals, 1, theta1);
  set_normal(normals, 3, theta3);
  std::vector<double> errors = zero_errors();
  accumulate_angle(normals, errors);

  EXPECT_NEAR(0.5*(theta1+theta3), error(errors, 0), tol);
  EXPECT_NEAR(theta1, error(errors, 1), tol);   // neighbors 0 (normal) and 2 (none)
  EXPECT_NEAR(theta3, error(errors, 3), tol);   // neighbor 0
  EXPECT_NEAR(0., error(errors, 2), tol);       // no normal
}

TEST_F(FourRightTrisErrorEstimate, facetAngle_elementWithoutNeighborNormalsAddsNothing)
{
  if (!is_valid_proc_size_for_test()) return;

  std::vector<stk::math::Vector3d> normals = no_normals();
  set_normal(normals, 2, 0.7); // its only neighbor, element 1, has no normal
  std::vector<double> errors = zero_errors();
  accumulate_angle(normals, errors);
  EXPECT_EQ(0., error(errors, 2));
  EXPECT_FALSE(std::isnan(error(errors, 2)));
}

TEST_F(FourRightTrisErrorEstimate, facetAngle_accumulatesOverCalls)
{
  if (!is_valid_proc_size_for_test()) return;

  const double theta = 0.4;
  std::vector<stk::math::Vector3d> normals = no_normals();
  set_normal(normals, 0, 0.);
  set_normal(normals, 1, theta);
  std::vector<double> errors = zero_errors();
  accumulate_angle(normals, errors);
  accumulate_angle(normals, errors); // e.g. a second surface
  EXPECT_NEAR(2.*theta, error(errors, 0), tol);
  EXPECT_NEAR(2.*theta, error(errors, 1), tol);
}

TEST_F(FourRightTrisErrorEstimate, facetAngle_unselectedElementsGetAndGiveNothing)
{
  if (!is_valid_proc_size_for_test()) return;

  std::vector<stk::math::Vector3d> normals = no_normals();
  set_normal(normals, 0, 0.);
  set_normal(normals, 1, 0.4);
  std::vector<double> errors = zero_errors();
  // Selecting nothing: no element is processed and none is a neighbor.
  accumulate_facet_angle_errors_from_interface_normals(mMesh, stk::mesh::Selector(), normals, errors);
  for (unsigned i=0; i<4; ++i)
    EXPECT_EQ(0., error(errors, i));
}

TEST_F(FourRightTrisErrorEstimate, facetChordal_isElementSizeTimesCoefficientTimesSqrtOfNormalJump)
{
  if (!is_valid_proc_size_for_test()) return;

  // Elements 0 and 1 with normals at angle theta, d = cos(theta).  In 2D err = h sqrt(1/32) sqrt(1 - d).
  const double theta = 0.3;
  std::vector<stk::math::Vector3d> normals = no_normals();
  set_normal(normals, 0, 0.);
  set_normal(normals, 1, theta);
  std::vector<double> errors = zero_errors();
  accumulate_chordal(normals, errors);

  const double jumpFactor = std::sqrt(1./32.) * std::sqrt(1. - std::cos(theta));
  EXPECT_NEAR(rms_edge_length(0) * jumpFactor, error(errors, 0), tol);
  EXPECT_NEAR(rms_edge_length(1) * jumpFactor, error(errors, 1), tol);
  EXPECT_NEAR(0., error(errors, 2), tol);
  EXPECT_NEAR(0., error(errors, 3), tol);
}

TEST_F(FourRightTrisErrorEstimate, facetChordal_flatInterfaceHasNoErrorAndIsolatedElementAddsNothing)
{
  if (!is_valid_proc_size_for_test()) return;

  std::vector<stk::math::Vector3d> normals = no_normals();
  for (unsigned i=0; i<4; ++i)
    set_normal(normals, i, 1.1);
  std::vector<double> errors = zero_errors();
  accumulate_chordal(normals, errors);
  for (unsigned i=0; i<4; ++i)
    EXPECT_NEAR(0., error(errors, i), tol);

  std::vector<stk::math::Vector3d> isolated = no_normals();
  set_normal(isolated, 2, 0.7);
  std::vector<double> isolatedErrors = zero_errors();
  accumulate_chordal(isolated, isolatedErrors);
  EXPECT_EQ(0., error(isolatedErrors, 2));
}

using TwoRightTetsErrorEstimate = ErrorEstimateFixture<TwoRightTets>;

TEST_F(TwoRightTetsErrorEstimate, facetChordal_3dUsesSqrtOneSixthCoefficient)
{
  if (!is_valid_proc_size_for_test()) return;

  if (is_element_present_on_proc(0))
  {
    ASSERT_EQ(1u, get_selected_side_attached_elements(mMesh, all(), elem(0)).size());
  }

  const double theta = 0.25;
  std::vector<stk::math::Vector3d> normals = no_normals();
  set_normal(normals, 0, stk::math::Vector3d(0., 0., 1.));
  set_normal(normals, 1, stk::math::Vector3d(std::sin(theta), 0., std::cos(theta)));
  std::vector<double> errors = zero_errors();
  accumulate_facet_chordal_errors_from_interface_normals(mMesh, FieldRef(get_coordinates_field()), all(), normals, errors);

  const double jumpFactor = std::sqrt(1./6.) * std::sqrt(1. - std::cos(theta));
  EXPECT_NEAR(rms_edge_length(0) * jumpFactor, error(errors, 0), tol);
  EXPECT_NEAR(rms_edge_length(1) * jumpFactor, error(errors, 1), tol);
}

} // namespace krino
