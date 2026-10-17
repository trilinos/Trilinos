// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MeshQuality.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_QualityMetric.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <cmath>

namespace krino {

// Mesh-driven coverage for the free functions in Akri_MeshQuality, exercised on
// single reference elements built through StkMeshFixture.  Gold values match the
// analytic element qualities verified in the quality_metric unit tests:
//   RegularTet / RegularTri (equilateral) -> 1.0 for both metrics.
//   RightTet (unit legs)  -> MeanRatio sqrt(2)/1.5^1.5, ScaledJacobian sqrt(2)/2.
//   RightTri (unit legs)  -> MeanRatio sqrt(3)/2,        ScaledJacobian 2/sqrt(6).
// Every element is on processor 0, so in parallel these also check that compute_mesh_quality's reduction gives the
// quality on processors without elements.

namespace {
const double tol = 1.e-6;

const double goldRightTetMeanRatio = std::sqrt(2.) / std::pow(1.5, 1.5);
const double goldRightTetScaledJacobian = std::sqrt(2.) / 2.;
const double goldRightTriMeanRatio = std::sqrt(3.) / 2.;
const double goldRightTriScaledJacobian = 2. / std::sqrt(6.);
}

//--------------------------------------------------------------------------
// compute_mesh_quality -- tets
//--------------------------------------------------------------------------

class MeshQualityTet : public StkMeshTetFixture
{
protected:
  template <typename MeshSpecType>
  void build(const MeshSpecType & meshSpec)
  {
    build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Selector active_selector() { return get_aux_meta().active_part(); }
  FieldRef coords_field() { return FieldRef(get_mesh().mesh_meta_data().coordinate_field()); }
};

class MeshQualityTri : public StkMeshTriFixture
{
protected:
  template <typename MeshSpecType>
  void build(const MeshSpecType & meshSpec)
  {
    build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
  stk::mesh::Selector active_selector() { return get_aux_meta().active_part(); }
  FieldRef coords_field() { return FieldRef(get_mesh().mesh_meta_data().coordinate_field()); }
};

TEST_F(MeshQualityTet, regularTet_meanRatioAndScaledJacobianAreBest)
{
  build(RegularTet{});

  EXPECT_NEAR(1., compute_mesh_quality(get_mesh(), active_selector(), coords_field(), MeanRatioQualityMetric{}), tol);
  EXPECT_NEAR(1., compute_mesh_quality(get_mesh(), active_selector(), coords_field(), ScaledJacobianQualityMetric{}), tol);
}

TEST_F(MeshQualityTet, rightTet_matchesAnalyticQuality)
{
  build(RightTet{});

  EXPECT_NEAR(goldRightTetMeanRatio, compute_mesh_quality(get_mesh(), active_selector(), coords_field(), MeanRatioQualityMetric{}), tol);
  EXPECT_NEAR(goldRightTetScaledJacobian, compute_mesh_quality(get_mesh(), active_selector(), coords_field(), ScaledJacobianQualityMetric{}), tol);
}

TEST_F(MeshQualityTet, emptySelector_returnsBestValue)
{
  build(RightTet{});

  // A selector that matches no elements (every element is in the active part,
  // so its complement is empty) leaves the running minimum at the metric's best
  // value.
  const MeanRatioQualityMetric metric;
  const stk::mesh::Selector nothing = !get_aux_meta().active_part();
  EXPECT_NEAR(metric.get_best_value_for_metric(), compute_mesh_quality(get_mesh(), nothing, coords_field(), metric), tol);
}

//--------------------------------------------------------------------------
// compute_mesh_quality -- 2d triangles
//--------------------------------------------------------------------------

TEST_F(MeshQualityTri, regularTri_meanRatioAndScaledJacobianAreBest)
{
  build(RegularTri{});

  EXPECT_NEAR(1., compute_mesh_quality(get_mesh(), active_selector(), coords_field(), MeanRatioQualityMetric{}), tol);
  EXPECT_NEAR(1., compute_mesh_quality(get_mesh(), active_selector(), coords_field(), ScaledJacobianQualityMetric{}), tol);
}

TEST_F(MeshQualityTri, rightTri_matchesAnalyticQuality)
{
  build(RightTri{});

  EXPECT_NEAR(goldRightTriMeanRatio, compute_mesh_quality(get_mesh(), active_selector(), coords_field(), MeanRatioQualityMetric{}), tol);
  EXPECT_NEAR(goldRightTriScaledJacobian, compute_mesh_quality(get_mesh(), active_selector(), coords_field(), ScaledJacobianQualityMetric{}), tol);
}

//--------------------------------------------------------------------------
// compute_minimum_element_quality (explicit element container)
//--------------------------------------------------------------------------

TEST_F(MeshQualityTet, computeMinimumElementQuality_overOwnedElementsMatchesMeshQuality)
{
  build(RightTet{});

  const ScaledJacobianQualityMetric metric;
  const CoordinatesFieldRef coordsField(get_mesh().mesh_meta_data().coordinate_field(), get_mesh().mesh_meta_data().spatial_dimension());
  const double localQuality = compute_minimum_element_quality(get_mesh(), coordsField, active_selector(), get_owned_elements(), metric);
  EXPECT_NEAR(compute_mesh_quality(get_mesh(), active_selector(), coords_field(), metric), stk::get_global_min(get_comm(), localQuality), tol);
  EXPECT_NEAR(goldRightTetScaledJacobian, stk::get_global_min(get_comm(), localQuality), tol);
}

//--------------------------------------------------------------------------
// compute_quality_if_node_is_moved_terminating_early_if_below_threshold
//--------------------------------------------------------------------------

TEST_F(MeshQualityTet, qualityIfNodeMoved_backToOriginalIsBestForwardToFlattenIsZero)
{
  const RegularTet meshSpec;
  build(meshSpec);
  if (!is_node_present_on_proc(3)) return; // the quality is computed from the elements of the node on this processor

  const ScaledJacobianQualityMetric metric;
  const stk::mesh::Entity apex = get_assigned_node_for_index(3);
  const stk::math::Vector3d apexLoc = get_node_coordinates(apex);

  // Moving the node back to its own location leaves the regular tet unchanged.
  const double qualityUnmoved = compute_quality_if_node_is_moved_terminating_early_if_below_threshold(
      get_mesh(), active_selector(), coords_field(), apex, apexLoc, metric, 0.);
  EXPECT_NEAR(1., qualityUnmoved, tol);

  // Collapsing the apex onto another node of the tet makes it degenerate
  // (zero volume), so the moved quality is zero.
  const stk::math::Vector3d collapsedLoc = get_node_coordinates(get_assigned_node_for_index(0));
  const double qualityFlattened = compute_quality_if_node_is_moved_terminating_early_if_below_threshold(
      get_mesh(), active_selector(), coords_field(), apex, collapsedLoc, metric, 0.);
  EXPECT_NEAR(0., qualityFlattened, tol);
  EXPECT_TRUE(metric.is_first_quality_metric_better_than_second(qualityUnmoved, qualityFlattened));
}

}
