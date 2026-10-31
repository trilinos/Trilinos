// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Edge.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_MeshSizeMetric.hpp>
#include <Akri_MetricEdgeMarking.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <array>
#include <cmath>
#include <functional>

namespace krino {

static constexpr double sqrt2 = 1.4142135623730951;

/*
 * FourRightTris (unit legs):
 *    4
 *    | \
 *    3---2       9 edges: 0-1, 0-3, 0-5, 1-2, 2-3, 3-4 have length 1;
 *    |\  |                1-3, 1-5, 2-4 have length sqrt(2).
 *    | \ |
 *    |  \|
 *    0---1
 *    | /
 *    5
 */
class FourRightTrisMetricRefinement : public RefinementFixture<FourRightTris>
{
public:
  FourRightTrisMetricRefinement()
  {
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    myMetricField = &meta.declare_field<double>(stk::topology::NODE_RANK, "nodal_metric");
    stk::mesh::put_field_on_mesh(*myMetricField, meta.universal_part(), 4, nullptr);

    set_valid_proc_sizes_for_test({1,2,4});
    this->build_mesh_from_spec(meshSpec);
  }

protected:
  FieldRef metric_field() const { return FieldRef(myMetricField); }
  FieldRef coords_field() const { return FieldRef(this->get_coordinates_field()); }

  // Diagonal metric with sizes hx, hy, as a function of position, at every node.
  void set_metric(const std::function<std::array<double,2>(const stk::math::Vector3d &)> & sizes)
  {
    auto data = metric_field().data<double, stk::mesh::ReadWrite>();
    for (auto * bucket : mMesh.buckets(stk::topology::NODE_RANK))
      for (auto node : *bucket)
      {
        const auto [hx, hy] = sizes(get_vector_field(mMesh, coords_field(), node, 2));
        auto m = data.entity_values(node);
        m(stk::mesh::ComponentIdx(0)) = 1./(hx*hx);
        m(stk::mesh::ComponentIdx(1)) = 0.;
        m(stk::mesh::ComponentIdx(2)) = 0.;
        m(stk::mesh::ComponentIdx(3)) = 1./(hy*hy);
      }
  }

  void set_isotropic_metric(const double h) { set_metric([h](const stk::math::Vector3d &) { return std::array<double,2>{{h, h}}; }); }

  // Marks edges too long and refines them; returns the number of nodes added (one per refined edge).
  size_t mark_and_refine(const MeshSizeMetric & sizeMetric, const double maxRelativeEdgeLength)
  {
    clear_edge_marker();
    mark_edges_too_long(mMesh, myRefinement, coords_field(), sizeMetric, maxRelativeEdgeLength, myEdgeBasedEdgeMarker);
    const size_t numNodesBefore = get_global_num_entities(mMesh, stk::topology::NODE_RANK);
    do_edge_refinement();
    return get_global_num_entities(mMesh, stk::topology::NODE_RANK) - numNodesBefore;
  }

  // As above, with the nodal metric field.
  size_t mark_and_refine(const double maxRelativeEdgeLength)
  {
    clear_edge_marker();
    const NodalMetricTensorFieldMeshSizeMetric sizeMetric(mMesh, metric_field());
    mark_edges_too_long(mMesh, myRefinement, coords_field(), sizeMetric, maxRelativeEdgeLength, myEdgeBasedEdgeMarker);
    const size_t numNodesBefore = get_global_num_entities(mMesh, stk::topology::NODE_RANK);
    do_edge_refinement();
    return get_global_num_entities(mMesh, stk::topology::NODE_RANK) - numNodesBefore;
  }

  std::array<double,4> get_metric(const stk::mesh::Entity node) const
  {
    const auto m = metric_field().data<double>().entity_values(node);
    return {{m(stk::mesh::ComponentIdx(0)), m(stk::mesh::ComponentIdx(1)), m(stk::mesh::ComponentIdx(2)), m(stk::mesh::ComponentIdx(3))}};
  }

  // Largest relative length in sizeMetric of the edges of the owned leaf (non-parent) elements, over all processors.
  double max_relative_edge_length(const MeshSizeMetric & sizeMetric) const
  {
    double maxLength = 0.;
    std::vector<Edge> elemEdges;
    for (auto * bucket : mMesh.get_buckets(stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() & !myRefinement.parent_part()))
      for (auto elem : *bucket)
      {
        fill_entity_edges(mMesh, elem, elemEdges);
        for (auto && edge : elemEdges)
        {
          const std::array<stk::mesh::Entity,2> nodes = get_edge_nodes(edge);
          const stk::math::Vector3d edgeVec = get_vector_field(mMesh, coords_field(), nodes[1], 2) - get_vector_field(mMesh, coords_field(), nodes[0], 2);
          maxLength = std::max(maxLength, sizeMetric.compute_relative_edge_length(nodes[0], nodes[1], edgeVec));
        }
      }
    return stk::get_global_max(mComm, maxLength);
  }

  // Largest extent along each axis of the edges of the owned leaf elements, over all processors.
  std::array<double,2> max_leaf_edge_extents() const
  {
    std::array<double,2> maxExtents{{0., 0.}};
    std::vector<Edge> elemEdges;
    for (auto * bucket : mMesh.get_buckets(stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().locally_owned_part() & !myRefinement.parent_part()))
      for (auto elem : *bucket)
      {
        fill_entity_edges(mMesh, elem, elemEdges);
        for (auto && edge : elemEdges)
        {
          const std::array<stk::mesh::Entity,2> nodes = get_edge_nodes(edge);
          const stk::math::Vector3d edgeVec = get_vector_field(mMesh, coords_field(), nodes[1], 2) - get_vector_field(mMesh, coords_field(), nodes[0], 2);
          for (unsigned d=0; d<2; ++d)
            maxExtents[d] = std::max(maxExtents[d], std::abs(edgeVec[d]));
        }
      }
    return {{stk::get_global_max(mComm, maxExtents[0]), stk::get_global_max(mComm, maxExtents[1])}};
  }

  size_t refine_too_long(const MeshSizeMetric & sizeMetric, const unsigned maxNumPasses)
  {
    return refine_edges_too_long(mMesh, myRefinement, coords_field(), sizeMetric, sqrt2, maxNumPasses);
  }

  stk::mesh::Field<double> * myMetricField{nullptr};
};

TEST_F(FourRightTrisMetricRefinement, isotropicMetricLargerThanEveryEdge_nothingRefined)
{
  if (!is_valid_proc_size_for_test()) return;

  set_isotropic_metric(2.); // relative lengths at most sqrt(2)/2
  EXPECT_EQ(0u, mark_and_refine(sqrt2));
}

TEST_F(FourRightTrisMetricRefinement, isotropicMetric_onlyEdgesLongerThanThresholdAreRefined)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size 1: unit edges have relative length 1, the three diagonals sqrt(2).  With threshold 1.2 only the diagonals.
  set_isotropic_metric(1.);
  EXPECT_EQ(3u, mark_and_refine(1.2));
}

TEST_F(FourRightTrisMetricRefinement, edgeAtThreshold_isNotRefined)
{
  if (!is_valid_proc_size_for_test()) return;

  // Strictly longer than the threshold is refined; the diagonals are exactly sqrt(2) long here.
  set_isotropic_metric(1.);
  EXPECT_EQ(0u, mark_and_refine(sqrt2 + 1.e-12));
}

TEST_F(FourRightTrisMetricRefinement, anisotropicMetric_refinesEdgesAlongTheSmallSize)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size 0.5 along x and 1 along y: edges with an x extent of 1 have metric length 2 (0-1, 2-3) or sqrt(5)
  // (1-3, 1-5, 2-4) and are refined; the vertical edges 0-3, 0-5, 1-2, 3-4 have length 1 and are not.
  set_metric([](const stk::math::Vector3d &) { return std::array<double,2>{{0.5, 1.}}; });
  EXPECT_EQ(5u, mark_and_refine(1.5));
}

TEST_F(FourRightTrisMetricRefinement, isotropicSizeMetrics_refineTheSameEdgesAsTheEquivalentMetricTensor)
{
  if (!is_valid_proc_size_for_test()) return;

  // Isotropic size 1 is the metric I; with threshold 1.2 only the three diagonals (length sqrt(2)) are too long.
  const ConstantIsotropicMeshSizeMetric constantSize(1.);
  EXPECT_EQ(3u, mark_and_refine(constantSize, 1.2));
}

TEST_F(FourRightTrisMetricRefinement, metricTensorFunction_refinesTheSameEdgesAsTheMetricTensorField)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size 0.5 along x and 1 along y, as a function: the same 5 edges as the field version above.
  const MetricTensorFunctionMeshSizeMetric metricFunction(mMesh, coords_field(),
      [](const stk::math::Vector3d &, double * m) { m[0] = 4.; m[1] = 0.; m[2] = 0.; m[3] = 1.; });
  EXPECT_EQ(5u, mark_and_refine(metricFunction, 1.5));
}

TEST_F(FourRightTrisMetricRefinement, constantMetricTensor_refinesEdgesAlongTheSmallSize)
{
  if (!is_valid_proc_size_for_test()) return;

  const ConstantMetricTensorMeshSizeMetric constantMetric(2, {4., 0., 0., 1.});
  EXPECT_EQ(5u, mark_and_refine(constantMetric, 1.5));
}

TEST_F(FourRightTrisMetricRefinement, isotropicSizeFunction_refinesOnlyWhereTheSizeIsSmall)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size 0.25 at y <= 0 and 10 above.  The edges with both ends at y <= 0 (0-1, 0-5, 1-5) have relative lengths
  // 4, 4 and 4 sqrt(2).  Edges from y = 0 to y = 1 (0-3, 1-2, 1-3) mix the two sizes, and their log-mean relative
  // lengths (about 1.06, 1.06, 1.5) stay under 3.
  const IsotropicSizeFunctionMeshSizeMetric sizeFunction(mMesh, coords_field(),
      [](const stk::math::Vector3d & x) { return (x[1] <= 0.) ? 0.25 : 10.; });
  EXPECT_EQ(3u, mark_and_refine(sizeFunction, 3.));
}

TEST_F(FourRightTrisMetricRefinement, newNodeMetric_isLogEuclideanMeanOfParents)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size h = 2^(x+y), so a midpoint of nodes with sizes h0, h1 gets size sqrt(h0 h1) (metric 1/(h0 h1)).
  set_metric([](const stk::math::Vector3d & x) { const double h = std::pow(2., x[0]+x[1]); return std::array<double,2>{{h, h}}; });

  ASSERT_GT(mark_and_refine(0.1), 0u); // refine every edge
  update_size_for_new_nodes(mMesh, myRefinement, NodalMetricTensorFieldMeshSizeMetric(mMesh, metric_field()));

  size_t numChecked = 0;
  myRefinement.for_each_newly_refined_node(
      [&](const stk::mesh::Entity child, const stk::mesh::Entity * parentBeg, const stk::mesh::Entity * parentEnd)
      {
        ASSERT_EQ(2, parentEnd - parentBeg);
        const double h0 = std::pow(2., get_vector_field(mMesh, coords_field(), parentBeg[0], 2)[0] + get_vector_field(mMesh, coords_field(), parentBeg[0], 2)[1]);
        const double h1 = std::pow(2., get_vector_field(mMesh, coords_field(), parentBeg[1], 2)[0] + get_vector_field(mMesh, coords_field(), parentBeg[1], 2)[1]);
        const auto m = get_metric(child);
        EXPECT_NEAR(1./(h0*h1), m[0], 1.e-12);
        EXPECT_NEAR(0., m[1], 1.e-12);
        EXPECT_NEAR(0., m[2], 1.e-12);
        EXPECT_NEAR(1./(h0*h1), m[3], 1.e-12);
        ++numChecked;
      });
  EXPECT_GT(numChecked, 0u);
}

TEST_F(FourRightTrisMetricRefinement, refineEdgesTooLong_constantIsotropicSize_convergesWithEveryEdgeShortEnough)
{
  if (!is_valid_proc_size_for_test()) return;

  const ConstantIsotropicMeshSizeMetric sizeMetric(0.3);
  EXPECT_EQ(0u, refine_too_long(sizeMetric, 6));
  EXPECT_LE(max_relative_edge_length(sizeMetric), sqrt2);
}

TEST_F(FourRightTrisMetricRefinement, refineEdgesTooLong_anisotropicMetricFunction_convergesAndRefinesMoreAlongTheSmallSize)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size 0.2 along x and 0.6 along y.
  const MetricTensorFunctionMeshSizeMetric sizeMetric(mMesh, coords_field(),
      [](const stk::math::Vector3d &, double * m) { m[0] = 1./(0.2*0.2); m[1] = 0.; m[2] = 0.; m[3] = 1./(0.6*0.6); });
  EXPECT_EQ(0u, refine_too_long(sizeMetric, 8));
  EXPECT_LE(max_relative_edge_length(sizeMetric), sqrt2);

  // Edges may be up to sqrt(2) 0.2 long along x but sqrt(2) 0.6 along y, so the y edges stay longer.
  const std::array<double,2> maxExtents = max_leaf_edge_extents();
  EXPECT_LE(maxExtents[0], sqrt2*0.2);
  EXPECT_GT(maxExtents[1], maxExtents[0]);
}

TEST_F(FourRightTrisMetricRefinement, refineEdgesTooLong_nothingTooLong_leavesMeshUnchanged)
{
  if (!is_valid_proc_size_for_test()) return;

  const ConstantIsotropicMeshSizeMetric sizeMetric(2.);
  const size_t numNodesBefore = get_global_num_entities(mMesh, stk::topology::NODE_RANK);
  EXPECT_EQ(0u, refine_too_long(sizeMetric, 4));
  EXPECT_EQ(numNodesBefore, get_global_num_entities(mMesh, stk::topology::NODE_RANK));
}

TEST_F(FourRightTrisMetricRefinement, refineEdgesTooLong_tooFewPasses_returnsNumberOfEdgesStillTooLong)
{
  if (!is_valid_proc_size_for_test()) return;

  const ConstantIsotropicMeshSizeMetric sizeMetric(0.05);
  const size_t numStillTooLong = refine_too_long(sizeMetric, 1);
  EXPECT_LT(0u, numStillTooLong);

  clear_edge_marker();
  EXPECT_EQ(numStillTooLong, stk::get_global_sum(mComm, mark_edges_too_long(mMesh, myRefinement, coords_field(), sizeMetric, sqrt2, myEdgeBasedEdgeMarker)));
  EXPECT_GT(max_relative_edge_length(sizeMetric), sqrt2);
}

TEST_F(FourRightTrisMetricRefinement, refineEdgesTooLong_nodalMetricField_updatesNewNodesAndConverges)
{
  if (!is_valid_proc_size_for_test()) return;

  // Size from 0.25 at x+y = -1 to 1 at x+y = 3, at the original nodes; new nodes get the mean of their parents.
  set_metric([](const stk::math::Vector3d & x) { const double h = 0.25*std::pow(4., (x[0]+x[1]+1.)/4.); return std::array<double,2>{{h, h}}; });
  const NodalMetricTensorFieldMeshSizeMetric sizeMetric(mMesh, metric_field());
  EXPECT_EQ(0u, refine_too_long(sizeMetric, 8));
  EXPECT_LE(max_relative_edge_length(sizeMetric), sqrt2);
}

} // namespace krino
