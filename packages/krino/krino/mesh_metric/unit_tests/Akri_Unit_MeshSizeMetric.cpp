// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSizeMetric.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_MetricTensorUtils.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace krino {

// TwoRightTris: nodes 0 (0,0), 1 (1,0), 2 (0,1), 3 (-1,0).  On 2 processors every node is present on both.
class MeshSizeMetricFixture : public StkMeshTriFixture
{
protected:
  MeshSizeMetricFixture()
  {
    stk::mesh::MetaData & meta = mMesh.mesh_meta_data();
    mySizeField = &meta.declare_field<double>(stk::topology::NODE_RANK, "size");
    stk::mesh::put_field_on_mesh(*mySizeField, meta.universal_part(), 1, nullptr);
    myMetricField = &meta.declare_field<double>(stk::topology::NODE_RANK, "metric");
    stk::mesh::put_field_on_mesh(*myMetricField, meta.universal_part(), 4, nullptr);
    build_mesh_from_spec(meshSpec);
  }

  FieldRef coords() const { return FieldRef(get_coordinates_field()); }
  FieldRef size_field() const { return FieldRef(mySizeField); }
  FieldRef metric_field() const { return FieldRef(myMetricField); }
  stk::mesh::Entity node(const unsigned i) const { return get_assigned_node_for_index(i); }
  stk::math::Vector3d x(const unsigned i) const { return get_vector_field(mMesh, coords(), node(i), 2); }
  stk::math::Vector3d edge(const unsigned i, const unsigned j) const { return x(j) - x(i); }

  void set_size(const unsigned i, const double h) { size_field().data<double, stk::mesh::ReadWrite>().entity_values(node(i))(stk::mesh::ComponentIdx(0)) = h; }
  double get_size(const unsigned i) const { return size_field().data<double>().entity_values(node(i))(stk::mesh::ComponentIdx(0)); }
  void set_metric(const unsigned i, const std::array<double,4> & m)
  {
    auto values = metric_field().data<double, stk::mesh::ReadWrite>().entity_values(node(i));
    for (unsigned c=0; c<4; ++c) values(stk::mesh::ComponentIdx(c)) = m[c];
  }
  std::array<double,4> get_metric(const unsigned i) const
  {
    const auto values = metric_field().data<double>().entity_values(node(i));
    return {{values(stk::mesh::ComponentIdx(0)), values(stk::mesh::ComponentIdx(1)), values(stk::mesh::ComponentIdx(2)), values(stk::mesh::ComponentIdx(3))}};
  }

  TwoRightTris meshSpec;
  stk::mesh::Field<double> * mySizeField{nullptr};
  stk::mesh::Field<double> * myMetricField{nullptr};
};

static constexpr double tol = 1.e-12;
static const std::array<std::pair<unsigned,unsigned>,3> edgesOfFirstTri{{{0,1}, {1,2}, {0,2}}};

static std::array<double,4> isotropic_metric(const double h) { return {{1./(h*h), 0., 0., 1./(h*h)}}; }

TEST_F(MeshSizeMetricFixture, constantIsotropic_isLengthOverSize)
{
  if (!is_valid_proc_size_for_test()) return;
  const ConstantIsotropicMeshSizeMetric sizeMetric(0.5);
  EXPECT_NEAR(edge(0,1).length()/0.5, sizeMetric.compute_relative_edge_length(node(0), node(1), edge(0,1)), tol);
}

TEST_F(MeshSizeMetricFixture, nodalIsotropicField_isLogMeanOfNodalRelativeLengths)
{
  if (!is_valid_proc_size_for_test()) return;
  set_size(0, 0.5);
  set_size(1, 2.);
  const NodalIsotropicFieldMeshSizeMetric sizeMetric(mMesh, size_field());
  const double length = edge(0,1).length();
  EXPECT_NEAR(metric::logarithmic_mean(length/0.5, length/2.), sizeMetric.compute_relative_edge_length(node(0), node(1), edge(0,1)), tol);
}

TEST_F(MeshSizeMetricFixture, isotropicSizeFunction_matchesNodalFieldSampledAtNodes)
{
  if (!is_valid_proc_size_for_test()) return;
  const auto sizeOfX = [](const stk::math::Vector3d & p) { return 0.3 + p[0] + 2.*p[1]; };
  for (unsigned i=0; i<meshSpec.nodeLocs.size(); ++i)
    set_size(i, sizeOfX(x(i)));
  const IsotropicSizeFunctionMeshSizeMetric functionMetric(mMesh, coords(), sizeOfX);
  const NodalIsotropicFieldMeshSizeMetric fieldMetric(mMesh, size_field());
  for (auto [i,j] : edgesOfFirstTri)
    EXPECT_NEAR(fieldMetric.compute_relative_edge_length(node(i), node(j), edge(i,j)),
        functionMetric.compute_relative_edge_length(node(i), node(j), edge(i,j)), tol);
}

TEST_F(MeshSizeMetricFixture, isotropicSizeAndEquivalentMetricTensor_giveTheSameLengths)
{
  if (!is_valid_proc_size_for_test()) return;
  set_size(0, 0.5); set_metric(0, isotropic_metric(0.5));
  set_size(1, 2.);  set_metric(1, isotropic_metric(2.));
  const NodalIsotropicFieldMeshSizeMetric sizeMetric(mMesh, size_field());
  const NodalMetricTensorFieldMeshSizeMetric tensorMetric(mMesh, metric_field());
  EXPECT_NEAR(sizeMetric.compute_relative_edge_length(node(0), node(1), edge(0,1)),
      tensorMetric.compute_relative_edge_length(node(0), node(1), edge(0,1)), tol);

  const ConstantIsotropicMeshSizeMetric constantSize(0.5);
  const std::array<double,4> constantIsotropicMetric = isotropic_metric(0.5);
  const ConstantMetricTensorMeshSizeMetric constantTensor(2, std::vector<double>(constantIsotropicMetric.begin(), constantIsotropicMetric.end()));
  EXPECT_NEAR(constantSize.compute_relative_edge_length(node(0), node(1), edge(0,1)),
      constantTensor.compute_relative_edge_length(node(0), node(1), edge(0,1)), tol);
}

TEST_F(MeshSizeMetricFixture, nodalMetricTensorField_matchesComputeRelativeEdgeLength)
{
  if (!is_valid_proc_size_for_test()) return;
  const std::array<double,4> m0{{4., 0.5, 0.5, 1.}};
  const std::array<double,4> m1{{1., -0.2, -0.2, 9.}};
  set_metric(0, m0);
  set_metric(1, m1);
  const NodalMetricTensorFieldMeshSizeMetric tensorMetric(mMesh, metric_field());
  EXPECT_NEAR(metric::compute_relative_edge_length(2, x(0), m0.data(), x(1), m1.data()),
      tensorMetric.compute_relative_edge_length(node(0), node(1), edge(0,1)), tol);
}

TEST_F(MeshSizeMetricFixture, metricTensorFunction_matchesNodalFieldSampledAtNodes)
{
  if (!is_valid_proc_size_for_test()) return;
  const auto metricOfX = [](const stk::math::Vector3d & p, double * m) { m[0] = 1.+p[0]; m[1] = 0.1; m[2] = 0.1; m[3] = 2.+p[1]; };
  for (unsigned i=0; i<meshSpec.nodeLocs.size(); ++i)
  {
    std::array<double,4> m;
    metricOfX(x(i), m.data());
    set_metric(i, m);
  }
  const MetricTensorFunctionMeshSizeMetric functionMetric(mMesh, coords(), metricOfX);
  const NodalMetricTensorFieldMeshSizeMetric fieldMetric(mMesh, metric_field());
  for (auto [i,j] : edgesOfFirstTri)
    EXPECT_NEAR(fieldMetric.compute_relative_edge_length(node(i), node(j), edge(i,j)),
        functionMetric.compute_relative_edge_length(node(i), node(j), edge(i,j)), tol);
}

TEST_F(MeshSizeMetricFixture, constantMetricTensor_rejectsInvalidMetric)
{
  EXPECT_ANY_THROW(ConstantMetricTensorMeshSizeMetric(2, {1., 0., 0.}));          // wrong size
  EXPECT_ANY_THROW(ConstantMetricTensorMeshSizeMetric(2, {1., 0.5, 0., 1.}));     // not symmetric
  EXPECT_ANY_THROW(ConstantMetricTensorMeshSizeMetric(2, {1., 0., 0., -1.}));     // not positive definite
  EXPECT_ANY_THROW(ConstantIsotropicMeshSizeMetric(0.));
}

TEST_F(MeshSizeMetricFixture, nodalFields_newNodeGetsLogEuclideanMeanOfParents)
{
  if (!is_valid_proc_size_for_test()) return;
  // Node 2 plays the new node, with parents 0 and 1.
  set_size(0, 0.5); set_metric(0, isotropic_metric(0.5));
  set_size(1, 2.);  set_metric(1, isotropic_metric(2.));
  const std::array<stk::mesh::Entity,2> parents{{node(0), node(1)}};

  NodalIsotropicFieldMeshSizeMetric(mMesh, size_field()).update_size_at_new_node(node(2), parents.data(), parents.data()+2);
  EXPECT_NEAR(1., get_size(2), tol); // geometric mean of 0.5 and 2

  NodalMetricTensorFieldMeshSizeMetric(mMesh, metric_field()).update_size_at_new_node(node(2), parents.data(), parents.data()+2);
  const auto m = get_metric(2);
  const auto gold = isotropic_metric(1.);
  for (unsigned c=0; c<4; ++c)
    EXPECT_NEAR(gold[c], m[c], tol);
}

} // namespace krino
