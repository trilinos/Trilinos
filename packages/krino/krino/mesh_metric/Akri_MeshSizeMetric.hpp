// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_METRIC_AKRI_MESHSIZEMETRIC_HPP_
#define KRINO_KRINO_MESH_METRIC_AKRI_MESHSIZEMETRIC_HPP_

#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>

#include <array>
#include <functional>
#include <vector>

namespace krino {

// Target mesh size, isotropic or anisotropic (a metric tensor), at mesh nodes, as used to decide whether edges are
// too long (modelled on ngs interface::MeshSizeMetric).  An edge is measured with the size at both its nodes: the
// logarithmic mean of its relative lengths at the two nodes, as ngs does.
class MeshSizeMetric
{
public:
  virtual ~MeshSizeMetric() = default;

  // Length of the edge from node0 to node1 relative to the target size; edge = x1 - x0.
  virtual double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1,
      const stk::math::Vector3d & edge) const = 0;

  // Sets the size at a node created by refinement from its parent nodes (the ends of a refined edge, or the nodes
  // of a refined element or face).  No-op for constant or analytic sizes.
  virtual void update_size_at_new_node(const stk::mesh::Entity /*newNode*/, const stk::mesh::Entity * /*parentBegin*/,
      const stk::mesh::Entity * /*parentEnd*/) const {}

  // After update_size_at_new_node on the new nodes here, makes the sizes consistent on ghost nodes.  No-op for
  // constant or analytic sizes.
  virtual void communicate_updated_sizes(const stk::mesh::BulkData & /*mesh*/) const {}
};

// Uniform isotropic size h.
class ConstantIsotropicMeshSizeMetric : public MeshSizeMetric
{
public:
  explicit ConstantIsotropicMeshSizeMetric(const double size);
  double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const override;
private:
  double mySize;
};

// Isotropic size given by a scalar nodal field.  New nodes get the geometric mean of their parents' sizes (the
// isotropic case of the log-Euclidean mean).
class NodalIsotropicFieldMeshSizeMetric : public MeshSizeMetric
{
public:
  NodalIsotropicFieldMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef sizeField);
  double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const override;
  void update_size_at_new_node(const stk::mesh::Entity newNode, const stk::mesh::Entity * parentBegin, const stk::mesh::Entity * parentEnd) const override;
  void communicate_updated_sizes(const stk::mesh::BulkData & mesh) const override;
private:
  double get_size(const stk::mesh::Entity node) const;
  FieldRef mySizeField;
};

// Isotropic size given by a function of position.
class IsotropicSizeFunctionMeshSizeMetric : public MeshSizeMetric
{
public:
  using SizeFunction = std::function<double(const stk::math::Vector3d &)>;
  IsotropicSizeFunctionMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const SizeFunction & size);
  double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const override;
private:
  const stk::mesh::BulkData & myMesh;
  FieldRef myCoordsField;
  SizeFunction mySize;
};

// Uniform metric tensor M (spatialDim*spatialDim row-major, symmetric positive definite).
class ConstantMetricTensorMeshSizeMetric : public MeshSizeMetric
{
public:
  ConstantMetricTensorMeshSizeMetric(const unsigned spatialDim, const std::vector<double> & metric);
  double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const override;
private:
  unsigned mySpatialDim;
  std::vector<double> myMetric;
};

// Metric tensor given by a nodal field holding M at each node (spatialDim*spatialDim row-major).  New nodes get the
// log-Euclidean mean of their parents' metrics.
class NodalMetricTensorFieldMeshSizeMetric : public MeshSizeMetric
{
public:
  NodalMetricTensorFieldMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef metricField);
  double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const override;
  void update_size_at_new_node(const stk::mesh::Entity newNode, const stk::mesh::Entity * parentBegin, const stk::mesh::Entity * parentEnd) const override;
  void communicate_updated_sizes(const stk::mesh::BulkData & mesh) const override;
private:
  std::vector<double> get_metric(const stk::mesh::Entity node) const;
  unsigned mySpatialDim;
  FieldRef myMetricField;
};

// Metric tensor given by a function of position, filling M (spatialDim*spatialDim row-major) at a point.
class MetricTensorFunctionMeshSizeMetric : public MeshSizeMetric
{
public:
  using MetricFunction = std::function<void(const stk::math::Vector3d &, double * metric)>;
  MetricTensorFunctionMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const MetricFunction & metric);
  double compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const override;
private:
  std::vector<double> get_metric(const stk::mesh::Entity node) const;
  const stk::mesh::BulkData & myMesh;
  unsigned mySpatialDim;
  FieldRef myCoordsField;
  MetricFunction myMetric;
};

} // namespace krino

#endif /* KRINO_KRINO_MESH_METRIC_AKRI_MESHSIZEMETRIC_HPP_ */
