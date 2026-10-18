// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshSizeMetric.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MetricTensorUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <cmath>

namespace krino {

static void check_nodal_field(const FieldRef field, const char * description)
{
  STK_ThrowRequireMsg(field.valid() && field.entity_rank() == stk::topology::NODE_RANK,
      "The " << description << " for a mesh size metric must be a nodal field.");
}

static std::vector<double> get_nodal_values(const FieldRef field, const unsigned numValues, const stk::mesh::Entity node)
{
  STK_ThrowRequireMsg(field.length(node) == numValues,
      "The field " << field.name() << " must have " << numValues << " components at every node.");
  const auto values = field.data<double>().entity_values(node);
  std::vector<double> result(numValues);
  for (unsigned i=0; i<numValues; ++i)
    result[i] = values(stk::mesh::ComponentIdx(i));
  return result;
}

static double compute_relative_edge_length_for_isotropic_sizes(const double size0, const double size1, const stk::math::Vector3d & edge)
{
  const double length = edge.length();
  return metric::logarithmic_mean(length/size0, length/size1);
}

//--------------------------------------------------------------------------------

ConstantIsotropicMeshSizeMetric::ConstantIsotropicMeshSizeMetric(const double size)
: mySize(size)
{
  STK_ThrowRequireMsg(size > 0., "The mesh size must be positive.");
}

double ConstantIsotropicMeshSizeMetric::compute_relative_edge_length(const stk::mesh::Entity, const stk::mesh::Entity, const stk::math::Vector3d & edge) const
{
  return edge.length()/mySize;
}

//--------------------------------------------------------------------------------

NodalIsotropicFieldMeshSizeMetric::NodalIsotropicFieldMeshSizeMetric(const stk::mesh::BulkData & /*mesh*/, const FieldRef sizeField)
: mySizeField(sizeField)
{
  check_nodal_field(sizeField, "size field");
}

double NodalIsotropicFieldMeshSizeMetric::get_size(const stk::mesh::Entity node) const
{
  return get_nodal_values(mySizeField, 1, node)[0];
}

double NodalIsotropicFieldMeshSizeMetric::compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const
{
  return compute_relative_edge_length_for_isotropic_sizes(get_size(node0), get_size(node1), edge);
}

void NodalIsotropicFieldMeshSizeMetric::update_size_at_new_node(const stk::mesh::Entity newNode, const stk::mesh::Entity * parentBegin, const stk::mesh::Entity * parentEnd) const
{
  STK_ThrowRequire(parentEnd > parentBegin);
  double sumLogSize = 0.;
  for (const stk::mesh::Entity * parent = parentBegin; parent != parentEnd; ++parent)
    sumLogSize += std::log(get_size(*parent));
  mySizeField.data<double, stk::mesh::ReadWrite>().entity_values(newNode)(stk::mesh::ComponentIdx(0)) = std::exp(sumLogSize/(parentEnd-parentBegin));
}

void NodalIsotropicFieldMeshSizeMetric::communicate_updated_sizes(const stk::mesh::BulkData & mesh) const
{
  stk::mesh::communicate_field_data(mesh, {&mySizeField.field()});
}

//--------------------------------------------------------------------------------

IsotropicSizeFunctionMeshSizeMetric::IsotropicSizeFunctionMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const SizeFunction & size)
: myMesh(mesh), myCoordsField(coordsField), mySize(size)
{
}

double IsotropicSizeFunctionMeshSizeMetric::compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const
{
  const unsigned dim = myMesh.mesh_meta_data().spatial_dimension();
  return compute_relative_edge_length_for_isotropic_sizes(mySize(get_vector_field(myMesh, myCoordsField, node0, dim)),
      mySize(get_vector_field(myMesh, myCoordsField, node1, dim)), edge);
}

//--------------------------------------------------------------------------------

ConstantMetricTensorMeshSizeMetric::ConstantMetricTensorMeshSizeMetric(const unsigned spatialDim, const std::vector<double> & metric)
: mySpatialDim(spatialDim), myMetric(metric)
{
  STK_ThrowRequireMsg(metric.size() == size_t(spatialDim) * spatialDim,
      "A " << spatialDim << "D metric tensor must have " << size_t(spatialDim) * spatialDim
           << " components.");
  STK_ThrowRequireMsg(metric::is_symmetric_positive_definite_metric(spatialDim, metric.data()), "The metric tensor must be symmetric positive definite.");
}

double ConstantMetricTensorMeshSizeMetric::compute_relative_edge_length(const stk::mesh::Entity, const stk::mesh::Entity, const stk::math::Vector3d & edge) const
{
  return metric::compute_relative_length_in_metric(mySpatialDim, myMetric.data(), edge);
}

//--------------------------------------------------------------------------------

NodalMetricTensorFieldMeshSizeMetric::NodalMetricTensorFieldMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef metricField)
: mySpatialDim(mesh.mesh_meta_data().spatial_dimension()), myMetricField(metricField)
{
  check_nodal_field(metricField, "metric field");
}

std::vector<double> NodalMetricTensorFieldMeshSizeMetric::get_metric(const stk::mesh::Entity node) const
{
  return get_nodal_values(myMetricField, mySpatialDim*mySpatialDim, node);
}

double NodalMetricTensorFieldMeshSizeMetric::compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const
{
  const std::vector<double> metric0 = get_metric(node0);
  const std::vector<double> metric1 = get_metric(node1);
  return metric::logarithmic_mean(metric::compute_relative_length_in_metric(mySpatialDim, metric0.data(), edge),
      metric::compute_relative_length_in_metric(mySpatialDim, metric1.data(), edge));
}

void NodalMetricTensorFieldMeshSizeMetric::update_size_at_new_node(const stk::mesh::Entity newNode, const stk::mesh::Entity * parentBegin, const stk::mesh::Entity * parentEnd) const
{
  std::vector<std::vector<double>> parentMetricValues;
  for (const stk::mesh::Entity * parent = parentBegin; parent != parentEnd; ++parent)
    parentMetricValues.push_back(get_metric(*parent));
  std::vector<const double *> parentMetrics;
  for (const auto & values : parentMetricValues)
    parentMetrics.push_back(values.data());

  const unsigned metricSize = mySpatialDim*mySpatialDim;
  STK_ThrowRequireMsg(myMetricField.length(newNode) == metricSize,
      "The field " << myMetricField.name() << " must have " << metricSize << " components at every new node.");
  const std::vector<double> mean = metric::compute_log_euclidean_mean_of_metrics(mySpatialDim, parentMetrics);
  auto newMetric = myMetricField.data<double, stk::mesh::ReadWrite>().entity_values(newNode);
  for (unsigned i=0; i<metricSize; ++i)
    newMetric(stk::mesh::ComponentIdx(i)) = mean[i];
}

void NodalMetricTensorFieldMeshSizeMetric::communicate_updated_sizes(const stk::mesh::BulkData & mesh) const
{
  stk::mesh::communicate_field_data(mesh, {&myMetricField.field()});
}

//--------------------------------------------------------------------------------

MetricTensorFunctionMeshSizeMetric::MetricTensorFunctionMeshSizeMetric(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const MetricFunction & metric)
: myMesh(mesh), mySpatialDim(mesh.mesh_meta_data().spatial_dimension()), myCoordsField(coordsField), myMetric(metric)
{
}

std::vector<double> MetricTensorFunctionMeshSizeMetric::get_metric(const stk::mesh::Entity node) const
{
  std::vector<double> metric(size_t(mySpatialDim) * mySpatialDim);
  myMetric(get_vector_field(myMesh, myCoordsField, node, mySpatialDim), metric.data());
  return metric;
}

double MetricTensorFunctionMeshSizeMetric::compute_relative_edge_length(const stk::mesh::Entity node0, const stk::mesh::Entity node1, const stk::math::Vector3d & edge) const
{
  const std::vector<double> metric0 = get_metric(node0);
  const std::vector<double> metric1 = get_metric(node1);
  return metric::logarithmic_mean(metric::compute_relative_length_in_metric(mySpatialDim, metric0.data(), edge),
      metric::compute_relative_length_in_metric(mySpatialDim, metric1.data(), edge));
}

} // namespace krino
