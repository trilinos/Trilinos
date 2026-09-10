// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_Objective_Mesh.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <numeric>

namespace krino {

struct MeshSmoothingMetrics
{
  double worstMeanRatio;
  std::vector<double> elemQualities;
  unsigned invertedCount;
  double objectiveValue;
  double gradientNorm;
};

// Poor parallel performance due to parallel_vector_concat
std::vector<double> get_global_sorted_element_qualities(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & ownedElems,
    const QualityMetric & qualityMetric)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> nodeLocations;

  std::vector<double> ownedElemQualities;
  for (auto && element : ownedElems)
  {
    fill_element_node_coordinates(mesh, element, coordsField, nodeLocations);
    const double elementQuality = qualityMetric.get_element_quality_metric(dim, nodeLocations);
    ownedElemQualities.push_back(elementQuality);
  }

  std::vector<double> elemQualities;
  stk::parallel_vector_concat(mesh.parallel(), ownedElemQualities, elemQualities);
  std::sort(elemQualities.begin(), elemQualities.end());
  return elemQualities;
}

void fill_smoothing_metrics(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & ownedElems,
    const DistributedVector & nodeCoords,
    const MeshObjectiveBase & smoothingObj,
    MeshSmoothingMetrics & meshQualityMetrics)
{
  MeanRatioQualityMetric qualityMetric;
  meshQualityMetrics.elemQualities = get_global_sorted_element_qualities(mesh, coordsField, ownedElems, qualityMetric);
  meshQualityMetrics.worstMeanRatio = meshQualityMetrics.elemQualities.front();

  meshQualityMetrics.invertedCount = 0;
  for (double elemQual : meshQualityMetrics.elemQualities)
    if (elemQual <= 0.0)
      ++meshQualityMetrics.invertedCount;

  meshQualityMetrics.objectiveValue = smoothingObj.compute_value(nodeCoords);

  DistributedVector grad;
  smoothingObj.fill_gradient(nodeCoords, grad);
  meshQualityMetrics.gradientNorm = std::sqrt(Dot(grad, grad));
}

void print_smoothing_metrics(const MeshSmoothingMetrics & meshQualityMetrics)
{
  krinolog << " worst mean ratio = " << meshQualityMetrics.worstMeanRatio << "\n";
  krinolog << " objective value = " << meshQualityMetrics.objectiveValue << "\n";
  krinolog << " gradient norm = " << meshQualityMetrics.gradientNorm << "\n";

  const size_t n = meshQualityMetrics.elemQualities.size();
  const double meanQ =
      std::accumulate(meshQualityMetrics.elemQualities.begin(), meshQualityMetrics.elemQualities.end(), 0.0) / n;
  krinolog << " Per-element mean ratio statistics:\n"
           << "  count:     " << n << "\n"
           << "  inverted:  " << meshQualityMetrics.invertedCount << "\n"
           << "  min:       " << meshQualityMetrics.elemQualities.front() << "\n"
           << "  10%:       " << meshQualityMetrics.elemQualities[n / 10] << "\n"
           << "  25%:       " << meshQualityMetrics.elemQualities[n / 4] << "\n"
           << "  50%:       " << meshQualityMetrics.elemQualities[n / 2] << "\n"
           << "  75%:       " << meshQualityMetrics.elemQualities[(3 * n) / 4] << "\n"
           << "  90%:       " << meshQualityMetrics.elemQualities[(9 * n) / 10] << "\n"
           << "  max:       " << meshQualityMetrics.elemQualities.back() << "\n"
           << "  mean:      " << meanQ << "\n";
}

void compute_and_print_smoothing_metrics(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & ownedElems,
    const DistributedVector & nodeCoords,
    const MeshObjectiveBase & smoothingObj)
{
  MeshSmoothingMetrics meshQualityMetrics;
  fill_smoothing_metrics(mesh, coordsField, ownedElems, nodeCoords, smoothingObj, meshQualityMetrics);
  print_smoothing_metrics(meshQualityMetrics);
}

void smooth_mesh_by_optimization(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const OptimizerFunction & optimizer,
    const MeshObjectiveBase & smoothingObj,
    const bool doPrintDiagnostics)
{
  const auto & nodes = smoothingObj.get_nodes();
  const auto & ownedElems = smoothingObj.get_owned_elements();

  DistributedVector nodeCoords = gather_node_coordinates(mesh, coordsField, nodes);

  if (doPrintDiagnostics)
  {
    krinolog << "Before smoothing:\n";
    compute_and_print_smoothing_metrics(mesh, coordsField, ownedElems, nodeCoords, smoothingObj);
  }

  optimizer(smoothingObj, nodeCoords);
  set_node_coordinates(mesh, coordsField, nodes, nodeCoords);

  if (doPrintDiagnostics)
  {
    krinolog << "After smoothing:\n";
    compute_and_print_smoothing_metrics(mesh, coordsField, ownedElems, nodeCoords, smoothingObj);
  }
}

void smooth_mesh_by_optimization(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const OptimizerFunction & optimizer,
    const ElementObjective & elemSmoothingObj,
    const GradientFilter & gradientFilter,
    const SolutionProjector & solutionProjector,
    const bool doPrintDiagnostics)
{
  MeshObjective smoothingObj(mesh, coordsField, elemSelector, elemSmoothingObj);
  smoothingObj.set_gradient_filter(gradientFilter);
  smoothingObj.set_solution_projector(solutionProjector);

  smooth_mesh_by_optimization(mesh, coordsField, optimizer, smoothingObj, doPrintDiagnostics);
}

}
