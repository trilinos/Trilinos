// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_KinematicMeshUtils.hpp>

#include <Akri_KinematicUtils.hpp>
#include <Akri_MetricTensorUtils.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_SethHillConfig.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <functional>

namespace krino {

std::vector<double> compute_element_ref_sizes_from_reference(
    const stk::mesh::BulkData & mesh,
    const FieldRef referenceCoordsField,
    const std::vector<stk::mesh::Entity> & elems,
    const SmoothReference & reference)
{
  STK_ThrowRequireMsg(mesh.mesh_meta_data().spatial_dimension() == 3,
      "Element reference sizes from a reference configuration are only supported for tet4 elements in 3D.");

  std::vector<double> elemRefSizes;
  elemRefSizes.reserve(elems.size());
  for (auto elem : elems)
  {
    STK_ThrowRequireMsg(mesh.num_nodes(elem) == 4,
        "Element reference sizes from a reference configuration are only supported for tet4 elements in 3D.");
    const std::array<stk::math::Vector3d, 4> elemNodeCoords = get_tetrahedron_vector(mesh, referenceCoordsField, mesh.begin_nodes(elem));
    elemRefSizes.push_back(kinematicUtils::compute_ideal_tet4_h(elemNodeCoords, reference));
  }
  return elemRefSizes;
}

static void check_nodal_metric(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const unsigned spatialDim, const double * metric)
{
  STK_ThrowRequireMsg(metric::is_symmetric_positive_definite_metric(spatialDim, metric),
      "The metric at node " << mesh.identifier(node) << " is not symmetric positive definite.  The metric is "
      "a spatial_dimension*spatial_dimension matrix in row-major order; a non-symmetric metric may mean its "
      "components are in an unexpected order.");
}

static void check_supported_element(const unsigned spatialDim, const size_t numElemNodes)
{
  STK_ThrowRequireMsg((spatialDim == 2 && numElemNodes == 3) || (spatialDim == 3 && numElemNodes == 4),
      "Element metric mappings are only supported for tri3 elements in 2D and tet4 elements in 3D.");
}

// Metric mapping of each of numElems elements from the metrics of its nodes, which
// fill_elem_nodal_metrics(e, elemNodalMetrics) gives for element e.
static std::vector<double> compute_element_metric_mappings(
    const unsigned spatialDim,
    const size_t numElems,
    const std::function<void(size_t, std::vector<const double *> &)> & fill_elem_nodal_metrics)
{
  std::vector<const double *> elemNodalMetrics;
  std::vector<double> elemMetricMapData;
  elemMetricMapData.reserve(numElems*spatialDim*spatialDim);
  for (size_t e=0; e<numElems; ++e)
  {
    elemNodalMetrics.clear();
    fill_elem_nodal_metrics(e, elemNodalMetrics);
    const std::vector<double> metricMapping = metric::compute_metric_mapping_data_from_nodal_metrics(spatialDim, elemNodalMetrics);
    elemMetricMapData.insert(elemMetricMapData.end(), metricMapping.begin(), metricMapping.end());
  }
  return elemMetricMapData;
}

std::vector<double> compute_element_metric_mappings_from_nodal_metric_field(
    const stk::mesh::BulkData & mesh,
    const FieldRef nodalMetricField,
    const std::vector<stk::mesh::Entity> & elems)
{
  STK_ThrowRequireMsg(nodalMetricField.valid() && nodalMetricField.entity_rank() == stk::topology::NODE_RANK,
      "The metric field for an anisotropic kinematic mesh objective must be a nodal field.");

  const unsigned spatialDim = mesh.mesh_meta_data().spatial_dimension();
  const unsigned metricSize = spatialDim*spatialDim;
  const auto nodalMetricData = nodalMetricField.data<double>();
  std::vector<std::vector<double>> elemNodalMetricValues;

  return compute_element_metric_mappings(spatialDim, elems.size(),
      [&](const size_t e, std::vector<const double *> & elemNodalMetrics)
      {
        const stk::mesh::Entity elem = elems[e];
        check_supported_element(spatialDim, mesh.num_nodes(elem));
        elemNodalMetricValues.clear();
        for (auto node : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
        {
          STK_ThrowRequireMsg(nodalMetricField.length(node) == metricSize,
              "The metric field " << nodalMetricField.name() << " must have " << metricSize << " components at every node of the elements.");
          const auto nodalMetric = nodalMetricData.entity_values(node);
          std::vector<double> & metric = elemNodalMetricValues.emplace_back(metricSize);
          for (unsigned i=0; i<metricSize; ++i)
            metric[i] = nodalMetric(stk::mesh::ComponentIdx(i));
          check_nodal_metric(mesh, node, spatialDim, metric.data());
        }
        for (const auto & metric : elemNodalMetricValues)
          elemNodalMetrics.push_back(metric.data());
      });
}

std::vector<double> compute_element_metric_mappings_from_nodal_metrics(
    const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    const MeshElemNodeIndices & elemNodeIndices,
    const std::vector<double> & nodalMetrics)
{
  const unsigned spatialDim = mesh.mesh_meta_data().spatial_dimension();
  const unsigned metricSize = spatialDim*spatialDim;
  STK_ThrowRequireMsg(nodalMetrics.size() == nodes.size()*metricSize,
      "Expected " << nodes.size()*metricSize << " nodal metric values (" << metricSize << " for each of the "
      << nodes.size() << " nodes), but got " << nodalMetrics.size() << ".");

  for (size_t n=0; n<nodes.size(); ++n)
    check_nodal_metric(mesh, nodes[n], spatialDim, &nodalMetrics[n*metricSize]);

  return compute_element_metric_mappings(spatialDim, elemNodeIndices.num_elements(),
      [&](const size_t e, std::vector<const double *> & elemNodalMetrics)
      {
        const ElemNodeIndices elemNodes = elemNodeIndices.element_node_indices(e);
        check_supported_element(spatialDim, elemNodes.size());
        for (const unsigned nodeIndex : elemNodes)
          elemNodalMetrics.push_back(&nodalMetrics[nodeIndex*metricSize]);
      });
}

} // namespace krino
