#include <Akri_Quality.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_QualityMetric.hpp>

#include <vector>

namespace krino {

static bool is_supported_topology(stk::topology topology)
{
  return topology.base() == stk::topology::TETRAHEDRON_4 ||
      topology.base() == stk::topology::SHELL_TRIANGLE_3 ||
      topology.base() == stk::topology::TRIANGLE_3_2D;
}

double compute_mesh_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    const QualityMetric &qualityMetric)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> nodeLocations;

  double quality = qualityMetric.get_best_value_for_metric();
  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector))
  {
    if (is_supported_topology(bucket->topology()))
    {
      for (auto && element : *bucket)
      {
        fill_element_node_coordinates(mesh, element, coordsField, nodeLocations);
        const double elementQuality = qualityMetric.get_element_quality_metric(dim, nodeLocations);
        quality = std::min(quality, elementQuality);
      }
    }
  }

  const double localQuality = quality;
  stk::all_reduce_min(mesh.parallel(), &localQuality, &quality, 1);

  return quality;
}

double compute_quality_if_node_is_moved_terminating_early_if_below_threshold(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldRef coordsField,
    stk::mesh::Entity node,
    const stk::math::Vector3d & newLocation,
    const QualityMetric &qualityMetric,
    const double qualityThreshold)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  double minQualityAfterMove = qualityMetric.get_best_value_for_metric();
  std::vector<stk::math::Vector3d> nodeLocations;

  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (elementSelector(mesh.bucket(elem)))
    {
      nodeLocations.clear();
      for (auto elemNode : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
      {
        if (elemNode == node)
          nodeLocations.push_back(newLocation);
        else
          nodeLocations.emplace_back(field_data<double>(coordsField, elemNode), dim);
      }

      const double elemQualityAfterMove = qualityMetric.get_element_quality_metric(dim, nodeLocations);

      if (qualityMetric.is_first_quality_metric_better_than_second(minQualityAfterMove, elemQualityAfterMove))
      {
        minQualityAfterMove = elemQualityAfterMove;
        if (qualityMetric.is_first_quality_metric_better_than_second(qualityThreshold, minQualityAfterMove))
          return minQualityAfterMove;
      }
    }
  }
  return minQualityAfterMove;
}

template<typename ELEMCONTAINER>
double compute_minimum_element_quality(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const ELEMCONTAINER & elements,
    const QualityMetric &qualityMetric)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<stk::math::Vector3d> nodeLocations;

  double quality = qualityMetric.get_best_value_for_metric();
  for (auto & element : elements)
  {
    if (elementSelector(mesh.bucket(element)))
    {
      STK_ThrowAssert(is_supported_topology(mesh.bucket(element).topology()));
      fill_element_node_coordinates(mesh, element, coordsField, nodeLocations);
      const double elementQuality = qualityMetric.get_element_quality_metric(dim, nodeLocations);
      quality = std::min(quality, elementQuality);
    }
  }
  return quality;
}

// Explicit template instantiation
template double compute_minimum_element_quality(const stk::mesh::BulkData & mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Selector & elementSelector,  const StkMeshEntities & elements, const QualityMetric &qualityMetric);
template double compute_minimum_element_quality(const stk::mesh::BulkData & mesh, const CoordinatesFieldRef coordsField, const stk::mesh::Selector & elementSelector, const std::vector<stk::mesh::Entity> & elements, const QualityMetric &qualityMetric);

}
