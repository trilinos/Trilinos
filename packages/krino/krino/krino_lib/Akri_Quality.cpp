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
double compute_mesh_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const QualityMetric &qualityMetric)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  std::vector<stk::math::Vector3d> nodeLocations;

  double quality = qualityMetric.get_best_value_for_metric();
  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector))
  {
    if (bucket->topology().base() == stk::topology::TETRAHEDRON_4 || bucket->topology().base() == stk::topology::TRIANGLE_3_2D)
    {
      for (auto && element : *bucket)
      {
        fill_element_node_coordinates(mesh, element, coordsField, nodeLocations);
        const double elementQuality = qualityMetric.get_element_quality_metric(nodeLocations);
        quality = std::min(quality, elementQuality);
      }
    }
  }

  const double localQuality = quality;
  stk::all_reduce_min(mesh.parallel(), &localQuality, &quality, 1);

  return quality;
}
}
