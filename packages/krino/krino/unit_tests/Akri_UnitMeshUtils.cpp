#include <Akri_UnitMeshUtils.hpp>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_math/StkVector.hpp>
#include <limits>

namespace krino {

stk::mesh::Entity find_local_node_closest_to_location(const stk::mesh::BulkData& mesh, const stk::mesh::Selector& nodeSelector, const FieldRef coordsField, const stk::math::Vector3d& location)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  double minSqrDistance = std::numeric_limits<double>::max();
  stk::mesh::Entity closest;

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, nodeSelector))
  {
    for(const auto node : *bucketPtr)
    {
      const stk::math::Vector3d nodeLocation = get_vector_field(mesh, coordsField, node, dim);
      const double nodeSqrDist = (nodeLocation-location).length_squared();
      if(nodeSqrDist < minSqrDistance)
      {
        minSqrDistance = nodeSqrDist;
        closest = node;
      }
    }
  }

  return closest;
}

}


