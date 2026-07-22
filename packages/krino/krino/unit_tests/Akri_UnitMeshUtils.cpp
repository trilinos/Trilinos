#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_UnitMeshUtils.hpp>

#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_OutputUtils.hpp>
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

void populate_bounding_box_mesh_and_activate(BoundingBoxMesh & bboxMesh, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize)
{
  bboxMesh.set_domain(krino::BoundingBoxMesh::BoundingBoxType(minCorner, maxCorner), meshSize);
  bboxMesh.populate_mesh();
  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  activate_all_entities(mesh, AuxMetaData::get(bboxMesh.meta_data()).active_part());
}

void generate_and_write_bounding_box_mesh(const stk::topology elemTopology, const stk::math::Vector3d & minCorner, const stk::math::Vector3d & maxCorner, const double meshSize, const std::string & filename)
{
  BoundingBoxMesh bboxMesh(elemTopology, MPI_COMM_WORLD);
  populate_bounding_box_mesh_and_activate(bboxMesh, minCorner, maxCorner, meshSize);
  output_mesh_with_fields(bboxMesh.bulk_data(), AuxMetaData::get(bboxMesh.meta_data()).active_part(), filename, 1, 0.0);
}

}


