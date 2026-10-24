#include <Akri_BoundingBox.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

namespace krino {

BoundingBox
compute_nodal_bbox( const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector, const FieldRef coordsField )
{
  const int ndim = mesh.mesh_meta_data().spatial_dimension();
  BoundingBox bbox;

  auto coordsFieldData = coordsField.data<double>();
  stk::mesh::BucketVector const& buckets = mesh.get_buckets( stk::topology::NODE_RANK, selector);
  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;
    for (stk::mesh::EntityIdx i : b.entities())
    {
      stk::math::Vector3d x(coordsFieldData.entity_values(b[i]), ndim);
      bbox.accommodate( x );
    }
  }

  return bbox;
}

BoundingBox
compute_nodal_bbox( const stk::mesh::BulkData & mesh, const FieldRef coordsField, const std::vector<stk::mesh::Entity> & nodes )
{
  const int ndim = mesh.mesh_meta_data().spatial_dimension();
  BoundingBox bbox;

  auto coordsFieldData = coordsField.data<double>();
  for ( auto && node : nodes )
  {
    stk::math::Vector3d x(coordsFieldData.entity_values(node), ndim);
    bbox.accommodate( x );
  }

  return bbox;
}

}
