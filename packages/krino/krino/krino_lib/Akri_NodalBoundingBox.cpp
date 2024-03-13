#include <Akri_BoundingBox.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

BoundingBox
compute_nodal_bbox( const stk::mesh::BulkData & mesh, const stk::mesh::Selector & selector, const FieldRef coordsField )
{
  const int ndim = mesh.mesh_meta_data().spatial_dimension();
  BoundingBox bbox;

  stk::mesh::BucketVector const& buckets = mesh.get_buckets( stk::topology::NODE_RANK, selector);
  for ( auto && bucket : buckets )
  {
    const stk::mesh::Bucket & b = *bucket;

    const size_t length = b.size();

    double *coordsData = field_data<double>(coordsField, b);

    for (size_t i = 0; i < length; ++i)
    {
      stk::math::Vector3d x(coordsData + i*ndim, ndim);
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

  for ( auto && node : nodes )
  {
    stk::math::Vector3d x(field_data<double>(coordsField, node), ndim);
    bbox.accommodate( x );
  }

  return bbox;
}

}
