#include <Akri_AnalyticSurf.hpp>
#include <Akri_NodalSurfaceDistance.hpp>

#include <Akri_NodalBoundingBox.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_CSG_Surface.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>


namespace krino {

void compute_nodal_surface_distance(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    Surface & surface,
    const double time,
    const double narrowBandSize,
    const double farFieldValue)
{
  const double actualFarFieldValue = (narrowBandSize != 0 && farFieldValue == 0) ? narrowBandSize : farFieldValue;
  const unsigned nDim = mesh.mesh_meta_data().spatial_dimension();
  const BoundingBox nodeBbox = compute_nodal_bbox(mesh, mesh.mesh_meta_data().universal_part(), coordsField);
  surface.prepare_to_compute(time, nodeBbox, narrowBandSize);

  stk::mesh::BucketVector const& buckets = mesh.get_buckets( stk::topology::NODE_RANK, stk::mesh::selectField(distanceField) );

  auto distData = distanceField.data<double, stk::mesh::ReadWrite>();
  auto coordData = coordsField.data<double>();

  for ( auto && bucket_ptr : buckets )
  {
    const stk::mesh::Bucket & b = *bucket_ptr;
    auto dist = distData.bucket_values(b);

    for (stk::mesh::EntityIdx n : b.entities())
    {
      const stk::math::Vector3d x(coordData.entity_values(b[n]), nDim);

      dist(n) = surface.possibly_truncated_point_signed_distance(x, narrowBandSize, actualFarFieldValue);
    }
  }

  stk::mesh::communicate_field_data(mesh, {&distanceField.field()});
}

void compute_nodal_distance_from_spheres(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    const std::vector<std::pair<stk::math::Vector3d,double>> & sphereData,
    const int sign)
{
  std::vector<std::shared_ptr<Surface>> spheres;
  for (auto & sphere : sphereData)
    spheres.emplace_back(new Sphere(sphere.first, sphere.second, sign));
  CSGTree tree = CSGTree::unite(spheres);
  CSG_Surface initializationSurfaces(tree);

  if (sign < 0)
    initializationSurfaces.set_inside_to_be_positive();

  compute_nodal_surface_distance(mesh, coordsField, distanceField, initializationSurfaces);
}

void compute_nodal_distance_from_plane(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    const stk::math::Vector3d & normal,
    const double offset)
{
  Plane plane(normal, offset, 1.0);
  compute_nodal_surface_distance(mesh, coordsField, distanceField, plane);
}

}
