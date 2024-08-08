#include <Akri_AnalyticSurf.hpp>
#include <Akri_NodalSurfaceDistance.hpp>

#include <Akri_Composite_Surface.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>


namespace krino {

void compute_nodal_surface_distance(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef distanceField, Composite_Surface & surfaces, const double time, const double narrowBandSize)
{
  const unsigned nDim = mesh.mesh_meta_data().spatial_dimension();
  const BoundingBox nodeBbox = compute_nodal_bbox(mesh, mesh.mesh_meta_data().universal_part(), coordsField);
  surfaces.prepare_to_compute(time, nodeBbox, narrowBandSize);

  stk::mesh::BucketVector const& buckets = mesh.get_buckets( stk::topology::NODE_RANK, stk::mesh::selectField(distanceField) );

  for ( auto && bucket_ptr : buckets )
  {
    const stk::mesh::Bucket & b = *bucket_ptr;
    const int length = b.size();
    double *dist = field_data<double>(distanceField, b);
    double * coord = field_data<double>(coordsField, b);

    for (int n = 0; n < length; ++n)
    {
      STK_ThrowAssert(&(dist[n]) != NULL);

      const stk::math::Vector3d x(&coord[nDim*n], nDim);

      dist[n] = surfaces.point_signed_distance_with_narrow_band(x, narrowBandSize);
    }
  }

  stk::mesh::communicate_field_data(mesh, {&distanceField.field()});
}

void compute_nodal_distance_from_spheres(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    const std::vector<std::pair<stk::math::Vector3d,double>> & spheres,
    const int sign)
{
  Composite_Surface initializationSurfaces("initialization surfaces");
  if (sign < 0)
    initializationSurfaces.set_composition_method(Composite_Surface::MAXIMUM_SIGNED_DISTANCE);
  for (auto & sphere : spheres)
    initializationSurfaces.add(new Sphere(sphere.first, sphere.second, sign));
  compute_nodal_surface_distance(mesh, coordsField, distanceField, initializationSurfaces);
}

void compute_nodal_distance_from_plane(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef distanceField,
    const stk::math::Vector3d & normal,
    const double offset)
{
  Composite_Surface initializationSurfaces("initialization surfaces");
  initializationSurfaces.add(new Plane(normal, offset, 1.0));
  compute_nodal_surface_distance(mesh, coordsField, distanceField, initializationSurfaces);
}

}
