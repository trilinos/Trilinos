// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Compute_Surface_Distance.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_NodalBoundingBox.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino{

void
Compute_Surface_Distance::calculate(
  const stk::mesh::BulkData & mesh,
  const stk::diag::Timer &parent_timer,
  const stk::mesh::Field<double>& coordinates,
  const stk::mesh::Field<double>& distance,
  const stk::mesh::Selector & surface_selector,
  const double narrowBandSize,
  const double farFieldValue)
{
  calculate(mesh, parent_timer, coordinates, distance, mesh.mesh_meta_data().universal_part(), surface_selector, narrowBandSize, farFieldValue);
}

static void compute_distance_to_facets(const stk::mesh::BulkData & mesh,
    const FacetedSurfaceBase & facets,
    const stk::mesh::Field<double>& coordinates,
    const stk::mesh::Field<double>& distance,
    const stk::mesh::Selector & volume_selector,
    const double narrowBandSize,
    const double userFarFieldValue)
{
  const int spatial_dimension = mesh.mesh_meta_data().spatial_dimension();
  const double farFieldValue = (userFarFieldValue > 0.0) ? userFarFieldValue : narrowBandSize; // Use farFieldValue, if specified, otherwise, use narrowBandSize
  for ( auto && bucket : mesh.get_buckets(stk::topology::NODE_RANK, volume_selector & stk::mesh::selectField(distance)) )
  {
    const size_t length = bucket->size();
    double *dist = stk::mesh::field_data(distance, *bucket);
    double * coord = stk::mesh::field_data(coordinates, *bucket);

    for (size_t n = 0; n < length; ++n)
    {
      STK_ThrowAssert(&(dist[n]) != NULL);

      const stk::math::Vector3d xvec(coord+n*spatial_dimension, spatial_dimension);
      dist[n] = facets.point_unsigned_distance(xvec, narrowBandSize, farFieldValue);
    }
  }
}

void print_facet_info(const FacetedSurfaceBase & facets, stk::ParallelMachine parallel)
{
  constexpr int vec_size = 3;
  std::array <unsigned,vec_size> local_sizes, global_min, global_max;
  local_sizes[0] = facets.storage_size();
  local_sizes[1] = facets.size();
  local_sizes[2] = facets.size()+facets.nonlocal_size();

  stk::all_reduce_min( parallel, local_sizes.data(), global_min.data(), vec_size );
  stk::all_reduce_max( parallel, local_sizes.data(), global_max.data(), vec_size );

  krinolog << "Compute Surface Distance: "<< stk::diag::dendl;
  krinolog << "   Local facet count:  min=" << global_min[1] << ", max=" << global_max[1] << stk::diag::dendl;
  krinolog << "   Total facet count:  min=" << global_min[2] << ", max=" << global_max[2] << stk::diag::dendl;
  krinolog << "   Memory usage (mb):  min=" << global_min[0]/(1024.0*1024.0) << ", max=" << global_max[0]/(1024.0*1024.0) << stk::diag::dendl;
}

void
Compute_Surface_Distance::calculate(
  const stk::mesh::BulkData & mesh,
  const stk::diag::Timer &parent_timer,
  const stk::mesh::Field<double>& coordinates,
  const stk::mesh::Field<double>& distance,
  const stk::mesh::Selector & volume_selector,
  const stk::mesh::Selector & surfaceSelector,
  const double narrowBandSize,
  const double farFieldValue)
{ /* %TRACE[ON]% */ Trace trace__("krino::Compute_Surface_Distance::compute_surface_distance(void)"); /* %TRACE% */

  stk::diag::Timer timer( "Compute Surface Distance", parent_timer );
  stk::diag::TimeBlock timer_(timer);

  std::unique_ptr<FacetedSurfaceBase> facets = build_mesh_surface(mesh.mesh_meta_data(), coordinates, surfaceSelector, +1);

  FieldRef coordsField(coordinates);
  const BoundingBox nodeBbox = compute_nodal_bbox(mesh, volume_selector & stk::mesh::selectField(distance), coordsField);
  facets->prepare_to_compute(0.0, nodeBbox, narrowBandSize); // Setup including communication of facets that are within this processors narrow band

  print_facet_info(*facets, mesh.parallel());

  compute_distance_to_facets(mesh, *facets, coordinates, distance, volume_selector, narrowBandSize, farFieldValue);
}

} // namespace krino
