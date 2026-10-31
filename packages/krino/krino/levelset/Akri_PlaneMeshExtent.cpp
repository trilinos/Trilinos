// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_PlaneMeshExtent.hpp>

#include <Akri_AllReduce.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingSurface.hpp>
#include <Akri_CSG_Surface.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_IC_Alg.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Surface_Manager.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>

namespace krino {

std::array<double,2> compute_min_max_of_nodes_along_direction(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::math::Vector3d & direction)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  double localMin = std::numeric_limits<double>::max();
  double localMax = std::numeric_limits<double>::lowest();

  const stk::mesh::Selector ownedSelector = mesh.mesh_meta_data().locally_owned_part();
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, ownedSelector))
  {
    for (auto node : *bucketPtr)
    {
      const double dot = Dot(direction, get_vector_field(mesh, coordsField, node, dim));
      localMin = std::min(localMin, dot);
      localMax = std::max(localMax, dot);
    }
  }

  const std::array<double,1> globalMin = all_reduce_min_array(mesh.parallel(), std::array<double,1>{localMin});
  const std::array<double,1> globalMax = all_reduce_max_array(mesh.parallel(), std::array<double,1>{localMax});
  return {globalMin[0], globalMax[0]};
}

namespace {

void append_planes_needing_offset_from_extent(Surface & surface, std::vector<Plane*> & planes)
{
  if (CSG_Surface * csgSurface = dynamic_cast<CSG_Surface *>(&surface))
  {
    for (auto && subsurface : csgSurface->get_subsurfaces())
      append_planes_needing_offset_from_extent(*subsurface, planes);
    return;
  }

  Plane * plane = dynamic_cast<Plane *>(&surface);
  if (plane && plane->needs_offset_from_extent())
    planes.push_back(plane);
}

bool have_same_components(const stk::math::Vector3d & a, const stk::math::Vector3d & b)
{
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

struct PlanesWithNormal
{
  stk::math::Vector3d normal;
  std::vector<Plane*> planes;
};

// The planes grouped by (unit) normal, in the order each normal is first found, which is the same on
// all processors.
std::vector<PlanesWithNormal> group_planes_by_normal(const std::vector<Plane*> & planes)
{
  std::vector<PlanesWithNormal> groups;
  for (Plane * plane : planes)
  {
    auto iter = std::find_if(groups.begin(), groups.end(),
        [&](const PlanesWithNormal & group) { return have_same_components(group.normal, plane->normal()); });
    if (iter == groups.end())
      groups.push_back(PlanesWithNormal{plane->normal(), {plane}});
    else
      iter->planes.push_back(plane);
  }
  return groups;
}

void log_offset_from_extent(const Plane & plane, const double minDot, const double maxDot)
{
  std::ostringstream os;
  os << std::setprecision(std::numeric_limits<double>::max_digits10)
     << "Plane with normal " << plane.normal() << " at fraction " << plane.fraction_of_extent()
     << " of the mesh extent [" << minDot << ", " << maxDot << "] along its normal has offset " << plane.offset();
  krinolog << os.str() << stk::diag::dendl;
}

}

void set_plane_offsets_from_mesh_extent(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const std::vector<Surface*> & surfaces)
{
  std::vector<Plane*> planes;
  for (Surface * surface : surfaces)
    append_planes_needing_offset_from_extent(*surface, planes);
  if (planes.empty())
    return;

  for (auto && group : group_planes_by_normal(planes))
  {
    const auto [minDot, maxDot] = compute_min_max_of_nodes_along_direction(mesh, coordsField, group.normal);
    for (Plane * plane : group.planes)
    {
      plane->set_offset_from_extent(minDot, maxDot);
      log_offset_from_extent(*plane, minDot, maxDot);
    }
  }
}

void set_plane_offsets_from_mesh_extent(const stk::mesh::BulkData & mesh)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);

  std::vector<Surface*> surfaces;
  for (auto && boundingSurface : surfaceManager.get_bounding_surfaces())
    surfaces.push_back(&boundingSurface->surface());
  for (auto && ls : surfaceManager.get_levelsets())
    if (ls->has_IC_surfaces())
      for (auto && subsurface : ls->get_IC_alg().get_initialization_surface().get_subsurfaces())
        surfaces.push_back(subsurface.get());

  set_plane_offsets_from_mesh_extent(mesh, AuxMetaData::get(meta).get_current_coordinates(), surfaces);
}

}
