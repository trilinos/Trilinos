// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Phase_Support::get_negative_levelset_interface_selector of a level set selects the interface sides on
// the boundary of its domain, touching the inside, whatever signs the phases are given for the other
// level sets: the same decomposition gives the same interface for complete phase tags (one level set
// per interface), one level set per phase, and bounding-surface phases (each phase negative for its own
// level set, with an elsewhere phase).

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AllReduce.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshDiagnostics.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_UnitMeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino
{

namespace
{

enum class InterfaceSelectorPolicy { PER_INTERFACE, PER_PHASE, PER_BOUNDING_SURFACE };

// The phases of two level sets that meet at the plane x = 0.13 below y = 0.37, with the elsewhere (or
// background) phase above it, on a tet4 mesh, cut only: the areas of the negative level set interface
// selectors of the two level sets.
std::array<double, 2> negative_level_set_interface_areas(const InterfaceSelectorPolicy policy)
{
  constexpr double xSplit = 0.13;
  constexpr double yTop = 0.37;
  std::vector<std::function<double(const stk::math::Vector3d &)>> levelSets{
      [](const stk::math::Vector3d & x) { return std::max(x[0] - xSplit, x[1] - yTop); },
      [](const stk::math::Vector3d & x) { return std::max(xSplit - x[0], x[1] - yTop); }};
  if (policy == InterfaceSelectorPolicy::PER_PHASE) levelSets.push_back([](const stk::math::Vector3d &) { return 0.; });

  BoundingBoxMesh bboxMesh(stk::topology::TETRAHEDRON_4, MPI_COMM_WORLD);
  stk::mesh::MetaData & meta = bboxMesh.meta_data();
  std::vector<LS_Field> lsFields;
  if (policy == InterfaceSelectorPolicy::PER_INTERFACE)
    lsFields = LSPerInterfacePolicy::setup_levelsets_on_all_blocks(meta, levelSets.size());
  else if (policy == InterfaceSelectorPolicy::PER_PHASE)
    lsFields = LSPerPhasePolicy::setup_levelsets_on_all_blocks(meta, levelSets.size());
  else
    lsFields = LSPerBoundingSurfacePolicy::setup_levelsets_on_all_blocks(meta, levelSets.size());

  CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
  cdfemSupport.set_coords_field(meta.coordinate_field());
  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_NODE);
  cdfemSupport.register_cdfem_snap_displacements_field();
  cdfemSupport.register_parent_node_ids_field();
  cdfemSupport.finalize_fields();
  bboxMesh.set_mesh_structure_type(FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  populate_bounding_box_mesh_and_activate(bboxMesh, stk::math::Vector3d(-1., -1., -1.), stk::math::Vector3d(1., 1., 1.), 0.5);

  stk::mesh::BulkData & mesh = bboxMesh.bulk_data();
  const FieldRef coordsField(meta.coordinate_field());
  for (const stk::mesh::Bucket * bucket : mesh.buckets(stk::topology::NODE_RANK))
    for (const stk::mesh::Entity node : *bucket)
      for (size_t i = 0; i < levelSets.size(); ++i)
        get_scalar_field(mesh, lsFields[i].isovar, node) = levelSets[i](get_vector_field(mesh, coordsField, node, 3));

  const AuxMetaData & auxMeta = AuxMetaData::get(meta);
  const Phase_Support & phaseSupport = Phase_Support::get(meta);
  const std::unique_ptr<InterfaceGeometry> geometry =
      create_levelset_geometry(3, auxMeta.active_part(), cdfemSupport, phaseSupport, lsFields);
  CDMesh::decompose_mesh(mesh, *geometry);

  const stk::mesh::Selector ownedActive = meta.locally_owned_part() & auxMeta.active_part();
  auto global_area = [&](const stk::mesh::Selector & selector)
  {
    double area = compute_volume_or_surface_area(mesh, meta.side_rank(), selector & ownedActive);
    all_reduce_sum(mesh.parallel(), area);
    return area;
  };

  std::array<double, 2> areas;
  for (size_t i = 0; i < 2; ++i)
    areas[i] = global_area(phaseSupport.get_negative_levelset_interface_selector(lsFields[i].identifier));
  return areas;
}

} // namespace

TEST(NegativeLevelSetInterfaceSelector, phasesOfTwoLevelSetsThatMeet_sameInterfaceForEveryPhaseSetup)
{
  const std::array<double, 2> perInterface = negative_level_set_interface_areas(InterfaceSelectorPolicy::PER_INTERFACE);
  const std::array<double, 2> perPhase = negative_level_set_interface_areas(InterfaceSelectorPolicy::PER_PHASE);
  const std::array<double, 2> perBoundingSurface = negative_level_set_interface_areas(InterfaceSelectorPolicy::PER_BOUNDING_SURFACE);
  for (size_t i = 0; i < 2; ++i)
  {
    EXPECT_LT(4., perInterface[i]) << "level set " << i; // the faces with the other phase and the elsewhere phase
    EXPECT_NEAR(perInterface[i], perPhase[i], 1.e-12 * perInterface[i]) << "level set " << i;
    EXPECT_NEAR(perInterface[i], perBoundingSurface[i], 1.e-12 * perInterface[i]) << "level set " << i;
  }
}

} // namespace krino
