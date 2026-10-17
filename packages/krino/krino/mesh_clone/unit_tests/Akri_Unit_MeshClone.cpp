// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>   // for declare_element
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_io/IossBridge.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshClone.hpp>
#include <Akri_MeshHelpers.hpp>

namespace krino {

namespace {

// A shell block whose triangle lies exactly on a solid tet face ("coincident
// shell") produces two distinct side-rank faces with the *identical* node set.
// When that interface is also split across a processor boundary, the mesh has
// two shared faces with the same nodes straddling ranks. STK resolves shared
// side-rank sharing by node set, so without id-based disambiguation it pairs
// the wrong faces across processors and BulkData::modification_end() dies in
// internal_resolve_shared_membership (FACE_RANK key mismatch under asserts, or
// a segfault in BucketRepository::get_partition in a release build).
//
// Build two tets sharing the interface face {1,2,3}, plus a coincident
// SHELL_TRI_3 element on those same three nodes. P0 owns the left tet (and the
// shell); P1 owns the right tet. Nodes 1,2,3 are shared across the P0/P1
// boundary, so the shell face and the solid interface face are two coincident
// shared faces straddling the parallel split. Only ranks 0 and 1 hold entities;
// any additional ranks are empty but still participate in modification_end, so
// the reproducer works on any run with >= 2 processors.
void build_coincident_shell_on_split_interface(stk::mesh::BulkData & mesh,
    stk::mesh::Part & solidL_part,
    stk::mesh::Part & solidR_part,
    stk::mesh::Part & shell_part,
    FieldRef coords_field)
{
  const int parallel_rank = mesh.parallel_rank();

  const std::vector<stk::mesh::EntityId> tetL_nodes = {1, 3, 2, 4};
  const std::vector<stk::mesh::EntityId> tetR_nodes = {1, 2, 3, 5};
  const std::vector<stk::mesh::EntityId> shell_nodes = {1, 2, 3};

  stk::mesh::PartVector solidL_parts{&solidL_part};
  stk::mesh::PartVector solidR_parts{&solidR_part};
  stk::mesh::PartVector shell_parts{&shell_part};

  mesh.modification_begin();
  {
    if (parallel_rank == 0)
    {
      stk::mesh::declare_element(mesh, solidL_parts, 1, tetL_nodes);
      stk::mesh::declare_element(mesh, shell_parts, 3, shell_nodes);
    }
    if (parallel_rank == 1)
    {
      stk::mesh::declare_element(mesh, solidR_parts, 2, tetR_nodes);
    }

    // Nodes 1,2,3 are shared between the two ranks.
    const int opp_rank = parallel_rank == 0 ? 1 : 0;
    for (stk::mesh::EntityId i = 1; i <= 3; ++i)
    {
      const auto node = mesh.get_entity(stk::topology::NODE_RANK, i);
      if (mesh.is_valid(node))
        mesh.add_node_sharing(node, opp_rank);
    }
  }
  mesh.modification_end();

  // Set coordinates so the geometry is well-defined.
  auto set_coords = [&](stk::mesh::EntityId id, double x, double y, double z)
  {
    const auto node = mesh.get_entity(stk::topology::NODE_RANK, id);
    if (mesh.is_valid(node))
    {
      double * c = field_data<double>(coords_field.field(), node);
      c[0] = x; c[1] = y; c[2] = z;
    }
  };
  set_coords(1, 0., 0., 0.);
  set_coords(2, 0., 0., 1.);
  set_coords(3, 0., 1., 0.);
  set_coords(4, -1., 0., 0.);
  set_coords(5, 1., 0., 0.);

  // Materialize the side-rank faces: the interior block-boundary face on the
  // solid interface AND the coincident shell face, both on nodes {1,2,3}.
  auto & aux_meta = AuxMetaData::get(mesh.mesh_meta_data());
  stk::mesh::create_all_sides(mesh,
      mesh.mesh_meta_data().universal_part(),
      stk::mesh::PartVector{&aux_meta.block_boundary_part()},
      false);
}

} // namespace

TEST(MeshCloneCoincidentShell, StashWithCoincidentShellOnParallelSplitInterface)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  // Inherently a parallel bug: needs the coincident interface split across a
  // processor boundary. The mesh is placed on ranks 0 and 1, so any run with at
  // least 2 processors exercises it; a serial run cannot reproduce it.
  if (parallel_size < 2) return;

  std::shared_ptr<stk::mesh::MetaData> meta =
      stk::mesh::MeshBuilder(pm).set_spatial_dimension(3).create_meta_data();
  std::unique_ptr<stk::mesh::BulkData> bulkPtr =
      stk::mesh::MeshBuilder(pm).set_spatial_dimension(3).create(meta);
  stk::mesh::BulkData & mesh = *bulkPtr;

  AuxMetaData::create(*meta);
  auto & aux_meta = AuxMetaData::get(*meta);

  FieldRef coords_field = aux_meta.register_field("coordinates",
      FieldType::VECTOR_3D, stk::topology::NODE_RANK, 1u, 1u, meta->universal_part());
  meta->set_coordinate_field(&coords_field.field());

  auto & solidL_part = meta->declare_part_with_topology("block_solidL", stk::topology::TETRAHEDRON_4);
  auto & solidR_part = meta->declare_part_with_topology("block_solidR", stk::topology::TETRAHEDRON_4);
  auto & shell_part  = meta->declare_part_with_topology("block_leads", stk::topology::SHELL_TRIANGLE_3_ALL_FACE_SIDES);
  stk::io::put_io_part_attribute(solidL_part);
  stk::io::put_io_part_attribute(solidR_part);
  stk::io::put_io_part_attribute(shell_part);

  meta->commit();

  // Mirror the STK IO reader: the source mesh resolves coincident shared faces
  // by entity id, so the source mesh itself is valid. The clone is what used to
  // omit this and crash.
  mesh.set_use_entity_ids_for_resolving_sharing(true);

  build_coincident_shell_on_split_interface(mesh, solidL_part, solidR_part, shell_part, coords_field);

  // Before the fix in Akri_MeshClone.cpp this throws/aborts inside the clone's
  // modification_end (FACE_RANK key mismatch); after the fix it succeeds.
  EXPECT_NO_THROW(MeshClone::stash_or_restore_mesh(mesh, 0));
}

}
