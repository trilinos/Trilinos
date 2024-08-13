// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <memory>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>

#include <Akri_MeshHelpers.hpp>

namespace krino {

void build_mesh(stk::mesh::BulkData & mesh, const std::vector<stk::mesh::EntityIdVector> & elem_nodes, const std::vector<int> & elem_procs, std::vector<stk::mesh::PartVector> & elem_parts)
{
  if (elem_nodes.empty()) return;

  const size_t numElem = elem_nodes.size();
  EXPECT_TRUE(elem_procs.size() == numElem);
  EXPECT_TRUE(elem_parts.size() == numElem);
  const size_t nodesPerElem = elem_nodes[0].size();
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  mesh.modification_begin();

  // Create nodes and elements
  std::set<stk::mesh::EntityId> nodes_I_declared;
  for (size_t e=0; e<numElem; ++e)
  {
    if (parallel_size == 1 || elem_procs[e] == parallel_rank)
    {
      for(size_t n=0; n<nodesPerElem; ++n)
      {
        stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];

        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
        if (!mesh.is_valid(node))
        {
          mesh.declare_node(nodeGlobalId);
          nodes_I_declared.insert(nodeGlobalId);
        }
      }
      stk::mesh::EntityId elemGlobalId = e+1;
      stk::mesh::declare_element(mesh, elem_parts[e], elemGlobalId, elem_nodes[e]);
    }
  }

  // declare node sharing
  if (parallel_size > 1)
  {
    for (size_t e=0; e<numElem; ++e)
    {
      if (elem_procs[e] != parallel_rank)
      {
        for(size_t n=0; n<nodesPerElem; ++n)
        {
          stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];
          if (nodes_I_declared.find(nodeGlobalId) != nodes_I_declared.end())
          {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
            mesh.add_node_sharing(node, elem_procs[e]);
          }
        }
      }
    }
  }

  mesh.modification_end();
}

namespace {
void test_and_cleanup_internal_side(stk::mesh::BulkData & mesh, const stk::mesh::Part & block_1, const stk::mesh::Part & block_2)
{
  EXPECT_TRUE(check_face_and_edge_relations(mesh));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));
  EXPECT_TRUE(check_shared_entity_nodes(mesh));

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_entities( mesh, mesh.mesh_meta_data().side_rank(), sides );

  EXPECT_TRUE(1 == sides.size());

  for (auto && side : sides)
  {
    EXPECT_TRUE(mesh.bucket(side).member(block_1));
    EXPECT_TRUE(mesh.bucket(side).member(block_2));
  }

  // cleanup
  mesh.modification_begin();
  for (auto && side : sides)
  {
    EXPECT_TRUE(disconnect_and_destroy_entity(mesh, side));
  }
  mesh.modification_end();
}
}

auto create_2D_mesh(const stk::ParallelMachine & pm)
{
  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(pm).set_spatial_dimension(2).create();
  return bulk;
}

TEST(MeshHelpers, DeclareElementSide)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 2) return;

  // This test will create a two element mesh (quad4 elements) on 1 or 2 processors.

  /*  Mesh
   *    4---5---6  P0 owns nodes 1,2,4,5; P, elem 1
   *    | 1 | 2 |  P1 : 3,6, elem 2
   *    1---2---3
   */


  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
  stk::mesh::Part& surface_1 = meta.declare_part_with_topology("surface_1", stk::topology::LINE_2);

  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 5, 4}, {2, 3, 6, 5} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1}, {&block_2}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  // All elements and nodes appear everywhere via aura
  stk::mesh::Entity element1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  stk::mesh::Entity element2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);
  stk::mesh::Entity node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  stk::mesh::Entity node5 = mesh.get_entity(stk::topology::NODE_RANK, 5);

  std::vector<stk::mesh::Entity> side_nodes;
  side_nodes.push_back(node2);
  side_nodes.push_back(node5);

  mesh.initialize_face_adjacent_element_graph();

  {
    // Case 1: Declare side on owning processor for both elements
    mesh.modification_begin();

    if (mesh.bucket(element1).owned())
    {
      const unsigned elem1_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element1, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element1, elem1_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }
    if (mesh.bucket(element2).owned())
    {
      const unsigned elem2_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element2, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element2, elem2_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }

    mesh.modification_end();

    test_and_cleanup_internal_side(mesh, block_1, block_2);
  }

  {
    // Case 2: Declare side of element1 on processor that owns element1
    mesh.modification_begin();

    if (mesh.bucket(element1).owned())
    {
      const unsigned elem1_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element1, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element1, elem1_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }

    mesh.modification_end();

    test_and_cleanup_internal_side(mesh, block_1, block_2);
  }

  {
    // Case 3: Declare side of element2 on processor that owns element2
    mesh.modification_begin();

    if (mesh.bucket(element2).owned())
    {
      const unsigned elem2_local_side_id = stk::mesh::get_entity_subcell_id(mesh, element2, meta.side_rank(), surface_1.topology(), side_nodes);
      mesh.declare_element_side(element2, elem2_local_side_id, stk::mesh::ConstPartVector{&surface_1});
    }

    mesh.modification_end();

    test_and_cleanup_internal_side(mesh, block_1, block_2);
  }
}

TEST(MeshHelpers, FullyCoincidentVolumeElements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // This test will create a two element mesh (quad4 elements), on more than 2 processors only
  // ranks 0 and 1 will have any elements. We test larger number of processors to ensure that
  // we get a parallel-consistent result to avoid potential parallel hangs in the full app.

  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& active_part = meta.declare_part("active");
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 3, 4}, {2, 3, 4, 1} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1}, {&block_1}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const bool ok = check_coincident_elements(mesh, active_part);
  EXPECT_FALSE(ok);
}

TEST(MeshHelpers, PartiallyCoincidentActiveVolumeElements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 2) return;

  // This test will create a two element mesh (quad4 elements) on 1 or 2 processors.

  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& active_part = meta.declare_part("active");
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 3, 4}, {1, 2, 5, 6} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1, &active_part}, {&block_1, &active_part}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const bool ok = check_coincident_elements(mesh, active_part);
  EXPECT_FALSE(ok);
}

TEST(MeshHelpers, NotCoincidentActiveDegenerateVolumeElements)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 2) return;

  // This test will create a two element mesh (quad4 elements) on 1 or 2 processors.

  auto meshPtr = create_2D_mesh(pm);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
  stk::mesh::Part& active_part = meta.declare_part("active");
  meta.commit();

  const std::vector<stk::mesh::EntityIdVector> elem_nodes{ {1, 2, 2, 3}, {3, 2, 2, 4} };
  const std::vector<int> elem_procs{0, 1};
  std::vector<stk::mesh::PartVector> elem_parts{{&block_1, &active_part}, {&block_1, &active_part}};

  build_mesh(mesh, elem_nodes, elem_procs, elem_parts);

  const bool ok = check_coincident_elements(mesh, active_part);
  EXPECT_TRUE(ok);
}

}
