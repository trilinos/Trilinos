// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include "mpi.h"                        // for MPI_COMM_SELF, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace stk { namespace mesh { class Part; } }

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;
using stk::unit_test_util::build_mesh;

namespace {

// Set up a very simple mesh with one element, one side, one node:
// rels: E->(S1,S2)->N
// parts: E in element_rank_part, unranked_part
//        S1, S2 in side_rank_part
//        element_ranked_part subset of unranked_superset_part
// modification cycle is left uncompleted

class UnitTestInducedPart2D : public stk::unit_test_util::MeshFixture
{
protected:
  UnitTestInducedPart2D()
    : MeshFixture(2) {}

  void setup_serial_mesh()
  {
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    unranked_part = &get_meta().declare_part("unranked_part");
    element_rank_part = &get_meta().declare_part_with_topology("element_rank_part", stk::topology::TRI_3_2D);
    element_rank_superset_part = &get_meta().declare_part("element_rank_superset_part", stk::topology::ELEMENT_RANK);
    side_rank_part = &get_meta().declare_part_with_topology("side_rank_part", stk::topology::LINE_2);
    unranked_superset_part = &get_meta().declare_part("unranked_superset_part");
    get_meta().declare_part_subset(*unranked_superset_part, *element_rank_part);
    get_meta().declare_part_subset(*element_rank_superset_part, *element_rank_part);
    get_meta().commit();

    get_bulk().modification_begin();

    stk::mesh::PartVector parts;
    parts.push_back(unranked_part);
    parts.push_back(element_rank_part);
    elem = get_bulk().declare_element(1 /*id*/, parts);

    parts.clear();
    node  = get_bulk().declare_node(1 /*id*/, parts);
    node2 = get_bulk().declare_node(2 /*id*/, parts);
    node3 = get_bulk().declare_node(3 /*id*/, parts);

    get_bulk().declare_relation(elem, node,   0 /*rel id*/);
    get_bulk().declare_relation(elem, node2,  1 /*rel id*/);
    get_bulk().declare_relation(elem, node3,  2 /*rel id*/);

    parts.clear();
    parts.push_back(side_rank_part);
    side1 = get_bulk().declare_element_side(elem, 0, parts);
    side2 = get_bulk().declare_element_side(elem, 2, parts);

  }

  stk::mesh::Part * unranked_part;
  stk::mesh::Part * element_rank_part;
  stk::mesh::Part * element_rank_superset_part;
  stk::mesh::Part * side_rank_part;
  stk::mesh::Part * unranked_superset_part;
  Entity elem;
  Entity side1;
  Entity side2;
  Entity node;
  Entity node2;
  Entity node3;

};

class UnitTestInducedPart3D : public stk::unit_test_util::MeshFixture
{
protected:
  UnitTestInducedPart3D()
    : MeshFixture(3) {}

  void setup_serial_mesh()
  {
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    unranked_part = &get_meta().declare_part("unranked_part");
    element_rank_part = &get_meta().declare_part_with_topology("element_rank_part", stk::topology::TET_4);
    element_rank_superset_part = &get_meta().declare_part("element_rank_superset_part", stk::topology::ELEMENT_RANK);
    side_rank_part = &get_meta().declare_part_with_topology("side_rank_part", stk::topology::TRI_3);
    edge_rank_part = &get_meta().declare_part_with_topology("edge_rank_part", stk::topology::LINE_2);
    unranked_superset_part = &get_meta().declare_part("unranked_superset_part");
    get_meta().declare_part_subset(*unranked_superset_part, *element_rank_part);
    get_meta().declare_part_subset(*element_rank_superset_part, *element_rank_part);
    get_meta().commit();

    get_bulk().modification_begin();

    stk::mesh::PartVector parts;
    parts.push_back(unranked_part);
    parts.push_back(element_rank_part);
    elem = get_bulk().declare_element(1 /*id*/, parts);

    parts.clear();
    parts.push_back(edge_rank_part);
    edge = get_bulk().declare_edge(1 /*id*/, parts);

    parts.clear();
    node                    = get_bulk().declare_node(1 /*id*/, parts);
    stk::mesh::Entity node1 = get_bulk().declare_node(2 /*id*/, parts);
    stk::mesh::Entity node2 = get_bulk().declare_node(3 /*id*/, parts);
    stk::mesh::Entity node3 = get_bulk().declare_node(4 /*id*/, parts);

    get_bulk().declare_relation(elem, node,  0 /*rel id*/);
    get_bulk().declare_relation(elem, node1, 1 /*rel id*/);
    get_bulk().declare_relation(elem, node2, 2 /*rel id*/);
    get_bulk().declare_relation(elem, node3, 3 /*rel id*/);

    get_bulk().declare_relation(edge, node,  0 /*rel id*/);
    get_bulk().declare_relation(edge, node1, 1 /*rel id*/);

    parts.clear();
    parts.push_back(side_rank_part);
    side = get_bulk().declare_element_side(elem, 0, parts);
    get_bulk().declare_relation(side, edge, 0);
    get_bulk().modification_end();
  }

  stk::mesh::Part * unranked_part;
  stk::mesh::Part * element_rank_part;
  stk::mesh::Part * element_rank_superset_part;
  stk::mesh::Part * side_rank_part;
  stk::mesh::Part * edge_rank_part;
  stk::mesh::Part * unranked_superset_part;
  Entity elem;
  Entity side;
  Entity edge;
  Entity node;
};

TEST_F( UnitTestInducedPart2D , verifyBasicInducedPart )
{
  if (stk::parallel_machine_size(get_comm()) == 1) {
    setup_serial_mesh();

    // Check that directly-induced parts are induced upon relation creation
    // before modification end.
    EXPECT_TRUE(get_bulk().bucket(node).member( *side_rank_part));
    EXPECT_TRUE(get_bulk().bucket(node).member( *element_rank_part));
    EXPECT_TRUE(get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_FALSE(get_bulk().bucket(node).member( *unranked_superset_part));
    EXPECT_TRUE(get_bulk().bucket(side1).member( *element_rank_part));
    EXPECT_TRUE(get_bulk().bucket(side1).member( *element_rank_superset_part));
    EXPECT_FALSE(get_bulk().bucket(side1).member( *unranked_superset_part));
    EXPECT_TRUE(get_bulk().bucket(side2).member( *element_rank_part));
    EXPECT_FALSE(get_bulk().bucket(side2).member( *unranked_superset_part));
    EXPECT_TRUE(get_bulk().bucket(side2).member( *element_rank_superset_part));

    get_bulk().modification_end();

    // Modification-end should not have changed induced parts
    EXPECT_TRUE(get_bulk().bucket(node).member( *side_rank_part));
    EXPECT_TRUE(get_bulk().bucket(node).member( *element_rank_part));
    EXPECT_TRUE(get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_FALSE(get_bulk().bucket(node).member( *unranked_superset_part));
    EXPECT_TRUE(get_bulk().bucket(side1).member( *element_rank_part));
    EXPECT_TRUE(get_bulk().bucket(side1).member( *element_rank_superset_part));
    EXPECT_FALSE(get_bulk().bucket(side1).member( *unranked_superset_part));
    EXPECT_TRUE(get_bulk().bucket(side2).member( *element_rank_part));
    EXPECT_FALSE(get_bulk().bucket(side2).member( *unranked_superset_part));
    EXPECT_TRUE(get_bulk().bucket(side2).member( *element_rank_superset_part));
  }
}

TEST_F( UnitTestInducedPart2D, verifyInducedPartCorrectnessWhenRelationsRemoved )
{
  if (stk::parallel_machine_size(get_comm()) == 1) {
    setup_serial_mesh();

    EXPECT_TRUE( get_bulk().bucket(side1).member( *element_rank_part));
    EXPECT_TRUE( get_bulk().bucket(side1).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(side1).member( *unranked_superset_part));

    EXPECT_TRUE( get_bulk().bucket(side2).member( *element_rank_part));
    EXPECT_TRUE( get_bulk().bucket(side2).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(side2).member( *unranked_superset_part));

    EXPECT_TRUE( get_bulk().bucket(node).member( *element_rank_part));
    EXPECT_TRUE( get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(node).member( *unranked_superset_part));

    // Destroy one relation from element to a side, confirm that the side that lost
    // the relation no longer has the element part.
    get_bulk().destroy_relation(elem, side1, 0 /*rel id*/);
    EXPECT_TRUE(!get_bulk().bucket(side1).member( *element_rank_part));
    EXPECT_TRUE(!get_bulk().bucket(side1).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(side1).member( *unranked_superset_part));

    get_bulk().destroy_relation(elem, side2, 2 /*rel id*/);
    EXPECT_TRUE(!get_bulk().bucket(side2).member( *element_rank_part));
    EXPECT_TRUE(!get_bulk().bucket(side2).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(side2).member( *unranked_superset_part));

    // Destroy one of the relations from side to node. Confirm that node still has
    // side part due to its remaining relation.
    get_bulk().destroy_relation(side1, node, 0 /*rel id*/);
    EXPECT_TRUE( get_bulk().bucket(node).member( *side_rank_part));
    EXPECT_TRUE( get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(node).member( *unranked_superset_part));

    // Destroy the other relations from side to node. Confirm that node no longer
    // has any induced parts from the side.
    get_bulk().destroy_relation(side2, node, 1 /*rel id*/);
    EXPECT_TRUE(!get_bulk().bucket(node).member( *side_rank_part));
    EXPECT_TRUE( get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(node).member( *unranked_superset_part));

    // Destroy relation from element to node.  Confirm that the node no longer has
    // any induced parts.
    get_bulk().destroy_relation(elem,node,0 /*rel id*/);
    EXPECT_TRUE(!get_bulk().bucket(node).member( *side_rank_part));
    EXPECT_TRUE(!get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(node).member( *unranked_superset_part));
  }
}


TEST_F( UnitTestInducedPart3D, verifyNotTransitive )
{
  if (stk::parallel_machine_size(get_comm()) == 1) {
    setup_serial_mesh();

    EXPECT_TRUE( get_bulk().bucket(elem).member( *element_rank_part));
    EXPECT_TRUE( get_bulk().bucket(elem).member( *element_rank_superset_part));
    EXPECT_TRUE( get_bulk().bucket(elem).member( *unranked_superset_part));

    EXPECT_TRUE( get_bulk().bucket(side).member( *side_rank_part ));
    EXPECT_TRUE( get_bulk().bucket(side).member( *element_rank_part));
    EXPECT_TRUE( get_bulk().bucket(side).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(side).member( *unranked_superset_part));

    EXPECT_TRUE( get_bulk().bucket(edge).member( *edge_rank_part ));
    EXPECT_TRUE( get_bulk().bucket(edge).member( *side_rank_part ));
    EXPECT_TRUE(!get_bulk().bucket(edge).member( *element_rank_part)); // See!  part induction is not transitive!
    EXPECT_TRUE(!get_bulk().bucket(edge).member( *element_rank_superset_part)); // See!  part induction is not transitive!
    EXPECT_TRUE(!get_bulk().bucket(edge).member( *unranked_superset_part));

    EXPECT_TRUE( get_bulk().bucket(node).member( *edge_rank_part ));
    EXPECT_TRUE( get_bulk().bucket(node).member( *side_rank_part ));
    EXPECT_TRUE( get_bulk().bucket(node).member( *element_rank_part));
    EXPECT_TRUE( get_bulk().bucket(node).member( *element_rank_superset_part));
    EXPECT_TRUE(!get_bulk().bucket(node).member( *unranked_superset_part));
  }
}

TEST ( UnitTestInducedPart, verifyForceNoInduce )
{
  stk::ParallelMachine pm = MPI_COMM_SELF;

  const unsigned spatial_dim = 2;
  const bool force_no_induce = true;

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dim, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& mesh = *bulkPtr;
  Part& unranked_part = meta_data.declare_part("unranked_part");
  Part& element_rank_part = meta_data.declare_part("element_rank_part", stk::topology::ELEMENT_RANK);
  Part& element_rank_part_no_induce = meta_data.declare_part("element_rank_part_no_induce", stk::topology::ELEMENT_RANK, force_no_induce);
  Part& element_rank_part_change_to_no_induce = meta_data.declare_part("element_rank_part_change_to_no_induce", stk::topology::ELEMENT_RANK);

  meta_data.force_no_induce(element_rank_part_change_to_no_induce);

  meta_data.commit();

  mesh.modification_begin();

  stk::mesh::PartVector parts;
  const stk::mesh::EntityId element_id = 1, node_id = 1;
  parts.push_back(&unranked_part);
  parts.push_back(&element_rank_part);
  parts.push_back(&element_rank_part_no_induce);
  parts.push_back(&element_rank_part_change_to_no_induce);
  Entity elem = mesh.declare_element(element_id, parts);

  parts.clear();
  Entity node = mesh.declare_node(node_id, parts);

  const stk::mesh::RelationIdentifier rel_0_id = 0, rel_1_id = 1;
  mesh.declare_relation(elem, node,  rel_0_id);
  mesh.declare_relation(elem, node,  rel_1_id);

  EXPECT_TRUE( mesh.bucket(node).member(element_rank_part) );
  EXPECT_FALSE( mesh.bucket(node).member(element_rank_part_no_induce) );
  EXPECT_FALSE( mesh.bucket(node).member(element_rank_part_change_to_no_induce) );
}

}
