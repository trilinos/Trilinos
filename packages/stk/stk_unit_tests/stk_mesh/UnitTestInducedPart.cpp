/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }





using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;

namespace {

// Set up a very simple mesh with one element, one side, one node:
// rels: E->(S1,S2)->N
// parts: E in element_rank_part, unranked_part
//        S1, S2 in side_rank_part
//        element_ranked_part subset of unranked_superset_part
// modification cycle is left uncompleted
#define SETUP_MESH()                            \
  stk::ParallelMachine pm = MPI_COMM_SELF;              \
                                                  \
  const unsigned spatial_dim = 2;                                       \
                                                                        \
  MetaData meta_data(spatial_dim);                                      \
  Part& unranked_part = meta_data.declare_part("unranked_part");        \
  Part& element_rank_part =                                             \
     meta_data.declare_part_with_topology("element_rank_part", stk::topology::TRI_3);        \
  Part& element_rank_superset_part =                                    \
    meta_data.declare_part("element_rank_superset_part", stk::topology::ELEMENT_RANK); \
  Part& side_rank_part =                                                \
    meta_data.declare_part_with_topology("side_rank_part", stk::topology::LINE_2); \
  Part& unranked_superset_part = meta_data.declare_part("unranked_superset_part"); \
  meta_data.declare_part_subset(unranked_superset_part, element_rank_part); \
  meta_data.declare_part_subset(element_rank_superset_part, element_rank_part); \
                                                                        \
  meta_data.commit();                                                   \
  BulkData mesh(meta_data, pm);                                         \
                                                                        \
  mesh.modification_begin();                                            \
                                                                        \
  stk::mesh::PartVector parts;                                          \
  parts.push_back(&unranked_part);                                       \
  parts.push_back(&element_rank_part);                                   \
  Entity elem = mesh.declare_entity(stk::topology::ELEMENT_RANK, 1 /*id*/, parts); \
                                                                        \
  parts.clear();                                                        \
  parts.push_back(&side_rank_part);                                      \
  Entity side1 = mesh.declare_entity(meta_data.side_rank(), 1 /*id*/, parts); \
  Entity side2 = mesh.declare_entity(meta_data.side_rank(), 2 /*id*/, parts); \
                                                                        \
  parts.clear();                                                        \
  Entity node  = mesh.declare_entity(stk::topology::NODE_RANK, 1 /*id*/, parts);      \
  Entity node2 = mesh.declare_entity(stk::topology::NODE_RANK, 2 /*id*/, parts);      \
  Entity node3 = mesh.declare_entity(stk::topology::NODE_RANK, 3 /*id*/, parts);      \
                                                                        \
  mesh.declare_relation(elem, side1,  0 /*rel id*/);                    \
  mesh.declare_relation(elem, node,   0 /*rel id*/);                    \
  mesh.declare_relation(elem, node2,  1 /*rel id*/);                    \
  mesh.declare_relation(elem, node3,  2 /*rel id*/);                    \
  mesh.declare_relation(elem, side2,  1 /*rel id*/);                    \
  mesh.declare_relation(side1, node,  0 /*rel id*/);                    \
  mesh.declare_relation(side1, node2,  1 /*rel id*/);                   \
  mesh.declare_relation(side2, node,  0 /*rel id*/);                    \
  mesh.declare_relation(side2, node3,  1 /*rel id*/);

TEST ( UnitTestInducedPart , verifyBasicInducedPart )
{
  SETUP_MESH();

  // Check that directly-induced parts are induced upon relation creation
  // before modification end.
  EXPECT_TRUE(mesh.bucket(node).member( side_rank_part));
  EXPECT_TRUE(mesh.bucket(side1).member( element_rank_part));
  EXPECT_TRUE(mesh.bucket(side2).member( element_rank_part));

  mesh.modification_end();

  // Modification-end should not have changed induced parts
  EXPECT_TRUE(mesh.bucket(node).member( side_rank_part));
  EXPECT_TRUE(mesh.bucket(side1).member( element_rank_part));
  EXPECT_TRUE(mesh.bucket(side2).member( element_rank_part));
}

TEST ( UnitTestInducedPart, verifyInducedPartCorrectnessWhenRelationsRemoved )
{
  SETUP_MESH();

  // Destroy one relation from element to a side, confirm that the side that lost
  // the relation no longer has the element part.
  mesh.destroy_relation(elem, side1, 0 /*rel id*/);
  EXPECT_TRUE(!mesh.bucket(side1).member( element_rank_part));

  // Destroy one of the relations from side to node. Confirm that node still has
  // side part due to its remaining relation.
  mesh.destroy_relation(side1, node, 0 /*rel id*/);
  EXPECT_TRUE(mesh.bucket(node).member( side_rank_part));

  // Destroy the other relations from side to node. Confirm that node no longer
  // has any induced parts.
  mesh.destroy_relation(side2, node, 0 /*rel id*/);
  EXPECT_TRUE(!mesh.bucket(node).member( side_rank_part));
}

TEST ( UnitTestInducedPart , verifySupersetsOfInducedPart )
{
  SETUP_MESH();

  // Check for superset/subset consistency in induced parts. If an entity is
  // induced into a part, it should also be induced into the supersets of
  // that part even if the superset parts are unranked.
  EXPECT_TRUE(mesh.bucket(side1).member( element_rank_superset_part));
  EXPECT_TRUE(mesh.bucket(side1).member( unranked_superset_part));
}

TEST ( UnitTestInducedPart, verifyForceNoInduce )
{
  stk::ParallelMachine pm = MPI_COMM_SELF;

  const unsigned spatial_dim = 2;
  const bool force_no_induce = true;

  MetaData meta_data(spatial_dim);
  Part& unranked_part = meta_data.declare_part("unranked_part");
  Part& element_rank_part = meta_data.declare_part("element_rank_part", stk::topology::ELEMENT_RANK);
  Part& element_rank_part_no_induce = meta_data.declare_part("element_rank_part_no_induce", stk::topology::ELEMENT_RANK, force_no_induce);
  Part& element_rank_part_change_to_no_induce = meta_data.declare_part("element_rank_part_change_to_no_induce", stk::topology::ELEMENT_RANK);

  meta_data.force_no_induce(element_rank_part_change_to_no_induce);

  meta_data.commit();
  BulkData mesh(meta_data, pm);

  mesh.modification_begin();

  stk::mesh::PartVector parts;
  const stk::mesh::EntityId element_id = 1, node_id = 1;
  parts.push_back(&unranked_part);
  parts.push_back(&element_rank_part);
  parts.push_back(&element_rank_part_no_induce);
  parts.push_back(&element_rank_part_change_to_no_induce);
  Entity elem = mesh.declare_entity(stk::topology::ELEMENT_RANK, element_id, parts);

  parts.clear();
  Entity node = mesh.declare_entity(stk::topology::NODE_RANK, node_id, parts);

  const stk::mesh::RelationIdentifier rel_0_id = 0, rel_1_id = 1;
  mesh.declare_relation(elem, node,  rel_0_id);
  mesh.declare_relation(elem, node,  rel_1_id);

  EXPECT_TRUE( mesh.bucket(node).member(element_rank_part) );
  EXPECT_FALSE( mesh.bucket(node).member(element_rank_part_no_induce) );
  EXPECT_FALSE( mesh.bucket(node).member(element_rank_part_change_to_no_induce) );
}

}
