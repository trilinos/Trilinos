/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>

using stk_classic::mesh::Part;
using stk_classic::mesh::fem::FEMMetaData;
using stk_classic::mesh::BulkData;
using stk_classic::mesh::Entity;

namespace {

bool has_part(const Entity& entity, const Part& part)
{
  return entity.bucket().member(part);
}

// Set up a very simple mesh with one element, one side, one node:
// rels: E->(S1,S2)->N
// parts: E in element_rank_part, unranked_part
//        S1, S2 in side_rank_part
//        element_ranked_part subset of unranked_superset_part
// modification cycle is left uncompleted
#define SETUP_MESH()                            \
  stk_classic::ParallelMachine pm = MPI_COMM_SELF;              \
                                                  \
  const unsigned spatial_dim = 2;                                       \
                                                                        \
  std::vector<std::string> entity_rank_names = stk_classic::mesh::fem::entity_rank_names(spatial_dim); \
  FEMMetaData meta_data(spatial_dim, entity_rank_names);                \
  Part& unranked_part = meta_data.declare_part("unranked_part");        \
  Part& element_rank_part = meta_data.declare_part("element_rank_part", meta_data.element_rank()); \
  Part& element_rank_superset_part = meta_data.declare_part("element_rank_superset_part", meta_data.element_rank()); \
  Part& side_rank_part = meta_data.declare_part("side_rank_part", meta_data.side_rank()); \
  Part& unranked_superset_part = meta_data.declare_part("unranked_superset_part"); \
  meta_data.declare_part_subset(unranked_superset_part, element_rank_part); \
  meta_data.declare_part_subset(element_rank_superset_part, element_rank_part); \
                                                                        \
  meta_data.commit();                                                   \
  BulkData mesh(FEMMetaData::get_meta_data(meta_data), pm);             \
                                                                        \
  mesh.modification_begin();                                            \
                                                                        \
  stk_classic::mesh::PartVector parts;                                          \
  parts.push_back(&unranked_part);                                       \
  parts.push_back(&element_rank_part);                                   \
  Entity& elem = mesh.declare_entity(meta_data.element_rank(), 1 /*id*/, parts); \
                                                                        \
  parts.clear();                                                        \
  parts.push_back(&side_rank_part);                                      \
  Entity& side1 = mesh.declare_entity(meta_data.side_rank(), 1 /*id*/, parts); \
  Entity& side2 = mesh.declare_entity(meta_data.side_rank(), 2 /*id*/, parts); \
                                                                        \
  parts.clear();                                                        \
  Entity& node = mesh.declare_entity(meta_data.node_rank(), 1 /*id*/, parts); \
                                                                        \
  mesh.declare_relation(elem, side1,  0 /*rel id*/);                    \
  mesh.declare_relation(elem, side2,  1 /*rel id*/);                    \
  mesh.declare_relation(side1, node,  0 /*rel id*/);                    \
  mesh.declare_relation(side2, node,  0 /*rel id*/);

STKUNIT_UNIT_TEST ( UnitTestInducedPart , verifyBasicInducedPart )
{
  SETUP_MESH();

  // Check that directly-induced parts are induced upon relation creation
  // before modification end.
  STKUNIT_EXPECT_TRUE(has_part(node, side_rank_part));
  STKUNIT_EXPECT_TRUE(has_part(side1, element_rank_part));
  STKUNIT_EXPECT_TRUE(has_part(side2, element_rank_part));

  mesh.modification_end();

  // Modification-end should not have changed induced parts
  STKUNIT_EXPECT_TRUE(has_part(node, side_rank_part));
  STKUNIT_EXPECT_TRUE(has_part(side1, element_rank_part));
  STKUNIT_EXPECT_TRUE(has_part(side2, element_rank_part));
}

STKUNIT_UNIT_TEST ( UnitTestInducedPart , verifyNotTransitiveInducedPart )
{
  SETUP_MESH();

  // Node should not have picked-up the element_rank_part indirectly through
  // it's relation to the sides because induced-parts are not supposed to be
  // transitive according to the STK_Mesh domain design.
  // TODO: Are we sure we don't want induced parts to be transitive??
  STKUNIT_EXPECT_TRUE(!has_part(node, element_rank_part));
}

STKUNIT_UNIT_TEST ( UnitTestInducedPart, verifyInducedPartCorrectnessWhenRelationsRemoved )
{
  SETUP_MESH();

  // Destroy one relation from element to a side, confirm that the side that lost
  // the relation no longer has the element part.
  mesh.destroy_relation(elem, side1, 0 /*rel id*/);
  STKUNIT_EXPECT_TRUE(!has_part(side1, element_rank_part));

  // Destroy one of the relations from side to node. Confirm that node still has
  // side part due to its remaining relation.
  mesh.destroy_relation(side1, node, 0 /*rel id*/);
  STKUNIT_EXPECT_TRUE(has_part(node, side_rank_part));

  // Destroy the other relations from side to node. Confirm that node no longer
  // has any induced parts.
  mesh.destroy_relation(side2, node, 0 /*rel id*/);
  STKUNIT_EXPECT_TRUE(!has_part(node, side_rank_part));
}

STKUNIT_UNIT_TEST ( UnitTestInducedPart , verifySupersetsOfInducedPart )
{
  SETUP_MESH();

  // Check for superset/subset consistency in induced parts. If an entity is
  // induced into a part, it should also be induced into the supersets of
  // that part even if the superset parts are unranked.
  STKUNIT_EXPECT_TRUE(has_part(side1, element_rank_superset_part));
  STKUNIT_EXPECT_TRUE(has_part(side1, unranked_superset_part));
}

}
