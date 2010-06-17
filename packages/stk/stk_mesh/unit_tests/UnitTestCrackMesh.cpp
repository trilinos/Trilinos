/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <unit_tests/UnitTestMesh.hpp>
#include <unit_tests/UnitTestBoxMeshFixture.hpp>

namespace {

void copy_node_and_break_relations( stk::mesh::BulkData  & mesh,
                                    stk::mesh::Entity    & old_node,
                                    stk::mesh::Entity    & new_node )
{


  std::vector<std::pair<stk::mesh::Entity *, unsigned> > relation_and_side_ordinal;

  stk::mesh::PairIterRelation rel = old_node.relations();

  for (; rel.first != rel.second; ++rel.first) {
    stk::mesh::Entity * entity = (rel.first->entity());
    unsigned side_ordinal = rel.first->identifier();

    relation_and_side_ordinal.push_back(
        std::pair<stk::mesh::Entity *, unsigned>(entity,side_ordinal)
        );
  }

  for ( std::vector<std::pair<stk::mesh::Entity *, unsigned > >::iterator
      itr = relation_and_side_ordinal.begin();
      itr != relation_and_side_ordinal.end();
      ++itr
      )
  {
    stk::mesh::Entity & entity = *(itr->first);
    unsigned side_ordinal = itr->second;

    mesh.destroy_relation( entity, old_node);
    mesh.declare_relation( entity , new_node, side_ordinal);
  }


}
}

STKUNIT_UNIT_TEST ( UnitTestCrackMesh , verifyBoxGhosting )
{
  const unsigned p_size = stk::parallel_machine_size( MPI_COMM_WORLD );
  if ( 1 < p_size ) { return ; }

  BoxMeshFixture fixture( MPI_COMM_WORLD );

  stk::mesh::BulkData & mesh = fixture.m_bulk_data;

  fixture.fill_mesh();

  stk::mesh::Entity & old_node = * fixture.m_nodes[0][0][1];
  unsigned old_node_id = old_node.identifier();

  unsigned new_node_id = 28;
  const stk::mesh::PartVector no_parts;

  mesh.modification_begin();

  stk::mesh::Entity & new_node = mesh.declare_entity(stk::mesh::Node, new_node_id, no_parts);

  copy_node_and_break_relations(mesh, old_node, new_node);

  mesh.modification_end();

  stk::mesh::Entity * old_node_ptr = mesh.get_entity(stk::mesh::Node, old_node_id);
  STKUNIT_EXPECT_TRUE( NULL == old_node_ptr || 0 == old_node_ptr->relations().size());

  STKUNIT_EXPECT_TRUE( 2 == new_node.relations().size());

}


