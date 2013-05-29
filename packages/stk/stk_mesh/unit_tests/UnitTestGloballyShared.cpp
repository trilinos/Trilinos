/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/fixtures/QuadFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;


STKUNIT_UNIT_TEST( UnitTestGloballyShared, keyhole_3x1 )
{
  // layout:
  // [ e_1, e_2, e_3 ] elements
  // [ p_0, p_1, p_0 ] processors
  //
  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 3;
  const unsigned NY = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign 1,5 to proc 0, all the rest to proc 1. The other
  // procs get nothing.
  std::map<unsigned,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    element_ids.push_back(3);
    parallel_distribution[0] = element_ids;
    element_ids.clear();
    element_ids.push_back(2);
    parallel_distribution[1] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::QuadFixture qf(MPI_COMM_WORLD,NX,NY);
  qf.m_meta.commit();
  if (p_rank <= 1) {
    qf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    qf.generate_mesh( empty_vector ) ;
  }

  BulkData & mesh = qf.m_bulk_data;

  stk::mesh::create_edges(mesh);

  // Quad edge ordinals:
  //            2
  //          -----
  //         |     |
  //       3 |     | 1
  //         |     |
  //          -----
  //            0

  // Verify that the entities and known and owned by the appropriate procs
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;
  Entity element_1 = mesh.get_entity(element_rank, 1);
  Entity element_2 = mesh.get_entity(element_rank, 2);
  Entity element_3 = mesh.get_entity(element_rank, 3);
  if (p_rank == 0) {
    STKUNIT_ASSERT_TRUE( element_1.is_valid() );
    STKUNIT_ASSERT_TRUE( element_2.is_valid() );
    STKUNIT_ASSERT_TRUE( element_3.is_valid() );
    STKUNIT_EXPECT_EQUAL( 0, element_1.owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1, element_2.owner_rank() );
    STKUNIT_EXPECT_EQUAL( 0, element_3.owner_rank() );
    // Verify global sharing of edges on element_1 and element_3
    // element_1:  edge_1 should be globally shared
    // element_3:  edge_3 should be globally shared
    stk::mesh::Entity const* element_1_edge_relations = mesh.begin_edges(element_1);
    const int num_element_1_edges = mesh.num_edges(element_1);
    STKUNIT_ASSERT_EQUAL( 4, num_element_1_edges );
    Entity element_1_edge_1 = element_1_edge_relations[1];
    STKUNIT_ASSERT_TRUE( element_1_edge_1.is_valid() );
    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_1_edge_1.key()) );

    stk::mesh::Entity const* element_3_edge_relations = mesh.begin_edges(element_3);
    const int num_element_3_edges = mesh.num_edges(element_3);
    STKUNIT_ASSERT_EQUAL( 4, num_element_3_edges );
    Entity element_3_edge_3 = element_3_edge_relations[3];
    STKUNIT_ASSERT_TRUE( element_3_edge_3.is_valid() );
    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_3_edge_3.key()) );

    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_1.key()) );
    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_2.key()) );
    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_3.key()) );
  }
  else if (p_rank == 1) {
    STKUNIT_ASSERT_TRUE( element_1.is_valid() );
    STKUNIT_ASSERT_TRUE( element_2.is_valid() );
    STKUNIT_ASSERT_TRUE( element_3.is_valid() );
    STKUNIT_EXPECT_EQUAL( 0, element_1.owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1, element_2.owner_rank() );
    STKUNIT_EXPECT_EQUAL( 0, element_3.owner_rank() );
    // Verify global sharing of edges on element_2
    // element_2:  edge_0 and edge_2 should _not be_ globally shared
    //             edge_1 and edge_3 should _be_ globally shared
    stk::mesh::Entity const* element_2_edge_relations = mesh.begin_edges(element_2);
    const int num_element_2_edges = mesh.num_edges(element_2);
    STKUNIT_ASSERT_EQUAL( 4, num_element_2_edges );

    Entity element_2_edge_0 = element_2_edge_relations[0];
    STKUNIT_ASSERT_TRUE( element_2_edge_0.is_valid() );

    Entity element_2_edge_2 = element_2_edge_relations[2];
    STKUNIT_ASSERT_TRUE( element_2_edge_2.is_valid() );

    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_2_edge_0.key()) );
    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_2_edge_2.key()) );

    Entity element_2_edge_1 = element_2_edge_relations[1];
    STKUNIT_ASSERT_TRUE( element_2_edge_1.is_valid() );
    Entity element_2_edge_3 = element_2_edge_relations[3];
    STKUNIT_ASSERT_TRUE( element_2_edge_3.is_valid() );

    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_2_edge_1.key()) );
    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_2_edge_3.key()) );

    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_1.key()) );
    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_2.key()) );
    STKUNIT_EXPECT_FALSE( mesh.in_shared(element_3.key()) );
  }
  else {
    STKUNIT_EXPECT_TRUE( !element_1.is_valid() );
    STKUNIT_EXPECT_TRUE( !element_2.is_valid() );
    STKUNIT_EXPECT_TRUE( !element_3.is_valid() );
  }
}


// Note:  The behavior of this unit test is not entirely the desired behavior.
// For Finite Element Models, we expect to disjointly distribute elements to the
// processors, so elements should never end up being globally shared. Stk_mesh
// does not know about this and so is happy to globally share elements if they
// are constructed identically on each involved processor and there is a
// downward relation to them across processors, as this test demonstrates.
STKUNIT_UNIT_TEST( UnitTestGloballyShared, constraints_element_2x1 )
{
  // layout:
  // [ c_1, c_2 ] constraints
  // [ p_0, p_1 ] processors
  // both constraints have downward relations to both elements
  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  std::vector< std::string > rank_names;
  rank_names.reserve( 5 );
  rank_names.push_back(std::string("NODE"));
  rank_names.push_back(std::string("EDGE"));
  rank_names.push_back(std::string("FACE"));
  rank_names.push_back(std::string("ELEMENT"));
  rank_names.push_back(std::string("CONSTRAINT"));

  const int spatial_dimension = 2;
  MetaData meta(spatial_dimension,rank_names);
  BulkData mesh(meta,MPI_COMM_WORLD);

  mesh.modification_begin();
  if (p_rank == 0) {
      Entity node_1 = mesh.declare_entity(stk::topology::NODE_RANK, 1);
      Entity node_2 = mesh.declare_entity(stk::topology::NODE_RANK, 2);
      Entity node_3 = mesh.declare_entity(stk::topology::NODE_RANK, 3);
      Entity node_4 = mesh.declare_entity(stk::topology::NODE_RANK, 4);
      Entity node_5 = mesh.declare_entity(stk::topology::NODE_RANK, 5);
      Entity node_6 = mesh.declare_entity(stk::topology::NODE_RANK, 6);

      Entity element_1 = mesh.declare_entity(stk::topology::ELEMENT_RANK, 1);
      Entity element_2 = mesh.declare_entity(stk::topology::ELEMENT_RANK, 2);

      mesh.declare_relation(element_1,node_1,0);
      mesh.declare_relation(element_1,node_2,1);
      mesh.declare_relation(element_1,node_3,2);
      mesh.declare_relation(element_1,node_4,3);

      mesh.declare_relation(element_2,node_2,0);
      mesh.declare_relation(element_2,node_5,1);
      mesh.declare_relation(element_2,node_6,2);
      mesh.declare_relation(element_2,node_3,3);

      Entity constraint_1 = mesh.declare_entity(stk::topology::CONSTRAINT_RANK, 1);
      mesh.declare_relation(constraint_1,element_1,0);
      mesh.declare_relation(constraint_1,element_2,1);
  }
  else if (p_rank == 1) {
      Entity node_1 = mesh.declare_entity(stk::topology::NODE_RANK, 1);
      Entity node_2 = mesh.declare_entity(stk::topology::NODE_RANK, 2);
      Entity node_3 = mesh.declare_entity(stk::topology::NODE_RANK, 3);
      Entity node_4 = mesh.declare_entity(stk::topology::NODE_RANK, 4);
      Entity node_5 = mesh.declare_entity(stk::topology::NODE_RANK, 5);
      Entity node_6 = mesh.declare_entity(stk::topology::NODE_RANK, 6);

      Entity element_1 = mesh.declare_entity(stk::topology::ELEMENT_RANK, 1);
      Entity element_2 = mesh.declare_entity(stk::topology::ELEMENT_RANK, 2);

      mesh.declare_relation(element_1,node_1,0);
      mesh.declare_relation(element_1,node_2,1);
      mesh.declare_relation(element_1,node_3,2);
      mesh.declare_relation(element_1,node_4,3);

      mesh.declare_relation(element_2,node_2,0);
      mesh.declare_relation(element_2,node_5,1);
      mesh.declare_relation(element_2,node_6,2);
      mesh.declare_relation(element_2,node_3,3);

      Entity constraint_2 = mesh.declare_entity(stk::topology::CONSTRAINT_RANK, 2);
      mesh.declare_relation(constraint_2,element_1,0);
      mesh.declare_relation(constraint_2,element_2,1);
  }
  else {
      // Do nothing on extra processors
  }
  mesh.modification_end();


  Entity element_1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1);
  Entity element_2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 2);

  Entity constraint_1 = mesh.get_entity(stk::topology::CONSTRAINT_RANK, 1);
  Entity constraint_2 = mesh.get_entity(stk::topology::CONSTRAINT_RANK, 2);

  // Verify that the entities and known and owned by the appropriate procs
  if (p_rank == 0) {
    STKUNIT_ASSERT_TRUE( element_1.is_valid() );
    STKUNIT_ASSERT_TRUE( element_2.is_valid() );
    STKUNIT_EXPECT_EQUAL( 0, constraint_1.owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1, constraint_2.owner_rank() );

    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_1.key()) );
    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_2.key()) );
  }
  else if (p_rank == 1) {
    STKUNIT_ASSERT_TRUE( element_1.is_valid() );
    STKUNIT_ASSERT_TRUE( element_2.is_valid() );
    STKUNIT_EXPECT_EQUAL( 0, constraint_1.owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1, constraint_2.owner_rank() );

    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_1.key()) );
    STKUNIT_EXPECT_TRUE( mesh.in_shared(element_2.key()) );
  }
  else {
    STKUNIT_EXPECT_TRUE( !element_1.is_valid() );
    STKUNIT_EXPECT_TRUE( !element_2.is_valid() );
  }
}
