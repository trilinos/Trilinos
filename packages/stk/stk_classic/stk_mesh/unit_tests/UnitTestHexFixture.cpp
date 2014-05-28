/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>

namespace {

using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::EntityId;

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_x )
{
  // Test 3x1x1 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(1,0,0), 2u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(2,0,0), 3u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_x )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(1,0,0), 2u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(2,0,0), 3u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_y )
{
  // Test 1x3x1 HexFixture structure
  const unsigned NX = 1;
  const unsigned NY = 3;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,1,0), 2u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,2,0), 3u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_y )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,1,0), 4u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,2,0), 7u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_z )
{
  // Test 1x1x3 HexFixture structure
  const unsigned NX = 1;
  const unsigned NY = 1;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,1), 2u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,2), 3u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_z )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,1), 10u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,2), 19u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_diag )
{
  // Test 3x3x3 HexFixture structure
  const unsigned NX = 3;
  const unsigned NY = 3;
  const unsigned NZ = 3;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();
  STKUNIT_EXPECT_EQUAL( hf.elem_id(0,0,0), 1u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(1,1,1), 14u );
  STKUNIT_EXPECT_EQUAL( hf.elem_id(2,2,2), 27u );
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, trivial_parallel_2 )
{
  // Test a customized element distribution with one element on proc 0 and 1
  // and none on the other procs

  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 2;
  const unsigned NY = 1;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. One element will go on rank 0, the other on rank 1, any
  // other ranks get nothing.
  std::map<unsigned,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    parallel_distribution[0] = element_ids;
    element_ids[0] = 2;
    parallel_distribution[1] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  if (p_rank <= 1) {
    hf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    hf.generate_mesh( empty_vector ) ;
  }

  stk::mesh::BulkData & mesh = hf.m_bulk_data;

  // Verify element_id 1 is owned by proc 0
  // Verify element_id 2 is owned by proc 1
  const EntityRank element_rank = hf.m_fem_meta.element_rank();
  Entity * entity_1 = mesh.get_entity(element_rank, 1);
  Entity * entity_2 = mesh.get_entity(element_rank ,2);
  if (p_rank <= 1) {
    STKUNIT_ASSERT_TRUE( entity_1 != NULL );
    STKUNIT_ASSERT_TRUE( entity_2 != NULL );
    STKUNIT_EXPECT_EQUAL( 0u, entity_1->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_2->owner_rank() );
  }
  else {
    STKUNIT_EXPECT_TRUE( entity_1 == NULL );
    STKUNIT_EXPECT_TRUE( entity_2 == NULL );
  }
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, disjoint_parallel_psizex1x1 )
{
  // Test a customized element distribution with one element on each proc

  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  const unsigned NX = p_size;
  const unsigned NY = 1;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign each processor an element such that rank+1 = elem_id
  std::map<unsigned,std::vector<EntityId> > parallel_distribution;
  for (unsigned p=0 ; p < p_size ; ++p) {
    std::vector< EntityId> element_ids;
    element_ids.push_back(p+1); // element id's start at 1
    parallel_distribution[p] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh(parallel_distribution[p_rank]);
  stk::mesh::BulkData & mesh = hf.m_bulk_data;
  const EntityRank element_rank = hf.m_fem_meta.element_rank();

  // We should always know about, and own, the element assigned to us
  Entity * my_entity    = mesh.get_entity(element_rank, p_rank + 1);
  STKUNIT_ASSERT_TRUE( my_entity != NULL );
  STKUNIT_EXPECT_EQUAL( p_rank, my_entity->owner_rank() );

  // If applicable, we know about the element on adjacent lower rank
  if (p_rank > 0) {
    Entity * prior_entity = mesh.get_entity(element_rank, p_rank);
    STKUNIT_ASSERT_TRUE( prior_entity != NULL );
    STKUNIT_EXPECT_EQUAL( p_rank - 1, prior_entity->owner_rank() );
  }

  // If applicable, we know about the element on adjacent higher rank
  if (p_rank < p_size - 1) {
    Entity * next_entity   = mesh.get_entity(element_rank, p_rank + 2);
    STKUNIT_ASSERT_TRUE( next_entity != NULL );
    STKUNIT_EXPECT_EQUAL( p_rank + 1, next_entity->owner_rank() );
  }
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, disjoint_parallel_4x2x1 )
{
  // layout: 4x2x1 hex mesh
  // elements:
  // [ e_1, e_2, e_3, e_4 ]
  // [ e_5, e_6, e_7, e_8 ]
  // processors:
  // [ p_0, p_1, p_1, p_1 ]
  // [ p_1, p_0, p_1, p_1 ]
  //

  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 4;
  const unsigned NY = 2;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign 1,6 to proc 0, all the rest to proc 1. The other
  // procs get nothing.
  std::map<unsigned,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    element_ids.push_back(6);
    parallel_distribution[0] = element_ids; // proc 0
    element_ids.clear();
    element_ids.push_back(2);
    element_ids.push_back(3);
    element_ids.push_back(4);
    element_ids.push_back(5);
    element_ids.push_back(7);
    element_ids.push_back(8);
    parallel_distribution[1] = element_ids; // proc 1
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  if (p_rank <= 1) {
    hf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    hf.generate_mesh( empty_vector ) ;
  }

  stk::mesh::BulkData & mesh = hf.m_bulk_data;

  // Verify that the entities and known and owned by the appropriate procs
  const EntityRank element_rank = hf.m_fem_meta.element_rank();
  Entity * entity_1 = mesh.get_entity(element_rank, 1);
  Entity * entity_2 = mesh.get_entity(element_rank, 2);
  Entity * entity_3 = mesh.get_entity(element_rank, 3);
  Entity * entity_4 = mesh.get_entity(element_rank, 4);
  Entity * entity_5 = mesh.get_entity(element_rank, 5);
  Entity * entity_6 = mesh.get_entity(element_rank, 6);
  Entity * entity_7 = mesh.get_entity(element_rank, 7);
  Entity * entity_8 = mesh.get_entity(element_rank, 8);
  if (p_rank == 0) {
    STKUNIT_ASSERT_TRUE( entity_1 != NULL );
    STKUNIT_ASSERT_TRUE( entity_2 != NULL );
    STKUNIT_ASSERT_TRUE( entity_3 != NULL );
    STKUNIT_ASSERT_TRUE( entity_4 == NULL );
    STKUNIT_ASSERT_TRUE( entity_5 != NULL );
    STKUNIT_ASSERT_TRUE( entity_6 != NULL );
    STKUNIT_ASSERT_TRUE( entity_7 != NULL );
    STKUNIT_ASSERT_TRUE( entity_8 == NULL );
    STKUNIT_EXPECT_EQUAL( 0u, entity_1->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_2->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_3->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_5->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 0u, entity_6->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_7->owner_rank() );
  }
  else if (p_rank == 1) {
    STKUNIT_ASSERT_TRUE( entity_1 != NULL );
    STKUNIT_ASSERT_TRUE( entity_2 != NULL );
    STKUNIT_ASSERT_TRUE( entity_3 != NULL );
    STKUNIT_ASSERT_TRUE( entity_4 != NULL );
    STKUNIT_ASSERT_TRUE( entity_5 != NULL );
    STKUNIT_ASSERT_TRUE( entity_6 != NULL );
    STKUNIT_ASSERT_TRUE( entity_7 != NULL );
    STKUNIT_ASSERT_TRUE( entity_8 != NULL );
    STKUNIT_EXPECT_EQUAL( 0u, entity_1->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_2->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_3->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_4->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_5->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 0u, entity_6->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_7->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_8->owner_rank() );
  }
  else {
    STKUNIT_EXPECT_TRUE( entity_1 == NULL );
    STKUNIT_EXPECT_TRUE( entity_2 == NULL );
    STKUNIT_EXPECT_TRUE( entity_3 == NULL );
    STKUNIT_EXPECT_TRUE( entity_4 == NULL );
    STKUNIT_EXPECT_TRUE( entity_5 == NULL );
    STKUNIT_EXPECT_TRUE( entity_6 == NULL );
    STKUNIT_EXPECT_TRUE( entity_7 == NULL );
    STKUNIT_EXPECT_TRUE( entity_8 == NULL );
  }
}

STKUNIT_UNIT_TEST( UnitTestHexFixture, disjoint_parallel_5x1x1 )
{
  // layout:
  // [ e_1, e_2, e_3, e_4, e_5 ] elements
  // [ p_0, p_1, p_1, p_1, p_0 ] processors
  //
  const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  // Skip unless p_size is at least 2
  if (p_size < 2)
    return;

  const unsigned NX = 5;
  const unsigned NY = 1;
  const unsigned NZ = 1;

  // map< processor, vector of element ids >, this is our custom parallel
  // distribution. Assign 1,5 to proc 0, all the rest to proc 1. The other
  // procs get nothing.
  std::map<unsigned,std::vector<EntityId> > parallel_distribution;
  {
    std::vector< EntityId> element_ids;
    element_ids.push_back(1);
    element_ids.push_back(5);
    parallel_distribution[0] = element_ids;
    element_ids.clear();
    element_ids.push_back(2);
    element_ids.push_back(3);
    element_ids.push_back(4);
    parallel_distribution[1] = element_ids;
  }

  // Create the fixture
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
  hf.m_fem_meta.commit();
  if (p_rank <= 1) {
    hf.generate_mesh(parallel_distribution[p_rank]);
  }
  else {
    std::vector<EntityId> empty_vector;
    hf.generate_mesh( empty_vector ) ;
  }

  stk::mesh::BulkData & mesh = hf.m_bulk_data;

  // Verify that the entities and known and owned by the appropriate procs
  const EntityRank element_rank = hf.m_fem_meta.element_rank();
  Entity * entity_1 = mesh.get_entity(element_rank, 1);
  Entity * entity_2 = mesh.get_entity(element_rank, 2);
  Entity * entity_3 = mesh.get_entity(element_rank, 3);
  Entity * entity_4 = mesh.get_entity(element_rank, 4);
  Entity * entity_5 = mesh.get_entity(element_rank, 5);
  if (p_rank == 0) {
    STKUNIT_ASSERT_TRUE( entity_1 != NULL );
    STKUNIT_ASSERT_TRUE( entity_2 != NULL );
    STKUNIT_ASSERT_TRUE( entity_3 == NULL );
    STKUNIT_ASSERT_TRUE( entity_4 != NULL );
    STKUNIT_ASSERT_TRUE( entity_5 != NULL );
    STKUNIT_EXPECT_EQUAL( 0u, entity_1->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_2->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_4->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 0u, entity_5->owner_rank() );
  }
  else if (p_rank == 1) {
    STKUNIT_ASSERT_TRUE( entity_1 != NULL );
    STKUNIT_ASSERT_TRUE( entity_2 != NULL );
    STKUNIT_ASSERT_TRUE( entity_3 != NULL );
    STKUNIT_ASSERT_TRUE( entity_4 != NULL );
    STKUNIT_ASSERT_TRUE( entity_5 != NULL );
    STKUNIT_EXPECT_EQUAL( 0u, entity_1->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_2->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_3->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 1u, entity_4->owner_rank() );
    STKUNIT_EXPECT_EQUAL( 0u, entity_5->owner_rank() );
  }
  else {
    STKUNIT_EXPECT_TRUE( entity_1 == NULL );
    STKUNIT_EXPECT_TRUE( entity_2 == NULL );
    STKUNIT_EXPECT_TRUE( entity_3 == NULL );
    STKUNIT_EXPECT_TRUE( entity_4 == NULL );
    STKUNIT_EXPECT_TRUE( entity_5 == NULL );
  }
}

} // end namespace
