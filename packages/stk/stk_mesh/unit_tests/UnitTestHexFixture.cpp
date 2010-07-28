/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <use_cases/HexFixture.hpp>

namespace {

  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_x )
  {
    const unsigned NX = 3;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][1] == 2 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][2] == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_x )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][1] == 2 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][2] == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_y )
  {
    const unsigned NX = 1;
    const unsigned NY = 3;
    const unsigned NZ = 1;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][1][0] == 2 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][2][0] == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_y )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][1][0] == 4 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][2][0] == 7 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_z )
  {
    const unsigned NX = 1;
    const unsigned NY = 1;
    const unsigned NZ = 3;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[1][0][0] == 2 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[2][0][0] == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_z )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[1][0][0] == 10 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[2][0][0] == 19 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_diag )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD);
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[0][0][0] == 1 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[1][1][1] == 14 );
    STKUNIT_EXPECT_TRUE( hf.m_elems_id[2][2][2] == 27 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, trivial_parallel_2 )
  {
    const unsigned NX = 2;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    // map< processor, vector of element ids >
    std::map<unsigned,std::vector<unsigned> > parallel_distribution;
    {
      std::vector<unsigned> element_ids;
      element_ids.push_back(1);
      parallel_distribution[0] = element_ids;
      element_ids[0] = 2;
      parallel_distribution[1] = element_ids;
    }
    const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (2 <= p_size) {
      HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD,&parallel_distribution[p_rank]);
      stk::mesh::BulkData & mesh = hf.m_bulk_data;
      // Verify element_id 1 is owned by proc 0
      // Verify element_id 2 is owned by proc 1
      stk::mesh::Entity * entity_1 = mesh.get_entity(stk::mesh::Element,1);
      stk::mesh::Entity * entity_2 = mesh.get_entity(stk::mesh::Element,2);
      if (p_rank <= 1) {
        STKUNIT_ASSERT_TRUE( entity_1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_2 != NULL );
        STKUNIT_EXPECT_TRUE( 0 == entity_1->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_2->owner_rank() );
      }
      else {
        STKUNIT_EXPECT_TRUE( entity_1 == NULL );
        STKUNIT_EXPECT_TRUE( entity_2 == NULL );
      }
    }
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, disjoint_parallel_5x1x1 )
  {
    // layout:
    // [ e_1, e_2, e_3, e_4, e_5 ] elements
    // [ p_0, p_1, p_1, p_1, p_0 ] processors
    //
    const unsigned NX = 5;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    // map< processor, vector of element ids >
    std::map<unsigned,std::vector<unsigned> > parallel_distribution;
    {
      std::vector<unsigned> element_ids;
      element_ids.push_back(1);
      element_ids.push_back(5);
      parallel_distribution[0] = element_ids;
      element_ids.clear();
      element_ids.push_back(2);
      element_ids.push_back(3);
      element_ids.push_back(4);
      parallel_distribution[1] = element_ids;
    }
    const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (2 <= p_size) {
      HexFixture<NX,NY,NZ> hf(MPI_COMM_WORLD,&parallel_distribution[p_rank]);
      stk::mesh::BulkData & mesh = hf.m_bulk_data;
      // Verify element_id 1 is owned by proc 0
      // Verify element_id 2 is owned by proc 1
      stk::mesh::Entity * entity_1 = mesh.get_entity(stk::mesh::Element,1);
      stk::mesh::Entity * entity_2 = mesh.get_entity(stk::mesh::Element,2);
      stk::mesh::Entity * entity_3 = mesh.get_entity(stk::mesh::Element,3);
      stk::mesh::Entity * entity_4 = mesh.get_entity(stk::mesh::Element,4);
      stk::mesh::Entity * entity_5 = mesh.get_entity(stk::mesh::Element,5);
      if (p_rank == 0) {
        STKUNIT_ASSERT_TRUE( entity_1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_2 != NULL );
        STKUNIT_ASSERT_TRUE( entity_3 == NULL );
        STKUNIT_ASSERT_TRUE( entity_4 != NULL );
        STKUNIT_ASSERT_TRUE( entity_5 != NULL );
        STKUNIT_EXPECT_TRUE( 0 == entity_1->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_2->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_4->owner_rank() );
        STKUNIT_EXPECT_TRUE( 0 == entity_5->owner_rank() );
      }
      else if (p_rank == 1) {
        STKUNIT_ASSERT_TRUE( entity_1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_2 != NULL );
        STKUNIT_ASSERT_TRUE( entity_3 != NULL );
        STKUNIT_ASSERT_TRUE( entity_4 != NULL );
        STKUNIT_ASSERT_TRUE( entity_5 != NULL );
        STKUNIT_EXPECT_TRUE( 0 == entity_1->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_2->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_3->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_4->owner_rank() );
        STKUNIT_EXPECT_TRUE( 0 == entity_5->owner_rank() );
      }
      else {
        STKUNIT_EXPECT_TRUE( entity_1 == NULL );
        STKUNIT_EXPECT_TRUE( entity_2 == NULL );
        STKUNIT_EXPECT_TRUE( entity_3 == NULL );
        STKUNIT_EXPECT_TRUE( entity_4 == NULL );
        STKUNIT_EXPECT_TRUE( entity_5 == NULL );
      }
    }
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, disjoint_parallel_5x1x1_skin )
  {
  }

} // end namespace

