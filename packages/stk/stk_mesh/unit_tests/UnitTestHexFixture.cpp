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

  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_x )
  {
    const unsigned NX = 3;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(1,0,0) == 2 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(2,0,0) == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_x )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(1,0,0) == 2 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(2,0,0) == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_y )
  {
    const unsigned NX = 1;
    const unsigned NY = 3;
    const unsigned NZ = 1;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,1,0) == 2 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,2,0) == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_y )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,1,0) == 4 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,2,0) == 7 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_1d_z )
  {
    const unsigned NX = 1;
    const unsigned NY = 1;
    const unsigned NZ = 3;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,1) == 2 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,2) == 3 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_z )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,1) == 10 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,2) == 19 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, elem_ids_3d_diag )
  {
    const unsigned NX = 3;
    const unsigned NY = 3;
    const unsigned NZ = 3;
    stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
    hf.meta_data.commit();
    hf.generate_mesh();
    STKUNIT_EXPECT_TRUE( hf.elem_id(0,0,0) == 1 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(1,1,1) == 14 );
    STKUNIT_EXPECT_TRUE( hf.elem_id(2,2,2) == 27 );
  }
  STKUNIT_UNIT_TEST( UnitTestHexFixture, trivial_parallel_2 )
  {
    const unsigned NX = 2;
    const unsigned NY = 1;
    const unsigned NZ = 1;
    // map< processor, vector of element ids >
    std::map<unsigned,std::vector<stk::mesh::EntityId> > parallel_distribution;
    {
      std::vector< stk::mesh::EntityId> element_ids;
      element_ids.push_back(1);
      parallel_distribution[0] = element_ids;
      element_ids[0] = 2;
      parallel_distribution[1] = element_ids;
    }
    const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (2 <= p_size) {
      stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
      hf.meta_data.commit();
      if (p_rank <= 1) {
        hf.generate_mesh(parallel_distribution[p_rank]);
      }
      else {
        std::vector<stk::mesh::EntityId> empty_vector;
        hf.generate_mesh( empty_vector ) ;
      }
      stk::mesh::BulkData & mesh = hf.bulk_data;
      stk::mesh::TopologicalMetaData & top = hf.top_data;
      // Verify element_id 1 is owned by proc 0
      // Verify element_id 2 is owned by proc 1
      stk::mesh::Entity * entity_1 = mesh.get_entity(top.element_rank,1);
      stk::mesh::Entity * entity_2 = mesh.get_entity(top.element_rank,2);
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
  STKUNIT_UNIT_TEST( UnitTestHexFixture, disjoint_parallel_psizex1x1 )
  {
    // layout: p_size x 1 x 1 hex mesh
    // elements:
    // [ e_1, e_2, e_3, ..., e_n ]
    // processors:
    // [ p_0, p_1, p_2, ..., p_{n-1} ]
    //
    const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (p_size == 4) {
      const unsigned NX = 4;
      //const unsigned NX = p_size;
      const unsigned NY = 1;
      const unsigned NZ = 1;
      // map< processor, vector of element ids >
      std::map<unsigned,std::vector<stk::mesh::EntityId> > parallel_distribution;
      for (unsigned p=0 ; p < p_size ; ++p) {
        std::vector< stk::mesh::EntityId> element_ids;
        element_ids.push_back(p+1); // element id's start at 1
        parallel_distribution[p] = element_ids;
      }
      stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);
      hf.meta_data.commit();
      hf.generate_mesh(parallel_distribution[p_rank]);
      stk::mesh::BulkData & mesh = hf.bulk_data;
      stk::mesh::TopologicalMetaData & top = hf.top_data;
      if (p_rank == 0) {
        stk::mesh::Entity * entity_1 = mesh.get_entity(top.element_rank,1);
        stk::mesh::Entity * entity_2 = mesh.get_entity(top.element_rank,2);
        STKUNIT_ASSERT_TRUE( entity_1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_2 != NULL );
        STKUNIT_EXPECT_TRUE( 0 == entity_1->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_2->owner_rank() );
      }
      else if (p_rank < p_size-1) {
        stk::mesh::Entity * entity_im1 = mesh.get_entity(top.element_rank,1+p_rank-1);
        stk::mesh::Entity * entity_i = mesh.get_entity(top.element_rank,1+p_rank);
        stk::mesh::Entity * entity_ip1 = mesh.get_entity(top.element_rank,1+p_rank+1);
        STKUNIT_ASSERT_TRUE( entity_im1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_i != NULL );
        STKUNIT_ASSERT_TRUE( entity_ip1 != NULL );
        STKUNIT_EXPECT_TRUE( p_rank-1 == entity_im1->owner_rank() );
        STKUNIT_EXPECT_TRUE( p_rank   == entity_i->owner_rank() );
        STKUNIT_EXPECT_TRUE( p_rank+1 == entity_ip1->owner_rank() );
      }
      else if (p_rank == p_size-1) {
        stk::mesh::Entity * entity_im1 = mesh.get_entity(top.element_rank,1+p_rank-1);
        stk::mesh::Entity * entity_i = mesh.get_entity(top.element_rank,1+p_rank);
        STKUNIT_ASSERT_TRUE( entity_im1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_i != NULL );
        STKUNIT_EXPECT_TRUE( (p_rank-1) == entity_im1->owner_rank() );
        STKUNIT_EXPECT_TRUE( p_rank == entity_i->owner_rank() );
      }
      else { // invalid
        STKUNIT_ASSERT_TRUE(false);
      }
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
    const unsigned NX = 4;
    const unsigned NY = 2;
    const unsigned NZ = 1;
    // map< processor, vector of element ids >
    std::map<unsigned,std::vector<stk::mesh::EntityId> > parallel_distribution;
    {
      std::vector< stk::mesh::EntityId> element_ids;
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
    const unsigned p_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    const unsigned p_size = stk::parallel_machine_size(MPI_COMM_WORLD);
    if (2 <= p_size) {
      stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);

      hf.meta_data.commit();
      if (p_rank <= 1) {
        hf.generate_mesh(parallel_distribution[p_rank]);
      }
      else {
        std::vector<stk::mesh::EntityId> empty_vector;
        hf.generate_mesh( empty_vector ) ;
      }

      stk::mesh::BulkData & mesh = hf.bulk_data;
      stk::mesh::TopologicalMetaData & top = hf.top_data;
      // Verify element_id 1 is owned by proc 0
      // Verify element_id 2 is owned by proc 1
      stk::mesh::Entity * entity_1 = mesh.get_entity(top.element_rank,1);
      stk::mesh::Entity * entity_2 = mesh.get_entity(top.element_rank,2);
      stk::mesh::Entity * entity_3 = mesh.get_entity(top.element_rank,3);
      stk::mesh::Entity * entity_4 = mesh.get_entity(top.element_rank,4);
      stk::mesh::Entity * entity_5 = mesh.get_entity(top.element_rank,5);
      stk::mesh::Entity * entity_6 = mesh.get_entity(top.element_rank,6);
      stk::mesh::Entity * entity_7 = mesh.get_entity(top.element_rank,7);
      stk::mesh::Entity * entity_8 = mesh.get_entity(top.element_rank,8);
      if (p_rank == 0) {
        STKUNIT_ASSERT_TRUE( entity_1 != NULL );
        STKUNIT_ASSERT_TRUE( entity_2 != NULL );
        STKUNIT_ASSERT_TRUE( entity_3 != NULL );
        STKUNIT_ASSERT_TRUE( entity_4 == NULL );
        STKUNIT_ASSERT_TRUE( entity_5 != NULL );
        STKUNIT_ASSERT_TRUE( entity_6 != NULL );
        STKUNIT_ASSERT_TRUE( entity_7 != NULL );
        STKUNIT_ASSERT_TRUE( entity_8 == NULL );
        STKUNIT_EXPECT_TRUE( 0 == entity_1->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_2->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_3->owner_rank() );
        //STKUNIT_EXPECT_TRUE( 1 == entity_4->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_5->owner_rank() );
        STKUNIT_EXPECT_TRUE( 0 == entity_6->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_7->owner_rank() );
        //STKUNIT_EXPECT_TRUE( 1 == entity_8->owner_rank() );
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
        STKUNIT_EXPECT_TRUE( 0 == entity_1->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_2->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_3->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_4->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_5->owner_rank() );
        STKUNIT_EXPECT_TRUE( 0 == entity_6->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_7->owner_rank() );
        STKUNIT_EXPECT_TRUE( 1 == entity_8->owner_rank() );
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
    std::map<unsigned,std::vector<stk::mesh::EntityId> > parallel_distribution;
    {
      std::vector< stk::mesh::EntityId> element_ids;
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
      stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);

      hf.meta_data.commit();
      if (p_rank <= 1) {
        hf.generate_mesh(parallel_distribution[p_rank]);
      }
      else {
        std::vector<stk::mesh::EntityId> empty_vector;
        hf.generate_mesh( empty_vector ) ;
      }

      stk::mesh::BulkData & mesh = hf.bulk_data;
      stk::mesh::TopologicalMetaData & top = hf.top_data;
      // Verify element_id 1 is owned by proc 0
      // Verify element_id 2 is owned by proc 1
      stk::mesh::Entity * entity_1 = mesh.get_entity(top.element_rank,1);
      stk::mesh::Entity * entity_2 = mesh.get_entity(top.element_rank,2);
      stk::mesh::Entity * entity_3 = mesh.get_entity(top.element_rank,3);
      stk::mesh::Entity * entity_4 = mesh.get_entity(top.element_rank,4);
      stk::mesh::Entity * entity_5 = mesh.get_entity(top.element_rank,5);
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
