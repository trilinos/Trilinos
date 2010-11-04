/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>


using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::TopologicalMetaData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::impl::PartRepository;
using stk::ParallelMachine;

using stk::mesh::EntityKey;
using stk::mesh::Entity;


namespace {

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestEntity,testEntityKey)
{
  
  EntityKey key_bad_zero = EntityKey();
  EntityKey key_good_0_1 = EntityKey( 0 , 1 );
  EntityKey key_good_1_1 = EntityKey( 1 , 1 );
  EntityKey key_good_2_10 = EntityKey( 2 , 10);
  EntityKey key_order_1_12 = EntityKey( 1 , 12 );
  EntityKey key_order_2_10 = EntityKey( 2 , 10 );


  STKUNIT_ASSERT( ! entity_key_valid(  key_bad_zero ) );
  STKUNIT_ASSERT(   entity_key_valid(  key_good_0_1 ) );
  STKUNIT_ASSERT(   entity_key_valid(  key_good_1_1 ) );
  STKUNIT_ASSERT(   entity_key_valid(  key_good_2_10 ) );

  STKUNIT_ASSERT( 0  == entity_rank( key_good_0_1));
  STKUNIT_ASSERT( 1  == entity_rank( key_good_1_1) );
  STKUNIT_ASSERT( 2  == entity_rank( key_good_2_10) );
  STKUNIT_ASSERT( 1  == entity_id( key_good_0_1) );
  STKUNIT_ASSERT( 1  == entity_id( key_good_1_1) );
  STKUNIT_ASSERT( 10 == entity_id( key_good_2_10) );


  STKUNIT_ASSERT(  key_order_1_12 <  key_order_2_10);
  STKUNIT_ASSERT( !( key_order_1_12 >  key_order_2_10));


  STKUNIT_ASSERT_THROW( EntityKey( ~0u , 1 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( EntityKey( 0 , ~stk::mesh::EntityKey::raw_key_type(0) ) , std::runtime_error );
}


STKUNIT_UNIT_TEST(UnitTestEntity,testEntityRepository)
{
  //Test Entity repository - covering EntityRepository.cpp/hpp
  const int spatial_dimension = 3;
  MetaData meta(TopologicalMetaData::entity_rank_names(spatial_dimension));
  Part & part = meta.declare_part( "another part"); 
  MPI_Barrier( MPI_COMM_WORLD );
  ParallelMachine pm = MPI_COMM_WORLD;
  BulkData bulk( meta , pm, 100 );
  const int rank = stk::parallel_machine_rank( pm );
  const int size = stk::parallel_machine_size( pm );
  std::vector<stk::mesh::Part *>  add_part;
  meta.commit();

  add_part.push_back ( & part );

  bulk.modification_begin();

  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id =  size * id_base +  rank;
     bulk.declare_entity( 0 , new_id+1 ,  add_part );
  }

  int new_id =  size * (++id_base) +  rank;
  stk::mesh::Entity & elem  =  bulk.declare_entity( 3 , new_id+1 ,  add_part );

  stk::mesh::impl::EntityRepository e;

  e.comm_clear( elem );

  e.comm_clear_ghosting( elem );

  const stk::mesh::Ghosting & ghost =  bulk.shared_aura();

  bulk.modification_end();

  STKUNIT_ASSERT_FALSE(e.erase_ghosting(elem, ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  STKUNIT_ASSERT_FALSE(e.erase_comm_info(elem, comm_info));

  STKUNIT_ASSERT(e.insert_comm_info(elem, comm_info));

  //Checking internal_create_entity

  e.internal_create_entity( stk::mesh::EntityKey( 3, 2 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, 5 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, 7 ));

  //Checking get_entity with invalid key - no rank or id
  {
    STKUNIT_ASSERT_THROW( 
        e.get_entity(stk::mesh::EntityKey()),
        std::runtime_error
        );
  }
}

//----------------------------------------------------------------------
}//namespace <anonymous>

