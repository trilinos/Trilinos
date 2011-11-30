/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>
#include <map>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/baseImpl/EntityRepository.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>

using stk::ParallelMachine;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::PartVector;
using stk::mesh::PartRelation;
using stk::mesh::EntityKey;
using stk::mesh::Entity;
using stk::mesh::Bucket;
using stk::mesh::impl::PartRepository;
using stk::mesh::impl::EntityRepository;

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

#ifndef NDEBUG
  STKUNIT_ASSERT_THROW( EntityKey( ~0u , 1 ) , std::logic_error );
  STKUNIT_ASSERT_THROW( EntityKey( 0 , ~stk::mesh::EntityKey::raw_key_type(0) ) , std::logic_error );
#endif // NDEBUG
}


STKUNIT_UNIT_TEST(UnitTestEntity,testEntityRepository)
{
  //Test Entity repository - covering EntityRepository.cpp/hpp
  const int spatial_dimension = 3;
  MetaData meta(stk::mesh::fem::entity_rank_names(spatial_dimension));
  Part & part = meta.declare_part( "another part");
  MPI_Barrier( MPI_COMM_WORLD );
  ParallelMachine pm = MPI_COMM_WORLD;
  BulkData bulk( meta , pm, 200 );
  const int rank = stk::parallel_machine_rank( pm );
  const int size = stk::parallel_machine_size( pm );
  std::vector<stk::mesh::Part *>  add_part;
  meta.commit();

  // Bail if not parallel. This test involves inducing errorneous conditions that
  // are only checked-for in parallel.
  if (size == 1) {
    return;
  }

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

  STKUNIT_ASSERT_FALSE(e.erase_ghosting(elem, ghost));

  const stk::mesh::EntityCommInfo comm_info( ghost.ordinal() , 0 );

  STKUNIT_ASSERT_FALSE(e.erase_comm_info(elem, comm_info));

  STKUNIT_ASSERT(e.insert_comm_info(elem, comm_info));

  //Coverage of verfify_parallel_attributes in BulkDataParallelVerify.cpp
  //for owned_closure = 1 AND recv_ghost = 1.
  //Also uses pack and unpack in DataTraits.cpp, DataTraitsClass.hpp and DataTraitsEnum.hpp
  STKUNIT_ASSERT_THROW(bulk.modification_end(), std::runtime_error);

  bulk.modification_begin();


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

  // stk::mesh::impl::EntityRepository::EntityMap eMap;
  stk::mesh::Entity & elem2  =  bulk.declare_entity( 3 , new_id+8 ,  add_part );
  stk::mesh::Entity & elem3  =  bulk.declare_entity( 3 , new_id+9 ,  add_part );
  stk::mesh::Entity & elem4  =  bulk.declare_entity( 3 , new_id+10 ,  add_part );

  e.internal_create_entity( stk::mesh::EntityKey( 3, new_id+8 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, new_id+9 ));
  e.internal_create_entity( stk::mesh::EntityKey( 3, new_id+10 ));

  typedef std::map<EntityKey,Entity*> EntityMap;
  EntityMap entity_map_array;

  entity_map_array[stk::mesh::EntityKey( 3, new_id+8 )] = &elem2;
  entity_map_array[stk::mesh::EntityKey( 3, new_id+9 )] = &elem3;
  entity_map_array[stk::mesh::EntityKey( 3, new_id+10 )] = &elem4;

  //Coverage of destroy_later in EntityRepository.cpp
  Bucket *nil_bucket =  bulk.buckets(3)[0];
  e.destroy_later(elem2, nil_bucket);
  //Call a second time for more coverage
  STKUNIT_ASSERT_THROW(e.destroy_later(elem2, nil_bucket), std::runtime_error);

  //Coverage of !comm_mesh_verify_parallel_consistency in BulkDataEndSync.cpp
  //in internal_modification_end function
  Bucket *nil_bucket2 =  bulk.buckets(0)[0];

  STKUNIT_ASSERT ( nil_bucket2 != NULL);

  e.destroy_later(elem3, nil_bucket2);

  STKUNIT_ASSERT_THROW(bulk.modification_end(), std::runtime_error);

  bulk.modification_begin();

  STKUNIT_ASSERT_THROW(e.destroy_later(elem2, nil_bucket), std::runtime_error);

  STKUNIT_ASSERT_THROW(bulk.modification_end(), std::runtime_error);

  bulk.modification_begin();
  STKUNIT_ASSERT_THROW(bulk.modification_end(), std::runtime_error);
}

//----------------------------------------------------------------------
}//namespace <anonymous>

