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
using stk::mesh::Field;


typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
enum { SpaceDim = 3 };

class UnitTestEntity {
public:
  UnitTestEntity(ParallelMachine pm);
  ~UnitTestEntity() {}

  const int spatial_dimension;
  MetaData meta;
  TopologicalMetaData top;
  Part & part; 
  Part & part_left;
  Part & part_right;
  BulkData bulk;
  std::vector<stk::mesh::Part *>  add_part;
  EntityKey key_bad_zero;
  EntityKey key_good_0_1;
  EntityKey key_good_1_1;
  EntityKey key_good_2_10;
  EntityKey key_order_1_12;
  EntityKey key_order_2_10;
  const int rank;
  const int size;
  VectorFieldType   & coordinates_field;
  ScalarFieldType   & temperature_field;
  ScalarFieldType   & volume_field;



};

UnitTestEntity::UnitTestEntity(ParallelMachine pm)
  :spatial_dimension(3)
  ,meta( TopologicalMetaData::entity_rank_names(spatial_dimension) )
  ,top( meta, spatial_dimension )
  ,part( meta.declare_part( "another part") )
  ,part_left (top.declare_part<shards::Hexahedron<8> >( "block_left" ))
  ,part_right (top.declare_part<shards::Hexahedron<8> >( "block_right" ))
  ,bulk( meta , pm, 100 )
  ,add_part()
  ,key_bad_zero (EntityKey())
  ,key_good_0_1 (EntityKey( 0 , 1 ))
  ,key_good_1_1 (EntityKey( 1 , 1 ))
  ,key_good_2_10 (EntityKey( 2 , 10 ))
  ,key_order_1_12 (EntityKey( 1 , 12 ))
  ,key_order_2_10 (EntityKey( 2 , 10 ))
  ,rank( bulk.parallel_rank() )
  ,size( bulk.parallel_size() )
  ,coordinates_field( meta.declare_field< VectorFieldType >( "coordinates" ))
  ,temperature_field( meta.declare_field< ScalarFieldType >( "temperature" ))
  ,volume_field( meta.declare_field< ScalarFieldType >( "volume" ))

{
  meta.commit();
}



namespace {

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST(UnitTestEntity,testEntityKey)
{
  UnitTestEntity uent(MPI_COMM_WORLD);

//   EntityId val = entity_id(key_good_2_10) ;

//   EntityKey good_key_gone_bad;
//   entity_good_key_gone_bad.key = val ;

  STKUNIT_ASSERT( ! entity_key_valid( uent.key_bad_zero ) );
  STKUNIT_ASSERT(   entity_key_valid( uent.key_good_0_1 ) );
  STKUNIT_ASSERT(   entity_key_valid( uent.key_good_1_1 ) );
  STKUNIT_ASSERT(   entity_key_valid( uent.key_good_2_10 ) );

  STKUNIT_ASSERT( 0  == entity_rank(uent.key_good_0_1));
  STKUNIT_ASSERT( 1  == entity_rank(uent.key_good_1_1) );
  STKUNIT_ASSERT( 2  == entity_rank(uent.key_good_2_10) );
  STKUNIT_ASSERT( 1  == entity_id(uent.key_good_0_1) );
  STKUNIT_ASSERT( 1  == entity_id(uent.key_good_1_1) );
  STKUNIT_ASSERT( 10 == entity_id(uent.key_good_2_10) );


  STKUNIT_ASSERT( uent.key_order_1_12 < uent.key_order_2_10);
  STKUNIT_ASSERT( !(uent.key_order_1_12 > uent.key_order_2_10));

//  STKUNIT_ASSERT( ! entity_key_valid( good_key_gone_bad ) );

//   std::cout.unsetf( std::ios::dec);
//   std::cout.setf( std::ios::hex);
//   std::cout << "TEST entity_key_type "
//             << ", key_good_2_10 = " << key_good_2_10.key
//             << ", good_key_gone_bad = " << good_key_gone_bad.key
//             << std::endl ;
//   std::cout.unsetf( std::ios::hex);
//   std::cout.setf( std::ios::dec);

  //uent.key01 = uent.key_default;

  STKUNIT_ASSERT_THROW( EntityKey( ~0u , 1 ) , std::runtime_error );
  STKUNIT_ASSERT_THROW( EntityKey( 0 , ~stk::mesh::EntityKey::raw_key_type(0) ) , std::runtime_error );
}


STKUNIT_UNIT_TEST(UnitTestEntity,testEntityRepository)
{
  //Test Entity repository - covering EntityRepository.cpp/hpp
  UnitTestEntity uent(MPI_COMM_WORLD);


  uent.add_part.push_back ( &uent.part );

  uent.bulk.modification_begin();

  int id_base = 0;
  for ( id_base = 0 ; id_base < 97 ; ++id_base )
  {
    int new_id = uent.size * id_base + uent.rank;
    uent.bulk.declare_entity( 0 , new_id+1 , uent.add_part );
  }

  int new_id = uent.size * (++id_base) + uent.rank;
  stk::mesh::Entity & elem  = uent.bulk.declare_entity( 3 , new_id+1 , uent.add_part );

  //new_id = uent.size * (++id_base) + uent.rank;
  // stk::mesh::Entity & elem2  = uent.bulk.declare_entity( 3 , new_id+1 , uent.add_part );

  stk::mesh::impl::EntityRepository e;

  e.comm_clear( elem );

  e.comm_clear_ghosting( elem );

  const stk::mesh::Ghosting & ghost = uent.bulk.shared_aura();

  uent.bulk.modification_end();

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

