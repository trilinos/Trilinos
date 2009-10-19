
#include <iostream>
#include <sstream>

#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

#include <unit_tests/UnitTestBulkData.hpp>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

void UnitTestBulkData::testChangeParts( ParallelMachine pm )
{
  static const char method[] =
    "stk::mesh::UnitTestBulkData::testChangeParts" ;

  std::cout << std::endl << method << std::endl ;

  const unsigned p_size = parallel_machine_size( pm );
  const unsigned p_rank = parallel_machine_rank( pm );

  // Meta data with entity ranks [0..9]
  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityRank_" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta( entity_names );

  Part & part_univ = meta.universal_part();
  Part & part_uses = meta.locally_used_part();
  Part & part_owns = meta.locally_owned_part();

  Part & part_A_0 = meta.declare_part( std::string("A_0") , 0 );
  Part & part_A_1 = meta.declare_part( std::string("A_1") , 1 );
  Part & part_A_2 = meta.declare_part( std::string("A_2") , 2 );
  Part & part_A_3 = meta.declare_part( std::string("A_3") , 3 );

  Part & part_B_0 = meta.declare_part( std::string("B_0") , 0 );
  // Part & part_B_1 = meta.declare_part( std::string("B_1") , 1 );
  Part & part_B_2 = meta.declare_part( std::string("B_2") , 2 );
  // Part & part_B_3 = meta.declare_part( std::string("B_3") , 3 );

  meta.commit();

  BulkData bulk( meta , pm , 100 );

  PartVector tmp(1);

  tmp[0] = & part_A_0 ;
  Entity & entity_0_1 = bulk.declare_entity(  0 , 1 , tmp );

  tmp[0] = & part_A_1 ;
  Entity & entity_1_1 = bulk.declare_entity(  1 , 1 , tmp );

  tmp[0] = & part_A_2 ;
  Entity & entity_2_1 = bulk.declare_entity(  2 , 1 , tmp );

  tmp[0] = & part_A_3 ;
  Entity & entity_3_1 = bulk.declare_entity( 3 , 1 , tmp );

  entity_0_1.bucket().supersets( tmp );
  STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
  STKUNIT_ASSERT( tmp[0] == & part_univ );
  STKUNIT_ASSERT( tmp[1] == & part_uses );
  STKUNIT_ASSERT( tmp[2] == & part_owns );
  STKUNIT_ASSERT( tmp[3] == & part_A_0 );

  entity_1_1.bucket().supersets( tmp );
  STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
  STKUNIT_ASSERT( tmp[0] == & part_univ );
  STKUNIT_ASSERT( tmp[1] == & part_uses );
  STKUNIT_ASSERT( tmp[2] == & part_owns );
  STKUNIT_ASSERT( tmp[3] == & part_A_1 );

  entity_2_1.bucket().supersets( tmp );
  STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
  STKUNIT_ASSERT( tmp[0] == & part_univ );
  STKUNIT_ASSERT( tmp[1] == & part_uses );
  STKUNIT_ASSERT( tmp[2] == & part_owns );
  STKUNIT_ASSERT( tmp[3] == & part_A_2 );

  entity_3_1.bucket().supersets( tmp );
  STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
  STKUNIT_ASSERT( tmp[0] == & part_univ );
  STKUNIT_ASSERT( tmp[1] == & part_uses );
  STKUNIT_ASSERT( tmp[2] == & part_owns );
  STKUNIT_ASSERT( tmp[3] == & part_A_3 );

  {
    tmp.resize(1);
    tmp[0] = & part_A_0 ;
    bulk.change_entity_parts( entity_0_1 , tmp );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
  }

  { // Add a new part:
    tmp.resize(1);
    tmp[0] = & part_B_0 ;
    bulk.change_entity_parts( entity_0_1 , tmp );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(5) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
    STKUNIT_ASSERT( tmp[4] == & part_B_0 );
  }

  { // Remove the part just added:
    tmp.resize(1);
    tmp[0] = & part_B_0 ;
    bulk.change_entity_parts( entity_0_1 , PartVector() , tmp );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
  }

  { // Relationship induced membership:
    bulk.declare_relation( entity_1_1 , entity_0_1 , 0 );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(5) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
    STKUNIT_ASSERT( tmp[4] == & part_A_1 );
  }

  { // Remove relationship induced membership:
    bulk.destroy_relation( entity_1_1 , entity_0_1 );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
  }

  { // Add a new part:
    tmp.resize(1);
    tmp[0] = & part_B_2 ;
    bulk.change_entity_parts( entity_2_1 , tmp );
    entity_2_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(5) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_2 );
    STKUNIT_ASSERT( tmp[4] == & part_B_2 );
  }

  { // Relationship induced membership:
    bulk.declare_relation( entity_2_1 , entity_0_1 , 0 );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(6) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
    STKUNIT_ASSERT( tmp[4] == & part_A_2 );
    STKUNIT_ASSERT( tmp[5] == & part_B_2 );
  }

  { // Remove relationship induced membership:
    bulk.destroy_relation( entity_2_1 , entity_0_1 );
    entity_0_1.bucket().supersets( tmp );
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
  }

  bulk.modification_end();

  //------------------------------
  // Now the parallel fun.  Existing entities should be shared
  // by all processes since they have the same identifiers.
  // They should also have the same parts.

  entity_0_1.bucket().supersets( tmp );
  if ( entity_0_1.owner_rank() == p_rank ) {
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
  }
  else {
    STKUNIT_ASSERT_EQUAL( size_t(3) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_A_0 );
  }

  entity_2_1.bucket().supersets( tmp );
  if ( entity_2_1.owner_rank() == p_rank ) {
    STKUNIT_ASSERT_EQUAL( size_t(5) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_2 );
    STKUNIT_ASSERT( tmp[4] == & part_B_2 );
  }
  else {
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_A_2 );
    STKUNIT_ASSERT( tmp[3] == & part_B_2 );
  }

  STKUNIT_ASSERT_EQUAL( size_t(p_size - 1) , entity_0_1.sharing().size() );
  STKUNIT_ASSERT_EQUAL( size_t(p_size - 1) , entity_1_1.sharing().size() );
  STKUNIT_ASSERT_EQUAL( size_t(p_size - 1) , entity_2_1.sharing().size() );
  STKUNIT_ASSERT_EQUAL( size_t(p_size - 1) , entity_3_1.sharing().size() );

  bulk.modification_begin();

  // Add a new part on the owning process:

  int ok_to_modify = entity_0_1.owner_rank() == p_rank ;

  try {
    tmp.resize(1);
    tmp[0] = & part_B_0 ;
    bulk.change_entity_parts( entity_0_1 , tmp );
    STKUNIT_ASSERT( ok_to_modify );
  }
  catch( const std::exception & x ) {
    STKUNIT_ASSERT( ! ok_to_modify );
  }

  entity_0_1.bucket().supersets( tmp );
  if ( entity_0_1.owner_rank() == p_rank ) {
    STKUNIT_ASSERT_EQUAL( size_t(5) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
    STKUNIT_ASSERT( tmp[4] == & part_B_0 );
  }
  else {
    STKUNIT_ASSERT_EQUAL( size_t(3) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_A_0 );
  }

  bulk.modification_end();

  entity_0_1.bucket().supersets( tmp );
  if ( entity_0_1.owner_rank() == p_rank ) {
    STKUNIT_ASSERT_EQUAL( size_t(5) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_owns );
    STKUNIT_ASSERT( tmp[3] == & part_A_0 );
    STKUNIT_ASSERT( tmp[4] == & part_B_0 );
  }
  else {
    STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
    STKUNIT_ASSERT( tmp[0] == & part_univ );
    STKUNIT_ASSERT( tmp[1] == & part_uses );
    STKUNIT_ASSERT( tmp[2] == & part_A_0 );
    STKUNIT_ASSERT( tmp[3] == & part_B_0 );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void UnitTestBulkData::testChangeParts_loop( ParallelMachine pm )
{
  enum { nPerProc = 10 };
  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  std::vector<EntityId> node_ids , edge_ids ;

  MetaData meta( fem_entity_type_names() );

  Part & part_univ = meta.universal_part();
  Part & part_uses = meta.locally_used_part();
  Part & part_owns = meta.locally_owned_part();

  PartVector edge_parts( nLocalEdge );
  for ( unsigned i = 0 ; i < nLocalEdge ; ++i ) {
    std::ostringstream name ;
    name << "EdgePart_" << i ;
    edge_parts[i] = & meta.declare_part( name.str() , Edge );
  }

  Part & edge_part_extra = meta.declare_part( "EdgePart_Extra" , Edge );

  meta.commit();

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used( meta.locally_used_part() );
  Selector select_all(  meta.universal_part() );

  std::vector<unsigned> local_count ;

  //------------------------------
  {
    BulkData bulk( meta , pm , 100 );

    generate_loop( bulk, edge_parts , false /* no aura */, nPerProc,
                   node_ids, edge_ids );

    PartVector tmp ;
    for ( unsigned i = 0 ; i < nLocalEdge ; ++i ) {
      const unsigned n = i + nPerProc * p_rank ;
      Entity * const edge = bulk.get_entity( 1 , edge_ids[n] );
      STKUNIT_ASSERT( edge != NULL );
      edge->bucket().supersets( tmp );
      STKUNIT_ASSERT( size_t(4) == tmp.size() );
      STKUNIT_ASSERT( tmp[0] == & part_univ );
      STKUNIT_ASSERT( tmp[1] == & part_uses );
      STKUNIT_ASSERT( tmp[2] == & part_owns );
      STKUNIT_ASSERT( tmp[3] == edge_parts[ n % edge_parts.size() ] );
    }

    for ( unsigned i = 0 ; i < nLocalNode ; ++i ) {
      const unsigned n = ( i + nPerProc * p_rank ) % node_ids.size();
      const unsigned e0 = n ;
      const unsigned e1 = ( n + edge_ids.size() - 1 ) % edge_ids.size();
      const unsigned ns = edge_parts.size();
      const unsigned n0 = e0 % ns ;
      const unsigned n1 = e1 % ns ;
      Part * const epart_0 = edge_parts[ n0 < n1 ? n0 : n1 ];
      Part * const epart_1 = edge_parts[ n0 < n1 ? n1 : n0 ];

      Entity * const node = bulk.get_entity( 0 , node_ids[n] );
      STKUNIT_ASSERT( node != NULL );
      node->bucket().supersets( tmp );
      if ( node->owner_rank() == p_rank ) {
        STKUNIT_ASSERT( size_t(5) == tmp.size() );
        STKUNIT_ASSERT( tmp[0] == & part_univ );
        STKUNIT_ASSERT( tmp[1] == & part_uses );
        STKUNIT_ASSERT( tmp[2] == & part_owns );
        STKUNIT_ASSERT( tmp[3] == epart_0 );
        STKUNIT_ASSERT( tmp[4] == epart_1 );
      }
      else {
        STKUNIT_ASSERT( size_t(4) == tmp.size() );
        STKUNIT_ASSERT( tmp[0] == & part_univ );
        STKUNIT_ASSERT( tmp[1] == & part_uses );
        STKUNIT_ASSERT( tmp[2] == epart_0 );
        STKUNIT_ASSERT( tmp[3] == epart_1 );
      }
    }

    bulk.modification_begin();

    if ( 0 == p_rank ) {

      for ( unsigned i = 0 ; i < nLocalEdge ; ++i ) {
        const unsigned n = i + nPerProc * p_rank ;

        PartVector add(1); add[0] = & edge_part_extra ;
        PartVector rem(1); rem[0] = edge_parts[ n % edge_parts.size() ];

        Entity * const edge = bulk.get_entity( 1 , edge_ids[n] );
        bulk.change_entity_parts( *edge , add , rem );
        edge->bucket().supersets( tmp );
        STKUNIT_ASSERT_EQUAL( size_t(4) , tmp.size() );
        STKUNIT_ASSERT( tmp[0] == & part_univ );
        STKUNIT_ASSERT( tmp[1] == & part_uses );
        STKUNIT_ASSERT( tmp[2] == & part_owns );
        STKUNIT_ASSERT( tmp[3] == & edge_part_extra );
      }
    }

    bulk.modification_end();

    for ( unsigned i = 0 ; i < nLocalNode ; ++i ) {
      const unsigned n = ( i + nPerProc * p_rank ) % node_ids.size();
      const unsigned e0 = n ;
      const unsigned e1 = ( n + edge_ids.size() - 1 ) % edge_ids.size();
      const unsigned ns = edge_parts.size();
      const unsigned n0 = e0 % ns ;
      const unsigned n1 = e1 % ns ;
      Part * ep_0 = e0 < nLocalEdge ? & edge_part_extra : edge_parts[n0] ;
      Part * ep_1 = e1 < nLocalEdge ? & edge_part_extra : edge_parts[n1] ;

      Part * epart_0 = ep_0->mesh_meta_data_ordinal() < ep_1->mesh_meta_data_ordinal() ? ep_0 : ep_1 ;
      Part * epart_1 = ep_0->mesh_meta_data_ordinal() < ep_1->mesh_meta_data_ordinal() ? ep_1 : ep_0 ;

      const size_t n_edge_part_count = epart_0 == epart_1 ? 1 : 2 ;

      Entity * const node = bulk.get_entity( 0 , node_ids[n] );
      STKUNIT_ASSERT( node != NULL );
      node->bucket().supersets( tmp );
      if ( node->owner_rank() == p_rank ) {
        STKUNIT_ASSERT_EQUAL( n_edge_part_count + 3 , tmp.size() );
        STKUNIT_ASSERT( tmp[0] == & part_univ );
        STKUNIT_ASSERT( tmp[1] == & part_uses );
        STKUNIT_ASSERT( tmp[2] == & part_owns );
        STKUNIT_ASSERT( tmp[3] == epart_0 );
        if ( 1 < n_edge_part_count ) STKUNIT_ASSERT( tmp[4] == epart_1 );
      }
      else {
        STKUNIT_ASSERT_EQUAL( n_edge_part_count + 2 , tmp.size() );
        STKUNIT_ASSERT( tmp[0] == & part_univ );
        STKUNIT_ASSERT( tmp[1] == & part_uses );
        STKUNIT_ASSERT( tmp[2] == epart_0 );
        if ( 1 < n_edge_part_count ) STKUNIT_ASSERT( tmp[3] == epart_1 );
      }
    }
  }
}

} // namespace mesh
} // namespace stk

