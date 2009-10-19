
#include <sstream>

#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

#include <unit_tests/UnitTestBulkData.hpp>

STKUNIT_UNIT_TEST(UnitTestingOfBulkData, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testBulkData( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testChangeOwner_nodes( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testChangeOwner_loop( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testChangeOwner_box( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testCreateMore_error( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testChangeParts( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testChangeParts_loop( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testDestroy_nodes( MPI_COMM_WORLD );
  stk::mesh::UnitTestBulkData::testDestroy_loop( MPI_COMM_WORLD );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

void UnitTestBulkData::modification_end(
  BulkData & mesh , bool aura )
{
  // Exactly match BulkData::modification_end
  // except for automatic regeneration of the ghosting aura.

  STKUNIT_ASSERT_EQUAL( false , mesh.m_sync_state );

  if ( 1 < mesh.m_parallel_size ) {
    mesh.internal_resolve_parallel_create_delete();
  }

  // Done with created and deleted entity lists.

  mesh.m_new_entities.clear();
   
  while ( ! mesh.m_del_entities.empty() ) {
    mesh.internal_destroy_entity( mesh.m_del_entities.back() );
    mesh.m_del_entities.pop_back();
  }  

  if ( 1 < mesh.m_parallel_size ) {
    if ( aura ) { mesh.internal_regenerate_shared_aura(); }
    mesh.internal_resolve_shared_membership();
  }

  mesh.internal_sort_bucket_entities();

  ++ mesh.m_sync_count ;

  mesh.m_sync_state = true ;
}

// Unit test the Part functionality in isolation:

void UnitTestBulkData::testBulkData( ParallelMachine pm )
{
  static const char method[] = "stk::mesh::UnitTestBulkData" ;

  std::cout << std::endl ;

  std::vector<std::string> entity_names(10);
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    std::ostringstream name ;
    name << "EntityType" << i ;
    entity_names[i] = name.str();
  }

  MetaData meta( entity_names );

  meta.commit();

  BulkData bulk( meta , pm , 100 );

  for ( size_t i = 0 ; i < 4 ; ++i ) {
    STKUNIT_ASSERT( i == bulk.synchronized_count() );
    STKUNIT_ASSERT( bulk.modification_end() );
    STKUNIT_ASSERT( bulk.modification_begin() );
  }

  std::vector<Part*> no_parts ;

  Entity * e[10] ;
  for ( size_t i = 0 ; i < 10 ; ++i ) {
    e[i] = & bulk.declare_entity(  i , 1 , no_parts );
  }

  STKUNIT_ASSERT( bulk.modification_end() );
  STKUNIT_ASSERT( bulk.modification_begin() );

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    STKUNIT_ASSERT( e[i] == bulk.get_entity( i , 1 ) );
  }

  bool ok = false ;
  try {
    bulk.declare_entity( 11 , 1 , no_parts );
  }
  catch( const std::exception & x ) {
    std::cout << method << " correctly caught: " << x.what() << std::endl ;
    ok = true ;
  }
  STKUNIT_ASSERT( ok );
  STKUNIT_ASSERT( bulk.modification_end() );

  ok = false ;
  try {
    bulk.declare_entity( 0 , 2 , no_parts );
  }
  catch( const std::exception & x ) {
    std::cout << method << " correctly caught: " << x.what() << std::endl ;
    ok = true ;
  }
  STKUNIT_ASSERT( ok );
}

//----------------------------------------------------------------------
// Testing for mesh entities without relations

void UnitTestBulkData::testChangeOwner_nodes( ParallelMachine pm )
{
  enum { nPerProc = 10 };
  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );

  MetaData meta( fem_entity_type_names() );

  const PartVector no_parts ;

  meta.commit();

  BulkData bulk( meta , pm , 100 );

  // Ids for all entities (all entities have type 0):

  std::vector<EntityId> ids( id_total );

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    ids[i] = i + 1;
  }

  // Declare just those entities in my range of ids:

  for ( unsigned i = id_begin ; i < id_end ; ++i ) {
    bulk.declare_entity( 0 , ids[i] , no_parts );
  }

  STKUNIT_ASSERT( bulk.modification_end() );

  // Verify that I only have entities in my range:

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    Entity * e = bulk.get_entity( 0 , ids[ i ] );
    if ( id_begin <= i && i < id_end ) {
      STKUNIT_ASSERT( NULL != e );
    }
    else {
      STKUNIT_ASSERT( NULL == e );
    }
  }

  // Test change owner no-op first:

  std::vector<EntityProc> change ;

  STKUNIT_ASSERT( bulk.modification_begin() );
  bulk.change_entity_owner( change );
  STKUNIT_ASSERT( bulk.modification_end() );

  for ( unsigned i = 0 ; i < id_total ; ++i ) {
    Entity * e = bulk.get_entity( 0 , ids[ i ] );
    if ( id_begin <= i && i < id_end ) {
      STKUNIT_ASSERT( NULL != e );
    }
    else {
      STKUNIT_ASSERT( NULL == e );
    }
  }

  // Can only test changing owner in parallel.

  if ( 1 < p_size ) {
    // Give my last two ids to the next process
    // Get the previous process' last two ids

    const unsigned p_give = ( p_rank + 1 ) % p_size ;
    const unsigned id_give = id_end - 2 ;
    const unsigned id_get  = ( id_begin + id_total - 2 ) % id_total ;

    STKUNIT_ASSERT( NULL != bulk.get_entity( 0 , ids[id_give] ) );
    STKUNIT_ASSERT( NULL != bulk.get_entity( 0 , ids[id_give+1] ) );
    STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[id_get] ) );
    STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[id_get+1] ) );
  
    change.resize(2);
    change[0].first = bulk.get_entity( 0 , ids[id_give] );
    change[0].second = p_give ;
    change[1].first = bulk.get_entity( 0 , ids[id_give+1] );
    change[1].second = p_give ;

    STKUNIT_ASSERT( bulk.modification_begin() );
    bulk.change_entity_owner( change );
    STKUNIT_ASSERT( bulk.modification_end() );

    STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[id_give] ) );
    STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[id_give+1] ) );
    STKUNIT_ASSERT( NULL != bulk.get_entity( 0 , ids[id_get] ) );
    STKUNIT_ASSERT( NULL != bulk.get_entity( 0 , ids[id_get+1] ) );
  }

  std::cout << std::endl
            << "P" << p_rank
            << ": UnitTestBulkData::testChangeOwner_nodes( NP = "
            << p_size << " ) SUCCESSFULL " << std::endl ;
}

//----------------------------------------------------------------------
// Testing for error using mesh entities without relations

void UnitTestBulkData::testCreateMore_error( ParallelMachine pm )
{
  enum { nPerProc = 10 };

  const unsigned p_size = parallel_machine_size( pm );
  const unsigned p_rank = parallel_machine_rank( pm );

  if ( 1 < p_size ) {

    std::cout << std::endl
              << "P" << p_rank
              << ": UnitTestBulkData::testCreateMore_error( NP = "
              << p_size << " ) TESTING FOR PARALLEL ERROR CATCHING"
              << std::endl ;
    std::cout.flush();

    const unsigned id_total = nPerProc * p_size ;
    const unsigned id_begin = nPerProc * p_rank ;
    const unsigned id_end   = nPerProc * ( p_rank + 1 );

    MetaData meta( fem_entity_type_names() );

    const PartVector no_parts ;

    meta.commit();

    BulkData bulk( meta , pm , 100 );

    // Ids for all entities (all entities have type 0):

    std::vector<EntityId> ids( id_total );

    for ( unsigned i = 0 ; i < id_total ; ++i ) { ids[i] = i + 1; }

    // Declare just those entities in my range of ids:

    for ( unsigned i = id_begin ; i < id_end ; ++i ) {
      bulk.declare_entity( 0 , ids[i] , no_parts );
    }

    STKUNIT_ASSERT( bulk.modification_end() );

    // Only one process create entities with previous process' last two ids

    const unsigned id_get  = ( id_begin + id_total - 2 ) % id_total ;

    STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[id_get] ) );
    STKUNIT_ASSERT( NULL == bulk.get_entity( 0 , ids[id_get+1] ) );

    STKUNIT_ASSERT( bulk.modification_begin() );

    if ( 1 == p_rank ) {
      bulk.declare_entity( 0 , ids[ id_get ] , no_parts );
      bulk.declare_entity( 0 , ids[ id_get + 1 ] , no_parts );
    }
  
    std::string error_msg ;
    bool exception_thrown = false ;
    try {
      bulk.modification_end();
    }
    catch( const std::exception & x ) {
      exception_thrown = true ;
      error_msg.append( x.what() );
      std::cerr.flush();
    }

    STKUNIT_ASSERT( exception_thrown );

    std::cout << std::endl
              << "P" << p_rank
              << ": UnitTestBulkData::testCreateMore_error( NP = "
              << p_size << " ) SUCCESSFULLY CAUGHT "
              << error_msg
              << std::endl ;
    std::cout.flush();
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void test_shift_loop( BulkData & mesh ,
                      const bool                     generate_aura ,
                      const unsigned                 nPerProc ,
                      const std::vector<EntityId> & node_ids ,
                      const std::vector<EntityId> & edge_ids )
{
  const unsigned p_rank = mesh.parallel_rank();
  const unsigned p_size = mesh.parallel_size();
  const unsigned id_total = nPerProc * p_size ;
  const unsigned id_begin = nPerProc * p_rank ;
  const unsigned id_end   = nPerProc * ( p_rank + 1 );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  const unsigned p_send  = ( p_rank + 1 ) % p_size ;
  const unsigned id_send = id_end - 2 ;
  const unsigned id_recv = ( id_begin + id_total - 2 ) % id_total ;

  stk::mesh::PartVector part_intersection;
  part_intersection.push_back(&( mesh.mesh_meta_data().locally_used_part() ));
  Selector select_used( part_intersection );

  std::vector<unsigned> local_count ;
  std::vector<EntityProc> change ;

  Entity * send_edge_1 = mesh.get_entity( 1 , edge_ids[ id_send ] );
  Entity * send_edge_2 = mesh.get_entity( 1 , edge_ids[ id_send + 1 ] );
  Entity * send_node_1 = send_edge_1->relations()[1].entity();
  Entity * send_node_2 = send_edge_2->relations()[1].entity();
  Entity * recv_edge_1 = mesh.get_entity( 1 , edge_ids[ id_recv ] );
  Entity * recv_edge_2 = mesh.get_entity( 1 , edge_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( NULL != send_edge_1 && p_rank == send_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL != send_edge_2 && p_rank == send_edge_2->owner_rank() );
  STKUNIT_ASSERT( NULL == recv_edge_1 || p_rank != recv_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL == recv_edge_2 || p_rank != recv_edge_2->owner_rank() );

  if ( p_rank == send_node_1->owner_rank() ) {
    EntityProc entry( send_node_1 , p_send );
    change.push_back( entry );
  }
  if ( p_rank == send_node_2->owner_rank() ) {
    EntityProc entry( send_node_2 , p_send );
    change.push_back( entry );
  }
  {
    EntityProc entry( send_edge_1 , p_send );
    change.push_back( entry );
  }
  {
    EntityProc entry( send_edge_2 , p_send );
    change.push_back( entry );
  }

  send_edge_1 = NULL ;
  send_edge_2 = NULL ;
  send_node_1 = NULL ;
  send_node_2 = NULL ;
  recv_edge_1 = NULL ;
  recv_edge_2 = NULL ;

  STKUNIT_ASSERT( mesh.modification_begin() );
  mesh.change_entity_owner( change );
  UnitTestBulkData::modification_end( mesh , generate_aura );

  send_edge_1 = mesh.get_entity( 1 , edge_ids[ id_send ] );
  send_edge_2 = mesh.get_entity( 1 , edge_ids[ id_send + 1 ] );
  recv_edge_1 = mesh.get_entity( 1 , edge_ids[ id_recv ] );
  recv_edge_2 = mesh.get_entity( 1 , edge_ids[ id_recv + 1 ] );

  STKUNIT_ASSERT( NULL == send_edge_1 || p_rank != send_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL == send_edge_2 || p_rank != send_edge_2->owner_rank() );
  STKUNIT_ASSERT( NULL != recv_edge_1 && p_rank == recv_edge_1->owner_rank() );
  STKUNIT_ASSERT( NULL != recv_edge_2 && p_rank == recv_edge_2->owner_rank() );

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT( local_count[0] == nLocalNode );
  STKUNIT_ASSERT( local_count[1] == nLocalEdge );

  const std::vector<EntityProc> & shared = mesh.shared_entities();

  STKUNIT_ASSERT( shared.size() == 2u );

  const unsigned n0 = id_recv ;
  const unsigned n1 = id_send ;

  if ( n0 < n1 ) {
    STKUNIT_ASSERT( shared[0].first->identifier() == node_ids[n0] );
    STKUNIT_ASSERT( shared[1].first->identifier() == node_ids[n1] );
  }
  else {
    STKUNIT_ASSERT( shared[1].first->identifier() == node_ids[n0] );
    STKUNIT_ASSERT( shared[0].first->identifier() == node_ids[n1] );
  }
}

}

void UnitTestBulkData::testChangeOwner_loop( ParallelMachine pm )
{
  enum { nPerProc = 10 };
  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalEdge = nPerProc ;

  std::vector<EntityId> node_ids , edge_ids ;

  MetaData meta( fem_entity_type_names() );

  meta.commit();

  Selector select_owned( meta.locally_owned_part() );

  //Create a select_used selector using a required-part and a part-union
  //just to get extra test coverage. The part-union will simply be
  //locally-owned and locally-used. This is redundant, since locally-owned
  //is a subset of locally-used.
  stk::mesh::PartVector part_union;
  part_union.push_back(&(meta.locally_owned_part()));
  part_union.push_back(&(meta.locally_used_part()));

  Selector select_used( meta.locally_used_part(), part_union );

  Selector select_all(  meta.universal_part() );

  PartVector no_parts ;

  std::vector<unsigned> local_count ;

  //------------------------------
  {
    BulkData bulk( meta , pm , 100 );

    generate_loop( bulk, no_parts , false /* no aura */, nPerProc, node_ids, edge_ids );

    count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    count_entities( select_all , bulk , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    // Shift loop by two nodes and edges.

    if ( 1 < p_size ) {
      test_shift_loop(bulk,false/* no aura */,nPerProc,node_ids,edge_ids);

      count_entities( select_used , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nLocalNode );
      STKUNIT_ASSERT( local_count[1] == nLocalEdge );

      count_entities( select_all , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nLocalNode );
      STKUNIT_ASSERT( local_count[1] == nLocalEdge );
    }
  }

  //------------------------------
  // Test shift starting with ghosting but not regenerated ghosting.
  {
    BulkData bulk( meta , pm , 100 );

    generate_loop( bulk, no_parts , true /* aura */, nPerProc, node_ids, edge_ids );

    count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    count_entities( select_all , bulk , local_count );
    const unsigned n_extra = 1 < p_size ? 2 : 0 ;
    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

    if ( 1 < p_size ) {
      test_shift_loop(bulk,false/* aura */,nPerProc,node_ids,edge_ids);

      count_entities( select_owned , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nPerProc );
      STKUNIT_ASSERT( local_count[1] == nPerProc );

      count_entities( select_used , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nLocalNode );
      STKUNIT_ASSERT( local_count[1] == nLocalEdge );

      // All of my ghosts were disrupted and therefore deleted:
      count_entities( select_all , bulk , local_count );
      STKUNIT_ASSERT_EQUAL( nLocalEdge , local_count[1] );
      STKUNIT_ASSERT_EQUAL( nLocalNode , local_count[0] );
    }
  }
  //------------------------------
  // Test shift starting with ghosting and regenerating ghosting.
  {
    BulkData bulk( meta , pm , 100 );

    generate_loop( bulk, no_parts , true /* aura */, nPerProc, node_ids, edge_ids );

    count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    count_entities( select_all , bulk , local_count );
    const unsigned n_extra = 1 < p_size ? 2 : 0 ;
    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

    if ( 1 < p_size ) {
      test_shift_loop(bulk,true/* aura */,nPerProc,node_ids,edge_ids);

      count_entities( select_owned , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nPerProc );
      STKUNIT_ASSERT( local_count[1] == nPerProc );

      count_entities( select_used , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nLocalNode );
      STKUNIT_ASSERT( local_count[1] == nLocalEdge );

      // All of my ghosts were regenerated:
      count_entities( select_all , bulk , local_count );
      STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
      STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );
    }
  }
  //------------------------------
  // Test bad owner change catching:
  if ( 1 < p_size ) {
    BulkData bulk( meta , pm , 100 );

    generate_loop( bulk, no_parts , true /* aura */, nPerProc, node_ids, edge_ids );

    count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT( local_count[0] == nLocalNode );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge );

    count_entities( select_all , bulk , local_count );
    const unsigned n_extra = 1 < p_size ? 2 : 0 ;
    STKUNIT_ASSERT( local_count[0] == nLocalNode + n_extra );
    STKUNIT_ASSERT( local_count[1] == nLocalEdge + n_extra );

    std::vector<EntityProc> change ;

    if ( 0 == p_rank ) {
      change.resize(4);
      // Error to change to bad owner:
      change[0].first = bulk.get_entity( 0 , node_ids[1] );
      change[0].second = p_size ;
      // Error to change a ghost:
      change[1].first = bulk.shared_aura().receive()[0] ;
      change[1].second = p_rank ;
      // Error to change to multiple owners:
      change[2].first = bulk.get_entity( 0 , node_ids[1] );
      change[2].second = ( p_rank + 1 ) % p_size ;
      change[3].first = change[2].first ;
      change[3].second = ( p_rank + 2 ) % p_size ;
    }

    STKUNIT_ASSERT( bulk.modification_begin() );

    std::string error_msg ;
    bool exception_thrown = false ;
    try {
      bulk.change_entity_owner( change );
    }
    catch( const std::exception & x ) {
      exception_thrown = true ;
      error_msg.assign( x.what() );
    }
    STKUNIT_ASSERT( exception_thrown );
    std::cout
      << std::endl
      << "  UnitTestBulkData::testChangeOwner_loop SUCCESSFULLY CAUGHT: "
      << error_msg << std::endl ;
    std::cout.flush();
  }
  //------------------------------
  // Test move one element with initial ghosting but not regenerated ghosting:
  // last processor give its shared node to P0
  if ( 1 < p_size ) {
    BulkData bulk( meta , pm , 100 );

    generate_loop( bulk , no_parts , true , nPerProc , node_ids , edge_ids );

    std::vector<EntityProc> change ;

    if ( p_rank + 1 == p_size ) {
      EntityProc entry ;
      entry.first = bulk.get_entity( 0 , node_ids[0] );
      entry.second = 0 ;
      STKUNIT_ASSERT_EQUAL( p_rank , entry.first->owner_rank() );
      change.push_back( entry );
    }

    STKUNIT_ASSERT( bulk.modification_begin() );
    bulk.change_entity_owner( change );
    UnitTestBulkData::modification_end( bulk , false );

    count_entities( select_owned , bulk , local_count );
    const unsigned n_node = p_rank == 0          ? nPerProc + 1 : (
                            p_rank + 1 == p_size ? nPerProc - 1 :
                                                   nPerProc );
              
    STKUNIT_ASSERT_EQUAL( n_node , local_count[0] );
    STKUNIT_ASSERT_EQUAL( (unsigned) nPerProc , local_count[1] );

    count_entities( select_used , bulk , local_count );
    STKUNIT_ASSERT_EQUAL( nLocalNode , local_count[0] );
    STKUNIT_ASSERT_EQUAL( nLocalEdge , local_count[1] );

    // Moving the node disrupted ghosting on first and last process
    count_entities( select_all , bulk , local_count );
    const unsigned n_extra = p_rank + 1 == p_size || p_rank == 0 ? 1 : 2 ;

    STKUNIT_ASSERT_EQUAL( nLocalNode + n_extra , local_count[0] );
    STKUNIT_ASSERT_EQUAL( nLocalEdge + n_extra , local_count[1] );
  }
  //------------------------------

  std::cout << std::endl
            << "P" << p_rank
            << ": UnitTestBulkData::testChangeOwner_loop( NP = "
            << p_size << " ) SUCCESSFULL" << std::endl ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Testing for collection of boxes

namespace {

void donate_one_element( BulkData & mesh , bool aura )
{
  const unsigned p_rank = mesh.parallel_rank();

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );

  std::vector<unsigned> before_count ;
  std::vector<unsigned> after_count ;

  count_entities( select_owned , mesh , before_count );

  // Change owner of an element on a process boundary
  // from P0 to P1, and then recount to confirm ownership change

  std::vector<EntityProc> change ;

  const std::vector<EntityProc> & shared = mesh.shared_entities(); 
  STKUNIT_ASSERT( ! shared.empty() );

  // A shared node:
  Entity * const node = shared[0].first ;
  Entity * elem = NULL ;

  for ( PairIterRelation rel = node->relations( 3 );
        ! rel.empty() && elem == NULL ; ++rel ) {
    elem = rel->entity();
    if ( elem->owner_rank() != p_rank ) { elem = NULL ; }
  }

  STKUNIT_ASSERT( elem );

  unsigned donated_nodes = 0 ;

  // Only process #0 donates an element and its owned nodes:
  if ( 0 == p_rank ) {
    EntityProc entry ;
    entry.first = elem ;
    entry.second = shared[0].second ;
    change.push_back( entry );
    for ( PairIterRelation
          rel = elem->relations(0) ; ! rel.empty() ; ++rel ) {
      if ( rel->entity()->owner_rank() == p_rank ) {
        entry.first = rel->entity();
        change.push_back( entry );
        ++donated_nodes ;
      }
    }
  }

  STKUNIT_ASSERT( mesh.modification_begin() );
  mesh.change_entity_owner( change );
  UnitTestBulkData::modification_end( mesh , aura );

  count_entities( select_owned , mesh , after_count );

  if ( 0 == p_rank ) {
    STKUNIT_ASSERT_EQUAL( before_count[3] - 1 , after_count[3] );
    STKUNIT_ASSERT_EQUAL( before_count[0] - donated_nodes, after_count[0] );
  }
}

void donate_all_shared_nodes( BulkData & mesh , bool aura )
{
  const unsigned p_rank = mesh.parallel_rank();

  Selector select_used( mesh.mesh_meta_data().locally_used_part() );

  std::vector<unsigned> before_count ;
  std::vector<unsigned> after_count ;

  count_entities( select_used , mesh , before_count );

  // Donate owned shared nodes to first sharing process.

  const std::vector<EntityProc> & shared = mesh.shared_entities(); 

  STKUNIT_ASSERT( ! shared.empty() );

  std::vector<EntityProc> change ;

  {
    Entity * node = NULL ;
    for ( std::vector<EntityProc>::const_iterator
          i =  shared.begin() ;
          i != shared.end() && i->first->entity_type() == 0 ; ++i ) {
      if ( i->first->owner_rank() == p_rank && node != i->first ) {
        change.push_back( *i );
        node = i->first ;
      }
    }
  }

  STKUNIT_ASSERT( mesh.modification_begin() );
  mesh.change_entity_owner( change );
  UnitTestBulkData::modification_end( mesh , aura );

  count_entities( select_used , mesh , after_count );

  STKUNIT_ASSERT( 3 <= after_count.size() );
  STKUNIT_ASSERT_EQUAL( before_count[0] , after_count[0] );
  STKUNIT_ASSERT_EQUAL( before_count[1] , after_count[1] );
  STKUNIT_ASSERT_EQUAL( before_count[2] , after_count[2] );
  STKUNIT_ASSERT_EQUAL( before_count[3] , after_count[3] );
}

}

void UnitTestBulkData::testChangeOwner_box( ParallelMachine pm )
{
  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );

  MetaData meta( fem_entity_type_names() );

  meta.commit();

  Selector select_owned( meta.locally_owned_part() );
  Selector select_used(  meta.locally_used_part() );
  Selector select_all(   meta.universal_part() );

  //------------------------------
  {
    BulkData bulk( meta , pm , 100 );

    generate_boxes( bulk , false /* no aura */ , root_box , local_box );

    if ( 1 < p_size ) {
      donate_one_element( bulk , false /* no aura */ );
    }
  }

  if ( 1 < p_size ) {
    BulkData bulk( meta , pm , 100 );

    generate_boxes( bulk , false /* no aura */ , root_box , local_box );

    donate_all_shared_nodes( bulk , false /* no aura */ );
  }
  //------------------------------
  if ( 1 < p_size ) {
    BulkData bulk( meta , pm , 100 );

    generate_boxes( bulk , false /* no aura */ , root_box , local_box );

    donate_one_element( bulk , false /* no aura */ );
  }
  //------------------------------
  // Introduce ghosts:
  if ( 1 < p_size ) {
    BulkData bulk( meta , pm , 100 );

    generate_boxes( bulk , true /* aura */ , root_box , local_box );

    std::vector<unsigned> used_count ;
    std::vector<unsigned> all_count ;

    count_entities( select_all , bulk , all_count );
    count_entities( select_used , bulk , used_count );

    STKUNIT_ASSERT( used_count[0] < all_count[0] );
    STKUNIT_ASSERT( used_count[3] < all_count[3] );

    donate_all_shared_nodes( bulk , false /* don't regenerate aura */ );

    count_entities( select_all , bulk , all_count );
    count_entities( select_used , bulk , used_count );

    STKUNIT_ASSERT_EQUAL( used_count[0] , all_count[0] );
    STKUNIT_ASSERT_EQUAL( used_count[3] , all_count[3] );
  }
  //------------------------------

  std::cout << std::endl
            << "P" << p_rank
            << ": UnitTestBulkData::testChangeOwner_box( NP = "
            << p_size << " ) SUCCESSFULL" << std::endl ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

