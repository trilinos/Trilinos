/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>

#include <unit_tests/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>

#include <unit_tests/UnitTestBulkData.hpp>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace {

/* Recursively split a box into ( up - ip ) sub-boxes */

typedef int BOX[3][2] ;

void box_partition( int ip , int up , int axis , 
                    const BOX box ,
                    BOX p_box[] )
{ 
  const int np = up - ip ;
  if ( 1 == np ) {
    p_box[ip][0][0] = box[0][0] ; p_box[ip][0][1] = box[0][1] ;
    p_box[ip][1][0] = box[1][0] ; p_box[ip][1][1] = box[1][1] ;
    p_box[ip][2][0] = box[2][0] ; p_box[ip][2][1] = box[2][1] ;
  } 
  else {  
    const int n = box[ axis ][1] - box[ axis ][0] ;
    const int np_low = np / 2 ;  /* Rounded down */
    const int np_upp = np - np_low ;
 
    const int n_upp = (int) (((double) n) * ( ((double)np_upp) / ((double)np)));
    const int n_low = n - n_upp ;
    const int next_axis = ( axis + 2 ) % 3 ;
 
    if ( np_low ) { /* P = [ip,ip+np_low) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;
 
      dbox[ axis ][1] = dbox[ axis ][0] + n_low ;
 
      box_partition( ip, ip + np_low, next_axis,
                     (const int (*)[2]) dbox, p_box );
    }
 
    if ( np_upp ) { /* P = [ip+np_low,ip+np_low+np_upp) */
      BOX dbox ;
      dbox[0][0] = box[0][0] ; dbox[0][1] = box[0][1] ;
      dbox[1][0] = box[1][0] ; dbox[1][1] = box[1][1] ;
      dbox[2][0] = box[2][0] ; dbox[2][1] = box[2][1] ;
 
      ip += np_low ;
      dbox[ axis ][0] += n_low ;
      dbox[ axis ][1]  = dbox[ axis ][0] + n_upp ;
 
      box_partition( ip, ip + np_upp, next_axis,
                     (const int (*)[2]) dbox, p_box );
    }
  }
}

}

void UnitTestBulkData::generate_boxes(
  BulkData  & mesh ,
  const bool  generate_aura ,
  const int   root_box[][2] ,
        int   local_box[][2] )
{
  const unsigned p_rank = mesh.parallel_rank();
  const unsigned p_size = mesh.parallel_size();
  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;
  const unsigned ngz = root_box[2][1] - root_box[2][0] ;
/*
  const unsigned e_global = ngx * ngy * ngz ;
  const unsigned n_global = ( ngx + 1 ) * ( ngy + 1 ) * ( ngz + 1 );
*/

  if ( 0 == p_rank ) {
    std::cout << "Global box = " << ngx << " x " << ngy << " x " << ngz 
              << std::endl ;
  }

  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  local_box[0][0] = p_box[ p_rank ][0][0] ;
  local_box[0][1] = p_box[ p_rank ][0][1] ;
  local_box[1][0] = p_box[ p_rank ][1][0] ;
  local_box[1][1] = p_box[ p_rank ][1][1] ;
  local_box[2][0] = p_box[ p_rank ][2][0] ;
  local_box[2][1] = p_box[ p_rank ][2][1] ;

  const unsigned nx = local_box[0][1] - local_box[0][0] ;
  const unsigned ny = local_box[1][1] - local_box[1][0] ;
  const unsigned nz = local_box[2][1] - local_box[2][0] ;

  const unsigned e_local = nx * ny * nz ;
  const unsigned n_local = ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 );

  // Create elements:

  std::vector<unsigned> local_count ;

  const PartVector no_parts ;

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    const EntityId n0 = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n1 = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n2 = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n3 = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const EntityId n4 = 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const EntityId n5 = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const EntityId n6 = 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const EntityId n7 = 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    const EntityId elem_id =  1 + i + j * ngx + k * ngx * ngy;

    Entity & node0 = mesh.declare_entity( 0 , n0 , no_parts );
    Entity & node1 = mesh.declare_entity( 0 , n1 , no_parts );
    Entity & node2 = mesh.declare_entity( 0 , n2 , no_parts );
    Entity & node3 = mesh.declare_entity( 0 , n3 , no_parts );
    Entity & node4 = mesh.declare_entity( 0 , n4 , no_parts );
    Entity & node5 = mesh.declare_entity( 0 , n5 , no_parts );
    Entity & node6 = mesh.declare_entity( 0 , n6 , no_parts );
    Entity & node7 = mesh.declare_entity( 0 , n7 , no_parts );
    Entity & elem  = mesh.declare_entity( 3 , elem_id , no_parts );

    mesh.declare_relation( elem , node0 , 0 );
    mesh.declare_relation( elem , node1 , 1 );
    mesh.declare_relation( elem , node2 , 2 );
    mesh.declare_relation( elem , node3 , 3 );
    mesh.declare_relation( elem , node4 , 4 );
    mesh.declare_relation( elem , node5 , 5 );
    mesh.declare_relation( elem , node6 , 6 );
    mesh.declare_relation( elem , node7 , 7 );
  }
  }
  }

  Selector select_owned( mesh.mesh_meta_data().locally_owned_part() );
  Selector select_used( mesh.mesh_meta_data().locally_used_part() );
  Selector select_all(  mesh.mesh_meta_data().universal_part() );

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT_EQUAL( e_local , local_count[3] );
  STKUNIT_ASSERT_EQUAL( 0u , local_count[2] );
  STKUNIT_ASSERT_EQUAL( 0u , local_count[1] );
  STKUNIT_ASSERT_EQUAL( n_local , local_count[0] );

  // Set up sharing:
  STKUNIT_ASSERT( mesh.internal_modification_end( generate_aura ) );

  // Verify declarations and sharing

  count_entities( select_used , mesh , local_count );
  STKUNIT_ASSERT( local_count[3] == e_local );
  STKUNIT_ASSERT( local_count[2] == 0 );
  STKUNIT_ASSERT( local_count[1] == 0 );
  STKUNIT_ASSERT( local_count[0] == n_local );

  for ( int k = local_box[2][0] ; k <= local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j <= local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i <= local_box[0][1] ; ++i ) {
    EntityType node_type = 0;
    EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
    Entity * const node = mesh.get_entity( node_type , node_id );
    STKUNIT_ASSERT( node != NULL );
    // Shared if on a processor boundary.
    const bool shared =
      ( k == local_box[2][0] && k != root_box[2][0] ) ||
      ( k == local_box[2][1] && k != root_box[2][1] ) ||
      ( j == local_box[1][0] && j != root_box[1][0] ) ||
      ( j == local_box[1][1] && j != root_box[1][1] ) ||
      ( i == local_box[0][0] && i != root_box[0][0] ) ||
      ( i == local_box[0][1] && i != root_box[0][1] );
    if (mesh.parallel_size() > 1) {
      STKUNIT_ASSERT_EQUAL( shared , ! node->sharing().empty() );
    }
  }
  }
  }

  size_t count_shared_node_pairs = 0 ;
  for ( unsigned p = 0 ; p < p_size ; ++p ) if ( p != p_rank ) {
    for ( int k = p_box[p][2][0] ; k <= p_box[p][2][1] ; ++k )
    if ( local_box[2][0] <= k && k <= local_box[2][1] ) {

      for ( int j = p_box[p][1][0] ; j <= p_box[p][1][1] ; ++j )
      if ( local_box[1][0] <= j && j <= local_box[1][1] ) {

        for ( int i = p_box[p][0][0] ; i <= p_box[p][0][1] ; ++i )
        if ( local_box[0][0] <= i && i <= local_box[0][1] ) {

          EntityType node_type = 0;
          EntityId node_id = 1 + i + j * (ngx+1) + k * (ngx+1) * (ngy+1);
          Entity * const node = mesh.get_entity( node_type , node_id );
          STKUNIT_ASSERT( node != NULL );
          // Must be shared with 'p'
          PairIterEntityComm iter = node->sharing();
          for ( ; ! iter.empty() && iter->proc != p ; ++iter );
          STKUNIT_ASSERT( ! iter.empty() );

          ++count_shared_node_pairs ;
        }
      }
    }
  }

  size_t count_shared_entities = 0 ;
  for ( std::vector<Entity*>::const_iterator
        i = mesh.entity_comm().begin() ; i != mesh.entity_comm().end() ; ++i ) {
    const PairIterEntityComm ec = (**i).sharing();
    count_shared_entities += ec.size();
  }
  STKUNIT_ASSERT_EQUAL( count_shared_entities , count_shared_node_pairs );

  delete[] p_box ;
}

} // namespace mesh
} // namespace stk

