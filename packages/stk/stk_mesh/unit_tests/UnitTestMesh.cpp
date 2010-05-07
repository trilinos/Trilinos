/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>

#include <unit_tests/UnitTestMesh.hpp>


namespace stk {
namespace unit_test {

/****************************************************************/

std::vector<std::string>  get_entity_rank_names ( int rank )
{
  std::vector<std::string>  ret_val;
  ret_val.push_back ( "Node" );
  if ( rank == 0 ) return ret_val;
  ret_val.push_back ( "Edge" );
  if ( rank == 1 ) return ret_val;
  ret_val.push_back ( "Face" );
  if ( rank == 2 ) return ret_val;
  ret_val.push_back ( "Element" );
  if ( rank == 3 ) return ret_val;
  for ( int i = 3 ; i != rank ; i++ )
  {
    std::stringstream  name;
    name << "Entity rank " << i;
    ret_val.push_back ( name.str() );
  }
  return ret_val;
}


UnitTestMesh::UnitTestMesh ( stk::ParallelMachine comm , unsigned block_size )
  : m_meta_data ( get_entity_rank_names ( MAX_RANK ) )
  , m_bulk_data ( m_meta_data , comm , block_size )
  , m_test_part ( m_meta_data.declare_part ( "Test Part" ) )
  , m_cell_part ( m_meta_data.declare_part ( "Cell list" , MAX_RANK ) )
  , m_part_A_0 ( m_meta_data.declare_part ( "Part A 0", 0 ) )
  , m_part_A_1 ( m_meta_data.declare_part ( "Part A 1", 1 ) )
  , m_part_A_2 ( m_meta_data.declare_part ( "Part A 2", 2 ) )
  , m_part_A_3 ( m_meta_data.declare_part ( "Part A 3", 3 ) )
  , m_part_A_superset ( m_meta_data.declare_part ( "Part A superset" ) )
  , m_part_B_0 ( m_meta_data.declare_part ( "Part B 0", 0 ) )
  , m_part_B_1 ( m_meta_data.declare_part ( "Part B 1", 1 ) )
  , m_part_B_2 ( m_meta_data.declare_part ( "Part B 2", 2 ) )
  , m_part_B_3 ( m_meta_data.declare_part ( "Part B 3", 3 ) )
  , m_part_B_superset ( m_meta_data.declare_part ( "Part B superset" ) )
  , m_comm_rank( stk::parallel_machine_rank( comm ) )
  , m_comm_size( stk::parallel_machine_size( comm ) )
  , m_previous_state ( stk::mesh::BulkData::MODIFIABLE )
{
  m_meta_data.declare_part_subset ( m_part_A_superset , m_part_A_0 );
  m_meta_data.declare_part_subset ( m_part_A_superset , m_part_A_1 );
  m_meta_data.declare_part_subset ( m_part_A_superset , m_part_A_2 );
  m_meta_data.declare_part_subset ( m_part_A_superset , m_part_A_3 );

  m_meta_data.declare_part_subset ( m_part_B_superset , m_part_B_0 );
  m_meta_data.declare_part_subset ( m_part_B_superset , m_part_B_1 );
  m_meta_data.declare_part_subset ( m_part_B_superset , m_part_B_2 );
  m_meta_data.declare_part_subset ( m_part_B_superset , m_part_B_3 );

  m_meta_data.commit ();
}


void UnitTestMesh::generate_boxes ( bool aura )
{
  enter_modification();
  const int root_box[3][2] = { { 0,4 } , { 0,5 } , { 0,6 } };
  int local_box[3][2] = { { 0,0 } , { 0,0 } , { 0,0 } };
  priv_generate_boxes( m_bulk_data , aura , root_box , local_box );
  exit_modification();
}



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

void UnitTestMesh::priv_generate_boxes(
  stk::mesh::BulkData  & mesh ,
  const bool  generate_aura ,
  const int   root_box[][2] ,
        int   local_box[][2] )
{
  const unsigned p_rank = mesh.parallel_rank();
  const unsigned p_size = mesh.parallel_size();
  const unsigned ngx = root_box[0][1] - root_box[0][0] ;
  const unsigned ngy = root_box[1][1] - root_box[1][0] ;
//  const unsigned ngz = root_box[2][1] - root_box[2][0] ;
/*
  const unsigned e_global = ngx * ngy * ngz ;
  const unsigned n_global = ( ngx + 1 ) * ( ngy + 1 ) * ( ngz + 1 );
*/


  BOX * const p_box = new BOX[ p_size ];

  box_partition( 0 , p_size , 2 , root_box , & p_box[0] );

  local_box[0][0] = p_box[ p_rank ][0][0] ;
  local_box[0][1] = p_box[ p_rank ][0][1] ;
  local_box[1][0] = p_box[ p_rank ][1][0] ;
  local_box[1][1] = p_box[ p_rank ][1][1] ;
  local_box[2][0] = p_box[ p_rank ][2][0] ;
  local_box[2][1] = p_box[ p_rank ][2][1] ;

  //const unsigned nx = local_box[0][1] - local_box[0][0] ;
  //const unsigned ny = local_box[1][1] - local_box[1][0] ;
  //const unsigned nz = local_box[2][1] - local_box[2][0] ;

  //const unsigned e_local = nx * ny * nz ;
  //const unsigned n_local = ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 );

  // Create elements:

  std::vector<unsigned> local_count ;

  const stk::mesh::PartVector no_parts ;

  for ( int k = local_box[2][0] ; k < local_box[2][1] ; ++k ) {
  for ( int j = local_box[1][0] ; j < local_box[1][1] ; ++j ) {
  for ( int i = local_box[0][0] ; i < local_box[0][1] ; ++i ) {
    const stk::mesh::EntityId n0 = 1 + (i+0) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n1 = 1 + (i+1) + (j+0) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n2 = 1 + (i+1) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n3 = 1 + (i+0) + (j+1) * (ngx+1) + (k+0) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n4 = 1 + (i+0) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n5 = 1 + (i+1) + (j+0) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n6 = 1 + (i+1) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);
    const stk::mesh::EntityId n7 = 1 + (i+0) + (j+1) * (ngx+1) + (k+1) * (ngx+1) * (ngy+1);

    const stk::mesh::EntityId elem_id =  1 + i + j * ngx + k * ngx * ngy;

    stk::mesh::Entity & node0 = mesh.declare_entity( 0 , n0 , no_parts );
    stk::mesh::Entity & node1 = mesh.declare_entity( 0 , n1 , no_parts );
    stk::mesh::Entity & node2 = mesh.declare_entity( 0 , n2 , no_parts );
    stk::mesh::Entity & node3 = mesh.declare_entity( 0 , n3 , no_parts );
    stk::mesh::Entity & node4 = mesh.declare_entity( 0 , n4 , no_parts );
    stk::mesh::Entity & node5 = mesh.declare_entity( 0 , n5 , no_parts );
    stk::mesh::Entity & node6 = mesh.declare_entity( 0 , n6 , no_parts );
    stk::mesh::Entity & node7 = mesh.declare_entity( 0 , n7 , no_parts );
    stk::mesh::Entity & elem  = mesh.declare_entity( 3 , elem_id , no_parts );

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

  // Set up sharing:
  mesh.modification_end();

  delete[] p_box ;
}

//----------------------------------------------------------------------


} // namespace unit_test
} // namespace stk



