/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <unit_tests/stk_utest_macros.hpp>

//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

//----------------------------------------------------------------------

namespace stk_unit_tests {

// Assume the following mesh of 4 hex8 elements.
//
// Global node and element numbering
//     3       7      11      15      19         
//     +-------+-------+-------+-------+        
//    /       /       /       /       /|       
//  4/      8/     12/     16/     20/ |       
//  +-------+-------+-------+-------+  |       
//  |       |       |       |       |  +18     
//  |  e1   |  e2   |  e3   |  e4   | /        
//  |       |       |       |       |/         
//  +-------+-------+-------+-------+          
//  1       5      9       13      17          
//
// Local node numbering
//     8       7
//     +-------+
//    /       /| 
//  5/      6/ | 
//  +-------+  |
//  |       |  +3 
//  |  e1   | /
//  |       |/ 
//  +-------+
//  1       2  
//----------------------------------------------------------------------

namespace {

void elem_node_ids( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] )
{
  if ( elem_id == 0 ) {
    std::cout << "use_case_1, elem_node_ids: ERROR, elem_id ("
        << elem_id << ") must be greater than 0." << std::endl;
    return;
  }

  const unsigned base = ( elem_id - 1 ) * 4 ;
  node_ids[0] = base + 1 ;
  node_ids[1] = base + 5 ;
  node_ids[2] = base + 6 ;
  node_ids[3] = base + 2 ;
  node_ids[4] = base + 4 ;
  node_ids[5] = base + 8 ;
  node_ids[6] = base + 7 ;
  node_ids[7] = base + 3 ;
}

}

#include <unit_tests/UnitTestUseCase_1.hpp>

//----------------------------------------------------------------------

void use_case_1_driver( stk::ParallelMachine comm )
{
  const unsigned field_data_chunk_size = 1000 ;
  const unsigned number_elements_left  = 1 ;
  const unsigned number_elements_right = 3 ;

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stkmesh use case 1, trivial two-block mesh" << std::endl
              << "  sizeof(Entity)   = " << sizeof(stk::mesh::Entity) << std::endl
              << "  sizeof(Relation) = " << sizeof(stk::mesh::Relation) << std::endl ;
    std::cout.flush();
  }

  use_case_1_verify_attributes();

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );
  stk::mesh::Part * left_block = NULL ;
  stk::mesh::Part * right_block = NULL ;

  use_case_1_declare_blocks( meta_data , & left_block , & right_block );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data , comm , field_data_chunk_size );

  use_case_1_generate_mesh( bulk_data ,
                                    *left_block , number_elements_left ,
                                    *right_block , number_elements_right );

  use_case_1_verify_mesh( bulk_data ,
                                  *left_block , number_elements_left ,
                                  *right_block , number_elements_right );
}

//----------------------------------------------------------------------

void use_case_1_declare_blocks(
  stk::mesh::MetaData & meta_data ,
  stk::mesh::Part ** left_block ,
  stk::mesh::Part ** right_block )
{
  stk::mesh::Part & left  = meta_data.declare_part( std::string("block_1") , stk::mesh::Element );
  stk::mesh::Part & right = meta_data.declare_part( std::string("block_2") , stk::mesh::Element );

  // Declare the intention of the parts to be for hex 8 elements.

  stk::mesh::set_cell_topology< shards::Hexahedron<8>  >( left );
  stk::mesh::set_cell_topology< shards::Hexahedron<8>  >( right );

  *left_block  = & left ;
  *right_block = & right ;
}

//----------------------------------------------------------------------

void use_case_1_generate_mesh(
  stk::mesh::BulkData & mesh ,
  stk::mesh::Part & left_block , unsigned number_left ,
  stk::mesh::Part & right_block , unsigned number_right )
{
  stk::mesh::EntityId elem_id = 1 ;
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];

  for ( unsigned j = 0 ; j < number_left ; ++j , ++elem_id ) {
    elem_node_ids( elem_id , node_ids );
    stk::mesh::declare_element( mesh , left_block , elem_id , node_ids );
  }

  for ( unsigned j = 0 ; j < number_right ; ++j , ++elem_id ) {
    elem_node_ids( elem_id , node_ids );
    stk::mesh::declare_element( mesh , right_block , elem_id , node_ids );
  }

  // No parallel stuff for now
}

//----------------------------------------------------------------------

void use_case_1_verify_attributes()
{
  for ( unsigned et = 0 ; et < 8 ; ++et ) {
    for ( unsigned k = 0 ; k < 3 ; ++k ) {
      for ( unsigned j = 0 ; j < 32 ; ++j ) {
        stk::mesh::relation_attr_type attr =
          stk::mesh::relation_attr( et , j , k );

        STKUNIT_ASSERT_EQUAL( et, stk::mesh::relation_entity_type( attr ) );
        STKUNIT_ASSERT_EQUAL( j , stk::mesh::relation_identifier( attr ) );
        STKUNIT_ASSERT_EQUAL( k , stk::mesh::relation_kind( attr ) );
      }
    }
  }
}

//----------------------------------------------------------------------

void use_case_1_verify_mesh(
  stk::mesh::BulkData & mesh ,
  stk::mesh::Part & left_block , unsigned number_left ,
  stk::mesh::Part & right_block , unsigned number_right )
{
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8> ::node_count ];

  std::vector<unsigned> entity_counts ;

  STKUNIT_ASSERT_EQUAL( stk::mesh::get_cell_topology( left_block ) ,
                        shards::getCellTopologyData< shards::Hexahedron<8>  >() );

  STKUNIT_ASSERT_EQUAL( stk::mesh::get_cell_topology( right_block ) ,
                        shards::getCellTopologyData< shards::Hexahedron<8> >() );

  stk::mesh::Selector selector_left(left_block);
  stk::mesh::count_entities( selector_left, mesh , entity_counts );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Node] ,
                        (unsigned)( number_left + 1 ) * 4 );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Edge] ,    0u );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Face] ,    0u );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Element] , (unsigned) number_left );

  stk::mesh::Selector selector_right(right_block);
  stk::mesh::count_entities( selector_right, mesh , entity_counts );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Node] ,
                        (unsigned)( number_right + 1 ) * 4 );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Edge] ,    0u );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Face] ,    0u );
  STKUNIT_ASSERT_EQUAL( entity_counts[stk::mesh::Element] , (unsigned) number_right );

  stk::mesh::EntityId elem_id = 1 ;

  for ( unsigned j = 0 ; j < number_left ; ++j , ++elem_id ) {
    elem_node_ids( elem_id , node_ids );

    stk::mesh::Entity * const elem = mesh.get_entity( stk::mesh::Element , elem_id );

    STKUNIT_ASSERT( elem != NULL );
    STKUNIT_ASSERT( has_superset( elem->bucket(), left_block ) );

    stk::mesh::PairIterRelation rel = elem->relations();

    STKUNIT_ASSERT( shards::Hexahedron<8> ::node_count <= rel.size() );

    for ( unsigned i = 0 ; i < shards::Hexahedron<8> ::node_count ; ++i ) {
      stk::mesh::Entity * const rel_node = rel[i].entity();
      STKUNIT_ASSERT_EQUAL( node_ids[i] , rel_node->identifier() );
      STKUNIT_ASSERT( has_superset( rel_node->bucket(), left_block ) );
    }
  }

  for ( unsigned j = 0 ; j < number_right ; ++j , ++elem_id ) {

    elem_node_ids( elem_id , node_ids );

    stk::mesh::Entity * const elem = mesh.get_entity( stk::mesh::Element , elem_id );

    STKUNIT_ASSERT( elem != NULL );
    STKUNIT_ASSERT( has_superset( elem->bucket(),right_block ) );

    stk::mesh::PairIterRelation rel = elem->relations();

    STKUNIT_ASSERT( shards::Hexahedron<8> ::node_count <= rel.size() );

    for ( unsigned i = 0 ; i < shards::Hexahedron<8> ::node_count ; ++i ) {
      stk::mesh::Entity * const rel_node = rel[i].entity();
      STKUNIT_ASSERT_EQUAL( node_ids[i] , rel_node->identifier() );
      STKUNIT_ASSERT( has_superset( rel_node->bucket(), right_block ) );
    }
  }
}

} // namespace stk_unit_tests


STKUNIT_UNIT_TEST(UnitTestUseCase_1, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  stk_unit_tests::use_case_1_driver(MPI_COMM_WORLD);
}

