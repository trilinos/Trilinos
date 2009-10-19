/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <iostream>

#include <unit_tests/stk_utest_macros.hpp>

//----------------------------------------------------------------------

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
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
//  |       |       |       |       |  +18    Z  Y
//  |  e1   |  e2   |  e3   |  e4   | /       | /
//  |       |       |       |       |/        |/
//  +-------+-------+-------+-------+         *--X
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

typedef shards::Hexahedron<8>  ElementTraits ;

enum { SpaceDim = 3 };

typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;

//----------------------------------------------------------------------

template< unsigned NType , stk::mesh::EntityType EType ,
          unsigned NRel , class field_type >
bool gather_field_data( const field_type & field ,
                        const stk::mesh::Entity     & entity ,
                        typename stk::mesh::FieldTraits< field_type >::data_type * dst )
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type T ;

  stk::mesh::PairIterRelation rel = entity.relations( EType );

  bool result = NRel == (unsigned) rel.size();

  if ( result ) {
    T * const dst_end = dst + NType * NRel ;
    for ( const T * src ;
          ( dst < dst_end ) &&
          ( src = field_data( field , * rel->entity() ) ) ;
          ++rel , dst += NType ) {
      stk::Copy<NType>( dst , src );
    }
    result = dst == dst_end ;
  }
  return result ;
}

//----------------------------------------------------------------------

void use_case_1_declare_blocks(
  stk::mesh::MetaData & meta_data ,
  stk::mesh::Part ** left_block ,
  stk::mesh::Part ** right_block );

void use_case_2_declare_fields(
  stk::mesh::MetaData & meta_data ,
  VectorFieldType ** coordinates_field ,
  ScalarFieldType ** temperature_field ,
  ScalarFieldType ** volume_field );

void use_case_1_generate_mesh(
  stk::mesh::BulkData & mesh ,
  stk::mesh::Part & left_block , unsigned number_left ,
  stk::mesh::Part & right_block , unsigned number_right );

void use_case_2_assign_field_data(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & coordinates_field ,
  const ScalarFieldType & temperature_field ,
  const ScalarFieldType & volume_field );

void use_case_4_assign_element_volume(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & coordinates_field ,
  const ScalarFieldType & volume_field );

void use_case_4_driver( MPI_Comm comm )
{
  const unsigned field_data_chunk_size = 1000 ;
  const unsigned number_elements_left  = 1 ;
  const unsigned number_elements_right = 3 ;

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stkmesh use case 4, elem-node-gather" << std::endl;
  }

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );
  stk::mesh::Part * left_block = NULL ;
  stk::mesh::Part * right_block = NULL ;
  VectorFieldType * coordinates_field = NULL ;
  ScalarFieldType * temperature_field = NULL ;
  ScalarFieldType * volume_field = NULL ;

  use_case_1_declare_blocks( meta_data, & left_block, & right_block );

  use_case_2_declare_fields( meta_data ,
                                     & coordinates_field ,
                                     & temperature_field ,
                                     & volume_field );

  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data , comm , field_data_chunk_size );

  use_case_1_generate_mesh( bulk_data ,
                                    *left_block , number_elements_left ,
                                    *right_block , number_elements_right );

  use_case_2_assign_field_data( bulk_data ,
                                        * coordinates_field ,
                                        * temperature_field ,
                                        * volume_field );

  use_case_4_assign_element_volume( bulk_data ,
                                            * coordinates_field ,
                                            * volume_field );
}

//----------------------------------------------------------------------

void use_case_4_assign_element_volume(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & coordinates_field ,
  const ScalarFieldType & volume_field )
{
  const std::vector<stk::mesh::Bucket*> & elem_buckets = mesh.buckets( stk::mesh::Element );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator
        i = elem_buckets.begin() ; i != elem_buckets.end() ; ++i ) {

    const stk::mesh::Bucket& bucket = **i;
    size_t n = bucket.size();

    for( size_t k = 0; k < n; ++k) {
      const stk::mesh::Entity & elem = bucket[k] ;
      double elem_coord[ ElementTraits::node_count ][ SpaceDim ];

      const bool gather_result =
        gather_field_data< SpaceDim , stk::mesh::Node , ElementTraits::node_count >
                         ( coordinates_field , elem , & elem_coord[0][0] );

      STKUNIT_ASSERT( gather_result );

      double base[3] ; 
      base[0] = elem_coord[0][0] ;
      base[1] = elem_coord[0][1] ;
      base[2] = elem_coord[0][2] ;

      for ( unsigned j = 0 ; j < ElementTraits::node_count ; ++j ) {
        elem_coord[j][0] -= base[0] ;
        elem_coord[j][1] -= base[1] ;
        elem_coord[j][2] -= base[2] ;
      }

      STKUNIT_ASSERT_EQUAL( elem_coord[0][0] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[0][1] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[0][2] , 0.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[1][0] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[1][1] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[1][2] , 0.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[2][0] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[2][1] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[2][2] , 0.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[3][0] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[3][1] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[3][2] , 0.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[4][0] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[4][1] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[4][2] , 1.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[5][0] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[5][1] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[5][2] , 1.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[6][0] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[6][1] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[6][2] , 1.0 );

      STKUNIT_ASSERT_EQUAL( elem_coord[7][0] , 0.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[7][1] , 1.0 );
      STKUNIT_ASSERT_EQUAL( elem_coord[7][2] , 1.0 );

      double * volume = field_data( volume_field , bucket.begin() );

      *volume = 1.0 ;
    }
  }
}

//----------------------------------------------------------------------

} // namespace stk_unit_tests

STKUNIT_UNIT_TEST(UnitTestUseCase_4, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  stk_unit_tests::use_case_4_driver(MPI_COMM_WORLD);
}

