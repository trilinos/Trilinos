/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <iostream>
#include <mpi.h>

#include <unit_tests/stk_utest_macros.hpp>

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

namespace {

void node_coordinates( unsigned node_id , double coord[] )
{
  const unsigned i_length = ( node_id - 1 ) / 4 ;
  const unsigned i_plane  = ( node_id - 1 ) % 4 ;

  coord[0] = i_length ;
  coord[1] = i_plane == 1 || i_plane == 2 ? 1.0 : 0.0 ;
  coord[2] = i_plane == 2 || i_plane == 3 ? 1.0 : 0.0 ;
}

}

typedef shards::Hexahedron<8>  ElementTraits ;

enum { SpaceDim = 3 };

typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;

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

void use_case_2_driver( MPI_Comm comm )
{
  const unsigned field_data_chunk_size = 1000 ;
  const unsigned number_elements_left  = 1 ;
  const unsigned number_elements_right = 3 ;

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stkmesh use case 2, trivial two-block mesh" << std::endl
              << "  sizeof(Entity)   = " << sizeof(stk::mesh::Entity) << std::endl
              << "  sizeof(Relation) = " << sizeof(stk::mesh::Relation) << std::endl ;
    std::cout.flush();
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
}

//----------------------------------------------------------------------

void use_case_2_declare_fields(
  stk::mesh::MetaData & meta_data ,
  VectorFieldType ** coordinates_field ,
  ScalarFieldType ** temperature_field ,
  ScalarFieldType ** volume_field )
{
  const std::string name_coordinates("coordinates");
  const std::string name_temperature("temperature");
  const std::string name_volume("volume");

  stk::mesh::Part & universal = meta_data.universal_part();

  // Field declarations:

  *coordinates_field =
    & meta_data.declare_field< VectorFieldType >( name_coordinates );

  *temperature_field =
    & meta_data.declare_field< ScalarFieldType >( name_temperature );

  *volume_field =
    & meta_data.declare_field< ScalarFieldType >( name_volume );

  // Field restrictions:

  stk::mesh::put_field(
    **coordinates_field , stk::mesh::Node , universal , SpaceDim );

  stk::mesh::put_field( **temperature_field, stk::mesh::Node, universal );
  stk::mesh::put_field( **volume_field , stk::mesh::Element , universal );
}

//----------------------------------------------------------------------

void use_case_2_assign_field_data(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & coordinates_field ,
  const ScalarFieldType & temperature_field ,
  const ScalarFieldType & volume_field )
{
  const std::vector<stk::mesh::Bucket*> & node_buckets = mesh.buckets( stk::mesh::Node );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator
        i = node_buckets.begin() ; i != node_buckets.end() ; ++i ) {

    const stk::mesh::Bucket & bucket = **i;
    size_t n = bucket.size();

    double * coord = stk::mesh::field_data( coordinates_field , bucket.begin() );
    double * temp  = stk::mesh::field_data( temperature_field , bucket.begin() );

    for ( size_t j = 0; j < n; ++j) {
      const unsigned node_id = bucket[j].identifier();
    
      node_coordinates( node_id , coord );
      *temp = 98.6 ;

      coord += 3;
      ++temp;
    }

  }

  const std::vector<stk::mesh::Bucket*> & elem_buckets = mesh.buckets( stk::mesh::Element );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator
        i = elem_buckets.begin() ; i != elem_buckets.end() ; ++i ) {
    const stk::mesh::Bucket & bucket = **i;

    double * volume = stk::mesh::field_data( volume_field , bucket.begin() );

    for( size_t j = 0; j < bucket.size(); ++j) {
      volume[j] = 0.0 ;
    }
  }
}

//----------------------------------------------------------------------

} // namespace stk_unit_tests


STKUNIT_UNIT_TEST(UnitTestUseCase_2, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  stk_unit_tests::use_case_2_driver(MPI_COMM_WORLD);
}

