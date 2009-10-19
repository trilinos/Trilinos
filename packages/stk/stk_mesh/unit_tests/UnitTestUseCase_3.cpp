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

void use_case_3_interpret_fields(
  const std::vector<const stk::mesh::FieldBase *> & output_fields );

void use_case_3_driver( MPI_Comm comm )
{
  const unsigned field_data_chunk_size = 1000 ;
  const unsigned number_elements_left  = 1 ;
  const unsigned number_elements_right = 3 ;

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stkmesh use case 3, anon field query" << std::endl;
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

  //--------------------------------

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

void use_case_3_interpret_fields(
  const std::vector<const stk::mesh::FieldBase *> & output_fields )
{
  std::cout << "stk_mesh use case 3, field interpretation {" << std::endl ;

  std::vector<const stk::mesh::FieldBase *>::const_iterator i ;

  for ( i = output_fields.begin() ; i != output_fields.end() ; ++i ) {
    const stk::mesh::FieldBase & f = **i ;

    std::cout << "  Field<" ;
    std::cout << f.data_traits().name ;
    for ( unsigned j = 0 ; j < f.rank() ; ++j ) {
      if ( j ) std::cout << " ," ;
      std::cout << " " << f.dimension_tags()[j]->name();
    }
    std::cout << "> " ;
    std::cout << f.name();
    std::cout << std::endl ;
  }
  std::cout <<  "}" << std::endl ;
  std::cout.flush();
}

} // namespace stk_unit_tests

STKUNIT_UNIT_TEST(UnitTestUseCase_3, testUnit)
{
  int mpi_rank = 0;
  int mpi_size = 0;
  
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
  
  stk_unit_tests::use_case_3_driver(MPI_COMM_WORLD);
}

