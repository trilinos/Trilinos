#ifndef unit_tests_UnitTestUseCase_1_hpp
#define unit_tests_UnitTestUseCase_1_hpp

#include <stk_util/parallel/Parallel.hpp>

void use_case_1_verify_attributes();

void use_case_1_declare_blocks(
  stk::mesh::MetaData & meta_data ,
  stk::mesh::Part ** left_block ,
  stk::mesh::Part ** right_block );

void use_case_1_generate_mesh(
  stk::mesh::BulkData & mesh ,
  stk::mesh::Part & left_block , unsigned number_left ,
  stk::mesh::Part & right_block , unsigned number_right );

void use_case_1_verify_mesh(
  stk::mesh::BulkData & mesh ,
  stk::mesh::Part & left_block , unsigned number_left ,
  stk::mesh::Part & right_block , unsigned number_right );

void use_case_1_driver( stk::ParallelMachine comm );

#endif

