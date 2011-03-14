/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>

#include <use_cases/UseCase_1.hpp>
#include <use_cases/UseCase_2.hpp>
#include <use_cases/UseCase_3.hpp>
#include <use_cases/UseCase_4.hpp>
#include <use_cases/UseCase_ElementDeath.hpp>
#include <use_cases/UseCase_Skinning.hpp>
#include <use_cases/UseCase_ChangeOwner.hpp>

#include <stk_mesh/base/Types.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

void printStatus(bool status)
{
  if (status) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }
}

int main ( int argc, char * argv[] )
{
  stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);
  const bool single_process =
    stk::parallel_machine_size( parallel_machine ) <= 1 ;

  bool status = true;

   if ( single_process ) {
    std::cout << "Use Case 1 ... ";
    bool local_status = true ;
    try {
      stk::mesh::use_cases::UseCase_1_Mesh mesh(parallel_machine);
      printStatus(local_status);
    }
    catch ( const std::exception & x ) {
      local_status = false ;
      printStatus(local_status);
      std::cout << x.what();
    }
    status = status && local_status;
   }

  if ( single_process ) {
    std::cout << "Use Case 2 ... ";
    bool local_status = true ;
    try {
      stk::mesh::use_cases::UseCase_2_Mesh mesh(parallel_machine);
      mesh.populate(1,3);
      local_status = stk::mesh::use_cases::verifyMesh(mesh,1,3);
      printStatus(local_status);
    }
    catch ( const std::exception & x ) {
      local_status = false ;
      printStatus(local_status);
      std::cout << x.what();
    }
    status = status && local_status;
  }

  if ( single_process ) {
    std::cout << "Use Case 3 ... ";
    bool local_status = true ;
    try {
      stk::mesh::use_cases::UseCase_3_Mesh mesh(parallel_machine);
      mesh.populate();
      local_status = stk::mesh::use_cases::verifyMesh(mesh);
      printStatus(local_status);
    }
    catch ( const std::exception & x ) {
      local_status = false ;
      printStatus(local_status);
      std::cout << x.what();
    }
    status = status && local_status;
  }

  if ( single_process ) {
    std::cout << "Use Case 4 ... ";
    stk::mesh::use_cases::UseCase_4_Mesh mesh(parallel_machine);
    mesh.populate();
    stk::mesh::use_cases::runAlgorithms(mesh);
    const bool local_status = stk::mesh::use_cases::verifyMesh(mesh);
    printStatus(local_status);
    status = status && local_status;
  }

  {
    std::cout << "Use Case Change Owner ... ";
    Grid2D_Fixture test( parallel_machine );
    const bool local_status = test.test_change_owner();
    printStatus(local_status);
    status = status && local_status;
  }

  {
    std::cout << "Use Case Change Owner with constraint ... ";
    const bool local_status = test_change_owner_with_constraint( parallel_machine );
    printStatus(local_status);
    status = status && local_status;
  }

  {
    std::cout << "Use Case Change Owner #2 ... ";
    const bool local_status = test_change_owner_2( parallel_machine );
    printStatus(local_status);
    status = status && local_status;
  }

  {
    std::cout << "Use Case Change Owner #3 ... ";
    const bool result = test_change_owner_3( parallel_machine );
    printStatus(result);
  }

  {
    std::cout << "Use Case Element Death 1 ... ";
    bool local_status = element_death_use_case_1(parallel_machine);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case Skinning 1 ... ";
    bool local_status = skinning_use_case_1(parallel_machine);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case Skinning 1b ... ";
    bool local_status = skinning_use_case_1b(parallel_machine);
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case Skinning 2 ... ";
    bool local_status = skinning_use_case_2(parallel_machine);
    printStatus(local_status);
    status = status && local_status;
  }

  bool collective_result = use_case::print_status(parallel_machine, status);
  int return_code = collective_result ? 0 : -1;

  stk::parallel_machine_finalize();

  return return_code;
}
