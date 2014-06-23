/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include <mesh/MeshUseCase_1.hpp>
#include <mesh/MeshUseCase_2.hpp>
#include <mesh/MeshUseCase_3.hpp>
#include <mesh/UseCase_ElementDeath.hpp>
#include <mesh/UseCase_Skinning.hpp>
#include <mesh/UseCase_ChangeOwner.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/DiagWriter.hpp>
#include <stk_mesh/base/EntityKey.hpp>

//----------------------------------------------------------------------

void printStatus(bool status)
{
  if (status) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }
}

int
main(
  int           argc,
  char **       argv)
{
  // Add my command line options to the option descriptions.
  boost::program_options::options_description desc("Use case options");
  desc.add_options()
    ("performance", "run performance test")
    ("mesh", boost::program_options::value<std::string>(), "run mesh file performance test");

  stk::get_options_description().add(desc);

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  stk::ParallelMachine parallel_machine = use_case_environment.m_comm;

  bool status = true;

  // Now call the use-case drivers based on command line options

  const bool single_process =
    stk::parallel_machine_size( parallel_machine ) <= 1 ;

  {

//
//    stk_use_cases::use_case_13_driver( parallel_machine );
//    stk::app::use_case_14_driver( parallel_machine, run_performance_case );
//    stk::app::use_case_23_driver( parallel_machine, run_performance_case );
//    stk::app::use_case_AD_driver( parallel_machine, run_performance_case );
//
//    // Use cases imported from stk_mesh/use_cases

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

  }

  bool collective_result = use_case::print_status(parallel_machine, status);
  return collective_result ? 0 : -1;
}
