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

// Prototypes for use-case driver functions
// (the functions live in UseCase_*.cpp).

namespace stk_use_cases {
void use_case_13_driver( stk::ParallelMachine );
}//namespace stk_use_cases

namespace stk {
namespace app {

void use_case_14_driver( stk::ParallelMachine, bool run_performance_test );
void use_case_23_driver( stk::ParallelMachine, bool run_performance_test );
void use_case_AD_driver( stk::ParallelMachine, bool run_performance_test );

} // namespace app
} // namespace stk

//----------------------------------------------------------------------

int
main(
  int           argc,
  char **       argv)
{
  // Add my command line options to the option descriptions.
  boost::program_options::options_description desc("Use case options");
  desc.add_options()
    ("performance", "run performance test")
    ( "use_case_13" , "use case 13" )
    ( "use_case_14" , "use case 14" )
    ( "use_case_23" , "use case 23" )
    ( "use_case_AD" , "use case AD" )
    ("mesh", boost::program_options::value<std::string>(), "run mesh file performance test");

  stk::get_options_description().add(desc);

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  boost::program_options::variables_map &vm = stk::get_variables_map();
  
  stk::ParallelMachine parallel_machine = use_case_environment.m_comm;

  // Now call the use-case drivers based on command line options

  bool run_performance_case = vm.count("performance") != 0;

  if (vm.count("use_case_13")) {
    stk_use_cases::use_case_13_driver( parallel_machine );
  }
  else if (vm.count("use_case_14")) {
    stk::app::use_case_14_driver( parallel_machine, run_performance_case );
  }
  else if (vm.count("use_case_23")) {
    stk::app::use_case_23_driver( parallel_machine, run_performance_case );
  }
  else if (vm.count("use_case_AD")) {
    stk::app::use_case_AD_driver( parallel_machine, run_performance_case );
  }
  else {
    stk_use_cases::use_case_13_driver( parallel_machine );
    stk::app::use_case_14_driver( parallel_machine, run_performance_case );
    stk::app::use_case_23_driver( parallel_machine, run_performance_case );
    stk::app::use_case_AD_driver( parallel_machine, run_performance_case );
  }

  // If we've made it this far, the use-case has passed
  use_case::print_status(parallel_machine, true);

  return 0;
}
