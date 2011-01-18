/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>
#include <mpi.h>

// the redstorm platform doesn't like Google test
#ifndef REDS
#include <gtest/gtest.h>
#endif

#include <stk_percept/fixtures/Fixture.hpp>

namespace stk { 
  namespace adapt { 
    namespace unit_tests {

  //OptionMaskParser dw_option_mask("use case diagnostic writer");

  void TEST_geom_volume(const stk::ParallelMachine comm);

int utest_main(int argc, char **argv) { 

  //dw_option_mask.mask("search", use_case::LOG_SEARCH, "log search diagnostics");

  // junk - FIXME
  //myMain3();
  //myMain2();
  boost::program_options::options_description desc("stk_adapt unit tests options");
    
  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
#if 0
  std::string range_mesh;
  desc.add_options()
    ("range_mesh",    boost::program_options::value<std::string>(&range_mesh), " range mesh")
    ("offset",       boost::program_options::value<double>()->default_value(0.1), "transfer use case 3 offset" )
    //    ("dw", boost::program_options::value<std::string>(), dw_option_mask.describe().c_str())
    ("scale",        boost::program_options::value<double>()->default_value(0.0), "transfer use case 3 scale." )
    ;

  stk::get_options_description().add(desc);
#endif

  //use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  percept::RunEnvironment run_environment(&argc, &argv);

  //boost::program_options::variables_map &vm = stk::get_variables_map();  

  std::cout << "Running main() from gtest_main.cc" << std::endl;
  bool result = true;

#ifndef REDS
  testing::InitGoogleTest(&argc, argv);  
  //  bool result = 0;
  try {
    //TEST_geom_volume(run_environment.m_comm);
    
    result = (RUN_ALL_TESTS() == 0); 
  }
  catch ( const std::exception * X ) {
    std::cout << "  unexpected exception POINTER: " << X->what() << std::endl;
    //exit(1);
  }
  catch ( const std::exception & X ) {
    std::cout << " stk_adapt::unit_tests::main unexpected exception: " << X.what() << std::endl;
    //exit(1);
  }
  catch( ... ) {
    std::cout << "  ... exception" << std::endl;
    //exit(1);
  }

#endif

#if doMPI
  MPI_Finalize(); 
#endif

  if (result) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  } else {
    std::cout << "End Result: TEST FAILED" << std::endl;
  }
  return result;
}

    }
  }
}

#include <stk_percept/pyencore.h>
#if !PY_PERCEPT
int main(int argc, char **argv) { 

  return stk::adapt::unit_tests::utest_main(argc, argv);
}
#endif
