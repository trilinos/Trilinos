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

namespace stk { namespace percept { namespace unit_tests {

    void test_shards_array();
    int testSweepMesher( stk::ParallelMachine parallel_machine );

  //OptionMaskParser dw_option_mask("use case diagnostic writer");

  void TEST_geom_volume(const stk::ParallelMachine comm);

int utest_main(int argc, char **argv) { 

  //dw_option_mask.mask("search", use_case::LOG_SEARCH, "log search diagnostics");

  // junk - FIXME
  //myMain3();
  //myMain2();
  bopt::options_description desc("stk_percept unit tests options");
    
  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
#if 0
  std::string range_mesh;
  desc.add_options()
    ("range_mesh",    bopt::value<std::string>(&range_mesh), " range mesh")
    ("offset",       bopt::value<double>()->default_value(0.1), "transfer use case 3 offset" )
    //    ("dw", boost::program_options::value<std::string>(), dw_option_mask.describe().c_str())
    ("scale",        bopt::value<double>()->default_value(0.0), "transfer use case 3 scale." )
    ;

  stk::get_options_description().add(desc);
#endif

  //use_case::UseCaseEnvironment use_case_environment(&argc, &argv);
  RunEnvironment run_environment(&argc, &argv);

  //bopt::variables_map &vm = stk::get_variables_map();  

  unsigned p_size = stk::parallel_machine_size(run_environment.m_comm);
  bool isParallel = p_size > 0; // FIXME
  if (!isParallel)
    {
      stk::percept::unit_tests::test_shards_array();

      // run a few tests of SweepMesher - these aren't unit tests, and just produce some output FIXME: move to use_cases directory?
      int res1 =  stk::percept::unit_tests::testSweepMesher(run_environment.m_comm);
      std::cout << "testSweepMesher result= " << res1 << std::endl;
    }

  std::cout << "Running main() from gtest_main.cc" << std::endl;
  bool result = true;

#ifndef REDS
  testing::InitGoogleTest(&argc, argv);  
  //  bool result = 0;
  try {
    //TEST_geom_volume(run_environment.m_comm);
    
    result = RUN_ALL_TESTS(); 
  }
  catch ( const std::exception * X ) {
    std::cout << "  unexpected exception POINTER: " << X->what() << std::endl;
    //exit(1);
  }
  catch ( const std::exception & X ) {
    std::cout << " stk_percept::unit_tests::main unexpected exception: " << X.what() << std::endl;
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

  return result;
}

}}}

#include <stk_percept/pyencore.h>
#if !PY_PERCEPT
int main(int argc, char **argv) { 

  return stk::percept::unit_tests::utest_main(argc, argv);
}
#endif
