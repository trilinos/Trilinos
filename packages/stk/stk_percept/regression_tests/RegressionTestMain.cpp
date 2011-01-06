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
  namespace percept { 
    namespace regression_tests {

    std::string s_working_directory = "./";

    int rtest_main(int argc, char **argv) { 

      //dw_option_mask.mask("search", use_case::LOG_SEARCH, "log search diagnostics");

      // junk - FIXME
      //myMain3();
      //myMain2();
      bopt::options_description desc("stk_percept regression tests options");
    
      // NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
#if 0
      desc.add_options()
        ("range_mesh",    bopt::value<std::string>(&range_mesh), " range mesh")
        ("offset",       bopt::value<double>()->default_value(0.1), "transfer use case 3 offset" )
        //    ("dw", boost::program_options::value<std::string>(), dw_option_mask.describe().c_str())
        ("scale",        bopt::value<double>()->default_value(0.0), "transfer use case 3 scale." )
        ;

#endif
      std::string range_mesh;
      desc.add_options()
        ("range_mesh",    bopt::value<std::string>(&range_mesh), " range mesh");

      
      stk::get_options_description().add(desc);

      RunEnvironment run_environment(&argc, &argv);

      bopt::variables_map &vm = stk::get_variables_map();  
      if (vm.count("directory"))
        {
          s_working_directory = vm["directory"].as<std::string>();
          //std::cout << "tmp 0 s_working_directory = " << s_working_directory << std::endl;
        }


      bool result = true;

#ifndef REDS
      testing::InitGoogleTest(&argc, argv);  
      //  bool result = 0;
      try {
        //TEST_geom_volume(run_environment.m_comm);
    
        result = RUN_ALL_TESTS(); 
      }
      catch ( const std::exception * X ) {
        std::cout << "RegressionTestMain::  unexpected exception POINTER: " << X->what() << std::endl;
        //exit(1);
      }
      catch ( const std::exception & X ) {
        std::cout << "RegressionTestMain:: stk_percept::regression_tests::main unexpected exception: " << X.what() << std::endl;
        //exit(1);
      }
      catch( ... ) {
        std::cout << "RegressionTestMain::  ... exception" << std::endl;
        //exit(1);
      }

#endif

#if doMPI
      MPI_Finalize(); 
#endif

      return result;
    }

    }
  }
}

//#include "pyencore.h"
//#if !PY_PERCEPT
int main(int argc, char **argv) { 

  //return 0;
  return stk::percept::regression_tests::rtest_main(argc, argv);
}
//#endif
