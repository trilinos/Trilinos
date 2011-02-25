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
#if !(defined(__PGI) || defined(REDS))
#include <gtest/gtest.h>
#endif

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/RunEnvironment.hpp>

namespace stk { 
  namespace adapt { 
    namespace regression_tests {

      std::string s_working_directory = "./";

      int rtest_main(int argc, char **argv) 
      { 
        bool debug_re = true;
        percept::RunEnvironment run_environment(&argc, &argv, debug_re);

        run_environment.clp.setDocString("stk_adapt regression tests options");

        // NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
        std::string range_mesh;
        run_environment.clp.setOption("range_mesh", &range_mesh, " range mesh");
      
        run_environment.processCommandLine(&argc, &argv);

        {
          s_working_directory = run_environment.directory_opt;
          if (debug_re) std::cout << "tmp 0 s_working_directory = " << s_working_directory << std::endl;
        }

        bool result = true;

#if !(defined(__PGI) || defined(REDS))
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
          std::cout << "RegressionTestMain:: stk_adapt::regression_tests::main unexpected exception: " << X.what() << std::endl;
          //exit(1);
        }
        catch( ... ) {
          std::cout << "RegressionTestMain::  ... exception" << std::endl;
          //exit(1);
        }

#endif

        return result;
      }

    }
  }
}

//#include "pyadapt.h"
//#if !PY_ADAPT
int main(int argc, char **argv) { 

  //return 0;
  return stk::adapt::regression_tests::rtest_main(argc, argv);
}
//#endif
