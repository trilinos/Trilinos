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
#include <stk_percept/RunEnvironment.hpp>

namespace stk { 
  namespace percept { 
    namespace regression_tests {

      std::string s_working_directory = "./";

      int rtest_main(int argc, char **argv) 
      { 
        EXCEPTWATCH;

        RunEnvironment run_environment(&argc, &argv);
        run_environment.clp.setDocString("stk_percept regression tests options");

        // NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
        std::string range_mesh;

        run_environment.clp.setOption("range_mesh", &range_mesh,  " range mesh");
        run_environment.processCommandLine(&argc, &argv);
      
        if (true)
          {
            s_working_directory = run_environment.directory_opt;
            //std::cout << "tmp 0 s_working_directory = " << s_working_directory << std::endl;
          }

        bool result = true;

#if defined(__PGI) || defined(REDS)
#else
        testing::InitGoogleTest(&argc, argv);  
        //  bool result = 0;
        try {
          //TEST_geom_volume(run_environment.m_comm);
    
          result = RUN_ALL_TESTS(); 
          //std::cout << "tmp result = " << result << std::endl;
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
        //MPI_Finalize(); 
#endif
        //std::cout << "tmp 0 result = " << result << std::endl;

        return result;
      }

    }
  }
}

//#include "pyencore.h"
//#if !PY_PERCEPT
int main(int argc, char **argv) { 

  int res = stk::percept::regression_tests::rtest_main(argc, argv);
  //std::cout << "tmp res = " << res << std::endl;
  
  return res;
}
//#endif
