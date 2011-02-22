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

        bool debug_re=true;
        percept::RunEnvironment run_environment(&argc, &argv, debug_re);

        run_environment.clp.setDocString("stk_adapt unit tests options");
    
        run_environment.processCommandLine(&argc, &argv);

        std::cout << "Running main() from gtest_main.cc" << std::endl;
        int exitcode = 0;

#if defined(__PGI) || defined(REDS)
#else
        testing::InitGoogleTest(&argc, argv);  
        try {
          //TEST_geom_volume(run_environment.m_comm);
    
          exitcode = RUN_ALL_TESTS(); 
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

        if (exitcode == 0) {
          std::cout << "End Result: TEST PASSED" << std::endl;
        } else {
          std::cout << "End Result: TEST FAILED" << std::endl;
        }
        return exitcode;
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
