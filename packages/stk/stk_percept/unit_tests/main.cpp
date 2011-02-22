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
    namespace unit_tests {

      void test_shards_array();
      int testSweepMesher( stk::ParallelMachine parallel_machine );

      //OptionMaskParser dw_option_mask("use case diagnostic writer");

      void TEST_geom_volume(const stk::ParallelMachine comm);

      int utest_main(int argc, char **argv) { 

        for (int i = 0; i < argc; ++i) {
          const std::string s(argv[i]);
          std::cout << "tmp 0 argv["<<i<<"]= " << s << std::endl;
        }

        std::cout << "tmp 0 argv["<<0<<"]= " << argv[0] << std::endl;

        RunEnvironment run_environment(&argc, &argv, true);
        run_environment.clp.setDocString("stk_percept unit tests options");
        run_environment.processCommandLine(&argc, &argv);

        unsigned p_size = stk::parallel_machine_size(run_environment.m_comm);
        bool isParallel = p_size > 0; // FIXME
        std::cout << "isParallel = " << isParallel << std::endl;
        if (!isParallel)
          {
            std::cout << "running testSweepMesher ..." << std::endl;
            stk::percept::unit_tests::test_shards_array();

            // run a few tests of SweepMesher - these aren't unit tests, and just produce some output FIXME: move to use_cases directory?
            int res1 =  stk::percept::unit_tests::testSweepMesher(run_environment.m_comm);
            std::cout << "testSweepMesher result= " << res1 << std::endl;
          }

        std::cout << "Running main() from gtest_main.cc" << std::endl;
        bool result = true;

#if defined(__PGI) || defined(REDS)
#else
        testing::InitGoogleTest(&argc, argv);  
        //  bool result = 0;
        try {
          //TEST_geom_volume(run_environment.m_comm);
    
          result = RUN_ALL_TESTS(); 
          //std::cout << "stk::percept: result = " << result << std::endl;
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

  int success = stk::percept::unit_tests::utest_main(argc, argv);
  //std::cout << "stk::percept: success = " << success << std::endl;
  return success;
}
#endif
