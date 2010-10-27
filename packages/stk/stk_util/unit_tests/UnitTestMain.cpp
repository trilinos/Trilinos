/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

using std::abort;

// FIXME: We currently have both CppUnit and Gtest unit-tests that rely
// on this main. We need to convert everything to Gtest at some point;
// that would mean the only code needed below would be:
// STKUNIT_MAIN(argc,argv);

int
main(
  int           argc,
  char **       argv)
{
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) {
    std::cerr << "MPI_Init FAILED" << std::endl ;
    abort();
  }

  // cppunit block
  {
    std::string test_names;

    for (int i = 0; i < argc; ++i) {
      if (std::string(argv[i]) == "-test") {
        if (i + 1 < argc)
          test_names = argv[i + 1];
        else 
          std::cout << "Running all tests" << std::endl;
      }
    }
      
    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry & registry = CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest( registry.makeTest() );
    runner.run(test_names);
  }

  // gtest block
  bool result;
  {
    std::cout << "Running main() from gtest_main.cc\n";
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }

  MPI_Finalize();

  return result;
}
