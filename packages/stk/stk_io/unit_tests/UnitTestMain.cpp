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
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

STKUNIT_MAIN(argc,argv)

//#include <cppunit/extensions/TestFactoryRegistry.h>
//#include <cppunit/ui/text/TestRunner.h>
//
//#include <iostream>
//#include <string>
//#include <cstdlib>
//#include <mpi.h>
//
//int main(int argc, char ** argv)
//{
//  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) {
//    std::cerr << "MPI_Init FAILED" << std::endl ;
//    std::abort();
//  }
//
//  {
//    std::string test_names;
//
//    for (int i = 0; i < argc; ++i) {
//      if (std::string(argv[i]) == "-test") {
//        if (i + 1 < argc)
//          test_names = argv[i + 1];
//        else 
//          std::cout << "Running all tests" << std::endl;
//      }
//    }
//      
//    CppUnit::TextUi::TestRunner runner;
//    CppUnit::TestFactoryRegistry & registry = CppUnit::TestFactoryRegistry::getRegistry();
//    runner.addTest( registry.makeTest() );
//    runner.run(test_names);
//  }
//
//  MPI_Finalize();
//
//  return 0;
//}
//
