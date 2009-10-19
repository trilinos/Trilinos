#include <mpi.h>

#include <iostream>
#include <string>
#include <utility>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>


int
main(
  int           argc,
  char **       argv)
{
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) {
    std::cerr << "MPI_Init FAILED" << std::endl ;
    std::abort();
  }

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

  MPI_Finalize();

  return 0;
}
