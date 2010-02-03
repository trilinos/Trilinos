
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <utility>

#include <gtest/gtest.h>

#include <stk_util/parallel/Parallel.hpp>

int main(int argc, char **argv)
{
#if defined( STK_HAS_MPI )
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) {
    std::cerr << "MPI_Init FAILED" << std::endl ;
    std::abort();
  }
  std::cout << "Running main() from gtest_main.cc\n";
#endif

  testing::InitGoogleTest(&argc, argv);

  bool result = RUN_ALL_TESTS();

#if defined( STK_HAS_MPI )
  MPI_Finalize();
#endif

  return result;
}


