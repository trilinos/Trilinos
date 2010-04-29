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


