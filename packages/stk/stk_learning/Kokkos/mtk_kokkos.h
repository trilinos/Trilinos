#ifndef MTK_KOKKOS_H
#define MTK_KOKKOS_H

#include <gtest/gtest.h>                // for AssertHelper, etc
#include <iomanip>                      // for operator<<
#include <iostream>                     // for basic_ostream::operator<<, etc
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <cmath>

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_HAVE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif

class MTK_Kokkos : public ::testing::Test {

 protected:
  
  virtual void SetUp() {
    std::cout << "MTK_Kokkos: Setting up fixture." << std::endl;
    
#ifdef KOKKOS_HAVE_OPENMP
    Kokkos::InitArguments init_args;
    char *num_threads_string = std::getenv("OMP_NUM_THREADS");
    init_args.num_threads = std::atoi(num_threads_string);
    Kokkos::initialize(init_args);
    Kokkos::OpenMP::print_configuration( std::cout );
#else
    Kokkos::initialize();
#endif
  }
  
  virtual void TearDown() {
    Kokkos::finalize_all();
    
#if _OPENMP
    omp_set_num_threads(1);
    ASSERT_EQ( 1 , omp_get_max_threads() );
#endif
    
  }
  
};

#endif
