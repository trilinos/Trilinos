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

#ifdef KOKKOS_ENABLE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif


class MTK_Kokkos : public ::testing::Test {

 protected:
  
  virtual void SetUp() {
    std::cout << "MTK_Kokkos: Setting up fixture." << std::endl;
  }
  
  virtual void TearDown() {
    
  }
  
};

#endif
