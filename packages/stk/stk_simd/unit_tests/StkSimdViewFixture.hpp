#ifndef STK_KOKKOS_H
#define STK_KOKKOS_H

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

class StkSimdViewFixture : public ::testing::Test {

 protected:
  
  // virtual void SetUp() {
  static void SetUpTestCase() { 
    std::cout << "STK_Kokkos: Setting up fixtures." << std::endl;
    Kokkos::initialize();
  }
  
  //virtual void TearDown() {
  static void TearDownTestCase() {
    Kokkos::finalize_all();
    printf("STK_Kokkos: Tearing down fixture.\n");
  }  
};

// helpful verification function
template <typename T>
T running_sum(T n) {
  if (n==0.0) {
    return 0.0;
  }
  return running_sum(n-1)+n;
}


#endif
