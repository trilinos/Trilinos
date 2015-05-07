#include <Kokkos_Core.hpp>

#include "test_macrodef.hpp"
#include "test_suite.hpp"

/// \file test_serial.hpp
/// \brief Test serial execution space
/// \author Kyungjoo Kim (kyukim@sandia.gov)

typedef Kokkos::Serial space_type;

using namespace std;
using namespace Example;

int g_funct_counter = 0;

int main(int argc, char *argv[]) {
  int r_val = 0;

  Kokkos::initialize();

  //__TestSuiteDoUnitTests__(float,int,unsigned int,Kokkos::Serial,void);
  //__TestSuiteDoUnitTests__(float,long,unsigned long,Kokkos::Serial,void);

  __TestSuiteDoUnitTests__(double,int,unsigned int,Kokkos::Serial,void);
  // __TestSuiteDoUnitTests__(double,long,unsigned long,Kokkos::Serial,void);

  // __TestSuiteDoUnitTests__(complex<float>,int,unsigned int,Kokkos::Serial,void);
  // __TestSuiteDoUnitTests__(complex<float>,long,unsigned long,Kokkos::Serial,void);

  // __TestSuiteDoUnitTests__(complex<double>,int,unsigned int,Kokkos::Serial,void);
  // __TestSuiteDoUnitTests__(complex<double>,long,unsigned long,Kokkos::Serial,void);

  Kokkos::finalize();


  string eval;
  __EVAL_STRING__(r_val, eval);
  cout << "Testing Kokkos::Serial::" << eval << endl;

  return r_val;
}
