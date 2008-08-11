#ifndef _fei_unit_test_runner_hpp_
#define _fei_unit_test_runner_hpp_

#include "fei_mpi.h"

namespace fei {
namespace unit {

class test_runner {
 public:
  test_runner();
  ~test_runner();

  int run_tests(int numProcs, int localProc, MPI_Comm comm);

};//class test_runner

}//namespace unit
}//namespace fei

#endif

