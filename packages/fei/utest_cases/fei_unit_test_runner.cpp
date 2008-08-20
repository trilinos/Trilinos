
#include "fei_macros.hpp"
#include "fei_mpi.h"
#include "fei_iostream.hpp"

#include "fei_unit_test_runner.hpp"
#include "fei_unit_register_tests.hpp"

#include "fei_unit_testcase.hpp"
#include "fei_unit_testcontainer.hpp"

namespace fei {
namespace unit {

test_runner::test_runner()
{
  register_tests();
}

test_runner::~test_runner()
{
  destroy_tests(get_testcontainer());
}

int test_runner::run_tests(int numProcs, int localProc, MPI_Comm comm)
{
  std::vector<fei::unit::testcase*>& all_tests = get_testcontainer();

  std::vector<fei::unit::testcase*>::iterator
    iter = all_tests.begin(), iter_end = all_tests.end();

  int tests_failed = 0;

  for(; iter != iter_end; ++iter) {
    fei::unit::testcase& tstcase = *(*iter);

    tstcase.setup(comm);

    try {
      bool result = tstcase.run(comm);
      if (result != true) {
        tests_failed = 1;
        FEI_COUT << "test FAILED" << FEI_ENDL;
      }
    }
    catch(std::runtime_error& exc) {
      tests_failed = 1;
      FEI_COUT << "test failed with message: " << exc.what() << FEI_ENDL;
    }
  }

  int global_tests_failed = 0;
#ifdef FEI_SER
  global_tests_failed = tests_failed;
#else
  MPI_Allreduce(&tests_failed, &global_tests_failed, 1,
                MPI_INT, MPI_MAX, comm);
#endif

  if (localProc == 0) {
    if (global_tests_failed == 0) {
      FEI_COUT << "\n\nAll Tests passed.\n" << FEI_ENDL;
    }
    else {
      FEI_COUT << "\n\nAt least 1 test FAILED.\n" << FEI_ENDL;
    }
  }

  return global_tests_failed;
}

}//namespace unit
}//namespace fei

