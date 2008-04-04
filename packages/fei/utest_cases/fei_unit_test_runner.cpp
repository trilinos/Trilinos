
#include "fei_iostream.hpp"
#include "fei_Exception.hpp"

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

void test_runner::run_tests(int numProcs, int localProc, MPI_Comm comm)
{
  std::vector<fei::unit::testcase*>& all_tests = get_testcontainer();

  std::vector<fei::unit::testcase*>::iterator
    iter = all_tests.begin(), iter_end = all_tests.end();

  bool all_tests_passed = true;

  for(; iter != iter_end; ++iter) {
    fei::unit::testcase& tstcase = *(*iter);

    tstcase.setup(comm);

    try {
      bool result = tstcase.run(comm);
      if (result != true) {
        all_tests_passed = false;
        FEI_COUT << "test FAILED" << FEI_ENDL;
      }
    }
    catch(fei::Exception& exc) {
      all_tests_passed = false;
      FEI_COUT << "test failed with message: " << exc.what() << FEI_ENDL;
    }
  }

  if (localProc == 0) {
    if (all_tests_passed) {
      FEI_COUT << "\n\nAll Tests passed.\n" << FEI_ENDL;
    }
    else {
      FEI_COUT << "\n\nAt least 1 test FAILED.\n" << FEI_ENDL;
    }
  }
}

}//namespace unit
}//namespace fei

