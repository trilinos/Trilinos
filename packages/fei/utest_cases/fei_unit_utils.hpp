#ifndef _fei_unit_utils_hpp_
#define _fei_unit_utils_hpp_

#include "fei_unit_testcase.hpp"

class test_utils : public fei::unit::testcase {
  bool run(MPI_Comm comm);
};

#endif

