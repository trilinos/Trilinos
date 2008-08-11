#ifndef _fei_unit_utils_hpp_
#define _fei_unit_utils_hpp_

#include "fei_unit_testcase.hpp"

class test_utils : public fei::unit::testcase {
 public:
  test_utils(){}
  ~test_utils(){}

  bool run(MPI_Comm comm);
};

#endif

