#ifndef _fei_unit_impl_utils_hpp_
#define _fei_unit_impl_utils_hpp_

#include "fei_unit_testcase.hpp"

class test_impl_utils : public fei::unit::testcase {
 public:
  test_impl_utils(){}
  ~test_impl_utils(){}

  bool run(MPI_Comm comm);
};

#endif

