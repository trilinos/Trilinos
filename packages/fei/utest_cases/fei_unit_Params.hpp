#ifndef _fei_unit_Params_hpp_
#define _fei_unit_Params_hpp_

#include <fei_unit_testcase.hpp>

class test_Params : public fei::unit::testcase {
 public:
  test_Params(){}
  ~test_Params(){}

  bool run(MPI_Comm comm);
};

#endif

