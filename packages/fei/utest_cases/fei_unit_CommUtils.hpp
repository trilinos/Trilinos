#ifndef _fei_unit_CommUtils_hpp_
#define _fei_unit_CommUtils_hpp_

#include <fei_unit_testcase.hpp>

class test_CommUtils : public fei::unit::testcase {
 public:
  test_CommUtils(){}
  ~test_CommUtils(){}

  bool run(MPI_Comm comm);
};

#endif

