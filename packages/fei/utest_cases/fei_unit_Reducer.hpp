#ifndef _fei_unit_Reducer_hpp_
#define _fei_unit_Reducer_hpp_

#include <fei_unit_testcase.hpp>

class test_Reducer : public fei::unit::testcase {
 public:
  test_Reducer(){}
  ~test_Reducer(){}

  bool run(MPI_Comm comm);
};

#endif

