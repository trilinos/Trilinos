#ifndef _fei_unit_ReverseMapper_hpp_
#define _fei_unit_ReverseMapper_hpp_

#include <fei_unit_testcase.hpp>

class test_ReverseMapper : public fei::unit::testcase {
 public:
  test_ReverseMapper(){}
  ~test_ReverseMapper(){}

  bool run(MPI_Comm comm);
};

#endif

