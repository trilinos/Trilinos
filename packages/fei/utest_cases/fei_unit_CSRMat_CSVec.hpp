#ifndef _fei_unit_CSRMat_CSVec_hpp_
#define _fei_unit_CSRMat_CSVec_hpp_

#include "fei_unit_testcase.hpp"

class test_csvec : public fei::unit::testcase {
 public:
  test_csvec(){}
  ~test_csvec(){}

  bool run(MPI_Comm comm);
};

#endif

