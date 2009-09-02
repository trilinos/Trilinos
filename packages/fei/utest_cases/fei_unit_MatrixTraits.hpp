#ifndef _fei_unit_MatrixTraits_hpp_
#define _fei_unit_MatrixTraits_hpp_

#include "fei_unit_testcase.hpp"

class test_mtraits : public fei::unit::testcase {
 public:
  test_mtraits(){}
  ~test_mtraits(){}

  bool run(MPI_Comm comm);
};

#endif

