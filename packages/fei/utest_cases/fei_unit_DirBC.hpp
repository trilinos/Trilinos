#ifndef _fei_unit_DirBC_hpp_
#define _fei_unit_DirBC_hpp_

#include <fei_unit_testcase.hpp>

class test_DirBC : public fei::unit::testcase {
 public:
  test_DirBC(){}
  ~test_DirBC(){}

  bool run(MPI_Comm comm);
};

#endif

