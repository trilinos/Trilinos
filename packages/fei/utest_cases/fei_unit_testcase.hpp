#ifndef _fei_unit_testcase_hpp_
#define _fei_unit_testcase_hpp_

#include "fei_mpi.h"

namespace fei {
namespace unit {

class testcase {
 public:
  testcase();
  virtual ~testcase(){}

  virtual void setup(MPI_Comm /*comm*/){}

  //return true if successful, false otherwise
  virtual bool run(MPI_Comm /*comm*/) = 0;

  virtual void teardown(){}
};//class testcase

}//namespace unit
}//namespace fei

#endif

