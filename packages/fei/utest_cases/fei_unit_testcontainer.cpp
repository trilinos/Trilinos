
#include "fei_unit_testcontainer.hpp"

namespace fei {
namespace unit {

std::vector<testcase*>& get_testcontainer()
{
  static std::vector<testcase*> testcontainer;
  return testcontainer;
}

void destroy_tests(std::vector<testcase*>& tests)
{
  std::vector<testcase*>::iterator
    iter = tests.begin(), iter_end = tests.end();

  for(; iter!=iter_end; ++iter) {
    delete *iter;
  }
}

}//namespace unit
}//namespace fei

