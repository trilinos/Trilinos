
#include "fei_unit_testcontainer.hpp"

#include "fei_unit_CSRMat_CSVec.hpp"
#include "fei_unit_utils.hpp"

namespace fei {
namespace unit {

void register_tests()
{
  std::vector<testcase*>& all_tests = get_testcontainer();

  all_tests.push_back(new test_csvec);
  all_tests.push_back(new test_utils);

}//register_tests

}//namespace unit
}//namespace fei

