
#ifndef _fei_unit_testcontainer_hpp_
#define _fei_unit_testcontainer_hpp_

#include "fei_unit_testcase.hpp"
#include <vector>

namespace fei {
namespace unit {

std::vector<testcase*>& get_testcontainer();

void destroy_tests(std::vector<testcase*>& tests);

}//namespace unit
}//namespace fei

#endif

