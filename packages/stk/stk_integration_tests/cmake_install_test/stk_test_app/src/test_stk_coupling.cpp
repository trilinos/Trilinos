#include "test_stk_coupling.hpp"

#include <iostream>

#include <stk_util/Version.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include <stk_coupling/SplitComms.hpp>

namespace test_stk_lib {

void test_stk_coupling()
{
  std::cout << "stk_coupling installation test, STK version: " << stk::version_string()
            << ", Coupling-Version: "<< stk::util::get_local_max_coupling_version() << std::endl;
}

}

