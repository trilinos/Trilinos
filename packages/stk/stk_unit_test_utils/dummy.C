// header only lib, but need a library for the jamfile...
#include "ParallelGtestOutput.hpp"

namespace stk {
namespace unit_test_util {

int stk_unit_test_util_dummy_function()
{
  create_parallel_output(0);
  return 42;
}

}
}


