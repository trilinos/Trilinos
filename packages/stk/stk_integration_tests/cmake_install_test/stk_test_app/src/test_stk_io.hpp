#ifndef TEST_STK_IO_HPP
#define TEST_STK_IO_HPP

#include <string>
#include <stk_util/parallel/Parallel.hpp>

namespace test_stk_lib {

void test_stk_io(stk::ParallelMachine comm, const std::string& meshSource, bool useAutoDecomp);

}

#endif

