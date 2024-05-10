#include "test_stk_search.hpp"

#include <iostream>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>

namespace test_stk_lib {

void test_stk_search()
{
  if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
    std::cout << "stk_search installation test, SearchMethod: " << stk::search::SearchMethod::KDTREE << std::endl;
  }
}

}

