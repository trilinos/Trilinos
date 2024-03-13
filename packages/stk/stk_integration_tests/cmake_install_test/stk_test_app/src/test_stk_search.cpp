#include "test_stk_search.hpp"

#include <iostream>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>

namespace test_stk_lib {

void test_stk_search()
{
  std::cout << "stk_search installation test, SearchMethod: " << stk::search::SearchMethod::KDTREE << std::endl;
}

}

