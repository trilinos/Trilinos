// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Panzer_UtilityAlgs.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(utility_algs, reorder)
  {
    // test 1
    {
      std::vector<int> order = {2, 0, 1};
      std::vector<char> data = {'a', 'b', 'c'};

      // run reorder algorithm
      reorder(order,[&data](int a,int b) { 
          std::swap(data[a],data[b]);   
        });

      TEST_ASSERT(order==std::vector<int>({0,1,2}));
      TEST_ASSERT(data==std::vector<char>({'c','a','b'}));
    } 

    // test 2
    {
      std::vector<int> order = {0, 1};
      std::vector<char> data = {'b', 'c'};

      // run reorder algorithm
      reorder(order,[&data](int a,int b) { 
          std::swap(data[a],data[b]);   
        });

      TEST_ASSERT(order==std::vector<int>({0,1}));
      TEST_ASSERT(data==std::vector<char>({'b','c'}));
    } 

    // test 3
    {
      std::vector<int> order = {1, 0};
      std::vector<char> data = {'b', 'c'};

      // run reorder algorithm
      reorder(order,[&data](int a,int b) { 
          std::swap(data[a],data[b]);   
        });

      TEST_ASSERT(order==std::vector<int>({0,1}));
      TEST_ASSERT(data==std::vector<char>({'c','b'}));
    } 
  }

} // end namespace panzer
