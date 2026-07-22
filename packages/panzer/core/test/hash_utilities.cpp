// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include "Panzer_HashUtils.hpp"
#include <unordered_map>

namespace panzer {

  TEUCHOS_UNIT_TEST(hash_utilities, pair_hash)
  {
    std::unordered_map<std::pair<int,int>,int,panzer::pair_hash> my_map;

    std::pair<int,int> p = std::make_pair(5,10);

    my_map[p] = 3;

    TEST_EQUALITY(my_map[p], 3);
  }

}
