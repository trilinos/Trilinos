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
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_StlMap_Utilities.hpp"
#include <map>
#include <string>

namespace panzer {

  TEUCHOS_UNIT_TEST(stlmap_utilities, all)
  {
    
    std::map<std::string,int> my_map;
    my_map["a"] = 1;
    my_map["b"] = 2;
    my_map["c"] = 3;
    
    TEST_ASSERT(panzer::getEntry(my_map,"a") == my_map.find("a")->second);
  }

}
