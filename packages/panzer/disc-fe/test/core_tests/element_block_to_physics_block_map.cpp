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

#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(ElementBlockToPhsyicsBlockMap,basic)
  {
    
    Teuchos::ParameterList p;
    p.set("eblock-0_0", "fluid");
    p.set("eblock-0_1", "fluid");
    p.set("eblock-1_0", "solid");
    p.set("eblock-1_1", "solid");
    
    std::map<std::string,std::string> b_to_p;
    
    panzer::buildBlockIdToPhysicsIdMap(b_to_p, p);

    TEST_EQUALITY(b_to_p["eblock-0_0"], "fluid");
    TEST_EQUALITY(b_to_p["eblock-0_1"], "fluid");
    TEST_EQUALITY(b_to_p["eblock-1_0"], "solid");
    TEST_EQUALITY(b_to_p["eblock-1_1"], "solid");
  }

}
