// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Panzer_CheckBCConsistency.hpp"
#include "Panzer_BC.hpp"

namespace panzer_test {

  TEUCHOS_UNIT_TEST(bc, check_bc_consistency)
  {
    std::vector<std::string> element_block_names{"liquid","solid","gas"};
    std::vector<std::string> sideset_names{"left","right"};

    panzer::BC bc1(0,panzer::BCT_Dirichlet,"left","liquid","MyEquationSetName","MyStrategy");
    panzer::BC bc2(1,panzer::BCT_Neumann,"right","solid","MyEquationSetName","MyStrategy");
    panzer::BC bc_wrong_sideset_name(2,panzer::BCT_Neumann,"bogus","liquid","MyEquationSetName","MyStrategy");
    panzer::BC bc_wrong_element_block_name(3,panzer::BCT_Neumann,"left","bogus","MyEquationSetName","MyStrategy");

    std::vector<panzer::BC> bcs;
    bcs.push_back(bc1);
    bcs.push_back(bc2);

    panzer::checkBCConsistency(element_block_names,
                               sideset_names,
                               bcs);


    bcs.push_back(bc_wrong_sideset_name);
    TEST_THROW(panzer::checkBCConsistency(element_block_names,
                                          sideset_names,
                                          bcs),
               std::runtime_error);

    bcs.clear();
    bcs.push_back(bc1);
    bcs.push_back(bc2);
    bcs.push_back(bc_wrong_element_block_name);
    TEST_THROW(panzer::checkBCConsistency(element_block_names,
                                          sideset_names,
                                          bcs),
               std::runtime_error);
    

  }

}
