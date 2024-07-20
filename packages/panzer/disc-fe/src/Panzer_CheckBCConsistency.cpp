// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_CheckBCConsistency.hpp"
#include "Panzer_BC.hpp"

namespace panzer {

  void checkBCConsistency(const std::vector<std::string>& element_block_names,
                          const std::vector<std::string>& sideset_names,
                          const std::vector<panzer::BC>& bcs)
  {
    for (const auto& bc : bcs) {
      const auto eb_search = std::find(element_block_names.begin(),
                                       element_block_names.end(),
                                       bc.elementBlockID());

      TEUCHOS_TEST_FOR_EXCEPTION(eb_search == element_block_names.end(),
                                 std::runtime_error,
                                 "ERROR: the element block \"" << bc.elementBlockID() 
                                 << "\" for boundary condition \"" << bc.bcID() 
                                 << "\" does not exist in the mesh.");

      const auto ss_search = std::find(sideset_names.begin(),
                                       sideset_names.end(),
                                       bc.sidesetID());

      TEUCHOS_TEST_FOR_EXCEPTION(ss_search == sideset_names.end(),
                                 std::runtime_error,
                                 "ERROR: the element block \"" << bc.sidesetID()
                                 << "\" for boundary condition \"" << bc.bcID()
                                 << "\" does not exist in the mesh.");
    }
  }

}
