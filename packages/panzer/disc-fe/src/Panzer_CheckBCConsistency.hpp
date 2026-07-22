// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CHECK_BC_CONSISTENCY
#define PANZER_CHECK_BC_CONSISTENCY

#include <vector>
#include <string>

namespace panzer {

  class BC;

  void checkBCConsistency(const std::vector<std::string>& element_block_names,
                          const std::vector<std::string>& sideset_names,
                          const std::vector<panzer::BC>& bcs);

}

#endif
