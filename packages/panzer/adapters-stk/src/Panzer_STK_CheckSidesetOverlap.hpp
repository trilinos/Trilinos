// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CHECK_SIDESET_OVERLAP_HPP
#define PANZER_CHECK_SIDESET_OVERLAP_HPP

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include <string>

namespace panzer_stk {

  /// Returns true if the sidesets overlap.
  bool checkSidesetOverlap(const std::string& side_a_name,
                           const std::string& side_b_name,
                           const panzer_stk::STK_Interface& mesh);
}

#endif
