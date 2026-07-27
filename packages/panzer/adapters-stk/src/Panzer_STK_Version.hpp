// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_VERSION_HPP
#define PANZER_STK_VERSION_HPP

#include <string>
#include "PanzerAdaptersSTK_config.hpp"

namespace panzer_stk {

  /// \brief Returns a human-readable string identifying the panzer_stk (STK mesh adapter) version being used.
  std::string version();

}
#endif
