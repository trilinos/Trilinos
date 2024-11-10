// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  Intrepid_Version.cpp
    \brief Returns current version of Intrepid.
*/

#include "Intrepid2_ConfigDefs.hpp"

namespace Intrepid2 {

std::string Intrepid_Version() {
  return("Intrepid2 Version 1.0 / Trilinos 12.4 - September 2016");
}

} // namespace Intrepid2
