// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//@HEADER
// ************************************************************************
//         copyright
// ************************************************************************
//@HEADER

#include <Zoltan2_Version.hpp>

/*! \file Zoltan2_Version.cpp
    \brief Implementation of a Trilinos convention.
*/

namespace Zoltan2 {

  std::string Zoltan2_Version() { 
    return("Zoltan2 in Trilinos " TRILINOS_VERSION_STRING);
  }

} // namespace Zoltan2 
