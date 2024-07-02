// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_Utilities.hpp"

namespace Ifpack2 {
namespace Details {

  std::string canonicalize(const std::string& precType) {
    // precTypeUpper is the upper-case version of precType.
    std::string precTypeUpper (precType);
    if (precTypeUpper.size () > 0) {
      for (size_t k = 0; k < precTypeUpper.size (); ++k) {
        precTypeUpper[k] = ::toupper(precTypeUpper[k]);
      }
    }
    return precTypeUpper;
  }

} // namespace Details
} // namespace Ifpack2
