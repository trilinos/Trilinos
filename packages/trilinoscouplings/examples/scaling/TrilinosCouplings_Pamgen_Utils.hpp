// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_Pamgen_Utils_hpp_
#define __TrilinosCouplings_Pamgen_Utils_hpp_
#include <iostream>

namespace TrilinosCouplings {
  // Error checking routine after a call to Create_Pamgen_Mesh
  // Dumps out error output to os and then throws if error occurred.
  void pamgen_error_check(std::ostream & os, long long cr_result);
}
#endif
