// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_NOXSolver_Def.hpp"

namespace Piro {

// Explicit template instantiation
// NOX currently only supports Scalar = double

template class NOXSolver<double>;

} // namespace Piro
