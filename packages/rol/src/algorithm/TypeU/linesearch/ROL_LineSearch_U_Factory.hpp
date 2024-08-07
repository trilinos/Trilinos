// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINESEARCH_U_FACTORY_H
#define ROL_LINESEARCH_U_FACTORY_H

#include "ROL_IterationScaling_U.hpp"
#include "ROL_PathBasedTargetLevel_U.hpp"
#include "ROL_BackTracking_U.hpp"
#include "ROL_CubicInterp_U.hpp"
#include "ROL_ScalarMinimizationLineSearch_U.hpp"

namespace ROL {
template<typename Real>
inline Ptr<LineSearch_U<Real>> LineSearchUFactory(ParameterList &parlist) {
  ELineSearchU els = StringToELineSearchU(
    parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Type","Cubic Interpolation"));
  switch(els) {
    case LINESEARCH_U_ITERATIONSCALING:     return makePtr<IterationScaling_U<Real>>(parlist);
    case LINESEARCH_U_PATHBASEDTARGETLEVEL: return makePtr<PathBasedTargetLevel_U<Real>>(parlist);
    case LINESEARCH_U_BACKTRACKING:         return makePtr<BackTracking_U<Real>>(parlist);
    case LINESEARCH_U_CUBICINTERP:          return makePtr<CubicInterp_U<Real>>(parlist);
    case LINESEARCH_U_BRENTS:
    case LINESEARCH_U_GOLDENSECTION:
    case LINESEARCH_U_BISECTION:            return makePtr<ScalarMinimizationLineSearch_U<Real>>(parlist);
    default:                                return nullPtr;
  }
}
} // namespace ROL

#endif
