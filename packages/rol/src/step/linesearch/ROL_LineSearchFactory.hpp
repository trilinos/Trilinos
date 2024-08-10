// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINESEARCHFACTORY_H
#define ROL_LINESEARCHFACTORY_H

#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

#include "ROL_LineSearch.hpp"
#include "ROL_IterationScaling.hpp"
#include "ROL_PathBasedTargetLevel.hpp"
#include "ROL_BackTracking.hpp"
#include "ROL_CubicInterp.hpp"
#include "ROL_Bisection.hpp"
#include "ROL_GoldenSection.hpp"
#include "ROL_Brents.hpp"
#include "ROL_ScalarMinimizationLineSearch.hpp"

namespace ROL {
  template<class Real>
  inline ROL::Ptr<LineSearch<Real> > LineSearchFactory(ROL::ParameterList &parlist) {
    ELineSearch els = StringToELineSearch(
      parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Type","Cubic Interpolation"));
    switch(els) {
      case LINESEARCH_ITERATIONSCALING:     return ROL::makePtr<IterationScaling<Real>>(parlist);
      case LINESEARCH_PATHBASEDTARGETLEVEL: return ROL::makePtr<PathBasedTargetLevel<Real>>(parlist);
      case LINESEARCH_BACKTRACKING:         return ROL::makePtr<BackTracking<Real>>(parlist);
      case LINESEARCH_CUBICINTERP:          return ROL::makePtr<CubicInterp<Real>>(parlist);
//      case LINESEARCH_BISECTION:            return ROL::makePtr<Bisection<Real>>(parlist);
//      case LINESEARCH_BRENTS:               return ROL::makePtr<Brents<Real>>(parlist);
//      case LINESEARCH_GOLDENSECTION:        return ROL::makePtr<GoldenSection<Real>>(parlist);
      case LINESEARCH_BRENTS:
      case LINESEARCH_GOLDENSECTION:
      case LINESEARCH_BISECTION:            return ROL::makePtr<ScalarMinimizationLineSearch<Real>>(parlist);
      default:                              return ROL::nullPtr;
    }
  }
}

#endif
