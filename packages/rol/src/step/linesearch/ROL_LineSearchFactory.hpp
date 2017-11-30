// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_LINESEARCHFACTORY_H
#define ROL_LINESEARCHFACTORY_H

#include "ROL_Types.hpp"

#include "Teuchos_ParameterList.hpp"
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
  inline ROL::Ptr<LineSearch<Real> > LineSearchFactory(Teuchos::ParameterList &parlist) {
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
