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

#pragma once
#ifndef ROL2_TYPEU_MAKE_LINESEARCH_H
#define ROL2_TYPEU_MAKE_LINESEARCH_H

#include "ROL2_TypeU_IterationScaling.hpp"
#include "ROL2_TypeU_PathBasedTargetLevel.hpp"
#include "ROL2_TypeU_BackTracking.hpp"
#include "ROL2_TypeU_CubicInterp.hpp"
#include "ROL2_TypeU_ScalarMinimizationLineSearch.hpp"

namespace ROL {
namespace ROL2 {
namespace TypeU {

template<typename Real>
inline Ptr<LineSearch<Real>> 
make_LineSearch(ParameterList &parlist) {
  using Type = typename LineSearch<Real>::Type;
  Type type = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Type","Cubic Interpolation");
  switch(type) {
    case Type::IterationScaling:     return makePtr<IterationScaling<Real>>(parlist);
    case Type::PathBasedTargetLevel: return makePtr<PathBasedTargetLevel<Real>>(parlist);
    case Type::BackTracking:         return makePtr<BackTracking<Real>>(parlist);
    case Type::CubicInterp:          return makePtr<CubicInterp<Real>>(parlist);
    case Type::Brents:
    case Type::GoldenSection:
    case Type::Bisection:            return makePtr<ScalarMinimizationLineSearch<Real>>(parlist);
    default:                                return nullPtr;
  }
}
} // namespace TypeU
} // namespace ROL2
} // namespace ROL

namespace ROL2 {
namespace TypeU {
using ROL::ROL2::TypeU::make_LineSearch;
} // namespace TypeU
} // namespace ROL2
#endif
