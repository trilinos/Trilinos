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

#ifndef ROL_OPTIMIZATIONPROBLEMFACTORY_H
#define ROL_OPTIMIZATIONPROBLEMFACTORY_H

#include "ROL_Ptr.hpp"
#include "ROL_OptimizationProblem.hpp"

/** @ingroup func_group
    \class ROL::OptimizationProblemFactory
    \brief Defines the pebbl OptimizationProblemFactory interface.

    ROL's OptimizationProblemFactory constructs a new (identical)
    instance of an optimization problem for use in pebbl.

    ---
*/


namespace ROL {

template <class Real>
class OptimizationProblemFactory {
public:
  virtual ~OptimizationProblemFactory(void) {}

  virtual Ptr<Objective<Real>>       buildObjective(void) = 0;
  virtual Ptr<Vector<Real>>          buildSolutionVector(void) = 0;

  virtual Ptr<BoundConstraint<Real>> buildBoundConstraint(void) {
    return nullPtr;
  }

  virtual Ptr<Constraint<Real>>      buildEqualityConstraint(void) {
    return nullPtr;
  }
  virtual Ptr<Vector<Real>>          buildEqualityMultiplier(void) {
    return nullPtr;
  }

  virtual Ptr<Constraint<Real>>      buildInequalityConstraint(void) {
    return nullPtr;
  }
  virtual Ptr<Vector<Real>>          buildInequalityMultiplier(void) {
    return nullPtr;
  }
  virtual Ptr<BoundConstraint<Real>> buildInequalityBoundConstraint(void) {
    return nullPtr;
  }

  virtual void update(void) {}

  Ptr<OptimizationProblem<Real>> build(void) {
    update();
    return makePtr<OptimizationProblem<Real>>(buildObjective(),
                                              buildSolutionVector(),
                                              buildBoundConstraint(),
                                              buildEqualityConstraint(),
                                              buildEqualityMultiplier(),
                                              buildInequalityConstraint(),
                                              buildInequalityMultiplier(),
                                              buildInequalityBoundConstraint());
  }

}; // class OptimizationProblemFactory

} // namespace ROL

#endif
