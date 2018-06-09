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

/** \file
    \brief  Contains definitions of test objective functions.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_TESTPROBLEMS_HPP
#define ROL_TESTPROBLEMS_HPP

#include "ROL_OptimizationProblem.hpp"

namespace ROL {

template<class Real>
class TestProblem {
public:
  virtual ~TestProblem(void) {}
  TestProblem(void) {}
  virtual Ptr<Objective<Real>>           getObjective(void) const = 0;
  virtual Ptr<Vector<Real>>              getInitialGuess(void) const = 0;
  virtual Ptr<Vector<Real>>              getSolution(const int i = 0) const = 0;
  virtual int                            getNumSolutions(void) const {
    return 1;
  }
  virtual Ptr<BoundConstraint<Real>>     getBoundConstraint(void) const {
    return nullPtr;
  }
  virtual Ptr<Constraint<Real>>          getEqualityConstraint(void) const {
    return nullPtr;
  }
  virtual Ptr<Vector<Real>>              getEqualityMultiplier(void) const {
    return nullPtr;
  }
  virtual Ptr<Constraint<Real>>          getInequalityConstraint(void) const {
    return nullPtr;
  }
  virtual Ptr<Vector<Real>>              getInequalityMultiplier(void) const {
    return nullPtr;
  }
  virtual Ptr<BoundConstraint<Real>>     getSlackBoundConstraint(void) const {
    return nullPtr;
  }

  void get(Ptr<OptimizationProblem<Real>> &problem,
           Ptr<Vector<Real>>              &x0,
           std::vector<Ptr<Vector<Real>>> &x) const {
    x0 = getInitialGuess()->clone(); x0->set(*getInitialGuess());
    x.resize(getNumSolutions());
    for (int i = 0; i < getNumSolutions(); ++i) {
      x[i] = getSolution(i)->clone(); x[i]->set(*getSolution(i));
    }
    
    problem = makePtr<OptimizationProblem<Real>>(getObjective(),
                                                 x0,
                                                 getBoundConstraint(),
                                                 getEqualityConstraint(),
                                                 getEqualityMultiplier(),
                                                 getInequalityConstraint(),
                                                 getInequalityMultiplier(),
                                                 getSlackBoundConstraint());
  }
};
} // namespace ROL

#endif
