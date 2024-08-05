// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
