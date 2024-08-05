// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_034_HPP
#define HS_PROBLEM_034_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_034 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return -x[0];
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {

    c[0] = x[1] - std::exp(x[0]);
    c[1] = x[2] - std::exp(x[1]);

  }
};
} // namespace HS_034


template<class Real> 
class Problem_034 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

public:

  Problem_034() : NP( dimension_x() ) {
    NP::setLower(0,   0.0);
    NP::setUpper(0, 100.0);
    NP::setLower(1,   0.0);
    NP::setUpper(1, 100.0);
    NP::setLower(2,   0.0);
    NP::setUpper(2,  10.0);
  }

  int dimension_x()  { return 3; }
  int dimension_ci() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_034::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_034::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0.0,1.05,2.9};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(0.0);
  }
 
  Real getSolutionObjectiveValue() {
    return -std::log(std::log(10));
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = { std::log(std::log(10)), std::log(10), 10.0 };

    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_034_HPP
