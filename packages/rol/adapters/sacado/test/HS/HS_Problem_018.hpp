// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_018_HPP
#define HS_PROBLEM_018_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_018 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return 0.01*x[0]*x[0] + x[1]*x[1];
  }
};

template<class Real> 
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = x[0]*x[1] - 25;
    c[1] = x[0]*x[0] + x[1]*x[1] - 25;
  }
};
}


template<class Real> 
class Problem_018 : public ROL::NonlinearProgram<Real> {

  
  
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;
  typedef ROL::NonlinearProgram<Real>     NP;

public:

  Problem_018() : NP( dimension_x() ) {
    NP::setLower(0,-2.0);
    NP::setUpper(0,50.0);
    NP::setLower(1,0.0);
    NP::setUpper(1,50.0);
  }  

  int dimension_x()  { return 2; }
  int dimension_ci() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_018::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_018::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {2,2};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(4.04);
  }
 
  Real getSolutionObjectiveValue() {
    return 5.0;
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {std::sqrt(250),std::sqrt(2.5)};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_018_HPP
