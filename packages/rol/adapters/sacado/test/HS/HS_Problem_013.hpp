// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_013_HPP
#define HS_PROBLEM_013_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_013 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]-2.0;
    return a*a + x[1]*x[1];  
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    ScalarT a = 1-x[0];
    c[0] = a*a*a - x[1];
  }
};
}


template<class Real> 
class Problem_013 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;
  typedef ROL::NonlinearProgram<Real>     NP;

public:

  Problem_013() : NP( dimension_x() ) {
    NP::setLower(0,0);
    NP::setLower(1,0);
  }  

  int dimension_x()  { return 2; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_013::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_013::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-2,-2};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(20);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(1);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {1,0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_013_HPP
