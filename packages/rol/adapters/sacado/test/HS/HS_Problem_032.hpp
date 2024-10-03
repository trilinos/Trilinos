// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_032_HPP
#define HS_PROBLEM_032_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_032 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]+3*x[1]+x[2];
    ScalarT b = x[0]-x[1];
    return a*a + 4*b*b;
  }
};


template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {

    c[0] = 1.0 - x[0] - x[1] - x[2];

  }
};



template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {

    c[0] = 6*x[1] + 4*x[2] - x[0]*x[0]*x[0] - 3.0;

  }
};
} // namespace HS_032


template<class Real> 
class Problem_032 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

public:

  Problem_032() : NP( dimension_x() ) {
    NP::setLower(0, 0.0);
    NP::setLower(1, 0.0);
    NP::setLower(2, 0.0);
  }

  int dimension_x()  { return 3; }
  int dimension_ce() { return 1; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_032::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_032::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0.1,0.7,0.2};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(7.2);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(1.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = { 0., 0., 1.0};

    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_032_HPP
