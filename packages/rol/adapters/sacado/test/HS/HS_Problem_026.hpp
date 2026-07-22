// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_026_HPP
#define HS_PROBLEM_026_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_026 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]-x[1];
    ScalarT b = x[1]-x[2];
    return a*a + b*b*b*b;
  }
};

template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    ScalarT a = 1.0+x[1];
    c[0] = a*a*x[0] + std::pow(x[2],4)-3.0;    
  }
};
}


template<class Real> 
class Problem_026 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>   NP;
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;

public:

  Problem_026() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x()  { return 3; }
  int dimension_ce() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_026::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_026::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-2.6,2.0,2.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(21.16);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real alpha = std::sqrt(139.0/108.0);
    Real beta = 61.0/54.0;
    Real a = std::pow(alpha-beta,1.0/3.0) - std::pow(alpha+beta,1.0/3.0) - 2.0/3.0;
    Real x1[] = {1.0,1.0,1.0};
    Real x2[] = {a,a,a};
    return ROL::CreatePartitionedVector(NP::createOptVector(x1),
                                        NP::createOptVector(x2));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_026_HPP
