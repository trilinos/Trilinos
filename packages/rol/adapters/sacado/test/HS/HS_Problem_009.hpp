// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_009_HPP
#define HS_PROBLEM_009_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_009 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    Real pi = 3.14159265359; 
    return std::sin(pi*x[0]/12.0)*std::cos(pi*x[1]/16.0);
  }
};

template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = 4.0*x[0] - 3.0*x[1];    
  }
};
}


template<class Real> 
class Problem_009 : public ROL::NonlinearProgram<Real> {

  
  
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;
  typedef ROL::NonlinearProgram<Real>   NP;

public:
 
  Problem_009() : NP( dimension_x() ) {
    NP::noBound(); 
  }
 
  int dimension_x()  { return 2; }
  int dimension_ce() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_009::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_009::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0,0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(0.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-0.5);
  }

  // x* = (12k-3,16k-4) for all integer k
  ROL::Ptr<const V> getSolutionSet() {

    Real x1[] = {-15,-20};
    Real x2[] = {-3, -4};
    Real x3[] = { 9, 12};

    return ROL::CreatePartitionedVector(NP::createOptVector(x1),
                                        NP::createOptVector(x2),
                                        NP::createOptVector(x3));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_009_HPP
