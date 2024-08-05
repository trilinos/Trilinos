// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_010_HPP
#define HS_PROBLEM_010_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_010 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return x[0]-x[1];
  }
};

template<class Real>
class EqCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = -3*x[0]*x[0] + 2*x[0]*x[1] - x[1]*x[1] + 1; 
  }
};
}


template<class Real> 
class Problem_010 : public ROL::NonlinearProgram<Real> {

  
  
  typedef ROL::Vector<Real>             V;
  typedef ROL::Objective<Real>          OBJ;
  typedef ROL::Constraint<Real>         CON;
  typedef ROL::NonlinearProgram<Real>   NP;

public:
 
  Problem_010() : NP( dimension_x() ) {
    NP::noBound(); 
  }
 
  int dimension_x()  { return 2; }
  int dimension_ce() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_010::Obj>>();
  }

  const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_010::EqCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-10,10};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(-20.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-1);
  }

  ROL::Ptr<const V> getSolutionSet() {

    Real x[] = {0,1};

    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_010_HPP
