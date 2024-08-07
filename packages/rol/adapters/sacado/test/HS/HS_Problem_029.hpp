// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_029_HPP
#define HS_PROBLEM_029_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_029 {
template<class Real> 
class Obj {
public:
  template<class ScalarT>
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return -x[0]*x[1]*x[2];
  }
};


template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c,
              const std::vector<ScalarT> &x,
              Real &tol ) {
    c[0] = -x[0]*x[0] -2*x[1]*x[1] -4*x[2]*x[2] + 48;
  }
};
}


template<class Real> 
class Problem_029 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;

public:

  Problem_029() : NP( dimension_x() ) {
    NP::noBound();
  }

  int dimension_x()  { return 3; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_029::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_029::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {1.0,1.0,1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(-1.0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-16*std::sqrt(2));
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real a = 4; Real b = 2*std::sqrt(2); Real c = 2;
    
    Real x1[] = { a, b, c};
    Real x2[] = { a,-b,-c};
    Real x3[] = {-a, b,-c};
    Real x4[] = {-a,-b, c};

    return ROL::CreatePartitionedVector(NP::createOptVector(x1),
                                        NP::createOptVector(x2), 
                                        NP::createOptVector(x3),
                                        NP::createOptVector(x4));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_029_HPP
