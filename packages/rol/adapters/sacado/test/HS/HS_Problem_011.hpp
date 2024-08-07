// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_011_HPP
#define HS_PROBLEM_011_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_011 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]-5.0;
    return a*a+x[1]*x[1]-25.0;
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
    void value( std::vector<ScalarT> &c,
                const std::vector<ScalarT> &x,
                Real &tol ) {
    c[0] = -x[0]*x[0]+x[1];
  }
};
}



template<class Real> 
class Problem_011 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;
  typedef ROL::NonlinearProgram<Real>     NP;


public:

  Problem_011() : NP( dimension_x() ) {
    NP::noBound();
  }  

  int dimension_x()  { return 2; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_011::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_011::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {4.9,0.1};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(-24.98);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-8.498464223);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real a = 7.5*std::sqrt(6) + std::sqrt(338.5);
    Real x[] = {(a-1/a)/std::sqrt(6),
                (a*a-2+1/(a*a))/6.0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_011_HPP
