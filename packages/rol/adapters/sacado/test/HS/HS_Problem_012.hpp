// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_012_HPP
#define HS_PROBLEM_012_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS { 

namespace HS_012 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return 0.5*x[0]*x[0] + x[1]*x[1] - x[0]*x[1] - 7*x[0] - 7*x[1];
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT> 
  void value( std::vector<ScalarT> &c, const std::vector<ScalarT> &x, Real &tol ) {
    c[0] = 25 - 4*x[0]*x[0] - x[1]*x[1];
  }
};
}


template<class Real> 
class Problem_012 : public ROL::NonlinearProgram<Real> {

  
  
  typedef ROL::Vector<Real>               V;
  typedef ROL::PartitionedVector<Real>    PV;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON;
  typedef ROL::NonlinearProgram<Real>     NP;


public:

  Problem_012() : NP( dimension_x() ) {
    NP::noBound();
  }  

  int dimension_x()  { return 2; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_012::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_012::InCon>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {0,0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(0);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-30);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {2,3};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_012_HPP
