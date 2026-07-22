// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_004_HPP
#define HS_PROBLEM_004_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_004 {
template<class Real> 
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    ScalarT a = x[0]+1.0;
    return a*a*a/3.0 + x[1];
  }
};

}

template<class Real> 
class Problem_004 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::Vector<Real>            V;
  typedef ROL::Objective<Real>         OBJ;
  typedef ROL::NonlinearProgram<Real>  NP;


public:

  Problem_004() : NP( dimension_x() ) {
    NP::setLower(0,1);
    NP::setLower(1,0); 
  }

  int dimension_x() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_004::Obj>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {1.125,0.125};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(3.323568);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(8.0/3.0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {1,0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

}

#endif // HS_PROBLEM_004_HPP
