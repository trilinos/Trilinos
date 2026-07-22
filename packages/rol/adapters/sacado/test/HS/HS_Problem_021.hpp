// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_021_HPP
#define HS_PROBLEM_021_HPP

#include "ROL_NonlinearProgram.hpp"


namespace HS {

namespace HS_021 {
template<class Real>
class Obj {
public:
  template<class ScalarT> 
  ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
    return 0.01*x[0]*x[0]+x[1]*x[1]-100;
  }
};

template<class Real>
class InCon {
public:
  template<class ScalarT>
  void value( std::vector<ScalarT> &c, 
              const std::vector<ScalarT> &x,
              Real &told ) {
    c[0] = 10*x[0]-x[1]-10;
  }

};

}

template<class Real> 
class Problem_021 : public ROL::NonlinearProgram<Real> {
 
  

  typedef ROL::NonlinearProgram<Real>     NP;
  typedef ROL::Vector<Real>               V;
  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Constraint<Real>           CON; 
private:
public:

  Problem_021() : NP( dimension_x() ) {
    NP::setLower(0,  2.0);
    NP::setUpper(0, 50.0);
    NP::setLower(1,-50.0);
    NP::setUpper(1, 50.0);
  }

  int dimension_x() { return 2; }
  int dimension_ci() { return 1; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_021::Obj>>();
  }

  const ROL::Ptr<CON> getInequalityConstraint() { 
    return ROL::makePtr<ROL::Sacado_StdConstraint<Real,HS_021::InCon>>();
  }  

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {-1.0,-1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return false; }
  
  Real getInitialObjectiveValue() { 
    return Real(-98.99);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(-99.96);
  }

  ROL::Ptr<const V> getSolutionSet() {
    const Real x[] = {2.0,0.0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

}

#endif // HS_PROBLEM_021_HPP
