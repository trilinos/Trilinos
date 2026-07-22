// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEM_003_HPP
#define HS_PROBLEM_003_HPP

#include "ROL_NonlinearProgram.hpp"

namespace HS {

namespace HS_003 {
template<class Real>
  class Obj {
  public:
    template<class ScalarT> 
    ScalarT value( const std::vector<ScalarT> &x, Real &tol ) {
      ScalarT a = x[1]-x[0];
      return x[1]+a*a*1e-5;
    }
  };
}

template<class Real>
class Problem_003 : public ROL::NonlinearProgram<Real> {

  

  typedef ROL::Objective<Real>            OBJ;
  typedef ROL::Vector<Real>               V;
  typedef ROL::NonlinearProgram<Real>     NP;

public:

  Problem_003() : NP( dimension_x() ) {
    NP::setLower(1,0);
  }

  int dimension_x() { return 2; }

  const ROL::Ptr<OBJ> getObjective() { 
    return ROL::makePtr<ROL::Sacado_StdObjective<Real,HS_003::Obj>>();
  }

  const ROL::Ptr<const V> getInitialGuess() {
    Real x[] = {10.0,1.0};
    return NP::createOptVector(x);
  };
   
  bool initialGuessIsFeasible() { return true; }
  
  Real getInitialObjectiveValue() { 
    return Real(1.00081);
  }
 
  Real getSolutionObjectiveValue() {
    return Real(0);
  }

  ROL::Ptr<const V> getSolutionSet() {
    Real x[] = {0,0};
    return ROL::CreatePartitionedVector(NP::createOptVector(x));
  }
 
};

} // namespace HS

#endif // HS_PROBLEM_003_HPP
