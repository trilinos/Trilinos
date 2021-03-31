// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL2_TYPEU_ALGORITHM_DEF_HPP
#define ROL2_TYPEU_ALGORITHM_DEF_HPP

#include "ROL_Types.hpp"
#include "ROL_ReduceLinearConstraint.hpp"
#include "ROL_ValidParameters.hpp"

namespace ROL {
namespace ROL2 {
namespace TypeU {

template<typename Real>
Algorithm<Real>::Algorithm()
  : status_(makePtr<CombinedStatusTest<Real>>()),
    state_(makePtr<AlgorithmState_U<Real>>()) {
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>());
}

template<typename Real>
void Algorithm<Real>::initialize(const Vector<Real> &x, const Vector<Real> &g) {
  if (state_->iterateVec == nullPtr) {
    state_->iterateVec = x.clone();
  }
  state_->iterateVec->set(x);
  if (state_->stepVec == nullPtr) {
    state_->stepVec = x.clone();
  }
  state_->stepVec->zero();
  if (state_->gradientVec == nullPtr) {
    state_->gradientVec = g.clone();
  }
  state_->gradientVec->set(g);
  if (state_->minIterVec == nullPtr) {
    state_->minIterVec = x.clone();
  }
  state_->minIterVec->set(x);
  state_->minIter = state_->iter;
  state_->minValue = state_->value;
}

template<typename Real>
void Algorithm<Real>::setStatusTest(const Ptr<StatusTest<Real>> &status,
                                      bool combineStatus) {
  if (!combineStatus) { // Do not combine status tests
    status_->reset();
  }
  status_->add(status); // Add user-defined StatusTest
}

template<typename Real>
std::vector<std::string> Algorithm<Real>::run( NewOptimizationProblem<Real> &problem,
                                                 std::ostream                 &outStream ) {
  if (problem.getProblemType() == TYPE_U) {
    std::vector<std::string> output = run(*problem.getPrimalOptimizationVector(),
                                          *problem.getDualOptimizationVector(),
                                          *problem.getObjective(),
                                          outStream);
    problem.finalizeIteration();
    return output;
  }
  else {
    throw Exception::NotImplemented(">>> ROL::Algorithm::run : Optimization problem is not Type U!");
  }
}

template<typename Real>
std::vector<std::string> Algorithm<Real>::run( Vector<Real>    &x,
                                                 Objective<Real> &obj,
                                                 std::ostream    &outStream ) {
  return run(x,x.dual(),obj,outStream);
}

template<typename Real>
std::vector<std::string> Algorithm<Real>::run( Vector<Real>     &x,
                                                 Objective<Real>  &obj,
                                                 Constraint<Real> &linear_con,
                                                 Vector<Real>     &linear_mul,
                                                 std::ostream     &outStream ) {
  return run(x,x.dual(),obj,linear_con,linear_mul,linear_mul.dual(),outStream);
}

template<typename Real>
std::vector<std::string> Algorithm<Real>::run( Vector<Real>       &x,
                                                 const Vector<Real> &g,
                                                 Objective<Real>    &obj,
                                                 Constraint<Real>   &linear_con,
                                                 Vector<Real>       &linear_mul,
                                                 const Vector<Real> &linear_c,
                                                 std::ostream       &outStream ) {
  Ptr<Vector<Real>> xfeas = x.clone(); xfeas->set(x);
  ReduceLinearConstraint<Real> rlc(makePtrFromRef(linear_con),xfeas,makePtrFromRef(linear_c));
  Ptr<Vector<Real>> s = x.clone(); s->zero();
  std::vector<std::string>
    output = run(*s,g,*rlc.transform(makePtrFromRef(obj)),outStream);
  rlc.project(x,*s);
  x.plus(*rlc.getFeasibleVector());
  return output;
}

template<typename Real>
void Algorithm<Real>::writeHeader( std::ostream& os ) const {
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::endl;
}

template<typename Real>
void Algorithm<Real>::writeName( std::ostream& os ) const {
  throw Exception::NotImplemented(">>> ROL::Algorithm::printName() is not implemented!");
}

template<typename Real>
void Algorithm<Real>::writeOutput( std::ostream& os, bool print_header ) const {
  os << std::scientific << std::setprecision(6);
  if ( print_header ) {
    os << printHeader();
  }
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::endl;
  }
  else {
    os << "  "; 
    os << std::setw(6)  << std::left << state_->iter;  
    os << std::setw(15) << std::left << state_->value; 
    os << std::setw(15) << std::left << state_->gnorm; 
    os << std::setw(15) << std::left << state_->snorm; 
    os << std::setw(10) << std::left << state_->nfval;              
    os << std::setw(10) << std::left << state_->ngrad;              
    os << std::endl;
  }
}

template<typename Real>
void Algorithm<Real>::writeExitStatus( std::ostream& os ) const {
  os << "Optimization Terminated with Status: ";
  os << EExitStatusToString(state_->statusFlag);
  os << std::endl;
}

template<typename Real>
Ptr<const AlgorithmState<Real>>& Algorithm<Real>::getState() const {
  return state_;
}

template<typename Real>
void Algorithm<Real>::reset() {
  state_->reset();
}

} // namespace TypeU
} // namespace ROL2
} // namespace ROL

namespace ROL2 {
namespace TypeU {
using ROL::ROL2::TypeU::Algorithm;
} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_ALGORITHM_DEF_HPP
