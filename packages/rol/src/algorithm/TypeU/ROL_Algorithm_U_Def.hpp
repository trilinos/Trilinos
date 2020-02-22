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

#ifndef ROL_ALGORITHM_U_DEF_H
#define ROL_ALGORITHM_U_DEF_H

#include "ROL_Types.hpp"
#include "ROL_ReduceLinearConstraint.hpp"
#include "ROL_ValidParameters.hpp"

namespace ROL {

template<typename Real>
Algorithm_U<Real>::Algorithm_U( void )
  : status_(makePtr<CombinedStatusTest<Real>>()),
    state_(makePtr<AlgorithmState_U<Real>>()) {
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>());
}

template<typename Real>
void Algorithm_U<Real>::initialize(const Vector<Real> &x, const Vector<Real> &g) {
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
void Algorithm_U<Real>::setStatusTest(const Ptr<StatusTest<Real>> &status,
                                      const bool combineStatus) {
  if (!combineStatus) { // Do not combine status tests
    status_->reset();
  }
  status_->add(status); // Add user-defined StatusTest
}

template<typename Real>
std::vector<std::string> Algorithm_U<Real>::run( Vector<Real>    &x,
                                                 Objective<Real> &obj,
                                                 std::ostream    &outStream ) {
  return run(x,x.dual(),obj,outStream);
}

template<typename Real>
std::vector<std::string> Algorithm_U<Real>::run( Vector<Real>     &x,
                                                 Objective<Real>  &obj,
                                                 Constraint<Real> &linear_con,
                                                 Vector<Real>     &linear_mul,
                                                 std::ostream     &outStream ) {
  return run(x,x.dual(),obj,linear_con,linear_mul,linear_mul.dual(),outStream);
}

template<typename Real>
std::vector<std::string> Algorithm_U<Real>::run( Vector<Real>       &x,
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
std::string Algorithm_U<Real>::printHeader( void ) const {
  std::stringstream hist;
  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::endl;
  return hist.str();
}

template<typename Real>
std::string Algorithm_U<Real>::printName( void ) const {
  throw Exception::NotImplemented(">>> ROL::Algorithm_U::printName() is not implemented!");
}

template<typename Real>
std::string Algorithm_U<Real>::print( const bool print_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( print_header ) {
    hist << printHeader();
  }
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::endl;
  }
  else {
    hist << "  "; 
    hist << std::setw(6)  << std::left << state_->iter;  
    hist << std::setw(15) << std::left << state_->value; 
    hist << std::setw(15) << std::left << state_->gnorm; 
    hist << std::setw(15) << std::left << state_->snorm; 
    hist << std::setw(10) << std::left << state_->nfval;              
    hist << std::setw(10) << std::left << state_->ngrad;              
    hist << std::endl;
  }
  return hist.str();
}

template<typename Real>
std::string Algorithm_U<Real>::printExitStatus( void ) const {
  std::stringstream hist;
  hist << "Optimization Terminated with Status: ";
  hist << EExitStatusToString(state_->statusFlag);
  hist << std::endl;
  return hist.str();
}

template<typename Real>
Ptr<const AlgorithmState<Real>> Algorithm_U<Real>::getState(void) const {
  return state_;
}

template<typename Real>
void Algorithm_U<Real>::reset(void) {
  state_->reset();
}

} // namespace ROL

#endif
