// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPE_ALGORITHM_DEF_H
#define ROL_TYPE_ALGORITHM_DEF_H

#include "ROL_Types.hpp"
#include "ROL_ReduceLinearConstraint.hpp"
#include "ROL_ValidParameters.hpp"
#include "ROL_ConstraintStatusTest.hpp"

namespace ROL {
namespace TypeE {

template<typename Real>
Algorithm<Real>::Algorithm()
  : status_(makePtr<CombinedStatusTest<Real>>()),
    state_(makePtr<AlgorithmState<Real>>()) {
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>());
}

template<typename Real>
void Algorithm<Real>::initialize(const Vector<Real> &x,
                                 const Vector<Real> &g,
                                 const Vector<Real> &mul,
                                 const Vector<Real> &c) {
  if (state_->iterateVec == nullPtr) {
    state_->iterateVec = x.clone();
  }
  state_->iterateVec->set(x);
  if (state_->lagmultVec == nullPtr) {
    state_->lagmultVec = mul.clone();
  }
  state_->lagmultVec->set(mul);
  if (state_->stepVec == nullPtr) {
    state_->stepVec = x.clone();
  }
  state_->stepVec->zero();
  if (state_->gradientVec == nullPtr) {
    state_->gradientVec = g.clone();
  }
  state_->gradientVec->set(g);
  if (state_->constraintVec == nullPtr) {
    state_->constraintVec = c.clone();
  }
  state_->constraintVec->zero();
  if (state_->minIterVec == nullPtr) {
    state_->minIterVec = x.clone();
  }
  state_->minIterVec->set(x);
  state_->minIter = state_->iter;
  state_->minValue = state_->value;
}

template<typename Real>
void Algorithm<Real>::setStatusTest(const Ptr<StatusTest<Real>> &status,
                                      const bool combineStatus) {
  if (!combineStatus) { // Do not combine status tests
    status_->reset();
  }
  status_->add(status); // Add user-defined StatusTest
}

template<typename Real>
void Algorithm<Real>::run( Problem<Real> &problem,
                           std::ostream     &outStream ) {
  if (problem.getProblemType() == TYPE_E) {
    run(*problem.getPrimalOptimizationVector(),
        *problem.getDualOptimizationVector(),
        *problem.getObjective(),
        *problem.getConstraint(),
        *problem.getMultiplierVector(),
        *problem.getResidualVector(),
        outStream);
        problem.finalizeIteration();
  }
  else {
    throw Exception::NotImplemented(">>> ROL::Algorithm::run : Optimization problem is not Type E!");
  }
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>     &x,
                           Objective<Real>  &obj,
                           Constraint<Real> &econ,
                           Vector<Real>     &emul,
                           std::ostream     &outStream ) {
  Problem<Real> problem(makePtrFromRef(obj), makePtrFromRef(x));
  problem.addConstraint("NEC",makePtrFromRef(econ),makePtrFromRef(emul));
  problem.finalize(false,false,outStream);
  run(problem,outStream);
  //run(x,x.dual(),obj,econ,emul,emul.dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>     &x,
                           Objective<Real>  &obj,
                           Constraint<Real> &econ,
                           Vector<Real>     &emul,
                           Constraint<Real> &linear_econ,
                           Vector<Real>     &linear_emul,
                           std::ostream     &outStream ) {
  Problem<Real> problem(makePtrFromRef(obj), makePtrFromRef(x));
  problem.addConstraint("NEC",makePtrFromRef(econ),makePtrFromRef(emul));
  problem.addLinearConstraint("LEC",makePtrFromRef(linear_econ),makePtrFromRef(linear_emul));
  problem.finalize(false,false,outStream);
  run(problem,outStream);
  //run(x,x.dual(),obj,econ,emul,emul.dual(),linear_econ,linear_emul,linear_emul.dual(),outStream);
}

template<typename Real>
void Algorithm<Real>::run( Vector<Real>       &x,
                           const Vector<Real> &g,
                           Objective<Real>    &obj,
                           Constraint<Real>   &econ,
                           Vector<Real>       &emul,
                           const Vector<Real> &eres,
                           Constraint<Real>   &linear_econ,
                           Vector<Real>       &linear_emul,
                           const Vector<Real> &linear_eres,
                           std::ostream       &outStream ) {
  Ptr<Vector<Real>> gp = g.clone(), erp = eres.clone(), lerp = linear_eres.clone();
  Problem<Real> problem(makePtrFromRef(obj), makePtrFromRef(x), gp);
  problem.addConstraint("NEC",makePtrFromRef(econ),makePtrFromRef(emul),erp,false);
  problem.addLinearConstraint("LEC",makePtrFromRef(linear_econ),makePtrFromRef(linear_emul),lerp,false);
  problem.finalize(false,false,outStream);
  run(problem,outStream);
  //Ptr<Vector<Real>> xfeas = x.clone(); xfeas->set(x);
  //ReduceLinearConstraint<Real> rlc(makePtrFromRef(linear_econ),xfeas,makePtrFromRef(linear_eres));
  //Ptr<Vector<Real>> s = x.clone(); s->zero();
  //void output = run(*s,g,*rlc.transform(makePtrFromRef(obj)),
  //                                      *rlc.transform(makePtrFromRef(econ)),emul,eres,outStream);
  //rlc.project(x,*s);
  //x.plus(*rlc.getFeasibleVector());
  //return output;
}

template<typename Real>
void Algorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "cnorm";
  os << std::setw(15) << std::left << "gLnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void Algorithm<Real>::writeName( std::ostream& os ) const {
  throw Exception::NotImplemented(">>> ROL::TypeE::Algorithm::writeName() is not implemented!");
}

template<typename Real>
void Algorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( write_header ) writeHeader(os);
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::endl;
  }
  else {
    os << "  "; 
    os << std::setw(6)  << std::left << state_->iter;  
    os << std::setw(15) << std::left << state_->value; 
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << state_->gnorm; 
    os << std::setw(15) << std::left << state_->snorm; 
    os << std::setw(10) << std::left << state_->nfval;              
    os << std::setw(10) << std::left << state_->ngrad;              
    os << std::endl;
  }
  os.flags(osFlags);
}

template<typename Real>
void Algorithm<Real>::writeExitStatus( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << "Optimization Terminated with Status: ";
  os << EExitStatusToString(state_->statusFlag);
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
Ptr<const AlgorithmState<Real>> Algorithm<Real>::getState() const {
//Ptr<const AlgorithmState<Real>>& Algorithm<Real>::getState() const {
  return state_;
}

template<typename Real>
void Algorithm<Real>::reset() {
  state_->reset();
}

} // namespace TypeE

} // namespace ROL

#endif
