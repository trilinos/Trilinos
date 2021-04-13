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

#ifndef ROL_TYPEE_FLETCHERALGORITHM_DEF_H
#define ROL_TYPEE_FLETCHERALGORITHM_DEF_H

#include "ROL_TypeU_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeE {

template<typename Real>
FletcherAlgorithm<Real>::FletcherAlgorithm( ParameterList &list )
  : TypeE::Algorithm<Real>::Algorithm(), list_(list), subproblemIter_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>(list));

  ParameterList& sublist = list.sublist("Step").sublist("Fletcher");
  sigma_                 = sublist.get("Penalty Parameter",                        1.0);
  delta_                 = sublist.get("Regularization Parameter",                 0.0);
  minDelta_              = sublist.get("Minimum Regularization Parameter",         1e-8);
  deltaUpdate_           = sublist.get("Regularization Parameter Decrease Factor", 1e-1);
  sigmaUpdate_           = sublist.get("Penalty Parameter Growth Factor",          2.0);
  modifySigma_           = sublist.get("Modify Penalty Parameter",                 false);
  maxSigma_              = sublist.get("Maximum Penalty Parameter",                1e8);
  minSigma_              = sublist.get("Minimum Penalty Parameter",                1e-6);     
  subStep_               = sublist.get("Subproblem Step Type",                     "Trust Region");
  int subiter            = sublist.get("Subproblem Iteration Limit",               100);
  // Verbosity setting
  verbosity_             = list.sublist("General").get("Output Level", 0);
  printHeader_           = verbosity_ > 2;
  bool print             = verbosity_ >= 2;
  // Set parameter list for subproblem solve
  list_.sublist("General").set("Output Level",(print ? verbosity_-1 : 0));
  list_.sublist("Step").set("Type", subStep_);
  list_.sublist("Status Test").set("Iteration Limit", subiter);
}

template<typename Real>
void FletcherAlgorithm<Real>::initialize( Vector<Real>             &x,
                                          const Vector<Real>       &g,
                                          const Vector<Real>       &l,
                                          const Vector<Real>       &c,
                                          FletcherObjectiveE<Real> &fobj,
                                          Constraint<Real>         &con,
                                          std::ostream             &outStream ) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  TypeE::Algorithm<Real>::initialize(x,g,l,c);

  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ncval = 0;
  state_->ngrad = 0;

  // Compute objective value
  fobj.reset(sigma_,delta_);
  fobj.update(x,UpdateType::Initial,state_->iter);
  merit_ = fobj.value(x,tol);
  state_->value = fobj.getObjectiveValue(x);
  fobj.gradient(*state_->gradientVec,x,tol);
  gpnorm_ = state_->gradientVec->norm();
  state_->gradientVec->set(*fobj.getLagrangianGradient(x));
  state_->gnorm = state_->gradientVec->norm();

  // Compute constraint violation
  state_->constraintVec->set(*fobj.getConstraintVec(x));
  state_->cnorm = state_->constraintVec->norm();

  // Update evaluation counters
  state_->ncval += fobj.getNumberConstraintEvaluations();
  state_->nfval += fobj.getNumberFunctionEvaluations();
  state_->ngrad += fobj.getNumberGradientEvaluations();
}

template<typename Real>
void FletcherAlgorithm<Real>::run( Vector<Real>       &x,
                                   const Vector<Real> &g,
                                   Objective<Real>    &obj,
                                   Constraint<Real>   &econ,
                                   Vector<Real>       &emul,
                                   const Vector<Real> &eres,
                                   std::ostream       &outStream ) {
  // Initialize Fletcher penalty data
  const Real one(1);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  Ptr<Vector<Real>> dwa_ = g.clone();
  FletcherObjectiveE<Real> fobj(makePtrFromRef(obj),makePtrFromRef(econ),x,g,eres,emul,list_);
  initialize(x,g,emul,eres,fobj,econ,outStream);
  Ptr<TypeU::Algorithm<Real>> algo;

  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Minimize Fletcher penalty
    algo = TypeU::AlgorithmFactory<Real>(list_);
    algo->run(x,g,fobj,outStream);
    subproblemIter_ = algo->getState()->iter;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Update iteration information
    state_->iter++;
    state_->iterateVec->set(x);
    state_->value = fobj.getObjectiveValue(x);
    state_->constraintVec->set(*fobj.getConstraintVec(x));
    state_->cnorm = state_->constraintVec->norm();
    state_->gradientVec->set(*fobj.getLagrangianGradient(x));
    state_->gnorm = state_->gradientVec->norm();
    state_->lagmultVec->set(*fobj.getMultiplierVec(x));
    emul.set(*state_->lagmultVec);
    merit_  = algo->getState()->value;
    gpnorm_ = algo->getState()->gnorm;

    // Update evaluation counters
    state_->nfval += fobj.getNumberFunctionEvaluations();
    state_->ngrad += fobj.getNumberGradientEvaluations();
    state_->ncval += fobj.getNumberConstraintEvaluations();

    // Update penalty parameters
    bool too_infeasible = state_->cnorm > static_cast<Real>(100.)*gpnorm_;
    bool too_feasible   = state_->cnorm < static_cast<Real>(1e-2)*gpnorm_;
    bool modified = false;
    if( too_infeasible && !modified && modifySigma_
        && algo->getState()->statusFlag == EXITSTATUS_CONVERGED) {
      sigma_ = std::min(sigma_*sigmaUpdate_, maxSigma_);
      modified = true;
    }
    if( too_feasible && !modified && modifySigma_
        && algo->getState()->statusFlag == EXITSTATUS_CONVERGED) {
      sigma_ = std::max(sigma_/sigmaUpdate_, minSigma_);
      modified = true;      
    }
    if( delta_ > minDelta_ && !modified ) {
      Real deltaNext = delta_ * deltaUpdate_;
      if( gpnorm_ < deltaNext ) {
        delta_ = deltaNext;
        modified = true;
      }
    }
    if( modified ) {
      fobj.reset(sigma_,delta_);
      merit_ = fobj.value(x,tol);
      fobj.gradient(*dwa_,x,tol);
      gpnorm_ = dwa_->norm();

      state_->nfval++;
      state_->ngrad++;
      state_->ncval++;
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }

  if (verbosity_ > 0) TypeE::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void FletcherAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if(verbosity_>1) {
    hist << std::string(114,'-') << std::endl;
    hist << "Fletcher exact penalty status output definitions" << std::endl << std::endl;
    hist << "  iter    - Number of iterates (steps taken)"            << std::endl;
    hist << "  fval    - Objective function value"                    << std::endl;
    hist << "  cnorm   - Norm of the constraint violation"            << std::endl;
    hist << "  gLnorm  - Norm of the gradient of the Lagrangian"      << std::endl;
    hist << "  snorm   - Norm of the step"                            << std::endl;
    hist << "  merit   - Penalty function value"                      << std::endl;
    hist << "  gpnorm  - Norm of the gradient of the penalty"         << std::endl;
    hist << "  penalty - Penalty parameter"                           << std::endl;
    hist << "  delta   - Feasibility tolerance"                       << std::endl;
    hist << "  #fval   - Number of times the objective was computed"  << std::endl;
    hist << "  #grad   - Number of times the gradient was computed"   << std::endl;
    hist << "  #cval   - Number of times the constraint was computed" << std::endl;
    hist << "  subIter - Number of iterations to solve subproblem"    << std::endl;
    hist << std::string(114,'-') << std::endl;
  }
  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "fval";
  hist << std::setw(15) << std::left << "cnorm";
  hist << std::setw(15) << std::left << "gLnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "merit";
  hist << std::setw(15) << std::left << "gpnorm";
  hist << std::setw(10) << std::left << "penalty";
  hist << std::setw(10) << std::left << "delta";
  hist << std::setw(8)  << std::left << "#fval";
  hist << std::setw(8)  << std::left << "#grad";
  hist << std::setw(8)  << std::left << "#cval";
  hist << std::setw(8)  << std::left << "subIter";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void FletcherAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Fletcher Exact Penalty Solver (Type E, Equality Constraints)";
  hist << std::endl;
  hist << "Subproblem Solver: " << subStep_ << std::endl;
  os << hist.str();
}

template<typename Real>
void FletcherAlgorithm<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( print_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->cnorm;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(15) << std::left << merit_;
    hist << std::setw(15) << std::left << gpnorm_;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << sigma_;
    hist << std::setw(10) << std::left << delta_;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << state_->ncval;
    hist << std::setw(8) << std::left << "---";
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->cnorm;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::setw(15) << std::left << merit_;
    hist << std::setw(15) << std::left << gpnorm_;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << sigma_;
    hist << std::setw(10) << std::left << delta_;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << state_->ncval;
    hist << std::setw(8) << std::left << subproblemIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeE
} // namespace ROL

#endif
