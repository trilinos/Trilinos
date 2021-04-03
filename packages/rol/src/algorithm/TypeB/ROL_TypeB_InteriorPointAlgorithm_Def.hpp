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

#ifndef ROL_TYPEB_PRIMALINTERIORPOINTALGORITHM_DEF_HPP
#define ROL_TYPEB_PRIMALINTERIORPOINTALGORITHM_DEF_HPP

#include "ROL_TypeU_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeB {

template<typename Real>
InteriorPointAlgorithm<Real>::InteriorPointAlgorithm(ParameterList &list)
  : TypeB::Algorithm<Real>::Algorithm(),
    list_(list), subproblemIter_(0), print_(false) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameters
  ParameterList& steplist = list.sublist("Step").sublist("Interior Point");
  state_->searchSize = steplist.get("Initial Barrier Parameter",           1.0);
  mumin_             = steplist.get("Minimum Barrier Parameter",           1e-4);
  mumax_             = steplist.get("Maximum Barrier Parameter",           1e8);
  rho_               = steplist.get("Barrier Penalty Reduction Factor",    0.5);
  useLinearDamping_  = steplist.get("Use Linear Damping",                  true);
  kappaD_            = steplist.get("Linear Damping Coefficient",          1.e-4);
  print_             = steplist.sublist("Subproblem").get("Print History", false);
  // Set parameters for step subproblem
  gtol_ = steplist.sublist("Subproblem").get("Initial Optimality Tolerance",  1e-2);
  stol_ = static_cast<Real>(1e-6)*gtol_;
  int maxit = steplist.sublist("Subproblem").get("Iteration Limit",       1000);
  list_.sublist("Status Test").set("Iteration Limit",      maxit);
  // Subproblem tolerance update parameters
  gtolrate_ = steplist.sublist("Subproblem").get("Optimality Tolerance Reduction Factor",     0.1);
  mingtol_  = static_cast<Real>(1e-2)*list.sublist("Status Test").get("Gradient Tolerance",   1e-8);
  // Get step name from parameterlist
  stepname_ = steplist.sublist("Subproblem").get("Step Type","Augmented Lagrangian");

  // Output settings
  verbosity_   = list.sublist("General").get("Output Level", 0);
  writeHeader_ = verbosity_ > 2;
  print_       = (verbosity_ > 2 ? true : print_);
  list_.sublist("General").set("Output Level",(print_ ? verbosity_ : 0));
}

template<typename Real>
void InteriorPointAlgorithm<Real>::initialize(Vector<Real>                 &x,
                                              const Vector<Real>           &g,
                                              InteriorPointObjective<Real> &ipobj,
                                              BoundConstraint<Real>        &bnd,
                                              Vector<Real>                 &pwa,
                                              std::ostream                 &outStream) {
  hasPolyProj_ = true;
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasPolyProj_ = false;
  }
  proj_->project(x,outStream);
  bnd.projectInterior(x);
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ngrad = 0;
  updateState(x,ipobj,bnd,pwa);
}


template<typename Real>
void InteriorPointAlgorithm<Real>::updateState(const Vector<Real>           &x,
                                               InteriorPointObjective<Real> &ipobj,
                                               BoundConstraint<Real>        &bnd,
                                               Vector<Real>                 &pwa,
                                               std::ostream                 &outStream) {
  const Real one(1);
  Real zerotol = std::sqrt(ROL_EPSILON<Real>());
  // Update objective and constraint
  if (state_-> iter == 0) {
    ipobj.update(x,UpdateType::Initial,state_->iter);
  }
  //else {
  //  ipobj.update(x,UpdateType::Accept,state_->iter);
  //}
  // Compute norm of the gradient of the Lagrangian
  state_->value = ipobj.getObjectiveValue(x, zerotol);
  //state_->gradientVec->set(*ipobj.getObjectiveGradient(x, zerotol));
  ipobj.gradient(*state_->gradientVec, x, zerotol);
  //state_->gnorm = state_->gradientVec->norm();
  pwa.set(x);
  pwa.axpy(-one,state_->gradientVec->dual());
  proj_->project(pwa,outStream);
  pwa.axpy(-one,x);
  state_->gnorm = pwa.norm();
  // Update state
  state_->nfval++;
  state_->ngrad++;
}

template<typename Real>
void InteriorPointAlgorithm<Real>::run( Vector<Real>          &x,
                                        const Vector<Real>    &g, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        std::ostream          &outStream ) {
  const Real one(1);
  Ptr<Vector<Real>> pwa = x.clone();
  // Initialize interior point data
  InteriorPointObjective<Real> ipobj(makePtrFromRef(obj),makePtrFromRef(bnd),
                                     x,g,useLinearDamping_,kappaD_,
                                     state_->searchSize);
  initialize(x,g,ipobj,bnd,*pwa,outStream);
  Ptr<TypeU::Algorithm<Real>> algo;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Solve interior point subproblem
    list_.sublist("Status Test").set("Gradient Tolerance",   gtol_);
    list_.sublist("Status Test").set("Step Tolerance",       stol_);
    algo = TypeU::AlgorithmFactory<Real>(list_);
    if (hasPolyProj_) algo->run(x,g,ipobj,
                                *proj_->getLinearConstraint(),
                                *proj_->getMultiplier(),
                                *proj_->getResidual(),outStream);
    else              algo->run(x,g,ipobj,outStream);
    subproblemIter_ = algo->getState()->iter;
    state_->nfval += algo->getState()->nfval;
    state_->ngrad += algo->getState()->ngrad;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    state_->iterateVec->set(x);

    // Update objective and constraint
    state_->iter++;

    // Update barrier parameter and subproblem tolerances
    if (algo->getState()->statusFlag == EXITSTATUS_CONVERGED) {
      if( (rho_< one && state_->searchSize > mumin_) || (rho_ > one && state_->searchSize < mumax_) ) {
        state_->searchSize *= rho_;
        ipobj.updatePenalty(state_->searchSize);
      }
      gtol_ *= gtolrate_; gtol_ = std::max(gtol_,mingtol_);
      stol_  = static_cast<Real>(1e-6)*gtol_;
    }

    // Update state
    updateState(x,ipobj,bnd,*pwa,outStream);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(109,'-') << std::endl;
    hist << "Interior Point Solver";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  fval     - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  penalty  - Penalty parameter for bound constraints" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  optTol   - Subproblem optimality tolerance" << std::endl;
    hist << "  subiter  - Number of subproblem iterations" << std::endl;
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "fval";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(10) << std::left << "penalty";
  hist << std::setw(8) << std::left << "#fval";
  hist << std::setw(8) << std::left << "#grad";
  hist << std::setw(10) << std::left << "optTol";
  hist << std::setw(8) << std::left << "subIter";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Interior Point Solver (Type B, Bound Constraints)";
  hist << std::endl;
  hist << "Subproblem Solver: " << stepname_ << std::endl;
  os << hist.str();
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( write_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << "---";
    hist << std::setw(8) << std::left << "---";
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << gtol_;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << subproblemIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeB
} // namespace ROL

#endif
