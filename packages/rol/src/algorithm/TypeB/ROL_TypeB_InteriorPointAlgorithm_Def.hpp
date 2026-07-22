// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_PRIMALINTERIORPOINTALGORITHM_DEF_HPP
#define ROL_TYPEB_PRIMALINTERIORPOINTALGORITHM_DEF_HPP

#include "ROL_TypeU_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeB {

template<typename Real>
InteriorPointAlgorithm<Real>::InteriorPointAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant)
  : TypeB::Algorithm<Real>::Algorithm(), secant_(secant),
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
    algo = TypeU::AlgorithmFactory<Real>(list_,secant_);
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
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(109,'-') << std::endl;
    os << "Interior Point Solver";
    os << " status output definitions" << std::endl << std::endl;
    os << "  iter     - Number of iterates (steps taken)" << std::endl;
    os << "  fval     - Objective function value" << std::endl;
    os << "  gnorm    - Norm of the gradient" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  penalty  - Penalty parameter for bound constraints" << std::endl;
    os << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    os << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    os << "  optTol   - Subproblem optimality tolerance" << std::endl;
    os << "  subiter  - Number of subproblem iterations" << std::endl;
    os << std::string(109,'-') << std::endl;
  }

  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "fval";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "penalty";
  os << std::setw(8) << std::left << "#fval";
  os << std::setw(8) << std::left << "#grad";
  os << std::setw(10) << std::left << "optTol";
  os << std::setw(8) << std::left << "subIter";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Interior Point Solver (Type B, Bound Constraints)";
  os << std::endl;
  os << "Subproblem Solver: " << stepname_ << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( write_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << "---";
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << "---";
    os << std::setw(8) << std::left << "---";
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << state_->snorm;
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << gtol_;
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << subproblemIter_;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeB
} // namespace ROL

#endif
