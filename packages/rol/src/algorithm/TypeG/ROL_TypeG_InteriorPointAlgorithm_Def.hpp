// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_INTERIORPOINTALGORITHM_DEF_H
#define ROL_TYPEG_INTERIORPOINTALGORITHM_DEF_H

#include "ROL_TypeE_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeG {

template<typename Real>
InteriorPointAlgorithm<Real>::InteriorPointAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant)
  : TypeG::Algorithm<Real>::Algorithm(), secant_(secant),
    list_(list), subproblemIter_(0), print_(false) {
  // Set status test
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>(list));

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
  ctol_ = steplist.sublist("Subproblem").get("Initial Feasibility Tolerance", 1e-2);
  stol_ = static_cast<Real>(1e-6)*std::min(gtol_,ctol_);
  int maxit = steplist.sublist("Subproblem").get("Iteration Limit",       1000);
  list_.sublist("Status Test").set("Iteration Limit",      maxit);
  // Subproblem tolerance update parameters
  gtolrate_ = steplist.sublist("Subproblem").get("Optimality Tolerance Reduction Factor",     0.1);
  ctolrate_ = steplist.sublist("Subproblem").get("Feasibility Tolerance Reduction Factor",    0.1);
  mingtol_  = static_cast<Real>(1e-2)*list.sublist("Status Test").get("Gradient Tolerance",   1e-8);
  minctol_  = static_cast<Real>(1e-2)*list.sublist("Status Test").get("Constraint Tolerance", 1e-8);
  // Get step name from parameterlist
  stepname_ = steplist.sublist("Subproblem").get("Step Type","Augmented Lagrangian");

  // Output settings
  verbosity_   = list.sublist("General").get("Output Level", 0);
  printHeader_ = verbosity_ > 2;
  print_       = (verbosity_ > 2 ? true : print_);
  list_.sublist("General").set("Output Level",(print_ ? verbosity_ : 0));
}

template<typename Real>
void InteriorPointAlgorithm<Real>::initialize(Vector<Real>                 &x,
                                              const Vector<Real>           &g,
                                              const Vector<Real>           &l,
                                              const Vector<Real>           &c,
                                              InteriorPointObjective<Real> &ipobj,
                                              BoundConstraint<Real>        &bnd,
                                              Constraint<Real>             &con,
                                              Vector<Real>                 &pwa,
                                              Vector<Real>                 &dwa,
                                              std::ostream                 &outStream) {
  hasPolyProj_ = true;
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasPolyProj_ = false;
  }
  proj_->project(x,outStream);
  bnd.projectInterior(x);
  // Initialize data
  TypeG::Algorithm<Real>::initialize(x,g,l,c);
  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ngrad = 0;
  state_->ncval = 0;
  updateState(x,l,ipobj,bnd,con,pwa,dwa,outStream);
}


template<typename Real>
void InteriorPointAlgorithm<Real>::updateState(const Vector<Real>           &x,
                                               const Vector<Real>           &l,
                                               InteriorPointObjective<Real> &ipobj,
                                               BoundConstraint<Real>        &bnd,
                                               Constraint<Real>             &con,
                                               Vector<Real>                 &pwa,
                                               Vector<Real>                 &dwa,
                                               std::ostream                 &outStream) {
  const Real one(1);
  Real zerotol = std::sqrt(ROL_EPSILON<Real>());
  // Update objective and constraint
  if (state_->iter == 0) {
    ipobj.update(x,UpdateType::Initial,state_->iter);
    con.update(x,UpdateType::Initial,state_->iter);
  }
  //else {
  //  ipobj.update(x,UpdateType::Accept,state_->iter);
  //  con.update(x,UpdateType::Accept,state_->iter);
  //}
  // Compute norm of the gradient of the Lagrangian
  state_->value = ipobj.getObjectiveValue(x, zerotol);
  //state_->gradientVec->set(*ipobj.getObjectiveGradient(x, zerotol));
  ipobj.gradient(*state_->gradientVec, x, zerotol);
  con.applyAdjointJacobian(dwa, l, x, zerotol);
  state_->gradientVec->plus(dwa);
  //state_->gnorm = state_->gradientVec->norm();
  pwa.set(x);
  pwa.axpy(-one,state_->gradientVec->dual());
  proj_->project(pwa,outStream);
  pwa.axpy(-one,x);
  state_->gnorm = pwa.norm();
  // Compute constraint violation
  con.value(*state_->constraintVec, x, zerotol);
  state_->cnorm = state_->constraintVec->norm();
  // Update state
  state_->nfval++;
  state_->ngrad++;
  state_->ncval++;
}

template<typename Real>
void InteriorPointAlgorithm<Real>::run( Vector<Real>          &x,
                                        const Vector<Real>    &g, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        Constraint<Real>      &econ,
                                        Vector<Real>          &emul,
                                        const Vector<Real>    &eres,
                                        std::ostream          &outStream ) {
  const Real one(1);
  Ptr<Vector<Real>> pwa = x.clone(), dwa = g.clone();
  // Initialize interior point data
  InteriorPointObjective<Real> ipobj(makePtrFromRef(obj),makePtrFromRef(bnd),
                                           x,g,useLinearDamping_,kappaD_,
                                           state_->searchSize);
  initialize(x,g,emul,eres,ipobj,bnd,econ,*pwa,*dwa,outStream);
  Ptr<TypeE::Algorithm<Real>> algo;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Solve interior point subproblem
    list_.sublist("Status Test").set("Gradient Tolerance",   gtol_);
    list_.sublist("Status Test").set("Constraint Tolerance", ctol_);
    list_.sublist("Status Test").set("Step Tolerance",       stol_);
    algo = TypeE::AlgorithmFactory<Real>(list_,secant_);
    if (hasPolyProj_) algo->run(x,g,ipobj,econ,emul,eres,
                                *proj_->getLinearConstraint(),
                                *proj_->getMultiplier(),
                                *proj_->getResidual(),outStream);
    else              algo->run(x,g,ipobj,econ,emul,eres,outStream);
    subproblemIter_ = algo->getState()->iter;
    state_->nfval += algo->getState()->nfval;
    state_->ngrad += algo->getState()->ngrad;
    state_->ncval += algo->getState()->ncval;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->lagmultVec->axpy(-one,emul);
    state_->snorm = state_->stepVec->norm();
    state_->snorm += state_->lagmultVec->norm();

    // Update iterate and Lagrange multiplier
    state_->iterateVec->set(x);
    state_->lagmultVec->set(emul);

    // Update objective and constraint
    state_->iter++;

    // Update barrier parameter and subproblem tolerances
    if (algo->getState()->statusFlag == EXITSTATUS_CONVERGED) {
      if( (rho_< one && state_->searchSize > mumin_) || (rho_ > one && state_->searchSize < mumax_) ) {
        state_->searchSize *= rho_;
        ipobj.updatePenalty(state_->searchSize);
      }
      gtol_ *= gtolrate_; gtol_ = std::max(gtol_,mingtol_);
      ctol_ *= ctolrate_; ctol_ = std::max(ctol_,minctol_);
      stol_  = static_cast<Real>(1e-6)*std::min(gtol_,ctol_);
    }

    // Update state
    updateState(x,emul,ipobj,bnd,econ,*pwa,*dwa);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) TypeG::Algorithm<Real>::writeExitStatus(outStream);
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
    os << "  cnorm    - Norm of the constraint" << std::endl;
    os << "  gLnorm   - Norm of the gradient of the Lagrangian" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  penalty  - Penalty parameter for bound constraints" << std::endl;
    os << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    os << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    os << "  #cval    - Cumulative number of times the constraint was evaluated" << std::endl;
    os << "  optTol   - Subproblem optimality tolerance" << std::endl;
    os << "  feasTol  - Subproblem feasibility tolerance" << std::endl;
    os << "  subiter  - Number of subproblem iterations" << std::endl;
    os << std::string(109,'-') << std::endl;
  }

  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "fval";
  os << std::setw(15) << std::left << "cnorm";
  os << std::setw(15) << std::left << "gLnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "penalty";
  os << std::setw(8) << std::left << "#fval";
  os << std::setw(8) << std::left << "#grad";
  os << std::setw(8) << std::left << "#cval";
  os << std::setw(10) << std::left << "optTol";
  os << std::setw(10) << std::left << "feasTol";
  os << std::setw(8) << std::left << "subIter";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Interior Point Solver (Type G, General Constraints)";
  os << std::endl;
  os << "Subproblem Solver: " << stepname_ << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void InteriorPointAlgorithm<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( print_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << "---";
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::setw(8) << std::left << state_->ncval;
    os << std::setw(10) << std::left << "---";
    os << std::setw(10) << std::left << "---";
    os << std::setw(8) << std::left << "---";
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->cnorm;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << state_->snorm;
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << state_->searchSize;
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << state_->nfval;
    os << std::setw(8) << std::left << state_->ngrad;
    os << std::setw(8) << std::left << state_->ncval;
    os << std::scientific << std::setprecision(2);
    os << std::setw(10) << std::left << gtol_;
    os << std::setw(10) << std::left << ctol_;
    os << std::scientific << std::setprecision(6);
    os << std::setw(8) << std::left << subproblemIter_;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeG
} // namespace ROL

#endif
