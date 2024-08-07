// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_IPIANOALGORITHM_DEF_HPP
#define ROL_TYPEP_IPIANOALGORITHM_DEF_HPP

namespace ROL {
namespace TypeP {

template<typename Real>
iPianoAlgorithm<Real>::iPianoAlgorithm(ParameterList &list) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("iPiano");
  t0_           = list.sublist("Status Test").get("Gradient Scale",             1.0);
  maxit_        = lslist.get("Reduction Iteration Limit",                        20);
  useConstBeta_ = lslist.get("Use Constant Beta",                             false);
  beta_         = lslist.get("Momentum Parameter",                             0.25);
  rhodec_       = lslist.get("Backtracking Rate",                               0.5);
  rhoinc_       = lslist.get("Increase Rate",                                   2.0);
  c1_           = lslist.get("Upper Interpolation Factor",                     1e-5);
  c2_           = lslist.get("Lower Interpolation Factor",                     1e-6);
  L_            = lslist.get("Initial Lipschitz Constant Estimate",         0.5/t0_);
  initProx_     = lslist.get("Apply Prox to Initial Guess",                   false);
  verbosity_    = list.sublist("General").get("Output Level",                     0);
  writeHeader_  = verbosity_ > 2;
}

template<typename Real>
void iPianoAlgorithm<Real>::initialize(Vector<Real>       &x,
                                                 const Vector<Real> &g,
                                                 Objective<Real>    &sobj,
                                                 Objective<Real>    &nobj,
                                                 Vector<Real>       &px,
                                                 Vector<Real>       &dg,
                                                 std::ostream       &outStream) {
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  // Initialize data
  TypeP::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  if (initProx_) {
    nobj.prox(*state_->iterateVec,x,t0_,ftol); state_->nprox++;
    x.set(*state_->iterateVec);
  }
  sobj.update(x,UpdateType::Initial,state_->iter);
  state_->svalue = sobj.value(x,ftol); state_->nsval++;
  nobj.update(x,UpdateType::Initial,state_->iter);
  state_->nvalue = nobj.value(x,ftol); state_->nnval++;
  state_->value  = state_->svalue + state_->nvalue;
  sobj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  dg.set(state_->gradientVec->dual());
  pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, dg, t0_, ftol);
  state_->snorm = state_->stepVec->norm();
  state_->gnorm = state_->snorm / t0_;
}

template<typename Real>
void iPianoAlgorithm<Real>::run( Vector<Real>       &x,
                                           const Vector<Real> &g, 
                                           Objective<Real>    &sobj,
                                           Objective<Real>    &nobj,
                                           std::ostream       &outStream ) {
  const Real half(0.5), one(1), two(2);
  // Initialize trust-region data
  Ptr<Vector<Real>> sP = x.clone(), xP = x.clone(), dg = x.clone(), xold = x.clone();
  initialize(x,g,sobj,nobj,*sP,*dg,outStream);
  Real strial(0), strialP(0), snormP(0), LP(0), alphaP(0), betaP(0), gs(0), b(0);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  bool accept(true);

  xold->set(x);

  // Output
  if (verbosity_ > 0) writeOutput(outStream, true);

  // Iterate spectral projected gradient
  while (status_->check(*state_)) {
    // Compute parameters alpha and beta
    if (!useConstBeta_) {
      b     = (c1_ + half * L_) / (c2_ + half * L_);
      beta_ = (b - one) / (b - half);
    }
    alpha_  = two * (1 - beta_) / (two * c2_ + L_);
    // Compute inertial step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-alpha_, *dg);
    state_->stepVec->axpy(beta_, x);
    state_->stepVec->axpy(-beta_, *xold);
    nobj.prox(*state_->iterateVec, *state_->stepVec, alpha_, tol); state_->nprox++;
    state_->stepVec->set(*state_->iterateVec);
    state_->stepVec->axpy(-one,x);
    state_->snorm = state_->stepVec->norm();
    // Compute smooth objective value
    sobj.update(*state_->iterateVec,UpdateType::Trial);
    nobj.update(*state_->iterateVec,UpdateType::Trial);
    strial = sobj.value(*state_->iterateVec,tol); state_->nsval++;
    gs     = state_->gradientVec->apply(*state_->stepVec);
    // Estimate Lipschitz constant of sobj
    if (strial <= state_->svalue + gs + half * L_ * state_->snorm * state_->snorm) {
      accept = true;
      for (int i = 0; i < maxit_; ++i) {
        // Store previously computed information
        sobj.update(*state_->iterateVec,UpdateType::Accept);
        nobj.update(*state_->iterateVec,UpdateType::Accept);
        LP      = L_;
        alphaP  = alpha_;
        betaP   = beta_;
        strialP = strial;
        snormP  = state_->snorm;
        xP->set(*state_->iterateVec);
        sP->set(*state_->stepVec);
        // Update alpha and beta with new Lipschitz constant estimate
        L_ /= rhoinc_;
        if (!useConstBeta_) {
          b     = (c1_ + half * L_) / (c2_ + half * L_);
          beta_ = (b - one) / (b - half);
        }
        alpha_  = two * (one - beta_) / (two * c2_ + L_);
        // Compute updated inertial step
        state_->stepVec->set(x);
        state_->stepVec->axpy(-alpha_, *dg);
        state_->stepVec->axpy(beta_, x);
        state_->stepVec->axpy(-beta_, *xold);
        nobj.prox(*state_->iterateVec, *state_->stepVec, alpha_, tol); state_->nprox++;
        state_->stepVec->set(*state_->iterateVec);
        state_->stepVec->axpy(-one,x);
        state_->snorm = state_->stepVec->norm();
        // Compute smooth objective value
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        strial = sobj.value(*state_->iterateVec,tol); state_->nsval++;
        gs     = state_->gradientVec->apply(*state_->stepVec);
        if (strial > state_->svalue + gs + half * L_ * state_->snorm * state_->snorm) {
          accept = false;
          L_     = LP;
          alpha_ = alphaP;
          beta_  = betaP;
          strial = strialP;
          state_->snorm = snormP;
          state_->iterateVec->set(*xP);
          state_->stepVec->set(*sP);
          break;
        }
      }
      if (accept) {
        sobj.update(*state_->iterateVec,UpdateType::Accept);
        nobj.update(*state_->iterateVec,UpdateType::Accept);
      }
      else {
        sobj.update(*state_->iterateVec,UpdateType::Revert);
        nobj.update(*state_->iterateVec,UpdateType::Revert);
      }
    }
    else {
      while (strial > state_->svalue + gs + half * L_ * state_->snorm * state_->snorm) {
        // Update alpha and beta with new Lipschitz constant estimate
        L_ /= rhodec_;
        if (!useConstBeta_) {
          b     = (c1_ + half * L_) / (c2_ + half * L_);
          beta_ = (b - one) / (b - half);
        }
        alpha_  = two * (one - beta_) / (two * c2_ + L_);
        // Compute updated inertial step
        state_->stepVec->set(x);
        state_->stepVec->axpy(-alpha_, *dg);
        state_->stepVec->axpy(beta_, x);
        state_->stepVec->axpy(-beta_, *xold);
        nobj.prox(*state_->iterateVec, *state_->stepVec, alpha_, tol); state_->nprox++;
        state_->stepVec->set(*state_->iterateVec);
        state_->stepVec->axpy(-one,x);
        state_->snorm = state_->stepVec->norm();
        // Compute smooth objective value
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        strial = sobj.value(*state_->iterateVec,tol); state_->nsval++;
        gs     = state_->gradientVec->apply(*state_->stepVec);
      }
      sobj.update(*state_->iterateVec,UpdateType::Accept);
      nobj.update(*state_->iterateVec,UpdateType::Accept);
    }
    // Update iteration
    state_->iter++;
    xold->set(x);
    x.set(*state_->iterateVec);
    state_->svalue = strial;
    state_->nvalue = nobj.value(x,tol); state_->nnval++;
    state_->value  = state_->svalue + state_->nvalue;
    sobj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    dg->set(state_->gradientVec->dual());
    // Compute proximal gradient for status check
    pgstep(*xP,*sP,nobj,x,*dg,t0_,tol);
    state_->gnorm = sP->norm() / t0_;

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeP::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void iPianoAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(109,'-') << std::endl;
    os << "iPiano: Inertial proximal algorithm for nonconvex optimization";
    os << " status output definitions" << std::endl << std::endl;
    os << "  iter     - Number of iterates (steps taken)" << std::endl;
    os << "  value    - Objective function value" << std::endl;
    os << "  gnorm    - Norm of the proximal gradient with parameter lambda" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  alpha    - Inertial gradient parameter" << std::endl;
    os << "  beta     - Inertial step parameter" << std::endl;
    os << "  L        - Lipschitz constant estimate" << std::endl;
    os << "  #sval    - Cumulative number of times the smooth objective function was evaluated" << std::endl;
    os << "  #nval    - Cumulative number of times the nonsmooth objective function was evaluated" << std::endl;
    os << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    os << "  #prox    - Cumulative number of times the proximal operator was computed" << std::endl; 
    os << std::string(109,'-') << std::endl;
  }

  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(15) << std::left << "alpha";
  os << std::setw(15) << std::left << "beta";
  os << std::setw(15) << std::left << "L";
  os << std::setw(10) << std::left << "#sval";
  os << std::setw(10) << std::left << "#nval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "#nprox"; 
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void iPianoAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "iPiano: Inertial Proximal Algorithm for Nonconvex Optimization (Type P)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void iPianoAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << L_;
    os << std::setw(10) << std::left << state_->nsval;
    os << std::setw(10) << std::left << state_->nnval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << state_->nprox; 
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << state_->snorm;
    os << std::setw(15) << std::left << alpha_;
    os << std::setw(15) << std::left << beta_;
    os << std::setw(15) << std::left << L_;
    os << std::setw(10) << std::left << state_->nsval;
    os << std::setw(10) << std::left << state_->nnval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << state_->nprox; 
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeP
} // namespace ROL

#endif
