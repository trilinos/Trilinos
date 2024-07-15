// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_PROXGRADIENTALGORITHM_DEF_HPP
#define ROL_TYPEP_PROXGRADIENTALGORITHM_DEF_HPP

namespace ROL {
namespace TypeP {

template<typename Real>
ProxGradientAlgorithm<Real>::ProxGradientAlgorithm(ParameterList &list) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("Line Search");
  maxit_        = lslist.get("Function Evaluation Limit",                        20);
  alpha0_       = lslist.get("Initial Step Size",                               1.0);
  normAlpha_    = lslist.get("Normalize Initial Step Size",                   false);
  alpha0bnd_    = lslist.get("Lower Bound for Initial Step Size",              1e-4);
  useralpha_    = lslist.get("User Defined Initial Step Size",                false);
  usePrevAlpha_ = lslist.get("Use Previous Step Length as Initial Guess",     false);
  c1_           = lslist.get("Sufficient Decrease Tolerance",                  1e-4);
  maxAlpha_     = lslist.get("Maximum Step Size",                           alpha0_);
  useAdapt_     = lslist.get("Use Adaptive Step Size Selection",               true);
  initProx_     = lslist.get("Apply Prox to Initial Guess",                   false);
  rhodec_       = lslist.sublist("Line-Search Method").get("Backtracking Rate", 0.5);
  rhoinc_       = lslist.sublist("Line-Search Method").get("Increase Rate"    , 2.0);
  t0_           = list.sublist("Status Test").get("Gradient Scale"            , 1.0);
  verbosity_    = list.sublist("General").get("Output Level",                     0);
  writeHeader_  = verbosity_ > 2;
}

template<typename Real>
void ProxGradientAlgorithm<Real>::initialize(Vector<Real>       &x,
                                             const Vector<Real> &g,
                                             Objective<Real>    &sobj,
                                             Objective<Real>    &nobj,
                                             Vector<Real>       &px,
                                             Vector<Real>       &dg,
                                             std::ostream       &outStream) {
  const Real one(1);
  // Initialize data
  TypeP::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  if (initProx_) {
    nobj.prox(*state_->iterateVec,x,state_->searchSize,ftol);
    state_->nprox++;
    x.set(*state_->iterateVec);
  }
  // Evaluate objective function
  sobj.update(x,UpdateType::Initial,state_->iter);
  nobj.update(x,UpdateType::Initial,state_->iter); 
  state_->svalue = sobj.value(x,ftol); state_->nsval++;
  state_->nvalue = nobj.value(x,ftol); state_->nnval++; 
  state_->value  = state_->svalue + state_->nvalue;
  // Evaluate gradient of smooth part
  sobj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  dg.set(state_->gradientVec->dual());
  // Compute initial step size as 2/L, where L = 2|f(x+s)-f(x)-f'(x)s|/||s||^2
  // is a lower estimate of the Lipschitz constant of f
  if (!useralpha_) {
    bool flag = maxAlpha_ == alpha0_;
    // Evaluate objective at Prox(x - t0 dg)
    pgstep(px, *state_->stepVec, nobj, x, dg, t0_, ftol);
    state_->snorm = state_->stepVec->norm();
    sobj.update(px,UpdateType::Trial);
    Real snew = sobj.value(px,ftol); 
    sobj.update(x,UpdateType::Revert);
    state_->nsval++;
    Real gs = state_->gradientVec->apply(*state_->stepVec);
    alpha0_ = (state_->snorm * state_->snorm) / std::abs(snew - state_->svalue - gs);
    alpha0_ = ((alpha0_ > alpha0bnd_) ? alpha0_ : one);
    if (flag) maxAlpha_ = alpha0_;
  }
  // Normalize initial CP step length
  if (normAlpha_)
    alpha0_ /= state_->gradientVec->norm();
  state_->searchSize = alpha0_;
  // Evaluate proximal gradient
  pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, dg, state_->searchSize, ftol);
  state_->snorm = state_->stepVec->norm();
  state_->gnorm = state_->snorm / state_->searchSize;
}

template<typename Real>
void ProxGradientAlgorithm<Real>::run( Vector<Real>       &x,
                                       const Vector<Real> &g, 
                                       Objective<Real>    &sobj,
                                       Objective<Real>    &nobj,
                                       std::ostream       &outStream ) {
  const Real one(1);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  // Initialize trust-region data
  Ptr<Vector<Real>> px = x.clone(), pxP = x.clone(), dg = x.clone();
  initialize(x,g,sobj,nobj,*px,*dg,outStream);
  Real strial(0), ntrial(0), Ftrial(0), Qk(0);
  Real strialP(0), ntrialP(0), FtrialP(0), alphaP(0);
  Real snorm(state_->snorm), searchSize(state_->searchSize);
  int ls_nfval = 0;
  bool incAlpha = false, accept = true;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  // Compute steepest descent step
  while (status_->check(*state_)) {
    accept = true;
    // Perform backtracking line search 
    state_->searchSize = searchSize;
    // Compute objective function values
    sobj.update(*state_->iterateVec,UpdateType::Trial);
    strial = sobj.value(*state_->iterateVec,tol);
    nobj.update(*state_->iterateVec,UpdateType::Trial);
    ntrial = nobj.value(*state_->iterateVec,tol); 
    Ftrial = strial + ntrial;
    ls_nfval = 1;
    // Compute decrease indicator
    Qk = state_->gradientVec->apply(*state_->stepVec) + ntrial - state_->nvalue;
    incAlpha = (Ftrial - state_->value <= c1_*Qk);
    if (verbosity_ > 1) {
      outStream << "  In TypeP::GradientAlgorithm: Line Search"  << std::endl;
      outStream << "    Step size:                        " << state_->searchSize   << std::endl;
      outStream << "    Trial smooth value:               " << strial               << std::endl;
      outStream << "    Trial nonsmooth value:            " << ntrial               << std::endl; 
      outStream << "    Computed reduction:               " << state_->value-Ftrial << std::endl;
      outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
      outStream << "    Increase alpha?:                  " << incAlpha             << std::endl;
    }
    if (incAlpha && useAdapt_) {
      ntrialP = ROL_INF<Real>();
      strialP = ROL_INF<Real>(); 
      FtrialP = ntrialP + strialP;
      while ( Ftrial - state_->value <= c1_*Qk
           && Ftrial <= FtrialP
           && state_->searchSize < maxAlpha_
           && ls_nfval < maxit_ ) {
        // Previous value was acceptable
        sobj.update(*state_->iterateVec,UpdateType::Accept);
        nobj.update(*state_->iterateVec,UpdateType::Accept); 
        // Backup previous values to avoid recomputation
        pxP->set(*state_->iterateVec);
        alphaP  = state_->searchSize;
        strialP = strial;
        ntrialP = ntrial; 
        FtrialP = Ftrial;
        // Increase search size
        state_->searchSize *= rhoinc_;
        state_->searchSize  = std::min(state_->searchSize,maxAlpha_);
        // Compute proximal gradient step with new search size
        pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, *dg, state_->searchSize, tol);
        // Compute objective function values
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        strial = sobj.value(*state_->iterateVec,tol);
        nobj.update(*state_->iterateVec,UpdateType::Trial);
        ntrial = nobj.value(*state_->iterateVec,tol); 
        Ftrial = strial + ntrial;
        ls_nfval++;
        // Compute decrease indicator
        Qk = state_->gradientVec->apply(*state_->stepVec) + ntrial - state_->nvalue;
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "    Step size:                        " << state_->searchSize   << std::endl;
          outStream << "    Trial smooth value:               " << strial               << std::endl;
          outStream << "    Trial nonsmooth value:            " << ntrial               << std::endl; 
          outStream << "    Computed reduction:               " << state_->value-Ftrial << std::endl;
          outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
          outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
          outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
        }
      }
      if (Ftrial - state_->value > c1_*Qk || Ftrial > FtrialP) {
        state_->iterateVec->set(*pxP);
        strial = strialP;
        ntrial = ntrialP; 
        Ftrial = FtrialP;
        state_->searchSize = alphaP;
        state_->stepVec->set(*state_->iterateVec);
        state_->stepVec->axpy(-one,x);
        accept = false;
      }
    }
    else {
      while ( Ftrial - state_->value > c1_*Qk && ls_nfval < maxit_ ) {
        // Decrease search size
        state_->searchSize *= rhodec_;
        // Compute proximal gradient step with new search size
        pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, *dg, state_->searchSize, tol);
        // Compute objective function values
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        strial = sobj.value(*state_->iterateVec,tol);
        nobj.update(*state_->iterateVec,UpdateType::Trial);
        ntrial = nobj.value(*state_->iterateVec,tol); 
        Ftrial = strial + ntrial;
        ls_nfval++;
        // Compute decrease indicator
        Qk = state_->gradientVec->apply(*state_->stepVec) + ntrial - state_->nvalue;
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "    Step size:                        " << state_->searchSize   << std::endl;
          outStream << "    Trial smooth value:               " << strial               << std::endl;
          outStream << "    Trial nonsmooth value:            " << ntrial               << std::endl; 
          outStream << "    Computed reduction:               " << state_->value-Ftrial << std::endl;
          outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
          outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
          outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
        }
      }
    }
    state_->nsval += ls_nfval;
    state_->nnval += ls_nfval;

    // Compute norm of step
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    state_->iter++;
    x.set(*state_->iterateVec);

    // Compute new value and gradient
    state_->svalue = strial;
    state_->nvalue = ntrial;
    state_->value  = Ftrial;
    if (accept) {
      sobj.update(x,UpdateType::Accept,state_->iter);
      nobj.update(x,UpdateType::Accept,state_->iter);
    }
    else {       
      sobj.update(x,UpdateType::Revert,state_->iter);
      nobj.update(x,UpdateType::Revert,state_->iter); 
    }
    sobj.gradient(*state_->gradientVec,x,tol);
    state_->ngrad++;
    dg->set(state_->gradientVec->dual());

    // Compute proximal gradient step with initial search size
    searchSize = state_->searchSize;
    if (!usePrevAlpha_ && !useAdapt_) searchSize = alpha0_;
    pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, *dg, searchSize, tol);
    snorm = state_->stepVec->norm();
    state_->gnorm = snorm / searchSize;

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeP::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(109,'-') << std::endl;
    os << "Proximal gradient descent";
    os << " status output definitions" << std::endl << std::endl;
    os << "  iter     - Number of iterates (steps taken)" << std::endl;
    os << "  value    - Objective function value" << std::endl;
    os << "  gnorm    - Norm of the proximal gradient with parameter alpha" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  alpha    - Line search step length" << std::endl;
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
  os << std::setw(10) << std::left << "#sval";
  os << std::setw(10) << std::left << "#nval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "#nprox"; 
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Proximal Gradient Descent with Bidirectional Line Search (Type P)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(15) << std::left << state_->searchSize;
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
