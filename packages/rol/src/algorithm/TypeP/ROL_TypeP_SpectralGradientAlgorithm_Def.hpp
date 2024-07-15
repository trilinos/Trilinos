// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_SPECTRALGRADIENTALGORITHM_DEF_HPP
#define ROL_TYPEP_SPECTRALGRADIENTALGORITHM_DEF_HPP

namespace ROL {
namespace TypeP {

template<typename Real>
SpectralGradientAlgorithm<Real>::SpectralGradientAlgorithm(ParameterList &list) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("Spectral Gradient");
  maxit_        = lslist.get("Function Evaluation Limit",                        20);
  lambda_       = lslist.get("Initial Spectral Step Size",                     -1.0);
  lambdaMin_    = lslist.get("Minimum Spectral Step Size",                     1e-8); 
  lambdaMax_    = lslist.get("Maximum Spectral Step Size",                      1e8); 
  sigma1_       = lslist.get("Lower Step Size Safeguard",                       0.1);
  sigma2_       = lslist.get("Upper Step Size Safeguard",                       0.9);
  rhodec_       = lslist.get("Backtracking Rate",                              1e-1);
  gamma_        = lslist.get("Sufficient Decrease Tolerance",                  1e-4);
  maxSize_      = lslist.get("Maximum Storage Size",                             10);
  initProx_     = lslist.get("Apply Prox to Initial Guess",                   false);
  t0_           = list.sublist("Status Test").get("Gradient Scale"            , 1.0);
  verbosity_    = list.sublist("General").get("Output Level",                     0);
  writeHeader_  = verbosity_ > 2;
}

template<typename Real>
void SpectralGradientAlgorithm<Real>::initialize(Vector<Real>       &x,
                                                 const Vector<Real> &g,
                                                 Objective<Real>    &sobj,
                                                 Objective<Real>    &nobj,
                                                 Vector<Real>       &px,
                                                 Vector<Real>       &dg,
                                                 std::ostream       &outStream) {
  const Real zero(0);
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
  if (lambda_ <= zero && state_->gnorm != zero)
    lambda_ = std::max(lambdaMin_,std::min(t0_,lambdaMax_));
  pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, dg, lambda_, ftol);
  state_->snorm = state_->stepVec->norm();
  state_->gnorm = state_->snorm / lambda_;
}

template<typename Real>
void SpectralGradientAlgorithm<Real>::run( Vector<Real>       &x,
                                           const Vector<Real> &g, 
                                           Objective<Real>    &sobj,
                                           Objective<Real>    &nobj,
                                           std::ostream       &outStream ) {
  const Real half(0.5), one(1), eps(std::sqrt(ROL_EPSILON<Real>()));
  // Initialize trust-region data
  Ptr<Vector<Real>> s = x.clone(), px = x.clone(), dg = x.clone(), y = g.clone(), xmin = x.clone();
  initialize(x,g,sobj,nobj,*s,*dg,outStream);
  Real strial(0), ntrial(0), Ftrial(0), Fmin(0), Fmax(0), Qk(0), alpha(1), rhoTmp(1);
  Real gs(0), ys(0), snorm(state_->snorm), ss(0), tol(std::sqrt(ROL_EPSILON<Real>()));
  int ls_nfval = 0;
  std::deque<Real> Fqueue; Fqueue.push_back(state_->value);

  Fmin = state_->value;
  xmin->set(x);

  // Output
  if (verbosity_ > 0) writeOutput(outStream, true);

  // Iterate spectral projected gradient
  while (status_->check(*state_)) {
    // Nonmonotone Linesearch
    ls_nfval = 0;
    sobj.update(*state_->iterateVec,UpdateType::Trial);
    strial = sobj.value(*state_->iterateVec,tol);
    nobj.update(*state_->iterateVec,UpdateType::Trial);
    ntrial = nobj.value(*state_->iterateVec,tol);
    Ftrial = strial + ntrial;
    ls_nfval++;
    alpha  = one;
    Fmax   = *std::max_element(Fqueue.begin(),Fqueue.end());
    gs     = state_->gradientVec->apply(*state_->stepVec);
    Qk     = gs + ntrial - state_->nvalue;
    if (verbosity_ > 1) {
      outStream << "  In TypeP::SpectralGradientAlgorithm Line Search"  << std::endl;
      outStream << "    Step size:                        " << alpha                << std::endl;
      outStream << "    Trial objective value:            " << Ftrial               << std::endl;
      outStream << "    Max stored objective value:       " << Fmax                 << std::endl;
      outStream << "    Computed reduction:               " << Fmax-Ftrial          << std::endl;
      outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -Qk*gamma_           << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
    }
    while (Ftrial > Fmax + gamma_*Qk && ls_nfval < maxit_) {
      // Compute reduction factor by minimizing 1D quadratic model
      rhoTmp = std::min(one,-half*Qk/(strial-state_->svalue-alpha*gs));
      // Safeguard step size selection with back tracking
      alpha  = ((sigma1_ <= rhoTmp && rhoTmp <= sigma2_) ? rhoTmp : rhodec_)*alpha;
      // Update iterate vector
      state_->iterateVec->set(x);
      state_->iterateVec->axpy(alpha,*state_->stepVec);
      // Recompute objective function values
      sobj.update(*state_->iterateVec,UpdateType::Trial);
      strial = sobj.value(*state_->iterateVec,tol);
      nobj.update(*state_->iterateVec,UpdateType::Trial);
      ntrial = nobj.value(*state_->iterateVec,tol);
      Ftrial = strial + ntrial;
      ls_nfval++;
      Qk     = alpha * gs + ntrial - state_->nvalue;
      if (verbosity_ > 1) {
        outStream << "  In TypeP::SpectralGradientAlgorithm: Line Search"  << std::endl;
        outStream << "    Step size:                        " << alpha                << std::endl;
        outStream << "    Trial objective value:            " << Ftrial               << std::endl;
        outStream << "    Max stored objective value:       " << Fmax                 << std::endl;
        outStream << "    Computed reduction:               " << Fmax-Ftrial          << std::endl;
        outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
        outStream << "    Sufficient decrease bound:        " << -Qk*gamma_           << std::endl;
        outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
      }
    }
    state_->nsval += ls_nfval;
    state_->nnval += ls_nfval;
    if (static_cast<int>(Fqueue.size()) == maxSize_) Fqueue.pop_front();
    Fqueue.push_back(Ftrial);

    // Update state
    state_->iter++;
    state_->value      = Ftrial;
    state_->svalue     = strial;
    state_->nvalue     = ntrial;
    state_->searchSize = alpha;
    state_->snorm      = alpha * snorm;
    state_->stepVec->scale(alpha);
    x.set(*state_->iterateVec);
    sobj.update(x,UpdateType::Accept,state_->iter);
    nobj.update(x,UpdateType::Accept,state_->iter);

    // Store the best iterate
    if (state_->value <= Fmin) {
      Fmin = state_->value;
      xmin->set(x);
    }

    // Compute spectral step length
    y->set(*state_->gradientVec);
    y->scale(-one);
    sobj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    dg->set(state_->gradientVec->dual());
    y->plus(*state_->gradientVec);
    ys            = y->apply(*state_->stepVec);
    ss            = state_->snorm * state_->snorm;
    lambda_       = (ys<=eps*state_->snorm ? lambdaMax_ : std::max(lambdaMin_,std::min(ss/ys,lambdaMax_)));

    // Compute spectral proximal gradient step
    pgstep(*state_->iterateVec, *state_->stepVec, nobj, x, *dg, lambda_, tol);
    snorm         = state_->stepVec->norm();
    state_->gnorm = snorm / lambda_; 

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  x.set(*xmin);
  state_->value = Fmin;
  if (verbosity_ > 0) TypeP::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void SpectralGradientAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(109,'-') << std::endl;
    os << "Spectral proximal gradient with nonmonotone line search";
    os << " status output definitions" << std::endl << std::endl;
    os << "  iter     - Number of iterates (steps taken)" << std::endl;
    os << "  value    - Objective function value" << std::endl;
    os << "  gnorm    - Norm of the proximal gradient with parameter lambda" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  alpha    - Line search step length" << std::endl;
    os << "  lambda   - Spectral step length" << std::endl;
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
  os << std::setw(15) << std::left << "lambda";
  os << std::setw(10) << std::left << "#sval";
  os << std::setw(10) << std::left << "#nval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "#nprox"; 
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void SpectralGradientAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Spectral Proximal Gradient with Nonmonotone Line Search (Type P)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void SpectralGradientAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(15) << std::left << lambda_;
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
    os << std::setw(15) << std::left << lambda_;
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
