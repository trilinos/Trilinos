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

#ifndef ROL_SPECTRALGRADIENTALGORITHM_B_DEF_H
#define ROL_SPECTRALGRADIENTALGORITHM_B_DEF_H

#include <deque>

namespace ROL {

template<typename Real>
SpectralGradientAlgorithm_B<Real>::SpectralGradientAlgorithm_B(ParameterList &list) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("Spectral Gradient");
  maxit_        = lslist.get("Function Evaluation Limit",                        20);
  lambda_       = lslist.get("Initial Spectral Step Size",                     -1.0);
  lambdaMin_    = lslist.get("Minimum Spectral Step Size",                     1e-3); 
  lambdaMax_    = lslist.get("Maximum Spectral Step Size",                      1e3); 
  sigma1_       = lslist.get("Lower Step Size Safeguard",                       0.1);
  sigma2_       = lslist.get("Upper Step Size Safeguard",                       0.9);
  rhodec_       = lslist.get("Backtracking Rate",                               0.5);
  gamma_        = lslist.get("Sufficient Decrease Tolerance",                  1e-4);
  maxSize_      = lslist.get("Maximum Storage Size",                             10);
  verbosity_    = list.sublist("General").get("Output Level",                     0);
  printHeader_  = verbosity_ > 2;
}

template<typename Real>
void SpectralGradientAlgorithm_B<Real>::initialize(Vector<Real>          &x,
                                                   const Vector<Real>    &g,
                                                   Objective<Real>       &obj,
                                                   BoundConstraint<Real> &bnd,
                                                   std::ostream &outStream) {
  const Real zero(0), one(1);
  if (proj_ == nullPtr)
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
  // Initialize data
  Algorithm_B<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  proj_->project(x,outStream);
  obj.update(x,UPDATE_INITIAL,state_->iter);
  state_->value = obj.value(x,ftol); state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream);
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  if (lambda_ <= zero && state_->gnorm != zero)
    lambda_ = std::max(lambdaMin_,std::min(one/state_->gnorm,lambdaMax_));
}

template<typename Real>
std::vector<std::string> SpectralGradientAlgorithm_B<Real>::run( Vector<Real>          &x,
                                                                 const Vector<Real>    &g, 
                                                                 Objective<Real>       &obj,
                                                                 BoundConstraint<Real> &bnd,
                                                                 std::ostream          &outStream ) {
  const Real half(0.5), one(1), eps(std::sqrt(ROL_EPSILON<Real>()));
  // Initialize trust-region data
  std::vector<std::string> output;
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone(), y = g.clone();
  Real ftrial(0), fmax(0), gs(0), alpha(1), alphaTmp(1);
  Real ys(0), ss(0), tol(std::sqrt(ROL_EPSILON<Real>()));
  int ls_nfval = 0;
  std::deque<Real> fqueue; fqueue.push_back(state_->value);

  // Output
  output.push_back(print(true));
  if (verbosity_ > 0) outStream << print(true);

  // Iterate spectral projected gradient
  state_->stepVec->set(state_->gradientVec->dual());
  while (status_->check(*state_)) {
    // Compute projected spectral step
    state_->iterateVec->set(x);
    state_->iterateVec->axpy(-lambda_,*state_->stepVec);
    proj_->project(*state_->iterateVec,outStream);
    s->set(*state_->iterateVec);
    s->axpy(-one,x);

    // Nonmonotone Linesearch
    ls_nfval = 0;
    obj.update(*state_->iterateVec,UPDATE_TRIAL);
    ftrial = obj.value(*state_->iterateVec,tol); ls_nfval++;
    alpha  = one;
    fmax   = *std::max_element(fqueue.begin(),fqueue.end());
    gs     = state_->gradientVec->apply(*s);
    if (verbosity_ > 1) {
      outStream << "  In SpectralGradientAlgorithm_B: Line Search"  << std::endl;
      outStream << "    Step size:                        " << alpha                << std::endl;
      outStream << "    Trial objective value:            " << ftrial               << std::endl;
      outStream << "    Max stored objective value:       " << fmax                 << std::endl;
      outStream << "    Computed reduction:               " << fmax-ftrial          << std::endl;
      outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -gs*gamma_*alpha     << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
    }
    while (ftrial > fmax + gamma_*alpha*gs && ls_nfval < maxit_) {
      alphaTmp = -half*alpha*alpha*gs/(ftrial-state_->value-alpha*gs);
      alpha    = (sigma1_*alpha <= alphaTmp && alphaTmp <= sigma2_*alpha) ? alphaTmp : rhodec_*alpha;
      state_->iterateVec->set(x);
      state_->iterateVec->axpy(alpha,*s);
      obj.update(*state_->iterateVec,UPDATE_TRIAL);
      ftrial = obj.value(*state_->iterateVec,tol); ls_nfval++;
      if (verbosity_ > 1) {
        outStream << "  In SpectralGradientAlgorithm_B: Line Search"  << std::endl;
        outStream << "    Step size:                        " << alpha                << std::endl;
        outStream << "    Trial objective value:            " << ftrial               << std::endl;
        outStream << "    Max stored objective value:       " << fmax                 << std::endl;
        outStream << "    Computed reduction:               " << fmax-ftrial          << std::endl;
        outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
        outStream << "    Sufficient decrease bound:        " << -gs*gamma_*alpha     << std::endl;
        outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
      }
    }
    state_->nfval += ls_nfval;
    if (static_cast<int>(fqueue.size()) == maxSize_) fqueue.pop_front();
    fqueue.push_back(ftrial);

    // Update state
    state_->iter++;
    state_->value = ftrial;
    state_->searchSize = alpha;
    x.set(*state_->iterateVec);
    obj.update(x,UPDATE_ACCEPT,state_->iter);

    // Compute spectral step length
    s->scale(alpha);
    y->set(*state_->gradientVec);
    y->scale(-one);
    obj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    y->plus(*state_->gradientVec);
    ys      = y->apply(*s);
    ss      = s->dot(*s);
    lambda_ = (ys<=eps ? lambdaMax_ : std::max(lambdaMin_,std::min(ss/ys,lambdaMax_)));
    state_->snorm = std::sqrt(ss);

    // Compute gradient step
    state_->stepVec->set(state_->gradientVec->dual());

    // Compute projected gradient norm
    s->set(x); s->axpy(-one,*state_->stepVec);
    proj_->project(*s,outStream);
    s->axpy(-one,x);
    state_->gnorm = s->norm();

    // Update Output
    output.push_back(print(printHeader_));
    if (verbosity_ > 0) outStream << print(printHeader_);
  }
  output.push_back(Algorithm_B<Real>::printExitStatus());
  if (verbosity_ > 0) outStream << Algorithm_B<Real>::printExitStatus();
  return output;
}

template<typename Real>
std::string SpectralGradientAlgorithm_B<Real>::printHeader( void ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(109,'-') << std::endl;
    hist << "Spectral projected gradient descent";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  value    - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  lambda   - Spectral step length" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "alpha";
  hist << std::setw(15) << std::left << "lambda";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::endl;
  return hist.str();
}

template<typename Real>
std::string SpectralGradientAlgorithm_B<Real>::printName( void ) const {
  std::stringstream hist;
  hist << std::endl << "Projected Spectral Gradient Method (Type B, Bound Constraints)" << std::endl;
  return hist.str();
}

template<typename Real>
std::string SpectralGradientAlgorithm_B<Real>::print( const bool print_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) {
    hist << printName();
  }
  if ( print_header ) {
    hist << printHeader();
  }
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(15) << std::left << lambda_;
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::setw(15) << std::left << state_->searchSize;
    hist << std::setw(15) << std::left << lambda_;
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::endl;
  }
  return hist.str();
}

} // namespace ROL

#endif
