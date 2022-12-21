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
                                             std::ostream       &outStream) {
  const Real one(1);
  // Initialize data
  TypeP::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  nobj.prox(px,x,state_->searchSize,ftol);
  state_->nprox++;
  x.set(px); //revisit with option to do initial prox or not
  // Evaluate objective function
  sobj.update(x,UpdateType::Initial,state_->iter);
  state_->svalue = sobj.value(x,ftol); 
  state_->nsval++;
  nobj.update(x,UpdateType::Initial,state_->iter); 
  state_->nvalue = nobj.value(x,ftol); 
  state_->nnval++; 
  state_->value  = state_->svalue + state_->nvalue;
  // Evaluate gradient of smooth part
  sobj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  // Evaluate proximal gradient
  state_->stepVec->set(x);
  state_->stepVec->axpy(-t0_,state_->gradientVec->dual());
  nobj.prox(px,*state_->stepVec,t0_,ftol);
  state_->nprox++;
  state_->stepVec->set(px);
  // Compute initial step size
  Real fnew = state_->svalue;
  if (!useralpha_) {
    // Evaluate objective at P(x - g)
    sobj.update(*state_->stepVec,UpdateType::Trial);
    fnew = sobj.value(*state_->stepVec,ftol); 
    sobj.update(x,UpdateType::Revert);
    state_->nsval++;
  }
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm()/t0_;
  state_->snorm = ROL_INF<Real>();
  if (!useralpha_) {
    const Real half(0.5);
    // Minimize quadratic interpolate to compute new alpha
    //Real gs    = state_->stepVec->dot(state_->gradientVec->dual());
    Real gs    = state_->stepVec->apply(*state_->gradientVec);
    Real denom = (fnew - state_->svalue - gs);
    bool flag  = maxAlpha_ == alpha0_;
    alpha0_ = ((denom > ROL_EPSILON<Real>()) ? -half*gs/denom : alpha0bnd_);
    alpha0_ = ((alpha0_ > alpha0bnd_) ? alpha0_ : one);
    if (flag) maxAlpha_ = alpha0_;
  }
  // Normalize initial CP step length
  if (normAlpha_)
    alpha0_ /= state_->gradientVec->norm();
  state_->searchSize = alpha0_;
}

template<typename Real>
void ProxGradientAlgorithm<Real>::run( Vector<Real>       &x,
                                       const Vector<Real> &g, 
                                       Objective<Real>    &sobj,
                                       Objective<Real>    &nobj,
                                       std::ostream       &outStream ) {
  const Real one(1);
  // Initialize trust-region data
  Ptr<Vector<Real>> px = x.clone(), s = x.clone();
  initialize(x,g,sobj,nobj,*px,outStream);
  Real strial(0), ntrial(0), Ftrial(0), strialP(0), ntrialP(0), FtrialP(0);
  Real Qk(0), alphaP(0), tol(std::sqrt(ROL_EPSILON<Real>()));
  int ls_nfval = 0;
  bool incAlpha = false, accept = true;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  // Compute steepest descent step
  state_->stepVec->set(state_->gradientVec->dual());
  while (status_->check(*state_)) {
    accept = true;
    // Perform backtracking line search 
    if (!usePrevAlpha_ && !useAdapt_) state_->searchSize = alpha0_;
    // Compute proximal gradient step with initial search size
    px->set(x);
    px->axpy(-state_->searchSize,*state_->stepVec);
    nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
    state_->nprox++;
    s->set(*state_->iterateVec);
    s->axpy(-one,x);
    // Compute objective function values
    sobj.update(*state_->iterateVec,UpdateType::Trial);
    strial = sobj.value(*state_->iterateVec,tol);
    nobj.update(*state_->iterateVec,UpdateType::Trial);
    ntrial = nobj.value(*state_->iterateVec,tol); 
    Ftrial = strial + ntrial;
    ls_nfval = 1;
    // Compute decrease indicator
    Qk = s->dot(*state_->stepVec) + ntrial - state_->nvalue;
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
        // Increase search size
        alphaP  = state_->searchSize;
        strialP = strial;
        ntrialP = ntrial; 
        FtrialP = Ftrial;
        state_->searchSize *= rhoinc_;
        state_->searchSize  = std::min(state_->searchSize,maxAlpha_);
        // Compute proximal gradient step with new search size
        px->set(x);
        px->axpy(-state_->searchSize,*state_->stepVec);
        nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
        state_->nprox++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        // Compute objective function values
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        strial = sobj.value(*state_->iterateVec,tol);
        nobj.update(*state_->iterateVec,UpdateType::Trial);
        ntrial = nobj.value(*state_->iterateVec,tol); 
        Ftrial = strial + ntrial;
        ls_nfval++;
        // Compute decrease indicator
        Qk = s->dot(*state_->stepVec) + ntrial - state_->nvalue;
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
        strial = strialP;
        ntrial = ntrialP; 
        Ftrial = FtrialP;
        state_->searchSize = alphaP;
        // Recompute proximal gradient step
        px->set(x);
        px->axpy(-state_->searchSize,*state_->stepVec);
        nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
        state_->nprox++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        accept = false;
      }
    }
    else {
      while ( Ftrial - state_->value > c1_*Qk && ls_nfval < maxit_ ) {
        // Decrease search size
        state_->searchSize *= rhodec_;
        // Compute proximal gradient step with new search size
        px->set(x);
        px->axpy(-state_->searchSize,*state_->stepVec);
        nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
        state_->nprox++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        // Compute objective function values
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        strial = sobj.value(*state_->iterateVec,tol);
        nobj.update(*state_->iterateVec,UpdateType::Trial);
        ntrial = nobj.value(*state_->iterateVec,tol); 
        Ftrial = strial + ntrial;
        ls_nfval++;
        // Compute decrease indicator
        Qk = s->dot(*state_->stepVec) + ntrial - state_->nvalue;
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
    state_->stepVec->set(*s);
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    x.set(*state_->iterateVec);

    // Compute new value and gradient
    state_->iter++;
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

    // Compute steepest descent step
    state_->stepVec->set(state_->gradientVec->dual());

    // Compute projected gradient norm
    s->set(x); s->axpy(-t0_,*state_->stepVec);
    nobj.prox(*px,*s,t0_,tol);
    state_->nprox++;
    s->axpy(-one,x);
    state_->gnorm = s->norm() / t0_;

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeP::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(109,'-') << std::endl;
    hist << "Proximal gradient descent";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  obj      - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the proximal gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  #sval    - Cumulative number of times the smooth objective function was evaluated" << std::endl;
    hist << "  #nval    - Cumulative number of times the nonsmooth objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  #prox    - Cumulative number of times the proximal operator was computed" << std::endl; 
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "alpha";
  hist << std::setw(10) << std::left << "#sval";
  hist << std::setw(10) << std::left << "#nval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "#nprox"; 
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Proximal Gradient Descent with Bidirectional Line Search (Type P)" << std::endl;
  os << hist.str();
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(10) << std::left << state_->nsval;
    hist << std::setw(10) << std::left << state_->nnval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << state_->nprox; 
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::setw(15) << std::left << state_->searchSize;
    hist << std::setw(10) << std::left << state_->nsval;
    hist << std::setw(10) << std::left << state_->nnval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << state_->nprox; 
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeP
} // namespace ROL

#endif
