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
void ProxGradientAlgorithm<Real>::initialize(Vector<Real>          &x,
                                             const Vector<Real>    &g,
                                             Objective<Real>       &sobj,
                                             Objective<Real>       &nobj,
																						 Vector<Real>          &px,
                                             std::ostream &outStream) {
  const Real one(1);
  // Initialize data
  TypeP::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  nobj.prox(px,x,state_->searchSize,ftol);
	x.set(px); //revisit with option to do initial prox or not
  sobj.update(x,UpdateType::Initial,state_->iter);
  nobj.update(x,UpdateType::Initial,state_->iter); 
	state_->value  = sobj.value(x,ftol); 
  state_->phival = nobj.value(x,ftol); 
  state_->nfval++;
	state_->nphi++; 
  sobj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-t0_,state_->gradientVec->dual());

  nobj.prox(px,*state_->stepVec,t0_,ftol);
	state_->stepVec->set(px); 
  Real fnew = state_->value;
  if (!useralpha_) {
    // Evaluate objective at P(x - g)
    sobj.update(*state_->stepVec,UpdateType::Trial);
		fnew = sobj.value(*state_->stepVec,ftol); 
    sobj.update(x,UpdateType::Revert);
    state_->nfval++;
  }
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm()/t0_;
  state_->snorm = ROL_INF<Real>();
  if (!useralpha_) {
    const Real half(0.5);
    // Minimize quadratic interpolate to compute new alpha
    //Real gs    = state_->stepVec->dot(state_->gradientVec->dual());
    Real gs    = state_->stepVec->apply(*state_->gradientVec);
    Real denom = (fnew - state_->value - gs);
    bool flag  = maxAlpha_ == alpha0_;
    alpha0_ = ((denom > ROL_EPSILON<Real>()) ? -half*gs/denom : alpha0bnd_);
    alpha0_ = ((alpha0_ > alpha0bnd_) ? alpha0_ : one);
    if (flag) maxAlpha_ = alpha0_;
  }
  // Normalize initial CP step length
  if (normAlpha_) {
    alpha0_ /= state_->gradientVec->norm();
  }
  state_->searchSize = alpha0_;
}

template<typename Real>
void ProxGradientAlgorithm<Real>::run( Vector<Real>          &x,
                                       const Vector<Real>    &g, 
                                       Objective<Real>       &sobj,
                                       Objective<Real>       &nobj,
                                       std::ostream          &outStream ) {
  const Real one(1);
  // Initialize trust-region data
	Ptr<Vector<Real>> px = x.clone(); 
  Ptr<Vector<Real>> s  = x.clone();
	initialize(x,g,sobj,nobj,*px,outStream);
  Real ftrial(0), phitrial(0), Qk(0), ftrialP(0),phitrialP(0), alphaP(0), tol(std::sqrt(ROL_EPSILON<Real>()));
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
    px->set(x);
    px->axpy(-state_->searchSize,*state_->stepVec);
    nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
    nobj.update(*state_->iterateVec,UpdateType::Trial);
    sobj.update(*state_->iterateVec,UpdateType::Trial);
    ftrial = sobj.value(*state_->iterateVec,tol);
		phitrial = nobj.value(*state_->iterateVec,tol); 
    ls_nfval = 1;
    s->set(*state_->iterateVec);
    s->axpy(-one,x);
    Qk = s->dot(*state_->stepVec) + phitrial - state_->phival;
    incAlpha = (state_->value + state_->phival - ftrial - phitrial >= -c1_*Qk);
    if (verbosity_ > 1) {
      outStream << "  In TypeP::GradientAlgorithm: Line Search"  << std::endl;
      outStream << "    Step size:                        " << state_->searchSize   << std::endl;
      outStream << "    Trial sobj value:                 " << ftrial               << std::endl;
      outStream << "    Trial nobj value:                 " << phitrial             << std::endl; 
			outStream << "    Computed reduction:               " << state_->value + state_->phival - phitrial -ftrial << std::endl;
      outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
      outStream << "    Increase alpha?:                  " << incAlpha             << std::endl;
    }
    if (incAlpha && useAdapt_) {
      ftrialP = ROL_INF<Real>();
			phitrialP = ROL_INF<Real>(); 
      while ( state_->value + state_->phival - ftrial - phitrial >= -c1_*Qk
           && ftrial+phitrial <= ftrialP + phitrialP
           && state_->searchSize < maxAlpha_
           && ls_nfval < maxit_ ) {
        // Previous value was acceptable
        sobj.update(*state_->iterateVec,UpdateType::Accept);
        nobj.update(*state_->iterateVec,UpdateType::Accept); 
				alphaP  = state_->searchSize;
        ftrialP = ftrial;
				phitrialP = phitrial; 
        state_->searchSize *= rhoinc_;
        state_->searchSize  = std::min(state_->searchSize,maxAlpha_);
        px->set(x);
        px->axpy(-state_->searchSize,*state_->stepVec);
        nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        nobj.update(*state_->iterateVec,UpdateType::Trial);
        ftrial = sobj.value(*state_->iterateVec,tol);
        phitrial = nobj.value(*state_->iterateVec,tol); 
				ls_nfval++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        Qk = s->dot(*state_->stepVec) + phitrial - state_->phival;
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "    Step size:                        " << state_->searchSize   << std::endl;
          outStream << "    Trial sobj value:            			" << ftrial               << std::endl;
          outStream << "    Trial nobj value:            			" << phitrial             << std::endl;
          outStream << "    Computed reduction:               " << state_->value + state_->phival - phitrial -ftrial << std::endl;
          outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
          outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
          outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
        }
      }
      if (state_->value + state_->phival - ftrial - phitrial < -c1_*Qk || phitrial + ftrial > phitrialP + ftrialP) {
        ftrial = ftrialP;
				phitrial = phitrialP; 
        state_->searchSize = alphaP;
        px->set(x);
        px->axpy(-state_->searchSize,*state_->stepVec);
        nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        accept = false;
      }
    }
    else {
      while ( state_->value + state_->phival - phitrial - ftrial < -c1_*Qk && ls_nfval < maxit_ ) {
        state_->searchSize *= rhodec_;
        px->set(x);
        px->axpy(-state_->searchSize,*state_->stepVec);
        nobj.prox(*state_->iterateVec,*px,state_->searchSize,tol);
				state_->nprox++;
        sobj.update(*state_->iterateVec,UpdateType::Trial);
        nobj.update(*state_->iterateVec,UpdateType::Trial);
        ftrial = sobj.value(*state_->iterateVec,tol);
				phitrial = nobj.value(*state_->iterateVec,tol); 
        ls_nfval++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        Qk = s->dot(*state_->stepVec) + phitrial - state_->phival;
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "    Step size:                        " << state_->searchSize   << std::endl;
          outStream << "    Trial sobj value:            			" << ftrial               << std::endl;
          outStream << "    Trial nobj value:            			" << phitrial             << std::endl;
          outStream << "    Computed reduction:               " << state_->value + state_->phival - phitrial -ftrial << std::endl;
          outStream << "    Dot product of gradient and step: " << Qk                   << std::endl;
          outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
          outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
        }
      }
    }
    state_->nfval += ls_nfval;
    state_->nphi  += ls_nfval; 
    // Compute norm of step
    state_->stepVec->set(*s);
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    x.set(*state_->iterateVec);

    // Compute new value and gradient
    state_->iter++;
    state_->value = ftrial;
		state_->phival = phitrial; 
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
    s->set(x); s->axpy(-one,*state_->stepVec);
    nobj.prox(*px,*s,state_->searchSize,tol);
    s->axpy(-one,x);
    state_->gnorm = s->norm();

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
    hist << "Projected gradient descent";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  sobj     - Smooth Objective function value" << std::endl;
    hist << "  nobj     - Nonsmooth Objective function value" << std::endl; 
		hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
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
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
	hist << std::setw(10) << std::left << "#nprox"; 
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void ProxGradientAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Proximal Gradient Descent with Backtracking Line Search (Type P)" << std::endl;
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
    hist << std::setw(10) << std::left << state_->nfval;
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
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
		hist << std::setw(10) << std::left << state_->nprox; 
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeP
} // namespace ROL

#endif
