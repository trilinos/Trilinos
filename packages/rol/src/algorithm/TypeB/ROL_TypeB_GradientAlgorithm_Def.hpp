// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_GRADIENTALGORITHM_DEF_HPP
#define ROL_TYPEB_GRADIENTALGORITHM_DEF_HPP

namespace ROL {
namespace TypeB {

template<typename Real>
GradientAlgorithm<Real>::GradientAlgorithm(ParameterList &list) {
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
  verbosity_    = list.sublist("General").get("Output Level",                     0);
  writeHeader_  = verbosity_ > 2;
}

template<typename Real>
void GradientAlgorithm<Real>::initialize(Vector<Real>          &x,
                                         const Vector<Real>    &g,
                                         Objective<Real>       &obj,
                                         BoundConstraint<Real> &bnd,
                                         std::ostream &outStream) {
  const Real one(1);
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
  }
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  proj_->project(x,outStream);
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream);
  Real fnew = state_->value;
  if (!useralpha_) {
    // Evaluate objective at P(x - g)
    obj.update(*state_->stepVec,UpdateType::Trial);
    fnew = obj.value(*state_->stepVec,ftol);
    obj.update(x,UpdateType::Revert);
    state_->nfval++;
  }
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
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
void GradientAlgorithm<Real>::run( Vector<Real>          &x,
                                   const Vector<Real>    &g, 
                                   Objective<Real>       &obj,
                                   BoundConstraint<Real> &bnd,
                                   std::ostream          &outStream ) {
  const Real one(1);
  // Initialize trust-region data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone();
  Real ftrial(0), gs(0), ftrialP(0), alphaP(0), tol(std::sqrt(ROL_EPSILON<Real>()));
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
    state_->iterateVec->set(x);
    state_->iterateVec->axpy(-state_->searchSize,*state_->stepVec);
    proj_->project(*state_->iterateVec,outStream);
    obj.update(*state_->iterateVec,UpdateType::Trial);
    ftrial = obj.value(*state_->iterateVec,tol);
    ls_nfval = 1;
    s->set(*state_->iterateVec);
    s->axpy(-one,x);
    gs = s->dot(*state_->stepVec);
    incAlpha = (state_->value - ftrial >= -c1_*gs);
    if (verbosity_ > 1) {
      outStream << "  In TypeB::GradientAlgorithm: Line Search"  << std::endl;
      outStream << "    Step size:                        " << state_->searchSize   << std::endl;
      outStream << "    Trial objective value:            " << ftrial               << std::endl;
      outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
      outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -gs*c1_              << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
      outStream << "    Increase alpha?:                  " << incAlpha             << std::endl;
    }
    if (incAlpha && useAdapt_) {
      ftrialP = ROL_INF<Real>();
      while ( state_->value - ftrial >= -c1_*gs
           && ftrial <= ftrialP
           && state_->searchSize < maxAlpha_
           && ls_nfval < maxit_ ) {
        // Previous value was acceptable
        obj.update(*state_->iterateVec,UpdateType::Accept);
        alphaP  = state_->searchSize;
        ftrialP = ftrial;
        state_->searchSize *= rhoinc_;
        state_->searchSize  = std::min(state_->searchSize,maxAlpha_);
        state_->iterateVec->set(x);
        state_->iterateVec->axpy(-state_->searchSize,*state_->stepVec);
        proj_->project(*state_->iterateVec,outStream);
        obj.update(*state_->iterateVec,UpdateType::Trial);
        ftrial = obj.value(*state_->iterateVec,tol);
        ls_nfval++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        gs = s->dot(*state_->stepVec);
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "    Step size:                        " << state_->searchSize   << std::endl;
          outStream << "    Trial objective value:            " << ftrial               << std::endl;
          outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
          outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
          outStream << "    Sufficient decrease bound:        " << -gs*c1_              << std::endl;
          outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
        }
      }
      if (state_->value - ftrial < -c1_*gs || ftrial > ftrialP) {
        ftrial = ftrialP;
        state_->searchSize = alphaP;
        state_->iterateVec->set(x);
        state_->iterateVec->axpy(-state_->searchSize,*state_->stepVec);
        proj_->project(*state_->iterateVec,outStream);
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        accept = false;
      }
    }
    else {
      while ( state_->value - ftrial < -c1_*gs && ls_nfval < maxit_ ) {
        state_->searchSize *= rhodec_;
        state_->iterateVec->set(x);
        state_->iterateVec->axpy(-state_->searchSize,*state_->stepVec);
        proj_->project(*state_->iterateVec,outStream);
        obj.update(*state_->iterateVec,UpdateType::Trial);
        ftrial = obj.value(*state_->iterateVec,tol);
        ls_nfval++;
        s->set(*state_->iterateVec);
        s->axpy(-one,x);
        gs = s->dot(*state_->stepVec);
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "    Step size:                        " << state_->searchSize   << std::endl;
          outStream << "    Trial objective value:            " << ftrial               << std::endl;
          outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
          outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
          outStream << "    Sufficient decrease bound:        " << -gs*c1_              << std::endl;
          outStream << "    Number of function evaluations:   " << ls_nfval             << std::endl;
        }
      }
    }
    state_->nfval += ls_nfval;

    // Compute norm of step
    state_->stepVec->set(*s);
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    x.set(*state_->iterateVec);

    // Compute new value and gradient
    state_->iter++;
    state_->value = ftrial;
    if (accept) obj.update(x,UpdateType::Accept,state_->iter);
    else        obj.update(x,UpdateType::Revert,state_->iter);
    obj.gradient(*state_->gradientVec,x,tol);
    state_->ngrad++;

    // Compute steepest descent step
    state_->stepVec->set(state_->gradientVec->dual());

    // Compute projected gradient norm
    s->set(x); s->axpy(-one,*state_->stepVec);
    proj_->project(*s,outStream);
    s->axpy(-one,x);
    state_->gnorm = s->norm();

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void GradientAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(109,'-') << std::endl;
    os << "Projected gradient descent";
    os << " status output definitions" << std::endl << std::endl;
    os << "  iter     - Number of iterates (steps taken)" << std::endl;
    os << "  value    - Objective function value" << std::endl;
    os << "  gnorm    - Norm of the gradient" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  alpha    - Line search step length" << std::endl;
    os << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    os << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    os << std::string(109,'-') << std::endl;
  }

  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(15) << std::left << "alpha";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void GradientAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Projected Gradient Descent with Backtracking Line Search (Type B, Bound Constraints)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void GradientAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(10) << std::left << state_->nfval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << state_->snorm;
    os << std::setw(15) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << state_->nfval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeB
} // namespace ROL

#endif
