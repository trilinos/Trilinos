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

#ifndef ROL_NEWTONKRYLOVALGORITHM_B_DEF_H
#define ROL_NEWTONKRYLOVALGORITHM_B_DEF_H

namespace ROL {

template<typename Real>
NewtonKrylovAlgorithm_B<Real>::NewtonKrylovAlgorithm_B(ParameterList           &list,
                                                       const Ptr<Secant<Real>> &secant)
  : secant_(secant), esec_(SECANT_USERDEFINED) {
  parseParameterList(list);

  if ( secant_ == nullPtr ) {
    secantName_ = list.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS");
    esec_ = StringToESecant(secantName_);
    secant_ = SecantFactory<Real>(list);
  }
  else {
    secantName_ = list.sublist("General").sublist("Secant").get("User Defined Secant Name",
                                                                "Unspecified User Defined Secant Method");
  }

  krylovName_ = list.sublist("General").sublist("Krylov").get("Type","Conjugate Gradients");
  ekv_ = StringToEKrylov(krylovName_);
  krylov_ = KrylovFactory<Real>(list);
}

template<typename Real>
NewtonKrylovAlgorithm_B<Real>::NewtonKrylovAlgorithm_B(ParameterList           &list,
                                                       const Ptr<Krylov<Real>> &krylov,
                                                       const Ptr<Secant<Real>> &secant)
  : secant_(secant), esec_(SECANT_USERDEFINED), krylov_(krylov), ekv_(KRYLOV_USERDEFINED) {
  parseParameterList(list);

  if ( secant_ == nullPtr ) {
    secantName_ = list.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS");
    esec_ = StringToESecant(secantName_);
    secant_ = SecantFactory<Real>(list);
  }
  else {
    secantName_ = list.sublist("General").sublist("Secant").get("User Defined Secant Name",
                                                                "Unspecified User Defined Secant Method");
  }

  krylovName_ = list.sublist("General").sublist("Krylov").get("User Defined Krylov Name",
                                                              "Unspecified User Defined Krylov Method");
}

template<typename Real>
void NewtonKrylovAlgorithm_B<Real>::parseParameterList(ParameterList &list) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("Line Search");
  maxit_        = lslist.get("Function Evaluation Limit",                        20);
  alpha0_       = lslist.get("Initial Step Size",                               1.0);
  useralpha_    = lslist.get("User Defined Initial Step Size",                false);
  usePrevAlpha_ = lslist.get("Use Previous Step Length as Initial Guess",     false);
  c1_           = lslist.get("Sufficient Decrease Tolerance",                  1e-4);
  rhodec_       = lslist.sublist("Line-Search Method").get("Backtracking Rate", 0.5);

  useSecantHessVec_ = list.sublist("General").sublist("Secant").get("Use as Hessian",        false);
  useSecantPrecond_ = list.sublist("General").sublist("Secant").get("Use as Preconditioner", false);

  verbosity_    = list.sublist("General").get("Output Level",                     0);
  printHeader_  = verbosity_ > 2;
}

template<typename Real>
void NewtonKrylovAlgorithm_B<Real>::initialize(Vector<Real>          &x,
                                               const Vector<Real>    &g,
                                               Objective<Real>       &obj,
                                               BoundConstraint<Real> &bnd,
                                               std::ostream &outStream) {
  const Real one(1);
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
  }
  // Initialize data
  Algorithm_B<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  proj_->project(x,outStream);
  state_->iterateVec->set(x);
  obj.update(x,UPDATE_INITIAL,state_->iter);    
  state_->value = obj.value(x,ftol); state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream);
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  if (!useralpha_) {
    alpha0_ = one;
  }
  state_->searchSize = alpha0_;
}

template<typename Real>
std::vector<std::string> NewtonKrylovAlgorithm_B<Real>::run( Vector<Real>          &x,
                                                             const Vector<Real>    &g, 
                                                             Objective<Real>       &obj,
                                                             BoundConstraint<Real> &bnd,
                                                             std::ostream          &outStream ) {
  const Real one(1);
  // Initialize trust-region data
  std::vector<std::string> output;
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone(), gp = x.clone(), gold = g.clone();
  Ptr<Vector<Real>> pwa = x.clone(), dwa = g.clone();
  Real ftrial(0), gs(0), tol(std::sqrt(ROL_EPSILON<Real>()));

  Ptr<LinearOperator<Real>> hessian, precond;

  // Output
  output.push_back(print(true));
  if (verbosity_ > 0) outStream << print(true);

  // Compute steepest descent step
  gp->set(state_->gradientVec->dual());
  while (status_->check(*state_)) {
    // Compute step
    hessian = makePtr<HessianPNK>(makePtrFromRef(obj),makePtrFromRef(bnd),
                                  state_->iterateVec,state_->gradientVec,state_->gnorm,
                                  secant_,useSecantHessVec_,pwa);
    precond = makePtr<PrecondPNK>(makePtrFromRef(obj),makePtrFromRef(bnd),
                                  state_->iterateVec,state_->gradientVec,state_->gnorm,
                                  secant_,useSecantPrecond_,dwa);
    flagKrylov_ = 0;
    krylov_->run(*s,*hessian,*state_->gradientVec,*precond,iterKrylov_,flagKrylov_);
    if (flagKrylov_ == 2 && iterKrylov_ <= 1) {
      s->set(*gp);
    }
    // Perform backtracking line search 
    if (!usePrevAlpha_) state_->searchSize = alpha0_;
    x.set(*state_->iterateVec);
    x.axpy(-state_->searchSize,*s);
    proj_->project(x,outStream);
    obj.update(x,UPDATE_TRIAL);
    ftrial = obj.value(x,tol); ls_nfval_ = 1;
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    gs = state_->stepVec->dot(*gp);
    if (verbosity_ > 1) {
      outStream << "  In NewtonKrylovAlgorithm_B: Line Search"  << std::endl;
      outStream << "    Step size:                        " << state_->searchSize   << std::endl;
      outStream << "    Trial objective value:            " << ftrial               << std::endl;
      outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
      outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -gs*c1_              << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval_            << std::endl;
    }
    while ( state_->value - ftrial < -c1_*gs && ls_nfval_ < maxit_ ) {
      state_->searchSize *= rhodec_;
      x.set(*state_->iterateVec);
      x.axpy(-state_->searchSize,*s);
      proj_->project(x,outStream);
      obj.update(x,UPDATE_TRIAL);
      ftrial = obj.value(x,tol); ls_nfval_++;
      state_->stepVec->set(x);
      state_->stepVec->axpy(-one,*state_->iterateVec);
      gs = state_->stepVec->dot(*gp);
      if (verbosity_ > 1) {
        outStream << std::endl;
        outStream << "    Step size:                        " << state_->searchSize   << std::endl;
        outStream << "    Trial objective value:            " << ftrial               << std::endl;
        outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
        outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
        outStream << "    Sufficient decrease bound:        " << -gs*c1_              << std::endl;
        outStream << "    Number of function evaluations:   " << ls_nfval_            << std::endl;
      }
    }
    state_->nfval += ls_nfval_;

    // Compute norm of step
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    state_->iterateVec->set(x);

    // Compute new value and gradient
    state_->iter++;
    state_->value = ftrial;
    obj.update(x,UPDATE_ACCEPT,state_->iter);
    gold->set(*state_->gradientVec);
    obj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    gp->set(state_->gradientVec->dual());

    // Compute projected gradient norm
    s->set(x); s->axpy(-one,*gp);
    proj_->project(*s,outStream);
    s->axpy(-one,x);
    state_->gnorm = s->norm();

    // Update secant
    secant_->updateStorage(x,*state_->gradientVec,*gold,*state_->stepVec,state_->snorm,state_->iter);

    // Update Output
    output.push_back(print(printHeader_));
    if (verbosity_ > 0) outStream << print(printHeader_);
  }
  output.push_back(Algorithm_B<Real>::printExitStatus());
  if (verbosity_ > 0) outStream << Algorithm_B<Real>::printExitStatus();
  return output;
}

template<typename Real>
std::vector<std::string> NewtonKrylovAlgorithm_B<Real>::run( NewOptimizationProblem<Real> &problem,
                                                             std::ostream                 &outStream ) {
  if (problem.getPolyhedralProjection() == nullPtr) {
    return Algorithm_B<Real>::run(problem,outStream);
  }
  else {
    throw Exception::NotImplemented(">>> NewtonKrylovAlgorithm_B::run : This algorithm cannot solve problems with linear equality constraints!");
  }
}

template<typename Real>
std::vector<std::string> NewtonKrylovAlgorithm_B<Real>::run( Vector<Real>          &x,
                                                             const Vector<Real>    &g,
                                                             Objective<Real>       &obj,
                                                             BoundConstraint<Real> &bnd,
                                                             Constraint<Real>      &linear_econ,
                                                             Vector<Real>          &linear_emul,
                                                             const Vector<Real>    &linear_eres,
                                                             std::ostream          &outStream ) {
  throw Exception::NotImplemented(">>> NewtonKrylovAlgorithm_B::run : This algorithm cannot solve problems with linear equality constraints!");
}

template<typename Real>
std::string NewtonKrylovAlgorithm_B<Real>::printHeader( void ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    if (!useSecantHessVec_) {
      hist << "Line-Search Projected Newton";
    }
    else {
      hist << "Line-Search Projected Quasi-Newton with " << secantName_ << " Hessian approximation";
    }
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  value    - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  ls_#fval - Number of the times the objective function was evaluated during the line search" << std::endl;
    hist << "  iterCG   - Number of Krylov iterations" << std::endl << std::endl;
    hist << "  flagGC   - Krylov flag" << std::endl;
    for( int flag = CG_FLAG_SUCCESS; flag != CG_FLAG_UNDEFINED; ++flag ) {
      hist << "    " << NumberToString(flag) << " - "
           << ECGFlagToString(static_cast<ECGFlag>(flag)) << std::endl;
    }            
    hist << std::string(114,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "alpha";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "#ls_fval";
  hist << std::setw(10) << std::left << "iterCG";
  hist << std::setw(10) << std::left << "flagCG";
  hist << std::endl;
  return hist.str();
}

template<typename Real>
std::string NewtonKrylovAlgorithm_B<Real>::printName( void ) const {
  std::stringstream hist;
  if (!useSecantHessVec_) {
    hist << std::endl << "Line-Search Projected Newton (Type B, Bound Constraints)" << std::endl;
  }
  else {
    hist << std::endl << "Line-Search Projected Quasi-Newton with "
         << secantName_ << " Hessian approximation" << std::endl;
  }
  return hist.str();
}

template<typename Real>
std::string NewtonKrylovAlgorithm_B<Real>::print( const bool print_header ) const {
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
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << "---";
    hist << std::setw(10) << std::left << "---";
    hist << std::setw(10) << std::left << "---";
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
    hist << std::setw(10) << std::left << ls_nfval_;
    hist << std::setw(10) << std::left << iterKrylov_;
    hist << std::setw(10) << std::left << flagKrylov_;
    hist << std::endl;
  }
  return hist.str();
}

} // namespace ROL

#endif
