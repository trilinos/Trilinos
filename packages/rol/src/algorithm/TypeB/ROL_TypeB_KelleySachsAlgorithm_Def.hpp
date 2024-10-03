// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_KELLEYSACHSALGORITHM_DEF_HPP
#define ROL_TYPEB_KELLEYSACHSALGORITHM_DEF_HPP

namespace ROL {
namespace TypeB {

template<typename Real>
KelleySachsAlgorithm<Real>::KelleySachsAlgorithm(ParameterList &list,
                                                 const Ptr<Secant<Real>> &secant) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  ParameterList &trlist = list.sublist("Step").sublist("Trust Region");
  // Trust-Region Parameters
  state_->searchSize = trlist.get("Initial Radius",            -1.0);
  delMax_ = trlist.get("Maximum Radius",                       ROL_INF<Real>());
  eta0_   = trlist.get("Step Acceptance Threshold",            0.05);
  eta1_   = trlist.get("Radius Shrinking Threshold",           0.05);
  eta2_   = trlist.get("Radius Growing Threshold",             0.9);
  gamma0_ = trlist.get("Radius Shrinking Rate (Negative rho)", 0.0625);
  gamma1_ = trlist.get("Radius Shrinking Rate (Positive rho)", 0.25);
  gamma2_ = trlist.get("Radius Growing Rate",                  2.5);
  TRsafe_ = trlist.get("Safeguard Size",                       100.0);
  eps_    = TRsafe_*ROL_EPSILON<Real>();
  // Krylov Parameters
  maxit_ = list.sublist("General").sublist("Krylov").get("Iteration Limit",    20);
  tol1_  = list.sublist("General").sublist("Krylov").get("Absolute Tolerance", 1e-4);
  tol2_  = list.sublist("General").sublist("Krylov").get("Relative Tolerance", 1e-2);
  // Algorithm-Specific Parameters
  minit_  = trlist.sublist("Kelley-Sachs").get("Maximum Number of Smoothing Iterations", 20);
  mu0_    = trlist.sublist("Kelley-Sachs").get("Sufficient Decrease Parameter",          1e-4);
  mu1_    = trlist.sublist("Kelley-Sachs").get("Post-Smoothing Decrease Parameter",      0.9999);
  eps0_   = trlist.sublist("Kelley-Sachs").get("Binding Set Tolerance",                  1e-3);
  beta_   = trlist.sublist("Kelley-Sachs").get("Post-Smoothing Backtracking Rate",       1e-2);
  alpha0_ = trlist.sublist("Kelley-Sachs").get("Initial Post-Smoothing Step Size",       1.0);
  // Output Parameters
  verbosity_   = list.sublist("General").get("Output Level",0);
  writeHeader_ = verbosity_ > 2;
  // Secant Information
  useSecantPrecond_ = list.sublist("General").sublist("Secant").get("Use as Preconditioner", false);
  useSecantHessVec_ = list.sublist("General").sublist("Secant").get("Use as Hessian",        false);
  // Initialize trust region model
  model_ = makePtr<TrustRegionModel_U<Real>>(list,secant);
  useSecantPrecond_ = list.sublist("General").sublist("Secant").get("Use as Preconditioner", false);
  useSecantHessVec_ = list.sublist("General").sublist("Secant").get("Use as Hessian",        false);
  if (secant == nullPtr) {
    esec_ = StringToESecant(list.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS"));
  }
}

template<typename Real>
void KelleySachsAlgorithm<Real>::initialize(Vector<Real>          &x,
                                            const Vector<Real>    &g,
                                            Objective<Real>       &obj,
                                            BoundConstraint<Real> &bnd,
                                            std::ostream &outStream) {
  const Real one(1);
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  nhess_ = 0;
  // Update approximate gradient and approximate objective function.
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(); 
  bnd.project(x);
  state_->iterateVec->set(x);
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  bnd.project(*state_->stepVec);
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  // Compute initial trust region radius if desired.
  if ( state_->searchSize <= static_cast<Real>(0) ) {
    state_->searchSize = state_->gradientVec->norm();
  }
}

template<typename Real>
void KelleySachsAlgorithm<Real>::run(Vector<Real>          &x,
                                     const Vector<Real>    &g, 
                                     Objective<Real>       &obj,
                                     BoundConstraint<Real> &bnd,
                                     std::ostream          &outStream ) {
  const Real one(1), xeps(ROL_EPSILON<Real>());
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  Real gfnorm(0), gfnormf(0), tol(0), stol(0), eps(0);
  Real ftrial(0), fnew(0), pRed(0), rho(1), alpha(1), lambda(1);
  int cnt(0);
  // Initialize trust-region data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone(), gfree = g.clone();
  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone(), pwa3 = x.clone();
  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone(), dwa3 = g.clone();

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  // Compute initial free gradient
  gfree->set(*state_->gradientVec);
  //bnd.pruneActive(*gfree,*state_->gradientVec,x,xeps,eps);
  //gfnorm = gfree->norm();
  pwa1->set(gfree->dual());
  bnd.pruneActive(*pwa1,state_->gradientVec->dual(),x,xeps,eps);
  gfree->set(pwa1->dual());
  gfnorm = gfree->norm();

  // Initial tolerances
  eps  = std::min(eps0_,std::sqrt(state_->gnorm));
  tol  = tol2_*std::sqrt(state_->gnorm);
  stol = std::min(tol1_,tol*gfnorm);

  while (status_->check(*state_)) {
    // Build trust-region model
    model_->setData(obj,*state_->iterateVec,*state_->gradientVec,ftol);

    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
    pRed = trpcg(*s,SPflag_,SPiter_,*gfree,x,*state_->gradientVec,
                 state_->searchSize,*model_,bnd,eps,stol,maxit_,
                 *pwa1,*dwa1,*pwa2,*dwa2,*pwa3,*dwa3,outStream);
    x.plus(*s);
    bnd.project(x);

    // Compute trial objective value
    obj.update(x,UpdateType::Trial);
    ftrial = obj.value(x,ftol); state_->nfval++;

    // Compute ratio of acutal and predicted reduction
    TRflag_ = TRUtils::SUCCESS;
    TRUtils::analyzeRatio<Real>(rho,TRflag_,state_->value,ftrial,pRed,eps_,outStream,verbosity_>1);

    // Check sufficient decrease
    if ( rho >= eta0_ ) {
      lambda = std::min(one, state_->searchSize/gfnorm);
      // Compute Scaled Measure || x - P( x - lam * PI(g) ) ||
      pwa1->set(*state_->iterateVec);
      pwa1->axpy(-lambda,gfree->dual());
      bnd.project(*pwa1);
      pwa1->axpy(-one,*state_->iterateVec);
      gfnormf = pwa1->norm();
      // Sufficient decrease?
      if (state_->value-ftrial < mu0_*gfnormf*state_->gnorm) {
        TRflag_ = TRUtils::QMINSUFDEC;
      }
      if ( verbosity_ > 1 ) {
        outStream << "    Norm of projected free gradient:         " << gfnormf                    << std::endl;
        outStream << "    Decrease lower bound (constraints):      " << mu0_*gfnormf*state_->gnorm << std::endl;
        outStream << "    Trust-region flag (constraints):         " << TRflag_                    << std::endl;
        outStream << std::endl;
      }
    }

    // Update algorithm state
    state_->iter++;
    // Accept/reject step and update trust region radius
    if ((rho < eta0_ && TRflag_ == TRUtils::SUCCESS) || (TRflag_ >= 2)) {
      state_->stepVec->set(x);
      state_->stepVec->axpy(-one,*state_->iterateVec);
      state_->snorm = state_->stepVec->norm();
      x.set(*state_->iterateVec);
      obj.update(x,UpdateType::Revert,state_->iter);
      // Decrease trust-region radius
      state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
    }
    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
             || (TRflag_ == TRUtils::POSPREDNEG)) {
      if (rho >= eta0_ && rho < eta1_) {
        // Decrease trust-region radius
        state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
      }
      else if (rho >= eta2_) {
        // Increase trust-region radius
        state_->searchSize = std::min(delMax_,gamma2_*state_->searchSize);
      }
      // Projected search
      cnt = 0;
      alpha = one;
      obj.gradient(*dwa1,x,ftol); state_->ngrad++;
      pwa2->set(dwa1->dual());
      pwa1->set(x); pwa1->axpy(-alpha/alpha0_,*pwa2);
      bnd.project(*pwa1);
      obj.update(*pwa1,UpdateType::Trial);
      fnew = obj.value(*pwa1,ftol); state_->nfval++;
      while ((fnew-ftrial) >= mu1_*(state_->value-ftrial)) {
        alpha *= beta_;
        pwa1->set(x); pwa1->axpy(-alpha/alpha0_,*pwa2);
        bnd.project(*pwa1);
        obj.update(*pwa1,UpdateType::Trial);
        fnew = obj.value(*pwa1,ftol); state_->nfval++;
        if ( cnt >= minit_ ) break;
        cnt++;
      }
      state_->stepVec->set(*pwa1);
      state_->stepVec->axpy(-one,*state_->iterateVec);
      state_->snorm = state_->stepVec->norm();
      // Store objective function and iteration information
      if (std::isnan(fnew)) {
        TRflag_ = TRUtils::TRNAN;
        rho     = -one;
        x.set(*state_->iterateVec);
        obj.update(x,UpdateType::Revert,state_->iter);
        // Decrease trust-region radius
        state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
      }
      else {
        // Update objective function information
        x.set(*pwa1);
        state_->iterateVec->set(x);
        state_->value = fnew;
        dwa1->set(*state_->gradientVec);
        obj.update(x,UpdateType::Accept,state_->iter);
        obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
        // Compute free gradient
        gfree->set(*state_->gradientVec);
        //bnd.pruneActive(*gfree,*state_->gradientVec,x,xeps,eps);
        //gfnorm = gfree->norm();
        pwa1->set(gfree->dual());
        bnd.pruneActive(*pwa1,state_->gradientVec->dual(),x,xeps,eps);
        gfree->set(pwa1->dual());
        gfnorm = gfree->norm();
        // Compute criticality measure
        pwa1->set(x);
        pwa1->axpy(-one,state_->gradientVec->dual());
        bnd.project(*pwa1);
        pwa1->axpy(-one,x);
        state_->gnorm = pwa1->norm();
        // Update secant information in trust-region model
        model_->update(x,*state_->stepVec,*dwa1,*state_->gradientVec,
                       state_->snorm,state_->iter);
        // Update tolerances
        eps  = std::min(eps0_,std::sqrt(state_->gnorm));
        tol  = tol2_*std::sqrt(state_->gnorm);
        stol = std::min(tol1_,tol*gfnorm);
      }
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void KelleySachsAlgorithm<Real>::run( Problem<Real> &problem,
                                      std::ostream  &outStream ) {
  if (problem.getPolyhedralProjection() == nullPtr) {
    TypeB::Algorithm<Real>::run(problem,outStream);
  }
  else {
    throw Exception::NotImplemented(">>> TypeB::KelleySachsAlgorithm::run : This algorithm cannot solve problems with linear equality constraints!");
  }
}

template<typename Real>
void KelleySachsAlgorithm<Real>::run( Vector<Real>          &x,
                                      const Vector<Real>    &g,
                                      Objective<Real>       &obj,
                                      BoundConstraint<Real> &bnd,
                                      Constraint<Real>      &linear_econ,
                                      Vector<Real>          &linear_emul,
                                      const Vector<Real>    &linear_eres,
                                      std::ostream          &outStream ) {
  throw Exception::NotImplemented(">>> TypeB::KelleySachsAlgorithm::run : This algorithm cannot solve problems with linear equality constraints!");
}

template<typename Real>
Real KelleySachsAlgorithm<Real>::trqsol(const Real xtx,
                                        const Real ptp,
                                        const Real ptx,
                                        const Real del) const {
  const Real zero(0);
  Real dsq = del*del;
  Real rad = ptx*ptx + ptp*(dsq-xtx);
  rad = std::sqrt(std::max(rad,zero));
  Real sigma(0);
  if (ptx > zero) {
    sigma = (dsq-xtx)/(ptx+rad);
  }
  else if (rad > zero) {
    sigma = (rad-ptx)/ptp;
  }
  else {
    sigma = zero;
  }
  return sigma;
}

template<typename Real>
Real KelleySachsAlgorithm<Real>::trpcg(Vector<Real> &w, int &iflag, int &iter,
                                       const Vector<Real> &g, const Vector<Real> &x,
                                       const Vector<Real> &g0,
                                       const Real del, TrustRegionModel_U<Real> &model,
                                       BoundConstraint<Real> &bnd, const Real eps,
                                       const Real tol, const int itermax,
                                       Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
                                       Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa,
                                       std::ostream &outStream) const {
  // p = step (primal)
  // q = hessian applied to step p (dual)
  // t = gradient (dual)
  // r = preconditioned gradient (primal)
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  const Real zero(0), half(0.5), one(1), two(2);
  Real rho(0), kappa(0), beta(0), sigma(0), alpha(0), pRed(0);
  Real rtr(0), tnorm(0), rnorm0(0), sMs(0), pMp(0), sMp(0);
  iter = 0; iflag = 0;
  // Initialize step
  w.zero();
  // Compute residual
  t.set(g); t.scale(-one);
  // Preconditioned residual
  applyFreePrecond(r,t,x,g0,eps,model,bnd,ftol,pwa,dwa);
  //rho    = r.dot(t.dual());
  rho    = r.apply(t);
  rnorm0 = std::sqrt(rho);
  if ( rnorm0 == zero ) {
    return zero;
  }
  // Initialize direction
  p.set(r);
  pMp = rho;
  // Iterate CG
  for (iter = 0; iter < itermax; ++iter) {
    // Apply Hessian to direction dir
    applyFreeHessian(q,p,x,g0,eps,model,bnd,ftol,pwa,dwa);
    // Compute sigma such that ||s+sigma*dir|| = del
    //kappa = p.dot(q.dual());
    kappa = p.apply(q);
    alpha = (kappa>zero) ? rho/kappa : zero;
    sigma = trqsol(sMs,pMp,sMp,del);
    // Check for negative curvature or if iterate exceeds trust region
    if (kappa <= zero || alpha >= sigma) {
      w.axpy(sigma,p);
      sMs = del*del;
      iflag = (kappa<=zero) ? 2 : 3;
      break;
    }
    pRed += half*alpha*rho;
    // Update iterate and residuals
    w.axpy(alpha,p);
    t.axpy(-alpha,q);
    applyFreePrecond(r,t,x,g0,eps,model,bnd,ftol,pwa,dwa);
    // Exit if residual tolerance is met
    //rtr   = r.dot(t.dual());
    rtr   = r.apply(t);
    tnorm = t.norm();
    if (rtr <= tol*tol || tnorm <= tol) {
      sMs   = sMs + two*alpha*sMp + alpha*alpha*pMp;
      iflag = 0;
      break;
    }
    // Compute p = r + beta * p
    beta = rtr/rho;
    p.scale(beta); p.plus(r);
    rho  = rtr;
    // Update dot products
    //   sMs = <s, inv(M)s> or <s, s> if equality constraint present
    //   sMp = <s, inv(M)p> or <s, p> if equality constraint present
    //   pMp = <p, inv(M)p> or <p, p> if equality constraint present
    sMs = sMs + two*alpha*sMp + alpha*alpha*pMp;
    sMp = beta*(sMp + alpha*pMp);
    pMp = rho + beta*beta*pMp;
  }
  // Check iteration count
  if (iter == itermax) {
    iflag = 1;
  }
  if (iflag > 1) {
    pRed += sigma*(rho-half*sigma*kappa);
  }
  if (iflag != 1) { 
    iter++;
  }
  return pRed;
}

template<typename Real>
void KelleySachsAlgorithm<Real>::applyFreeHessian(Vector<Real> &hv,
                                                  const Vector<Real> &v,
                                                  const Vector<Real> &x,
                                                  const Vector<Real> &g,
                                                  Real eps,
                                                  TrustRegionModel_U<Real> &model,
                                                  BoundConstraint<Real> &bnd,
                                                  Real &tol,
                                                  Vector<Real> &pwa,
                                                  Vector<Real> &dwa) const {
  const Real xeps(ROL_EPSILON<Real>());
  pwa.set(v);
  bnd.pruneActive(pwa,g.dual(),x,xeps,eps);
  model.hessVec(hv,pwa,x,tol); nhess_++;
  pwa.set(hv.dual());
  bnd.pruneActive(pwa,g.dual(),x,xeps,eps);
  hv.set(pwa.dual());
  pwa.set(v);
  bnd.pruneInactive(pwa,g.dual(),x,xeps,eps);
  hv.plus(pwa.dual());
  //pwa.set(v);
  //bnd.pruneActive(pwa,g,x,xeps,eps);
  //model.hessVec(hv,pwa,x,tol); nhess_++;
  //bnd.pruneActive(hv,g,x,xeps,eps);
  //pwa.set(v);
  //bnd.pruneInactive(pwa,g,x,xeps,eps);
  //dwa.set(pwa.dual());
  //bnd.pruneInactive(dwa,g,x,xeps,eps);
  //hv.plus(dwa);
}

template<typename Real>
void KelleySachsAlgorithm<Real>::applyFreePrecond(Vector<Real> &hv,
                                                  const Vector<Real> &v,
                                                  const Vector<Real> &x,
                                                  const Vector<Real> &g,
                                                  Real eps,
                                                  TrustRegionModel_U<Real> &model,
                                                  BoundConstraint<Real> &bnd,
                                                  Real &tol,
                                                  Vector<Real> &pwa,
                                                  Vector<Real> &dwa) const {
  const Real xeps(ROL_EPSILON<Real>());
  pwa.set(v.dual());
  bnd.pruneActive(pwa,g.dual(),x,xeps,eps);
  dwa.set(pwa.dual());
  model.precond(hv,dwa,x,tol);
  bnd.pruneActive(hv,g.dual(),x,xeps,eps);
  pwa.set(v.dual());
  bnd.pruneInactive(pwa,g.dual(),x,xeps,eps);
  hv.plus(pwa);
  //dwa.set(v);
  //bnd.pruneActive(dwa,g,x,xeps,eps);
  //model.precond(hv,dwa,x,tol);
  //bnd.pruneActive(hv,g,x,xeps,eps);
  //dwa.set(v);
  //bnd.pruneInactive(dwa,g,x,xeps,eps);
  //pwa.set(dwa.dual());
  //bnd.pruneInactive(pwa,g,x,xeps,eps);
  //hv.plus(pwa);
}

template<typename Real>
void KelleySachsAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(114,'-') << std::endl;
    os << " Kelley-Sachs trust-region method status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)" << std::endl;
    os << "  value   - Objective function value" << std::endl; 
    os << "  gnorm   - Norm of the gradient" << std::endl;
    os << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    os << "  delta   - Trust-Region radius" << std::endl;
    os << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    os << "  #grad   - Number of times the gradient was computed" << std::endl;
    os << "  #hess   - Number of times the Hessian was applied" << std::endl;
    os << std::endl;
    os << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      os << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    os << std::endl;
    os << "  iterCG - Number of Truncated CG iterations" << std::endl << std::endl;
    os << "  flagGC - Trust-Region Truncated CG flag" << std::endl;
    for( int flag = CG_FLAG_SUCCESS; flag != CG_FLAG_UNDEFINED; ++flag ) {
      os << "    " << NumberToString(flag) << " - "
           << ECGFlagToString(static_cast<ECGFlag>(flag)) << std::endl;
    }            
    os << std::string(114,'-') << std::endl;
  }
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(15) << std::left << "delta";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "#hess";
  os << std::setw(10) << std::left << "tr_flag";
  os << std::setw(10) << std::left << "iterCG";
  os << std::setw(10) << std::left << "flagCG";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void KelleySachsAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Kelley-Sachs Trust-Region Method (Type B, Bound Constraints)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void KelleySachsAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(15) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << state_->nfval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << nhess_;
    os << std::setw(10) << std::left << "---";
    os << std::setw(10) << std::left << "---";
    os << std::setw(10) << std::left << "---";
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
    os << std::setw(10) << std::left << nhess_;
    os << std::setw(10) << std::left << TRflag_;
    os << std::setw(10) << std::left << SPiter_;
    os << std::setw(10) << std::left << SPflag_;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeB
} // namespace ROL

#endif

// ORIGINAL KELLEY-SACHS ALGORITHM
//template<typename Real>
//std::vector<std::string> KelleySachsAlgorithm<Real>::run(Vector<Real>          &x,
//                                                         const Vector<Real>    &g, 
//                                                         Objective<Real>       &obj,
//                                                         BoundConstraint<Real> &bnd,
//                                                         std::ostream          &outStream ) {
//  const Real one(1);
//  Real ftol = std::sqrt(ROL_EPSILON<Real>()), tol1(1e-2);
//  Real gfnorm(0), gfnormf(0), tol(0), stol(0), eps(0);
//  Real ftrial(0), fnew(0), pRed(0), rho(1), alpha(1), lambda(1);
//  int rflag(0), cnt(0);
//  // Initialize trust-region data
//  initialize(x,g,obj,bnd,outStream);
//  Ptr<Vector<Real>> s = x.clone(), gfree = g.clone();
//  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone(), pwa3 = x.clone();
//  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone(), dwa3 = g.clone();
//
//  // Output
//  if (verbosity_ > 0) writeOutput(outStream,true);
//
//  // Initial tolerances
//  rflag = 0;
//  eps   = std::min(eps0_,std::sqrt(state_->gnorm));
//  tol   = tol1*std::min(tol0_,std::sqrt(state_->gnorm));
//
//  // Compute initial free gradient
//  gfree->set(*state_->gradientVec);
//  bnd.pruneActive(*gfree,*state_->gradientVec,x,eps);
//  gfnorm = gfree->norm();
//  stol = tol*gfnorm;
//
//  while (status_->check(*state_)) {
//    // Build trust-region model
//    model_->setData(obj,*state_->iterateVec,*state_->gradientVec);
//
//    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
//    pRed = trpcg(*s,SPflag_,SPiter_,*gfree,x,*state_->gradientVec,
//                 state_->searchSize,*model_,bnd,eps,stol,maxit_,
//                 *pwa1,*dwa1,*pwa2,*dwa2,*pwa3,*dwa3,outStream);
//    x.plus(*s);
//    bnd.project(x);
//
//    // Compute trial objective value
//    obj.update(x,false);
//    ftrial = obj.value(x,ftol); state_->nfval++;
//
//    // Compute ratio of acutal and predicted reduction
//    TRflag_ = TRUtils::SUCCESS;
//    TRUtils::analyzeRatio<Real>(rho,TRflag_,state_->value,ftrial,pRed,eps_,outStream,verbosity_>1);
//
//    // Check sufficient decrease
//    if ( rho >= eta0_ ) {
//      lambda = std::min(one, state_->searchSize/gfnorm);
//      // Compute Scaled Measure || x - P( x - lam * PI(g) ) ||
//      pwa1->set(*state_->iterateVec);
//      pwa1->axpy(-lambda,gfree->dual());
//      bnd.project(*pwa1);
//      pwa1->axpy(-one,*state_->iterateVec);
//      gfnormf = pwa1->norm();
//      // Sufficient decrease?
//      if (state_->value-ftrial < mu0_*gfnormf*state_->gnorm) {
//        TRflag_ = TRUtils::QMINSUFDEC;
//      }
//      if ( verbosity_ > 1 ) {
//        outStream << "    Norm of projected free gradient:         " << gfnormf                    << std::endl;
//        outStream << "    Decrease lower bound (constraints):      " << mu0_*gfnormf*state_->gnorm << std::endl;
//        outStream << "    Trust-region flag (constraints):         " << TRflag_                    << std::endl;
//        outStream << std::endl;
//      }
//    }
//
//    // Update algorithm state
//    state_->iter++;
//    // Accept/reject step and update trust region radius
//    if ((rho < eta0_ && TRflag_ == TRUtils::SUCCESS) || (TRflag_ >= 2)) {
//      state_->stepVec->set(x);
//      state_->stepVec->axpy(-one,*state_->iterateVec);
//      state_->snorm = state_->stepVec->norm();
//      rflag = 1;
//      x.set(*state_->iterateVec);
//      obj.update(x,false,state_->iter);
//      // Decrease trust-region radius
//      state_->searchSize = gamma1_*state_->searchSize;
//    }
//    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
//             || (TRflag_ == TRUtils::POSPREDNEG)) {
//      if (rho >= eta2_ && rflag == 0 && state_->searchSize != delMax_) {
//        state_->stepVec->set(x);
//        state_->stepVec->axpy(-one,*state_->iterateVec);
//        state_->snorm = state_->stepVec->norm();
//        x.set(*state_->iterateVec);
//        obj.update(x,false,state_->iter);
//        // Increase trust-region radius
//        state_->searchSize = std::min(delMax_,gamma2_*state_->searchSize);
//      }
//      else {
//        if (rho >= eta0_ && rho < eta1_) {
//          // Decrease trust-region radius
//          state_->searchSize = gamma1_*state_->searchSize;
//        }
//        // Projected search
//        cnt = 0;
//        alpha = one;
//        obj.gradient(*dwa1,x,ftol); state_->ngrad++;
//        pwa2->set(dwa1->dual());
//        pwa1->set(x); pwa1->axpy(-alpha/alpha0_,*pwa2);
//        bnd.project(*pwa1);
//        obj.update(*pwa1,false);
//        fnew = obj.value(*pwa1,ftol); state_->nfval++;
//        while ((fnew-ftrial) >= mu1_*(state_->value-ftrial)) {
//          alpha *= beta_;
//          pwa1->set(x); pwa1->axpy(-alpha/alpha0_,*pwa2);
//          bnd.project(*pwa1);
//          obj.update(*pwa1,false);
//          fnew = obj.value(*pwa1,ftol); state_->nfval++;
//          if ( cnt >= minit_ ) break;
//          cnt++;
//        }
//        state_->stepVec->set(*pwa1);
//        state_->stepVec->axpy(-one,*state_->iterateVec);
//        state_->snorm = state_->stepVec->norm();
//        // Store objective function and iteration information
//        if (std::isnan(fnew)) {
//          TRflag_ = TRUtils::TRNAN;
//          rho     = -one;
//          x.set(*state_->iterateVec);
//          obj.update(x,false,state_->iter);
//          // Decrease trust-region radius
//          state_->searchSize = gamma1_*state_->searchSize;
//        }
//        else {
//          // Update objective function information
//          x.set(*pwa1);
//          state_->iterateVec->set(x);
//          state_->value = fnew;
//          dwa1->set(*state_->gradientVec);
//          obj.update(x,true,state_->iter);
//          obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
//          // Compute free gradient
//          gfree->set(*state_->gradientVec);
//          bnd.pruneActive(*gfree,*state_->gradientVec,x,eps);
//          gfnorm = gfree->norm();
//          stol = tol*gfnorm;
//          // Compute criticality measure
//          pwa1->set(x);
//          pwa1->axpy(-one,state_->gradientVec->dual());
//          bnd.project(*pwa1);
//          pwa1->axpy(-one,x);
//          state_->gnorm = pwa1->norm();
//          // Update secant information in trust-region model
//          model_->update(x,*state_->stepVec,*dwa1,*state_->gradientVec,
//                         state_->snorm,state_->iter);
//          // Update tolerances
//          rflag = 0;
//          eps   = std::min(eps0_,std::sqrt(state_->gnorm));
//          tol   = tol1*std::min(tol0_,std::sqrt(state_->gnorm));
//        }
//      }
//    }
//
//    // Update Output
//    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
//  }
//  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
//}
