// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_COLEMANLIALGORITHM_DEF_HPP
#define ROL_TYPEB_COLEMANLIALGORITHM_DEF_HPP

namespace ROL {
namespace TypeB {

template<typename Real>
ColemanLiAlgorithm<Real>::ColemanLiAlgorithm(ParameterList &list,
                                         const Ptr<Secant<Real>> &secant) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  ParameterList &trlist = list.sublist("Step").sublist("Trust Region");
  // Trust-Region Parameters
  state_->searchSize = trlist.get("Initial Radius",            -1.0);
  delMax_    = trlist.get("Maximum Radius",                       ROL_INF<Real>());
  eta0_      = trlist.get("Step Acceptance Threshold",            0.05);
  eta1_      = trlist.get("Radius Shrinking Threshold",           0.05);
  eta2_      = trlist.get("Radius Growing Threshold",             0.9);
  gamma0_    = trlist.get("Radius Shrinking Rate (Negative rho)", 0.0625);
  gamma1_    = trlist.get("Radius Shrinking Rate (Positive rho)", 0.25);
  gamma2_    = trlist.get("Radius Growing Rate",                  2.5);
  TRsafe_    = trlist.get("Safeguard Size",                       100.0);
  eps_       = TRsafe_*ROL_EPSILON<Real>();
  interpRad_ = trlist.get("Use Radius Interpolation",             false);
  // Nonmonotone Parameters
  storageNM_ = trlist.get("Nonmonotone Storage Size",             0);
  useNM_     = (storageNM_ <= 0 ? false : true);
  // Krylov Parameters
  maxit_ = list.sublist("General").sublist("Krylov").get("Iteration Limit",    20);
  tol1_  = list.sublist("General").sublist("Krylov").get("Absolute Tolerance", 1e-4);
  tol2_  = list.sublist("General").sublist("Krylov").get("Relative Tolerance", 1e-2);
  // Algorithm-Specific Parameters
  ROL::ParameterList &lmlist = trlist.sublist("Coleman-Li");
  mu0_       = lmlist.get("Sufficient Decrease Parameter",                             1e-2);
  spexp_     = lmlist.get("Relative Tolerance Exponent",                               1.0);
  spexp_     = std::max(static_cast<Real>(1),std::min(spexp_,static_cast<Real>(2)));
  alphaMax_  = lmlist.get("Relaxation Safeguard",                                      0.999);
  alphaMax_  = (alphaMax_ >= static_cast<Real>(1) ? static_cast<Real>(0.5) : alphaMax_);
  // Output Parameters
  verbosity_   = list.sublist("General").get("Output Level",0);
  writeHeader_ = verbosity_ > 2;
  // Secant Information
  useSecantPrecond_ = list.sublist("General").sublist("Secant").get("Use as Preconditioner", false);
  useSecantHessVec_ = list.sublist("General").sublist("Secant").get("Use as Hessian",        false);
  ESecantMode mode = SECANTMODE_BOTH;
  if (useSecantPrecond_ && !useSecantHessVec_)      mode = SECANTMODE_INVERSE;
  else if (useSecantHessVec_ && !useSecantPrecond_) mode = SECANTMODE_FORWARD;
  // Initialize trust region model
  model_ = makePtr<TrustRegionModel_U<Real>>(list,secant,mode);
  if (secant == nullPtr) {
    esec_ = StringToESecant(list.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS"));
  }
}

template<typename Real>
void ColemanLiAlgorithm<Real>::initialize(Vector<Real>          &x,
                                          const Vector<Real>    &g,
                                          Objective<Real>       &obj,
                                          BoundConstraint<Real> &bnd,
                                          std::ostream &outStream) {
  const Real one(1);
  hasEcon_ = true;
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasEcon_ = false;
  }
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  nhess_ = 0;
  // Update approximate gradient and approximate objective function.
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>();
  proj_->getBoundConstraint()->projectInterior(x); state_->nproj++;
  state_->iterateVec->set(x);
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,ftol);
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream); state_->nproj++;
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  // Compute initial trust region radius if desired.
  if ( state_->searchSize <= static_cast<Real>(0) ) {
    state_->searchSize = state_->gradientVec->norm();
  }
  // Initialize null space projection
  if (hasEcon_) {
    rcon_ = makePtr<ReducedLinearConstraint<Real>>(proj_->getLinearConstraint(),
                                                   makePtrFromRef(bnd),
                                                   makePtrFromRef(x));
    ns_   = makePtr<NullSpaceOperator<Real>>(rcon_,x,
                                             *proj_->getResidual());
  }
}

template<typename Real>
void ColemanLiAlgorithm<Real>::run(Vector<Real>          &x,
                                   const Vector<Real>    &g,
                                   Objective<Real>       &obj,
                                   BoundConstraint<Real> &bnd,
                                   std::ostream          &outStream ) {
  const Real zero(0), one(1), half(0.5);
  Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  Real tol(0), stol(0), snorm(0);
  Real ftrial(0), pRed(0), rho(1), alpha(1);
  // Initialize trust-region data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone(), pwa3 = x.clone();
  Ptr<Vector<Real>> pwa4 = x.clone(), pwa5 = x.clone();
  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone(), dwa3 = g.clone();
  // Initialize nonmonotone data
  Real rhoNM(0), sigmac(0), sigmar(0), sBs(0), gs(0);
  Real fr(state_->value), fc(state_->value), fmin(state_->value);
  TRUtils::ETRFlag TRflagNM;
  int L(0);

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Build trust-region model (use only to encapsulate Hessian/secant)
    model_->setData(obj,*state_->iterateVec,*state_->gradientVec,tol0);

    // Run Truncated CG
    // TODO: Model is: 1/2 (x-xk)' (B + Einv(D)) + g'(x-xk)
    //       applyHessian returns (B+Einv(D))v
    SPflag_ = 0;
    SPiter_ = 0;
    tol     = std::min(tol1_,tol2_*std::pow(state_->gnorm,spexp_));
    stol    = tol; //zero;
    pwa5->set(state_->gradientVec->dual());
    snorm   = dtrpcg(*state_->stepVec,SPflag_,SPiter_,*state_->gradientVec,x,*pwa5,
                     state_->searchSize,*model_,bnd,tol,stol,
                     *pwa1,*dwa1,*pwa2,*dwa2,*pwa3,*pwa4,*dwa3,outStream);
    if (verbosity_ > 1) {
      outStream << "  Computation of CG step"               << std::endl;
      outStream << "    CG step length:                   " << snorm   << std::endl;
      outStream << "    Number of CG iterations:          " << SPiter_ << std::endl;
      outStream << "    CG flag:                          " << SPflag_ << std::endl;
      outStream << std::endl;
    }

    // Relax CG step so that it is interior
    snorm = dgpstep(*pwa1,*state_->stepVec,x,one,outStream);
    alpha = std::max(alphaMax_, one-snorm);
    pwa1->scale(alpha);
    state_->stepVec->set(*pwa1);
    state_->snorm = alpha * snorm;
    x.plus(*state_->stepVec);

    // Compute predicted reduction
    model_->hessVec(*dwa1,*pwa1,x,tol); nhess_++;
    gs   = state_->gradientVec->apply(*state_->stepVec);
    sBs  = dwa1->apply(*state_->stepVec);
    pRed = - half * sBs - gs;

    // Compute trial objective value
    obj.update(x,UpdateType::Trial);
    ftrial = obj.value(x,tol0);
    state_->nfval++;

    // Compute ratio of acutal and predicted reduction
    TRflag_ = TRUtils::SUCCESS;
    TRUtils::analyzeRatio<Real>(rho,TRflag_,state_->value,ftrial,pRed,eps_,outStream,verbosity_>1);
    if (useNM_) {
      TRUtils::analyzeRatio<Real>(rhoNM,TRflagNM,fr,ftrial,pRed+sigmar,eps_,outStream,verbosity_>1);
      TRflag_ = (rho < rhoNM ? TRflagNM : TRflag_);
      rho     = (rho < rhoNM ?    rhoNM :    rho );
    }

    // Update algorithm state
    state_->iter++;
    // Accept/reject step and update trust region radius
    if ((rho < eta0_ && TRflag_ == TRUtils::SUCCESS) || (TRflag_ >= 2)) { // Step Rejected
      x.set(*state_->iterateVec);
      obj.update(x,UpdateType::Revert,state_->iter);
      if (interpRad_ && (rho < zero && TRflag_ != TRUtils::TRNAN)) {
        // Negative reduction, interpolate to find new trust-region radius
        state_->searchSize = TRUtils::interpolateRadius<Real>(*state_->gradientVec,*state_->stepVec,
          state_->snorm,pRed,state_->value,ftrial,state_->searchSize,gamma0_,gamma1_,eta2_,
          outStream,verbosity_>1);
      }
      else { // Shrink trust-region radius
        state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
      }
    }
    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
             || (TRflag_ == TRUtils::POSPREDNEG)) { // Step Accepted
      state_->value = ftrial;
      obj.update(x,UpdateType::Accept,state_->iter);
      if (useNM_) {
        sigmac += pRed; sigmar += pRed;
        if (ftrial < fmin) { fmin = ftrial; fc = fmin; sigmac = zero; L = 0; }
        else {
          L++;
          if (ftrial > fc)     { fc = ftrial; sigmac = zero;   }
          if (L == storageNM_) { fr = fc;     sigmar = sigmac; }
        }
      }
      // Increase trust-region radius
      if (rho >= eta2_) state_->searchSize = std::min(gamma2_*state_->searchSize, delMax_);
      // Compute gradient at new iterate
      dwa1->set(*state_->gradientVec);
      obj.gradient(*state_->gradientVec,x,tol0);
      state_->ngrad++;
      state_->gnorm = TypeB::Algorithm<Real>::optimalityCriterion(x,*state_->gradientVec,*pwa1,outStream);
      state_->iterateVec->set(x);
      // Update secant information in trust-region model
      model_->update(x,*state_->stepVec,*dwa1,*state_->gradientVec,
                     state_->snorm,state_->iter);
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
Real ColemanLiAlgorithm<Real>::dgpstep(Vector<Real> &s, const Vector<Real> &w,
                                 const Vector<Real> &x, const Real alpha,
                                 std::ostream &outStream) const {
  s.set(x); s.axpy(alpha,w);
  proj_->getBoundConstraint()->projectInterior(s); state_->nproj++;
  s.axpy(static_cast<Real>(-1),x);
  return s.norm();
}

template<typename Real>
Real ColemanLiAlgorithm<Real>::dtrqsol(const Real xtx,
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
Real ColemanLiAlgorithm<Real>::dtrpcg(Vector<Real> &w, int &iflag, int &iter,
                                const Vector<Real> &g, const Vector<Real> &x,
                                const Vector<Real> &gdual,
                                const Real del, TrustRegionModel_U<Real> &model,
                                BoundConstraint<Real> &bnd,
                                const Real tol, const Real stol,
                                Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
                                Vector<Real> &t, Vector<Real> &pwa1,
                                Vector<Real> &pwa2, Vector<Real> &dwa,
                                std::ostream &outStream) const {
  // p = step (primal)
  // q = hessian applied to step p (dual)
  // t = gradient (dual)
  // r = preconditioned gradient (primal)
  Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  const Real zero(0), one(1), two(2);
  Real rho(0), kappa(0), beta(0), sigma(0), alpha(0);
  Real rtr(0), tnorm(0), sMs(0), pMp(0), sMp(0);
  iter = 0; iflag = 0;
  // Initialize step
  w.zero();
  // Compute residual
  t.set(g); t.scale(-one);
  // Preconditioned residual
  applyPrecond(r,t,x,gdual,model,bnd,tol0,dwa,pwa1);
  //rho = r.dot(t.dual());
  rho = r.apply(t);
  // Initialize direction
  p.set(r);
  pMp = (!hasEcon_ ? rho : p.dot(p)); // If no equality constraint, used preconditioned norm
  // Iterate CG
  for (iter = 0; iter < maxit_; ++iter) {
    // Apply Hessian to direction dir
    applyHessian(q,p,x,gdual,model,bnd,tol0,pwa1,pwa2);
    // Compute sigma such that ||s+sigma*dir|| = del
    //kappa = p.dot(q.dual());
    kappa = p.apply(q);
    alpha = (kappa>zero) ? rho/kappa : zero;
    sigma = dtrqsol(sMs,pMp,sMp,del);
    // Check for negative curvature or if iterate exceeds trust region
    if (kappa <= zero || alpha >= sigma) {
      w.axpy(sigma,p);
      sMs = del*del;
      iflag = (kappa<=zero) ? 2 : 3;
      break;
    }
    // Update iterate and residuals
    w.axpy(alpha,p);
    t.axpy(-alpha,q);
    applyPrecond(r,t,x,gdual,model,bnd,tol0,dwa,pwa1);
    // Exit if residual tolerance is met
    //rtr   = r.dot(t.dual());
    rtr   = r.apply(t);
    tnorm = t.norm();
    if (rtr <= stol*stol || tnorm <= tol) {
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
    pMp = (!hasEcon_ ? rho : p.dot(p)) + beta*beta*pMp;
  }
  // Check iteration count
  if (iter == maxit_) {
    iflag = 1;
  }
  if (iflag != 1) {
    iter++;
  }
  return std::sqrt(sMs); // w.norm();
}

template<typename Real>
void ColemanLiAlgorithm<Real>::applyC(Vector<Real> &Cv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &x,
                                      const Vector<Real> &g,
                                      BoundConstraint<Real> &bnd,
                                      Vector<Real> &pwa) const {
  bnd.applyInverseScalingFunction(pwa,v,x,g);
  bnd.applyScalingFunctionJacobian(Cv,pwa,x,g);
}

template<typename Real>
void ColemanLiAlgorithm<Real>::applyHessian(Vector<Real> &hv,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x,
                                            const Vector<Real> &g,
                                            TrustRegionModel_U<Real> &model,
                                            BoundConstraint<Real> &bnd,
                                            Real &tol,
                                            Vector<Real> &pwa1,
                                            Vector<Real> &pwa2) const {
  model.hessVec(hv,v,x,tol); nhess_++;
  applyC(pwa2,v,x,g,bnd,pwa1);
  hv.plus(pwa2.dual());
}

template<typename Real>
void ColemanLiAlgorithm<Real>::applyPrecond(Vector<Real> &hv,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x,
                                            const Vector<Real> &g,
                                            TrustRegionModel_U<Real> &model,
                                            BoundConstraint<Real> &bnd,
                                            Real &tol,
                                            Vector<Real> &dwa,
                                            Vector<Real> &pwa) const {
  model.precond(hv,v,x,tol);
}

template<typename Real>
void ColemanLiAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(114,'-') << std::endl;
    os << " Coleman-Li affine-scaling trust-region method status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)" << std::endl;
    os << "  value   - Objective function value" << std::endl;
    os << "  gnorm   - Norm of the gradient" << std::endl;
    os << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    os << "  delta   - Trust-Region radius" << std::endl;
    os << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    os << "  #grad   - Number of times the gradient was computed" << std::endl;
    os << "  #hess   - Number of times the Hessian was applied" << std::endl;
    os << "  #proj   - Number of times the projection was applied" << std::endl;
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
  os << std::setw(10) << std::left << "#proj";
  os << std::setw(10) << std::left << "tr_flag";
  os << std::setw(10) << std::left << "iterCG";
  os << std::setw(10) << std::left << "flagCG";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void ColemanLiAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Coleman-Li Affine-Scaling Trust-Region Method (Type B, Bound Constraints)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void ColemanLiAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(10) << std::left << state_->nproj;
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
    os << std::setw(10) << std::left << state_->nproj;
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
