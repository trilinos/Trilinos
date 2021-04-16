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

#ifndef ROL_TYPEB_LINMOREALGORITHM_DEF_HPP
#define ROL_TYPEB_LINMOREALGORITHM_DEF_HPP

namespace ROL {
namespace TypeB {

template<typename Real>
LinMoreAlgorithm<Real>::LinMoreAlgorithm(ParameterList &list,
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
  // Krylov Parameters
  maxit_ = list.sublist("General").sublist("Krylov").get("Iteration Limit",    20);
  tol1_  = list.sublist("General").sublist("Krylov").get("Absolute Tolerance", 1e-4);
  tol2_  = list.sublist("General").sublist("Krylov").get("Relative Tolerance", 1e-2);
  // Algorithm-Specific Parameters
  ROL::ParameterList &lmlist = trlist.sublist("Lin-More");
  minit_     = lmlist.get("Maximum Number of Minor Iterations",                        10);
  mu0_       = lmlist.get("Sufficient Decrease Parameter",                             1e-2);
  spexp_     = lmlist.get("Relative Tolerance Exponent",                               1.0);
  spexp_     = std::max(static_cast<Real>(1),std::min(spexp_,static_cast<Real>(2)));
  redlim_    = lmlist.sublist("Cauchy Point").get("Maximum Number of Reduction Steps", 10);
  explim_    = lmlist.sublist("Cauchy Point").get("Maximum Number of Expansion Steps", 10);
  alpha_     = lmlist.sublist("Cauchy Point").get("Initial Step Size",                 1.0);
  normAlpha_ = lmlist.sublist("Cauchy Point").get("Normalize Initial Step Size",       false); 
  interpf_   = lmlist.sublist("Cauchy Point").get("Reduction Rate",                    0.1);
  extrapf_   = lmlist.sublist("Cauchy Point").get("Expansion Rate",                    10.0);
  qtol_      = lmlist.sublist("Cauchy Point").get("Decrease Tolerance",                1e-8);
  interpfPS_ = lmlist.sublist("Projected Search").get("Backtracking Rate",             0.5);
  pslim_     = lmlist.sublist("Projected Search").get("Maximum Number of Steps",       20);
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
void LinMoreAlgorithm<Real>::initialize(Vector<Real>          &x,
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
  proj_->project(x,outStream); state_->nproj++;
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
  // Normalize initial CP step length
  if (normAlpha_) {
    alpha_ /= state_->gradientVec->norm();
  }
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
void LinMoreAlgorithm<Real>::run(Vector<Real>          &x,
                                 const Vector<Real>    &g, 
                                 Objective<Real>       &obj,
                                 BoundConstraint<Real> &bnd,
                                 std::ostream          &outStream ) {
  const Real zero(0);
  Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  Real gfnorm(0), gfnormf(0), tol(0), stol(0), snorm(0);
  Real ftrial(0), pRed(0), rho(1), q(0);
  int flagCG(0), iterCG(0), maxit(0);
  // Initialize trust-region data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone();
  Ptr<Vector<Real>> gmod = g.clone(), gfree = g.clone();
  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone(), pwa3 = x.clone();
  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone(), dwa3 = g.clone();

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Build trust-region model
    model_->setData(obj,*state_->iterateVec,*state_->gradientVec);

    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
    // Compute Cauchy point (TRON notation: x = x[1])
    snorm = dcauchy(*state_->stepVec,alpha_,q,*state_->iterateVec,
                    state_->gradientVec->dual(),state_->searchSize,
                    *model_,*dwa1,*dwa2,outStream); // Solve 1D optimization problem for alpha
    x.plus(*state_->stepVec);                       // Set x = x[0] + alpha*g
    pRed = -q;

    // Model gradient at s = x[1] - x[0]
    gmod->set(*dwa1); // hessVec from Cauchy point computation
    gmod->plus(*state_->gradientVec);
    gfree->set(*gmod);
    //bnd.pruneActive(*gfree,x,zero);
    pwa1->set(gfree->dual());
    bnd.pruneActive(*pwa1,x,zero);
    gfree->set(pwa1->dual());
    if (hasEcon_) {
      applyFreePrecond(*pwa1,*gfree,x,*model_,bnd,tol0,*dwa1,*pwa2);
      gfnorm = pwa1->norm();
    }
    else {
      gfnorm = gfree->norm();
    }
    SPiter_ = 0; SPflag_ = 0;
    if (verbosity_ > 1) {
      outStream << "    Norm of free gradient components: " << gfnorm << std::endl;
      outStream << std::endl;
    }

    // Trust-region subproblem solve loop
    maxit = maxit_;
    for (int i = 0; i < minit_; ++i) {
      // Run Truncated CG
      flagCG = 0; iterCG = 0;
      gfnormf = zero;
      tol     = std::min(tol1_,tol2_*std::pow(gfnorm,spexp_));
      stol    = tol; //zero;
      if (gfnorm > zero) {
        snorm = dtrpcg(*s,flagCG,iterCG,*gfree,x,
                       state_->searchSize,*model_,bnd,tol,stol,maxit,
                       *pwa1,*dwa1,*pwa2,*dwa2,*pwa3,*dwa3,outStream);
        maxit   -= iterCG;
        SPiter_ += iterCG;
        if (verbosity_ > 1) {
          outStream << "  Computation of CG step"               << std::endl;
          outStream << "    Current face (i):                 " << i       << std::endl;
          outStream << "    CG step length:                   " << snorm   << std::endl;
          outStream << "    Number of CG iterations:          " << iterCG  << std::endl;
          outStream << "    CG flag:                          " << flagCG  << std::endl;
          outStream << "    Total number of iterations:       " << SPiter_ << std::endl;
          outStream << std::endl;
        }

        // Projected search
        snorm = dprsrch(x,*s,q,gmod->dual(),*model_,bnd,*pwa1,*dwa1,outStream);
        pRed += -q;

        // Model gradient at s = (x[i+1]-x[i]) - (x[i]-x[0])
        state_->stepVec->plus(*s);
        gmod->plus(*dwa1); // gmod = H(x[i+1]-x[i]) + H(x[i]-x[0]) + g
        gfree->set(*gmod);
        //bnd.pruneActive(*gfree,x,zero);
        pwa1->set(gfree->dual());
        bnd.pruneActive(*pwa1,x,zero);
        gfree->set(pwa1->dual());
        if (hasEcon_) {
          applyFreePrecond(*pwa1,*gfree,x,*model_,bnd,tol0,*dwa1,*pwa2);
          gfnormf = pwa1->norm();
        }
        else {
          gfnormf = gfree->norm();
        }
        if (verbosity_ > 1) {
          outStream << "    Norm of free gradient components: " << gfnormf       << std::endl;
          outStream << std::endl;
        }
      }

      // Termination check
      if (gfnormf <= tol) {
        SPflag_ = 0;
        break;
      }
      else if (SPiter_ >= maxit_) {
        SPflag_ = 1;
        break;
      }
      else if (flagCG == 2) {
        SPflag_ = 2;
        break;
      }
      else if (flagCG == 3) {
        SPflag_ = 3;
        break;
      }

      // Update free gradient norm
      gfnorm = gfnormf;
    }
    state_->snorm = state_->stepVec->norm();

    // Compute trial objective value
    obj.update(x,UpdateType::Trial);
    ftrial = obj.value(x,tol0);
    state_->nfval++;

    // Compute ratio of acutal and predicted reduction
    TRflag_ = TRUtils::SUCCESS;
    TRUtils::analyzeRatio<Real>(rho,TRflag_,state_->value,ftrial,pRed,eps_,outStream,verbosity_>1);

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
Real LinMoreAlgorithm<Real>::dgpstep(Vector<Real> &s, const Vector<Real> &w,
                                 const Vector<Real> &x, const Real alpha,
                                 std::ostream &outStream) const {
  s.set(x); s.axpy(alpha,w);
  proj_->project(s,outStream); state_->nproj++;
  s.axpy(static_cast<Real>(-1),x);
  return s.norm();
}

template<typename Real>
Real LinMoreAlgorithm<Real>::dcauchy(Vector<Real> &s,
                                     Real &alpha,
                                     Real &q,
                                     const Vector<Real> &x,
                                     const Vector<Real> &g,
                                     const Real del,
                                     TrustRegionModel_U<Real> &model,
                                     Vector<Real> &dwa, Vector<Real> &dwa1,
                                     std::ostream &outStream) {
  const Real half(0.5);
  // const Real zero(0); // Unused
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  bool interp = false;
  Real gs(0), snorm(0);
  // Compute s = P(x[0] - alpha g[0])
  snorm = dgpstep(s,g,x,-alpha,outStream);
  if (snorm > del) {
    interp = true;
  }
  else {
    model.hessVec(dwa,s,x,tol); nhess_++;
    gs = s.dot(g);
    //q  = half * s.dot(dwa.dual()) + gs;
    q  = half * s.apply(dwa) + gs;
    interp = (q > mu0_*gs);
  }
  // Either increase or decrease alpha to find approximate Cauchy point
  int cnt = 0;
  if (interp) {
    bool search = true;
    while (search) {
      alpha *= interpf_;
      snorm = dgpstep(s,g,x,-alpha,outStream);
      if (snorm <= del) {
        model.hessVec(dwa,s,x,tol); nhess_++;
        gs = s.dot(g);
        //q  = half * s.dot(dwa.dual()) + gs;
        q  = half * s.apply(dwa) + gs;
        search = (q > mu0_*gs) && (cnt < redlim_);
      }
      cnt++;
    }
  }
  else {
    bool search = true;
    Real alphas = alpha;
    Real qs = q;
    dwa1.set(dwa);
    while (search) {
      alpha *= extrapf_;
      snorm = dgpstep(s,g,x,-alpha,outStream);
      if (snorm <= del && cnt < explim_) {
        model.hessVec(dwa,s,x,tol); nhess_++;
        gs = s.dot(g);
        //q  = half * s.dot(dwa.dual()) + gs;
        q  = half * s.apply(dwa) + gs;
        if (q <= mu0_*gs && std::abs(q-qs) > qtol_*std::abs(qs)) {
          dwa1.set(dwa);
          search = true;
          alphas = alpha;
          qs     = q;
        }
        else {
          q = qs;
          dwa.set(dwa1);
          search = false;
        }
      }
      else {
        search = false;
      }
      cnt++;
    }
    alpha = alphas;
    snorm = dgpstep(s,g,x,-alpha,outStream);
  }
  if (verbosity_ > 1) {
    outStream << "  Cauchy point"                         << std::endl;
    outStream << "    Step length (alpha):              " << alpha << std::endl;
    outStream << "    Step length (alpha*g):            " << snorm << std::endl;
    outStream << "    Model decrease (pRed):            " << -q    << std::endl;
    if (!interp) {
      outStream << "    Number of extrapolation steps:    " << cnt << std::endl;
    }
  }
  return snorm;
}

template<typename Real>
Real LinMoreAlgorithm<Real>::dprsrch(Vector<Real> &x, Vector<Real> &s,
                                     Real &q, const Vector<Real> &g,
                                     TrustRegionModel_U<Real> &model,
                                     BoundConstraint<Real> &bnd,
                                     Vector<Real> &pwa, Vector<Real> &dwa,
                                     std::ostream &outStream) {
  const Real half(0.5);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  Real beta(1), snorm(0), gs(0);
  int nsteps = 0;
  // Reduce beta until sufficient decrease is satisfied
  bool search = true;
  while (search) {
    nsteps++;
    snorm = dgpstep(pwa,s,x,beta,outStream);
    model.hessVec(dwa,pwa,x,tol); nhess_++;
    gs = pwa.dot(g);
    //q  = half * pwa.dot(dwa.dual()) + gs;
    q  = half * pwa.apply(dwa) + gs;
    if (q <= mu0_*gs || nsteps > pslim_) {
      search = false;
    }
    else {
      beta *= interpfPS_;
    }
  }
  s.set(pwa);
  x.plus(s);
  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "  Projected search"                     << std::endl;
    outStream << "    Step length (beta):               " << beta   << std::endl;
    outStream << "    Step length (beta*s):             " << snorm  << std::endl;
    outStream << "    Model decrease (pRed):            " << -q     << std::endl;
    outStream << "    Number of steps:                  " << nsteps << std::endl;
  }
  return snorm;
}

template<typename Real>
Real LinMoreAlgorithm<Real>::dtrqsol(const Real xtx,
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
Real LinMoreAlgorithm<Real>::dtrpcg(Vector<Real> &w, int &iflag, int &iter,
                                    const Vector<Real> &g, const Vector<Real> &x,
                                    const Real del, TrustRegionModel_U<Real> &model,
                                    BoundConstraint<Real> &bnd,
                                    const Real tol, const Real stol, const int itermax,
                                    Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
                                    Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa,
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
  applyFreePrecond(r,t,x,model,bnd,tol0,dwa,pwa);
  //rho = r.dot(t.dual());
  rho = r.apply(t);
  // Initialize direction
  p.set(r);
  pMp = (!hasEcon_ ? rho : p.dot(p)); // If no equality constraint, used preconditioned norm
  // Iterate CG
  for (iter = 0; iter < itermax; ++iter) {
    // Apply Hessian to direction dir
    applyFreeHessian(q,p,x,model,bnd,tol0,pwa);
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
    applyFreePrecond(r,t,x,model,bnd,tol0,dwa,pwa);
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
  if (iter == itermax) {
    iflag = 1;
  }
  if (iflag != 1) { 
    iter++;
  }
  return std::sqrt(sMs); // w.norm();
}

template<typename Real>
void LinMoreAlgorithm<Real>::applyFreeHessian(Vector<Real> &hv,
                                              const Vector<Real> &v,
                                              const Vector<Real> &x,
                                              TrustRegionModel_U<Real> &model,
                                              BoundConstraint<Real> &bnd,
                                              Real &tol,
                                              Vector<Real> &pwa) const {
  const Real zero(0);
  pwa.set(v);
  bnd.pruneActive(pwa,x,zero);
  model.hessVec(hv,pwa,x,tol); nhess_++;
  pwa.set(hv.dual());
  bnd.pruneActive(pwa,x,zero);
  hv.set(pwa.dual());
  //pwa.set(v);
  //bnd.pruneActive(pwa,x,zero);
  //model.hessVec(hv,pwa,x,tol); nhess_++;
  //bnd.pruneActive(hv,x,zero);
}

template<typename Real>
void LinMoreAlgorithm<Real>::applyFreePrecond(Vector<Real> &hv,
                                              const Vector<Real> &v,
                                              const Vector<Real> &x,
                                              TrustRegionModel_U<Real> &model,
                                              BoundConstraint<Real> &bnd,
                                              Real &tol,
                                              Vector<Real> &dwa,
                                              Vector<Real> &pwa) const {
  if (!hasEcon_) {
    const Real zero(0);
    pwa.set(v.dual());
    bnd.pruneActive(pwa,x,zero);
    dwa.set(pwa.dual());
    model.precond(hv,dwa,x,tol);
    bnd.pruneActive(hv,x,zero);
    //dwa.set(v);
    //bnd.pruneActive(dwa,x,zero);
    //model.precond(hv,dwa,x,tol);
    //bnd.pruneActive(hv,x,zero);
  }
  else {
    // Perform null space projection
    rcon_->setX(makePtrFromRef(x));
    ns_->update(x);
    pwa.set(v.dual());
    ns_->apply(hv,pwa,tol);
  }
}

template<typename Real>
void LinMoreAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    hist << " Lin-More trust-region method status output definitions" << std::endl << std::endl;
    hist << "  iter    - Number of iterates (steps taken)" << std::endl;
    hist << "  value   - Objective function value" << std::endl; 
    hist << "  gnorm   - Norm of the gradient" << std::endl;
    hist << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  delta   - Trust-Region radius" << std::endl;
    hist << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    hist << "  #grad   - Number of times the gradient was computed" << std::endl;
    hist << "  #hess   - Number of times the Hessian was applied" << std::endl;
    hist << "  #proj   - Number of times the projection was applied" << std::endl;
    hist << std::endl;
    hist << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      hist << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    hist << std::endl;
    if (minit_ > 0) {
      hist << "  iterCG - Number of Truncated CG iterations" << std::endl << std::endl;
      hist << "  flagGC - Trust-Region Truncated CG flag" << std::endl;
      for( int flag = CG_FLAG_SUCCESS; flag != CG_FLAG_UNDEFINED; ++flag ) {
        hist << "    " << NumberToString(flag) << " - "
             << ECGFlagToString(static_cast<ECGFlag>(flag)) << std::endl;
      }
    }
    hist << std::string(114,'-') << std::endl;
  }
  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "delta";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "#hess";
  hist << std::setw(10) << std::left << "#proj";
  hist << std::setw(10) << std::left << "tr_flag";
  if (minit_ > 0) {
    hist << std::setw(10) << std::left << "iterCG";
    hist << std::setw(10) << std::left << "flagCG";
  }
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void LinMoreAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Lin-More Trust-Region Method (Type B, Bound Constraints)" << std::endl;
  os << hist.str();
}

template<typename Real>
void LinMoreAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    hist << std::setw(15) << std::left << state_->searchSize;
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << nhess_;
    hist << std::setw(10) << std::left << state_->nproj;
    hist << std::setw(10) << std::left << "---";
    if (minit_ > 0) {
      hist << std::setw(10) << std::left << "---";
      hist << std::setw(10) << std::left << "---";
    }
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
    hist << std::setw(10) << std::left << nhess_;
    hist << std::setw(10) << std::left << state_->nproj;
    hist << std::setw(10) << std::left << TRflag_;
    if (minit_ > 0) {
      hist << std::setw(10) << std::left << SPiter_;
      hist << std::setw(10) << std::left << SPflag_;
    }
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeB
} // namespace ROL

#endif
