// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_LSecantBALGORITHM_DEF_HPP
#define ROL_TYPEB_LSecantBALGORITHM_DEF_HPP

namespace ROL {
namespace TypeB {

template<typename Real>
LSecantBAlgorithm<Real>::LSecantBAlgorithm(ParameterList &list,
                                 const Ptr<Secant<Real>> &secant) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  ParameterList &trlist = list.sublist("Step").sublist("Line Search");
  state_->searchSize = static_cast<Real>(1);
  // Krylov Parameters
  maxit_ = list.sublist("General").sublist("Krylov").get("Iteration Limit",    20);
  tol1_  = list.sublist("General").sublist("Krylov").get("Absolute Tolerance", 1e-4);
  tol2_  = list.sublist("General").sublist("Krylov").get("Relative Tolerance", 1e-2);
  // Algorithm-Specific Parameters
  ROL::ParameterList &lmlist = trlist.sublist("Quasi-Newton").sublist("L-Secant-B");
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
  interpfPS_ = trlist.sublist("Line-Search Method").get("Backtracking Rate",           0.5);
  // Output Parameters
  verbosity_   = list.sublist("General").get("Output Level",0);
  writeHeader_ = verbosity_ > 2;
  // Secant Information
  useSecantPrecond_ = true;
  useSecantHessVec_ = true;
  ESecantMode mode = SECANTMODE_BOTH;
  if (secant == nullPtr) {
    esec_   = StringToESecant(list.sublist("General").sublist("Secant").get("Type","Limited-Memory Secant"));
    secant_ = SecantFactory<Real>(list,mode);
  }
}

template<typename Real>
void LSecantBAlgorithm<Real>::initialize(Vector<Real>          &x,
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
  // Update approximate gradient and approximate objective function.
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(); 
  proj_->project(x,outStream); state_->nproj++;
  state_->iterateVec->set(x);
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,ftol); state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
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
void LSecantBAlgorithm<Real>::run(Vector<Real>          &x,
                                  const Vector<Real>    &g, 
                                  Objective<Real>       &obj,
                                  BoundConstraint<Real> &bnd,
                                  std::ostream          &outStream ) {
  const Real zero(0), one(1);
  Real tol0 = ROL_EPSILON<Real>(), ftol = ROL_EPSILON<Real>();
  Real gfnorm(0), gs(0), ftrial(0), q(0), tol(0), stol(0), snorm(0);
  // Initialize trust-region data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone();
  Ptr<Vector<Real>> gmod = g.clone(), gfree = g.clone();
  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone(), pwa3 = x.clone();
  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone(), dwa3 = g.clone();

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    /**** SOLVE LINEARLY CONSTRAINED QUADRATIC SUBPROBLEM ****/
    // Compute Cauchy point
    snorm = dcauchy(*s,
                    alpha_,
                    q,
                    x,
                    state_->gradientVec->dual(),
                    *secant_,
                    *dwa1,
                    *dwa2,
                    outStream); // Approximately solve 1D optimization problem for alpha
    x.plus(*s);                 // Set x = proj(x[0] - alpha*g)

    // Model gradient at s = x[1] - x[0]
    gmod->set(*dwa1); // hessVec from Cauchy point computation
    gmod->plus(*state_->gradientVec);
    gfree->set(*gmod);
    //bnd.pruneActive(*gfree,x,zero);
    pwa1->set(gfree->dual());
    bnd.pruneActive(*pwa1,x,zero);
    gfree->set(pwa1->dual());
    if (hasEcon_) {
      applyFreePrecond(*pwa1,*gfree,x,*secant_,bnd,tol0,*dwa1,*pwa2);
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

    // Linearly constrained quadratic subproblem solve
    // Run CG for subspace minimization
    tol    = std::min(tol1_,tol2_*std::pow(gfnorm,spexp_));
    stol   = tol; //zero;
    if (gfnorm > zero) {
      snorm = dpcg(*s, SPflag_, SPiter_, *gfree, x, *secant_, bnd,
                   tol, stol, maxit_,
                   *pwa1, *dwa1, *pwa2, *dwa2, *pwa3, *dwa3,
                   outStream);
      if (verbosity_ > 1) {
        outStream << "  Computation of CG step"               << std::endl;
        outStream << "    CG step length:                   " << snorm   << std::endl;
        outStream << "    Number of CG iterations:          " << SPiter_ << std::endl;
        outStream << "    CG flag:                          " << SPflag_ << std::endl;
        outStream << std::endl;
      }
      x.plus(*s);
      // Project to ensure that step is feasible
      // May lose feasibility because of numerical errors
      proj_->project(x,outStream); state_->nproj++;
    }
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    gs = state_->gradientVec->apply(*state_->stepVec);

    // Projected search
    state_->snorm = dsrch(*state_->iterateVec,*state_->stepVec,ftrial,state_->searchSize,
                          state_->value,gs,obj,*pwa1,outStream);
    x.set(*state_->iterateVec);

    // Update algorithm state
    state_->iter++;
    state_->value = ftrial;
    obj.update(x,UpdateType::Accept,state_->iter);
    dwa1->set(*state_->gradientVec);
    obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
    state_->gnorm = TypeB::Algorithm<Real>::optimalityCriterion(x,*state_->gradientVec,*pwa1,outStream);

    // Update secant information
    secant_->updateStorage(x,*state_->gradientVec,*dwa1,*state_->stepVec,
                           state_->snorm,state_->iter);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
Real LSecantBAlgorithm<Real>::dgpstep(Vector<Real> &s, const Vector<Real> &w,
                                const Vector<Real> &x, const Real alpha,
                                std::ostream &outStream) const {
  const Real one(1);
  s.set(x); s.axpy(alpha,w);
  proj_->project(s,outStream); state_->nproj++;
  s.axpy(-one,x);
  return s.norm();
}

template<typename Real>
Real LSecantBAlgorithm<Real>::dcauchy(Vector<Real> &s,
                                      Real &alpha,
                                      Real &q,
                                      const Vector<Real> &x,
                                      const Vector<Real> &g,
                                      Secant<Real> &secant,
                                      Vector<Real> &dwa, Vector<Real> &dwa1,
                                      std::ostream &outStream) {
  const Real half(0.5);
  bool interp = false;
  Real gs(0), snorm(0);
  // Compute s = P(x[0] - alpha g[0])
  snorm  = dgpstep(s,g,x,-alpha,outStream);
  secant.applyB(dwa,s);
  gs     = s.dot(g);
  q      = half * s.apply(dwa) + gs;
  interp = (q > mu0_*gs);
  // Either increase or decrease alpha to find approximate Cauchy point
  int cnt = 0;
  if (interp) {
    bool search = true;
    while (search) {
      alpha *= interpf_;
      snorm  = dgpstep(s,g,x,-alpha,outStream);
      secant.applyB(dwa,s);
      gs     = s.dot(g);
      q      = half * s.apply(dwa) + gs;
      search = (q > mu0_*gs) && (cnt < redlim_);
      cnt++;
    }
  }
  else {
    bool search = true;
    Real alphas = alpha;
    Real qs     = q;
    dwa1.set(dwa);
    while (search) {
      alpha *= extrapf_;
      snorm  = dgpstep(s,g,x,-alpha,outStream);
      if (cnt < explim_) {
        secant.applyB(dwa,s);
        gs = s.dot(g);
        q  = half * s.apply(dwa) + gs;
        if (q <= mu0_*gs && std::abs(q-qs) > qtol_*std::abs(qs)) {
          dwa1.set(dwa);
          search = true;
          alphas = alpha;
          qs     = q;
        }
        else {
          q      = qs;
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
    outStream << "    Model decrease:                   " << -q    << std::endl;
    if (!interp) {
      outStream << "    Number of extrapolation steps:    " << cnt << std::endl;
    }
  }
  return snorm;
}

template<typename Real>
Real LSecantBAlgorithm<Real>::dsrch(Vector<Real> &x, Vector<Real> &s, Real &fnew, Real &beta,
                                    Real fold, Real gs, Objective<Real> &obj,
                                    Vector<Real> &pwa, std::ostream &outStream) {
  const Real one(1);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  Real snorm(0);
  int nsteps = 0;
  // Reduce beta until sufficient decrease is satisfied
  bool search = true;
  beta = one;
  while (search) {
    nsteps++;
    pwa.set(x); pwa.axpy(beta,s);
    obj.update(pwa,UpdateType::Trial);
    fnew  = obj.value(pwa,tol); state_->nfval++;
    if (fnew <= fold + mu0_*beta*gs) search = false;
    else                             beta  *= interpfPS_;
  }
  s.scale(beta);
  x.plus(s);
  snorm = s.norm();
  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "  Line search"                          << std::endl;
    outStream << "    Step length (beta):               " << beta   << std::endl;
    outStream << "    Step length (beta*s):             " << snorm  << std::endl;
    outStream << "    New objective value:              " << fnew   << std::endl;
    outStream << "    Old objective value:              " << fold   << std::endl;
    outStream << "    Descent verification (gs):        " << gs     << std::endl;
    outStream << "    Number of steps:                  " << nsteps << std::endl;
  }
  return snorm;
}

template<typename Real>
Real LSecantBAlgorithm<Real>::dpcg(Vector<Real> &w, int &iflag, int &iter,
                                  const Vector<Real> &g, const Vector<Real> &x,
                                  Secant<Real> &secant,
                                  BoundConstraint<Real> &bnd,
                                  const Real tol, const Real stol, const int itermax,
                                  Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
                                  Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa,
                                  std::ostream &outStream) const {
  // p = step (primal)
  // q = hessian applied to step p (dual)
  // t = gradient (dual)
  // r = preconditioned gradient (primal)
  Real tol0 = ROL_EPSILON<Real>(), tolB(0);
  const Real minIntSize = ROL_EPSILON<Real>();
  const Real zero(0), half(0.5), one(1), two(2);
  Real rho(0), kappa(0), beta(0), sigma(0), alpha(0), alphatmp(0);
  Real rtr(0), tnorm(0), sMs(0), pMp(0), sMp(0);
  iter = 0; iflag = 0;
  // Initialize step
  w.zero();
  // Compute residual
  t.set(g); t.scale(-one);
  // Preconditioned residual
  applyFreePrecond(r,t,x,secant,bnd,tol0,dwa,pwa);
  //rho = r.dot(t.dual());
  rho = r.apply(t);
  // Initialize direction
  p.set(r);
  pMp = (!hasEcon_ ? rho : p.dot(p)); // If no equality constraint, used preconditioned norm
  // Iterate CG
  for (iter = 0; iter < itermax; ++iter) {
    // Apply Hessian to direction dir
    applyFreeHessian(q,p,x,secant,bnd,tol0,pwa);
    // Compute step length
    kappa = p.apply(q);
    alpha = rho/kappa;
    // Check if iterate is feasible
    pwa.set(x); pwa.plus(w); pwa.axpy(alpha,p);
    if (!bnd.isFeasible(pwa)) { // Bisection to find the largest feasible step size
      sigma = zero;
      tolB  = minIntSize * alpha;
      while (alpha-sigma > tolB) {
        alphatmp = sigma + half * (alpha - sigma);
        pwa.set(x); pwa.plus(w); pwa.axpy(alphatmp,p);
        if (!bnd.isFeasible(pwa)) alpha = alphatmp;
        else                      sigma = alphatmp;
      }
      w.axpy(sigma,p);
      t.axpy(-sigma,q);
      sMs = sMs + two*sigma*sMp + sigma*sigma*pMp;
      iflag = 2;
      break;
    }
    // Update iterate and residual
    w.axpy(alpha,p);
    t.axpy(-alpha,q);
    applyFreePrecond(r,t,x,secant,bnd,tol0,dwa,pwa);
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
  if (iter == itermax) iflag = 1;
  if (iflag != 1) iter++;
  return std::sqrt(sMs); // w.norm();
}

template<typename Real>
void LSecantBAlgorithm<Real>::applyFreeHessian(Vector<Real> &hv,
                                              const Vector<Real> &v,
                                              const Vector<Real> &x,
                                              Secant<Real> &secant,
                                              BoundConstraint<Real> &bnd,
                                              Real &tol,
                                              Vector<Real> &pwa) const {
  const Real zero(0);
  pwa.set(v);
  bnd.pruneActive(pwa,x,zero);
  secant.applyB(hv,pwa);
  pwa.set(hv.dual());
  bnd.pruneActive(pwa,x,zero);
  hv.set(pwa.dual());
}

template<typename Real>
void LSecantBAlgorithm<Real>::applyFreePrecond(Vector<Real> &hv,
                                              const Vector<Real> &v,
                                              const Vector<Real> &x,
                                              Secant<Real> &secant,
                                              BoundConstraint<Real> &bnd,
                                              Real &tol,
                                              Vector<Real> &dwa,
                                              Vector<Real> &pwa) const {
  if (!hasEcon_) {
    const Real zero(0);
    pwa.set(v.dual());
    bnd.pruneActive(pwa,x,zero);
    dwa.set(pwa.dual());
    secant.applyH(hv,dwa);
    bnd.pruneActive(hv,x,zero);
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
void LSecantBAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(114,'-') << std::endl;
    os << " L-Secant-B line search method status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)" << std::endl;
    os << "  value   - Objective function value" << std::endl; 
    os << "  gnorm   - Norm of the gradient" << std::endl;
    os << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    os << "  LSpar   - Line-Search parameter" << std::endl;
    os << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    os << "  #grad   - Number of times the gradient was computed" << std::endl;
    os << "  #proj   - Number of times the projection was applied" << std::endl;
    os << "  iterCG - Number of Truncated CG iterations" << std::endl << std::endl;
    os << "  flagGC - Trust-Region Truncated CG flag" << std::endl;
    os << "       0 - Converged" << std::endl;
    os << "       1 - Iteration Limit Exceeded" << std::endl;
    os << "       2 - Bounds Exceeded" << std::endl;
    os << std::string(114,'-') << std::endl;
  }
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(15) << std::left << "LSpar";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "#proj";
  os << std::setw(10) << std::left << "iterCG";
  os << std::setw(10) << std::left << "flagCG";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void LSecantBAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "L-Secant-B Line-Search Method (Type B, Bound Constraints)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void LSecantBAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(10) << std::left << state_->nproj;
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
    os << std::setw(10) << std::left << state_->nproj;
    os << std::setw(10) << std::left << SPiter_;
    os << std::setw(10) << std::left << SPflag_;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeB
} // namespace ROL

#endif
