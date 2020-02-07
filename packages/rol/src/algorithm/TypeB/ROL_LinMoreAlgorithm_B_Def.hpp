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

#ifndef ROL_LINMOREALGORITHM_B_DEF_H
#define ROL_LINMOREALGORITHM_B_DEF_H

namespace ROL {

template<typename Real>
LinMoreAlgorithm_B<Real>::LinMoreAlgorithm_B(ParameterList &list) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  maxit_       = parlist.sublist("General").sublist("Krylov").get("Iteration Limit",20);
  tol1_        = parlist.sublist("General").sublist("Krylov").get("Absolute Tolerance",1e-4);
  tol2_        = parlist.sublist("General").sublist("Krylov").get("Relative Tolerance",1e-2;
  verbosity_   = list.sublist("General").get("Output Level",0);
  printHeader_ = verbosity_ > 2;
}

template<typename Real>
void LinMoreAlgorithm_B<Real>::initialize(Vector<Real>          &x,
                                          const Vector<Real>    &g,
                                          Objective<Real>       &obj,
                                          BoundConstraint<Real> &bnd,
                                          std::ostream &outStream) {
  const Real one(1);
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(x),makePtrFromRef(bnd));
    hasEcon_ = false;
  }
  if (hasEcon_) {
    rcon_ = makePtr<ReducedConstraint<Real>>(proj_->getLinearConstraint(),
                                             makePtrFromRef(bnd),
                                             makePtrFromRef(x));
    ran_  = proj_->getMultiplier()->dual().clone();
  }
  // Initialize data
  Algorithm_B<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  proj_->project(x);
  obj.update(x,true,state_->iter);    
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec);
  Real fnew = state_->value;
  if (!useralpha_) {
    // Evaluate objective at P(x - g)
    obj.update(*state_->stepVec,false);
    fnew = obj.value(*state_->stepVec,ftol);
    state_->nfval++;
  }
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  if (!useralpha_) {
    const Real half(0.5);
    // Minimize quadratic interpolate to compute new alpha
    Real gs    = state_->stepVec->dot(state_->gradientVec->dual());
    Real denom = (fnew - state_->value - gs);
    alpha0_ = ((denom > ROL_EPSILON<Real>()) ? -half*gs/denom : alpha0bnd_);
    alpha0_ = ((alpha0_ > alpha0bnd_) ? alpha0_ : 10.0*one);
  }
}

template<typename Real>
std::vector<std::string> LinMoreAlgorithm_B<Real>::run(Vector<Real>          &x,
                                                       const Vector<Real>    &g, 
                                                       Objective<Real>       &obj,
                                                       BoundConstraint<Real> &bnd,
                                                       std::ostream          &outStream ) {
  const Real one(1);
  // Initialize trust-region data
  std::vector<std::string> output;
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone();
  Real ftrial(0), gs(0), tol(std::sqrt(ROL_EPSILON<Real>()));
  int iflag = 0, iter = 0;

  // Output
  output.push_back(print(true));
  if (verbosity_ > 0) outStream << print(true);

  while (status_->check(*state_)) {
    // Solve trust-region subproblem
    solve(*state_->stepVec,x,*state_->gradientVec,state_->snorm,iflag,iter,
          state_->searchSize,obj,bnd,outStream);

    // Compute trial objective value
    x.plus(*state_->stepVec);
    obj.update(x,false);
    ftrial = obj.value(x,tol);

    // Compute ratio of acutal and predicted reduction
    TRflag_ = TR_FLAG_SUCCESS;
    analyzeRatio(rho,TRflag_,state_->value,ftrial,pRed,outStream);

    // Update algorithm state
    state_->iter++;
    // Accept/reject step and update trust region radius
    if ((rho < eta0_ && TRflag_ == TR_FLAG_SUCCESS)
        || (TRflag_ >= 2)) { // Step Rejected
      x.set(*state_->iterateVec);
      obj.update(x,false,state_->iter);
      if (rho < zero && TRflag_ != TR_FLAG_NAN) {
        // Negative reduction, interpolate to find new trust-region radius
        state_->searchSize = interpolateRadius(*state_->gradientVec,*state_->stepVec,
          state_->snorm,pRed,state_->value,ftrial,state_->searchSize,outStream);
      }
      else { // Shrink trust-region radius
        state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
      }
    }
    else if ((rho >= eta0_ && TRflag_ != TR_FLAG_NPOSPREDNEG)
             || (TRflag_ == TR_FLAG_POSPREDNEG)) { // Step Accepted
      state_->iterateVec->set(x);
      state_->value = ftrial;
      obj.update(x,true,state_->iter);
      // Increase trust-region radius
      if (rho >= eta2_) state_->searchSize *= gamma2_;
      // Compute gradient at new iterate
      gvec->set(*state_->gradientVec);
      obj.gradient(*state_->gradientVec,x,tol);
      // Update secant information in trust-region model
      model_->update(x,*state_->stepVec,*gvec,*state_->gradientVec,
                     state_->snorm,state_->iter);
    }

    // Update Output
    output.push_back(print(printHeader_));
    if (verbosity_ > 0) outStream << print(printHeader_);
  }
  output.push_back(Algorithm_B<Real>::printExitStatus());
  if (verbosity_ > 0) outStream << Algorithm_B<Real>::printExitStatus();
  return output;
}

template<typename Real>
void LinMoreAlgorithm_B<Real>::solve( Vector<Real>           &s,
                                      const Vector<Real>     &x,
                                      const Vector<Real>     &g,
                                      Real                   &snorm,
                                      int                    &iflag,
                                      int                    &iter,
                                      const Real              del,
                                      Objective<Real> &obj, BoundConstraint<Real> &bnd,
                                      std::ostream &outStream ) {
  const Real zero(0), half(0.5), one(1);
  Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  Real gfnorm(0), gfnormf(0), tol(0), stol(0);
  int dim = s.dimension();
  // Compute Cauchy point (TRON notation: x_ = x[1])
  snorm = dcauchy(*s_,alpha_,x,g.dual(),del,obj,
                  *pwa1_,*pwa2_,*dwa1_,outStream); // Solve 1D optimization problem for alpha_
  x_->set(x);                                      // TRON notation: model.getIterate() = x[0]
  x_->plus(*s_);                                   // Set x_ = x[0] + alpha_*g
  proj_->project(*x_);                             // Project x_ onto feasible set

  // Model gradient at s = x[1] - x[0]
  s.set(*x_); s.axpy(-one,x); // s_ = x[i+1]-x[0]
  obj.hessVec(*g_,s,x,tol0);
  g_->plus(g);
  bnd.pruneActive(*g_,*x_,zero);
  gfnorm = g_->norm();
  if (verbosity_ > 1) {
    std::cout << std::endl;
    std::cout << "  Computation of Cauchy point"          << std::endl;
    std::cout << "    Norm of Cauchy point:             " << x_->norm() << std::endl;
    std::cout << "    Norm of free gradient components: " << gfnorm     << std::endl;
  }

  // Main step computation loop
  // There are at most dim iterations since at least one
  // face becomes active at each iteration
  iter = 0;
  for (int i = 0; i < dim; ++i) {
    // Run Truncated CG
    int flagCG = 0, iterCG = 0;
    tol  = tol2_*gfnorm;
    stol = zero;
    snorm = dtrpcg(*s_,flagCG,iterCG,*g_,*x_,del,obj,bnd,
                   tol,stol,maxit_,
                   *pwa1_,*dwa1_,*pwa2_,*dwa2_,pwa3_,dwa3_);
    iter += iterCG;
    if (verbosity_ > 1) {
      std::cout << std::endl;
      std::cout << "  Computation of CG step"               << std::endl;
      std::cout << "    Number of faces:                  " << dim        << std::endl;
      std::cout << "    Current face (i):                 " << i          << std::endl;
      std::cout << "    CG step length:                   " << snorm      << std::endl;
      std::cout << "    Number of CG iterations:          " << iterCG     << std::endl;
      std::cout << "    CG flag:                          " << flagCG     << std::endl;
      std::cout << "    Total number of iterations:       " << iter       << std::endl;
    }

    // Projected search
    snorm = dprsrch(*x_,*s_,g_->dual(),obj,*pwa1_,*dwa1_);
    if (verbosity_ > 1) {
      std::cout << "    Step length (beta*s):             " << snorm      << std::endl;
      std::cout << "    Iterate length:                   " << x_->norm() << std::endl;
    }

    // Model gradient at s = x[i+1] - x[0]
    s.set(*x_); s.axpy(-one,x); // s_ = x[i+1]-x[0]
    obj.hessVec(*g_,s,x,tol0);
    g_->plus(g);
    bnd.pruneActive(*g_,*x_,zero);
    gfnormf = g_->norm();
    if (verbosity_ > 0) {
      std::cout << std::endl;
      std::cout << "  Update model gradient"                << std::endl;
      std::cout << "    Step length:                      " << s.norm()   << std::endl;
      std::cout << "    Norm of free gradient components: " << gfnormf    << std::endl;
      std::cout << std::endl;
    }

    // Termination check
    if (gfnormf <= tol2_*gfnorm) {
      iflag = 0;
      break;
    }
    else if (iter >= maxit_) {
      iflag = 1;
      break;
    }
    else if (flagCG == 2) {
      iflag = 2;
      break;
    }
    else if (flagCG == 3) {
      iflag = 3;
      break;
    }

    // Update free gradient norm
    gfnorm = gfnormf;
  }
  // Update norm of step and update model predicted reduction
  snorm = s.norm();
  Real gs(0);
  obj.hessVec(*dwa1_,s,x,tol0);
  gs = s.dot(g->dual());
  pRed_ = -half * s.dot(dwa1_->dual()) + gs;
}

template<typename Real>
Real LinMoreAlgorithm_B<Real>::dgpstep(Vector<Real> &s, const Vector<Real> &w,
                                 const Vector<Real> &x, const Real alpha) const {
  s.set(x); s.axpy(alpha,w);
  proj_->project(s);
  s.axpy(static_cast<Real>(-1),x);
  return s.norm();
}

template<typename Real>
Real LinMoreAlgorithm_B<Real>::dcauchy(Vector<Real> &s, Real &alpha,
                                       const Vector<Real> &x, const Vector<Real> &g,
                                       const Real del, Objective<Real> &obj,
                                       Vector<Real> &dwa) {
  const Real half(0.5), one(1), mu0(0.01), interpf(0.1), extrapf(10);
  // const Real zero(0); // Unused
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  bool interp = false;
  Real q(0), gs(0), snorm(0);
  // Compute s = P(x[0] - alpha g[0].dual())
  snorm = dgpstep(s,g,x,-alpha);
  if (snorm > del) {
    interp = true;
  }
  else {
    obj.hessVec(dwa,s,x,tol);
    gs = s.dot(g);
    q  = half * s.dot(dwa.dual()) + gs;
    interp = (q > mu0*gs);
  }
  // Either increase or decrease alpha to find approximate Cauchy point
  if (interp) {
    bool search = true;
    while (search) {
      alpha *= interpf;
      snorm = dgpstep(s,g,x,-alpha);
      if (snorm <= del) {
        obj.hessVec(dwa,s,x,tol);
        gs = s.dot(g);
        q  = half * s.dot(dwa.dual()) + gs;
        search = (q > mu0*gs);
      }
    }
  }
  else {
    bool search = true;
    Real alphas = alpha;
    while (search) {
      alpha *= extrapf;
      snorm = dgpstep(s,g,x,-alpha);
      if (snorm <= del) {
        obj.hessVec(dwa,s,x,tol);
        gs = s.dot(g);
        q  = half * s.dot(dwa.dual()) + gs;
        if (q < mu0*gs) {
          search = true;
          alphas = alpha;
        }
      }
      else {
        search = false;
      }
    }
    alpha = alphas;
    snorm = dgpstep(s,g,x,-alpha);
  }
  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "  Cauchy point"                         << std::endl;
    std::cout << "    Step length (alpha):              " << alpha_     << std::endl;
    std::cout << "    Step length (alpha*g):            " << snorm      << std::endl;
  }
  return snorm;
}

template<typename Real>
Real LinMoreAlgorithm_B<Real>::dprsrch(Vector<Real> &x, Vector<Real> &s,
                                       const Vector<Real> &g, Objective<Real> &obj,
                                       Vector<Real> &pwa, Vector<Real> &dwa,
                                       std::ostream &outStream) {
  const Real half(0.5), one(1), mu0(0.01), interpf(0.5);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  Real beta(1), snorm(0), q(0), gs(0);
  int nsteps = 0;
  // Reduce beta until sufficient decrease is satisfied
  bool search = true;
  while (search) {
    nsteps++;
    snorm = dgpstep(pwa,s,x,beta);
    applyFreeHessian(dwa,pwa,x,obj,tol);
    gs = pwa.dot(g);
    q  = half * s.dot(dwa.dual()) + gs;
    if (q <= mu0*gs || nsteps > psmax_) {
      search = false;
    }
    else {
      beta *= interpf;
    }
  }
  snorm = dgpstep(pwa,s,x,beta);
  s.set(pwa);
  x.plus(s);
  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "  Projected search"                     << std::endl;
    outStream << "    Step length (beta):               " << beta       << std::endl;
    outStream << "    Number of steps:                  " << nsteps     << std::endl;
  }
  return snorm;
}

template<typename Real>
Real LinMoreAlgorithm_B<Real>::dtrqsol(const Real xtx,
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
Real LinMoreAlgorithm_B<Real>::dtrpcg(Vector<Real> &w, int &iflag, int &iter,
                                      const Vector<Real> &g, const Vector<Real> &x,
                                      const Real del, Objective<Real> &obj,
                                      BoundConstraint<Real> &bnd,
                                      const Real tol, const Real stol, const int itermax,
                                      Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
                                      Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa) const {
  Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  const Real zero(0), one(1), two(2);
  Real rho(0), tnorm(0), rnorm(0), rnorm0(0), kappa(0), beta(0), sigma(0), alpha(0), rtr(0);
  Real sMs(0), pMp(0), sMp(0);
  iter = 0; iflag = 0;
  // Initialize step
  w.zero();
  // Compute residual
  t.set(g); t.scale(-one);
  // Preconditioned residual
  applyFreePrecond(r,t,x,obj,bnd,tol0,dwa);
  rho    = r.dot(t.dual());
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
    applyFreeHessian(q,p,x,obj,bnd,tol0,pwa);
    // Compute sigma such that ||s+sigma*dir|| = del
    kappa = p.dot(q.dual());
    alpha = (kappa>zero) ? rho/kappa : zero;
    sigma = dtrqsol(sMs,pMp,sMp,del);
    // Check for negative curvature or if iterate exceeds trust region
    if (kappa <= zero || alpha >= sigma) {
      w.axpy(sigma,p);
      iflag = (kappa<=zero) ? 2 : 3;
      break;
    }
    // Update iterate and residuals
    w.axpy(alpha,p);
    t.axpy(-alpha,q);
    applyFreePrecond(r,t,x,obj,bnd,tol0,dwa);
    // Exit if residual tolerance is met
    rtr   = r.dot(t.dual());
    rnorm = std::sqrt(rtr);
    tnorm = t.norm();
    if (rnorm <= stol || tnorm <= tol) {
      iflag = 0;
      break;
    }
    // Compute p = r + beta * p
    beta = rtr/rho;
    p.scale(beta); p.plus(r);
    rho  = rtr;
    // Update dot products
    //   sMs = <s, inv(M)s>
    //   sMp = <s, inv(M)p>
    //   pMp = <p, inv(M)p>
    sMs  = sMs + two*alpha*sMp + alpha*alpha*pMp;
    sMp  = beta*(sMp + alpha*pMp);
    pMp  = rho + beta*beta*pMp;
  }
  // Check iteration count
  if (iter == itermax) {
    iflag = 1;
  }
  if (iflag != 1) { 
    iter++;
  }
  return w.norm();
}

template<typename Real>
void LinMoreAlgorithm_B<Real>::applyFreeHessian(Vector<Real> &hv,
                                                const Vector<Real> &v,
                                                const Vector<Real> &x,
                                                Objective<Real> &obj,
                                                BoundConstraint<Real> &bnd,
                                                Real &tol,
                                                Vector<Real> &pwa) {
  const Real zero(0);
  pwa.set(v);
  bnd.pruneActive(pwa,x,zero);
  obj.hessVec(hv,pwa,x,tol);
  bnd.pruneActive(hv,x,zero);
}

template<typename Real>
void LinMoreAlgorithm_B<Real>::applyFreePrecond(Vector<Real> &hv,
                                                const Vector<Real> &v,
                                                const Vector<Real> &x,
                                                Objective<Real> &obj,
                                                BoundConstraint<Real> &bnd,
                                                Real &tol,
                                                Vector<Real> &dwa) {
  if (!hasEcon_) {
    const Real zero(0);
    dwa.set(v);
    bnd.pruneActive(dwa,x,zero);
    obj.precond(hv,dwa,x,tol);
    bnd.pruneActive(hv,x,zero);
  }
  else {
    // Perform null space projection
    rcon_->setX(makePtrFromRef(x));
    NullSpaceOperator<Real> ns(rcon_,makePtrFromRef(x),ran_);
    ns.apply(hv,v,tol);
  }
}

template<typename Real>
std::string LinMoreAlgorithm_B<Real>::printHeader( void ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(109,'-') << std::endl;
    hist << "Lin-More trust-region method";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  value    - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  del      - Trust-region radius" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "del";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::endl;
  return hist.str();
}

template<typename Real>
std::string LinMoreAlgorithm_B<Real>::printName( void ) const {
  std::stringstream hist;
  hist << std::endl << "Lin-More Trust-Region Method" << std::endl;
  return hist.str();
}

template<typename Real>
std::string LinMoreAlgorithm_B<Real>::print( const bool print_header ) const {
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
    hist << std::endl;
  }
  return hist.str();
}

} // namespace ROL

#endif
