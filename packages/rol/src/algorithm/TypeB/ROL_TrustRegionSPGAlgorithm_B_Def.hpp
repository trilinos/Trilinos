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

#ifndef ROL_TRUSTREGIONSPGALGORITHM_B_DEF_H
#define ROL_TRUSTREGIONSPGALGORITHM_B_DEF_H

#include <deque>

namespace ROL {

template<typename Real>
TrustRegionSPGAlgorithm_B<Real>::TrustRegionSPGAlgorithm_B(ParameterList &list,
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
  verbosity_    = trlist.sublist("General").get("Output Level",                     0);
  // Algorithm-Specific Parameters
  ROL::ParameterList &lmlist = trlist.sublist("SPG");
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
  // Spectral projected gradient parameters
  lsmax_     = lmlist.sublist("SPG").get("Function Evaluation Limit",                  20);
  lambdaMin_ = lmlist.sublist("SPG").get("Minimum Spectral Step Size",                 1e-8); 
  lambdaMax_ = lmlist.sublist("SPG").get("Maximum Spectral Step Size",                 5e2); 
  sigma1_    = lmlist.sublist("SPG").get("Lower Step Size Safeguard",                  0.1);
  sigma2_    = lmlist.sublist("SPG").get("Upper Step Size Safeguard",                  0.9);
  rhodec_    = lmlist.sublist("SPG").get("Backtracking Rate",                          0.5);
  gamma_     = lmlist.sublist("SPG").get("Sufficient Decrease Tolerance",              1e-4);
  maxSize_   = lmlist.sublist("SPG").get("Maximum Storage Size",                       10);
  maxit_     = lmlist.sublist("SPG").get("Iteration Limit",                            25);
  tol1_      = lmlist.sublist("SPG").get("Absolute Tolerance",                         1e-4);
  tol2_      = lmlist.sublist("SPG").get("Relative Tolerance",                         1e-2);
  // Output Parameters
  verbosity_   = list.sublist("General").get("Output Level",0);
  printHeader_ = verbosity_ > 2;
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
void TrustRegionSPGAlgorithm_B<Real>::initialize(Vector<Real>          &x,
                                                 const Vector<Real>    &g,
                                                 Objective<Real>       &obj,
                                                 BoundConstraint<Real> &bnd,
                                                 std::ostream &outStream) {
  const Real one(1);
  if (proj_ == nullPtr)
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
  // Initialize data
  Algorithm_B<Real>::initialize(x,g);
  nhess_ = 0;
  // Update approximate gradient and approximate objective function.
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(); 
  proj_->project(x,outStream);
  state_->iterateVec->set(x);
  obj.update(x,UPDATE_INITIAL,state_->iter);
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream);
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  // Normalize initial CP step length
  if (normAlpha_) alpha_ /= state_->gradientVec->norm();
  // Compute initial trust region radius if desired.
  if ( state_->searchSize <= static_cast<Real>(0) )
    state_->searchSize = state_->gradientVec->norm();
}

template<typename Real>
std::vector<std::string> TrustRegionSPGAlgorithm_B<Real>::run(Vector<Real>          &x,
                                                              const Vector<Real>    &g, 
                                                              Objective<Real>       &obj,
                                                              BoundConstraint<Real> &bnd,
                                                              std::ostream          &outStream ) {
  const Real zero(0), one(1);
  Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  Real ftrial(0), pRed(0), rho(1), q(0);
  // Initialize trust-region data
  std::vector<std::string> output;
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> gmod = g.clone();
  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone();
  Ptr<Vector<Real>> pwa3 = x.clone(), pwa4 = x.clone();
  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone();

  // Output
  output.push_back(print(true));
  if (verbosity_ > 0) outStream << print(true);

  while (status_->check(*state_)) {
    // Build trust-region model
    model_->setData(obj,*state_->iterateVec,*state_->gradientVec);

    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
    // Compute Cauchy point (TRON notation: x = x[1])
    dcauchy(*state_->stepVec,alpha_,q,*state_->iterateVec,
            state_->gradientVec->dual(),state_->searchSize,
            *model_,*dwa1,*dwa2,outStream); // Solve 1D optimization problem for alpha
    x.plus(*state_->stepVec);               // Set x = x[0] + alpha*g

    // Model gradient at s = x[1] - x[0]
    gmod->set(*dwa1); // hessVec from Cauchy point computation
    gmod->plus(*state_->gradientVec);

    // Apply SPG starting from the Cauchy point
    dpsg(x,q,*gmod,*state_->iterateVec,state_->searchSize,*model_,
         *pwa1,*pwa2,*pwa3,*pwa4,*dwa1,outStream);
    pRed = -q;
    state_->stepVec->set(x); state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Compute trial objective value
    obj.update(x,UPDATE_TRIAL);
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
      obj.update(x,UPDATE_REVERT,state_->iter);
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
      obj.update(x,UPDATE_ACCEPT,state_->iter);
      // Increase trust-region radius
      if (rho >= eta2_) state_->searchSize = std::min(gamma2_*state_->searchSize, delMax_);
      // Compute gradient at new iterate
      dwa1->set(*state_->gradientVec);
      obj.gradient(*state_->gradientVec,x,tol0);
      state_->ngrad++;
      state_->gnorm = Algorithm_B<Real>::optimalityCriterion(x,*state_->gradientVec,*pwa1,outStream);
      state_->iterateVec->set(x);
      // Update secant information in trust-region model
      model_->update(x,*state_->stepVec,*dwa1,*state_->gradientVec,
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
Real TrustRegionSPGAlgorithm_B<Real>::dgpstep(Vector<Real> &s, const Vector<Real> &w,
                                              const Vector<Real> &x, const Real alpha,
                                              std::ostream &outStream) const {
  s.set(x); s.axpy(alpha,w);
  proj_->project(s,outStream);
  s.axpy(static_cast<Real>(-1),x);
  return s.norm();
}

template<typename Real>
Real TrustRegionSPGAlgorithm_B<Real>::dcauchy(Vector<Real> &s,
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
void TrustRegionSPGAlgorithm_B<Real>::dpsg(Vector<Real> &y,
                                           Real         &q,
                                           Vector<Real> &gmod,
                                           const Vector<Real> &x,
                                           Real del,
                                           TrustRegionModel_U<Real> &model,
                                           Vector<Real> &pwa,
                                           Vector<Real> &pwa1,
                                           Vector<Real> &pwa2,
                                           Vector<Real> &pwa3,
                                           Vector<Real> &dwa,
                                           std::ostream &outStream) {
  // Use SPG to approximately solve TR subproblem:
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)>  subject to y\in C, ||y|| \le del
  //
  //   Inpute:
  //       y = Cauchy step
  //       x = Current iterate
  //       g = Current gradient
  const Real half(0.5), one(1), eps(std::sqrt(ROL_EPSILON<Real>()));
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  Real alpha(1), mtrial(0), sHs(0), alphaTmp(1), mmax(0);
  int ls_nfval(0);
  std::deque<Real> mqueue; mqueue.push_back(q);

  // Compute initial projected gradient norm
  pwa1.set(gmod.dual());
  pwa.set(y);
  pwa.axpy(-one,pwa1);
  dproj(pwa,x,del,pwa2,pwa3,outStream);
  pwa.axpy(-one,y);
  Real gnorm = pwa.norm();
  const Real gtol = std::min(tol1_,tol2_*gnorm);

  // Compute initial step
  Real lambda = std::max(lambdaMin_,std::min(one/gmod.norm(),lambdaMax_));
  pwa.set(y); pwa.axpy(-lambda,pwa1);        // pwa = y - lambda gmod.dual()
  dproj(pwa,x,del,pwa2,pwa3,outStream); // pwa = P(y - lambda gmod.dual())
  pwa.axpy(-one,y);                          // pwa = P(y - lambda gmod.dual()) - y = step
  Real gs = gmod.apply(pwa);                 // gs  = <step, model gradient>
  Real ss = pwa.dot(pwa);                    // Norm squared of step

  if (verbosity_ > 1)
    outStream << "  Spectral Projected Gradient"          << std::endl;

  SPiter_ = 0;
  while (SPiter_ < maxit_) {
    SPiter_++;

    // Evaluate model
    model.hessVec(dwa,pwa,x,tol); nhess_++; // dwa = H step
    sHs    = dwa.apply(pwa);                // sHs = <step, H step>
    mtrial = q + gs + half * sHs;           // Current model value

    // Perform nonmonotone line search
    ls_nfval = 0;
    alpha = one;
    mmax  = *std::max_element(mqueue.begin(),mqueue.end());
    while (mtrial > mmax + gamma_*alpha*gs && ls_nfval < lsmax_) {
      alphaTmp = -half*alpha*alpha*gs/(mtrial-q-alpha*gs);
      alpha    = (sigma1_*alpha <= alphaTmp && alphaTmp <= sigma2_*alpha) ? alphaTmp : rhodec_*alpha;
      mtrial   = q + alpha * (gs + half * alpha * sHs); ls_nfval++;
    }
    y.axpy(alpha,pwa);    // y = y + alpha pwa
    q = mtrial;           // New model value at y
    gmod.axpy(alpha,dwa); // New model gradient at y
    if (static_cast<int>(mqueue.size())==maxSize_) mqueue.pop_front();
    mqueue.push_back(q);

    // Compute projected gradient norm
    pwa1.set(gmod.dual());
    pwa.set(y);
    pwa.axpy(-one,pwa1);
    dproj(pwa,x,del,pwa2,pwa3,outStream);
    pwa.axpy(-one,y);
    gnorm = pwa.norm();

    if (verbosity_ > 1) {
      outStream << "    Iterate:                          " << SPiter_ << std::endl;
      outStream << "    Spectral step length (lambda):    " << lambda  << std::endl;
      outStream << "    Step length (alpha):              " << alpha   << std::endl;
      outStream << "    Model decrease (pRed):            " << -q      << std::endl;
      outStream << "    Optimality criterion:             " << gnorm   << std::endl;
    }
    if (gnorm < gtol) break;

    // Compute new spectral step
    lambda = (sHs<=eps ? lambdaMax_ : std::max(lambdaMin_,std::min(ss/sHs,lambdaMax_)));
    pwa.set(y); pwa.axpy(-lambda,pwa1);
    dproj(pwa,x,del,pwa2,pwa3,outStream);
    pwa.axpy(-one,y);
    gs = gmod.apply(pwa);
    ss = pwa.dot(pwa);
  }
  SPflag_ = (SPiter_==maxit_) ? 1 : 0;
}

template<typename Real>
void TrustRegionSPGAlgorithm_B<Real>::dproj(Vector<Real> &x,
                                            const Vector<Real> &x0,
                                            Real del,
                                            Vector<Real> &y,
                                            Vector<Real> &p,
                                            std::ostream &outStream) const {
  const Real half(0.5), one(1), tol(1e-2*std::sqrt(ROL_EPSILON<Real>()));
  const Real fudge(1.0-1e-6);
  Real e(0), a0(0), a1(1), a(0.5);
  y.set(x);
  proj_->project(y,outStream);
  p.set(y); p.axpy(-one,x0);
  e = p.norm();
  if (e <= del) {
    x.set(y);
    return;
  }
  while (a1-a0 > tol) {
    y.set(x); y.scale(a); y.axpy(one-a,x0);
    proj_->project(y,outStream);
    p.set(y); p.axpy(-one,x0);
    e = p.norm();
    if (e < fudge*del) a0 = a; // Move lower end point
    else if (e > del)  a1 = a; // Move upper end point
    else               break;  // Exit if (1-1e-6)*del <= snorm <= del
    a = a0 + half*(a1-a0);
  }
  x.set(y);
}

template<typename Real>
std::string TrustRegionSPGAlgorithm_B<Real>::printHeader( void ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    hist << " SPG trust-region method status output definitions" << std::endl << std::endl;
    hist << "  iter    - Number of iterates (steps taken)" << std::endl;
    hist << "  value   - Objective function value" << std::endl; 
    hist << "  gnorm   - Norm of the gradient" << std::endl;
    hist << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  delta   - Trust-Region radius" << std::endl;
    hist << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    hist << "  #grad   - Number of times the gradient was computed" << std::endl;
    hist << "  #hess   - Number of times the Hessian was applied" << std::endl;
    hist << std::endl;
    hist << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      hist << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    hist << std::endl;
    hist << "  iterSPG - Number of Spectral Projected Gradient iterations" << std::endl << std::endl;
    hist << "  flagSPG - Trust-Region Truncated CG flag" << std::endl;
    hist << "    0 - Converged" << std::endl;
    hist << "    1 - Iteration Limit Exceeded" << std::endl;
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
  hist << std::setw(10) << std::left << "tr_flag";
  hist << std::setw(10) << std::left << "iterSPG";
  hist << std::setw(10) << std::left << "flagSPG";
  hist << std::endl;
  return hist.str();
}

template<typename Real>
std::string TrustRegionSPGAlgorithm_B<Real>::printName( void ) const {
  std::stringstream hist;
  hist << std::endl << "SPG Trust-Region Method (Type B, Bound Constraints)" << std::endl;
  return hist.str();
}

template<typename Real>
std::string TrustRegionSPGAlgorithm_B<Real>::print( const bool print_header ) const {
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
    hist << std::setw(15) << std::left << state_->searchSize;
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << nhess_;
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
    hist << std::setw(10) << std::left << nhess_;
    hist << std::setw(10) << std::left << TRflag_;
    hist << std::setw(10) << std::left << SPiter_;
    hist << std::setw(10) << std::left << SPflag_;
    hist << std::endl;
  }
  return hist.str();
}

} // namespace ROL

#endif
