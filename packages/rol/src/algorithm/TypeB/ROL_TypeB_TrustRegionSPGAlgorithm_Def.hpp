// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_TRUSTREGIONSPGALGORITHM_DEF_HPP
#define ROL_TYPEB_TRUSTREGIONSPGALGORITHM_DEF_HPP

#include <deque>

namespace ROL {
namespace TypeB {

template<typename Real>
TrustRegionSPGAlgorithm<Real>::TrustRegionSPGAlgorithm(ParameterList &list,
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
  verbosity_    = trlist.sublist("General").get("Output Level",   0);
  // Nonmonotone Parameters
  storageNM_ = trlist.get("Nonmonotone Storage Size",             0);
  useNM_     = (storageNM_ <= 0 ? false : true);
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
  lambdaMin_ = lmlist.sublist("Solver").get("Minimum Spectral Step Size",          1e-8);
  lambdaMax_ = lmlist.sublist("Solver").get("Maximum Spectral Step Size",          1e8);
  gamma_     = lmlist.sublist("Solver").get("Sufficient Decrease Tolerance",       1e-4);
  maxSize_   = lmlist.sublist("Solver").get("Maximum Storage Size",                10);
  maxit_     = lmlist.sublist("Solver").get("Iteration Limit",                     25);
  tol1_      = lmlist.sublist("Solver").get("Absolute Tolerance",                  1e-4);
  tol2_      = lmlist.sublist("Solver").get("Relative Tolerance",                  1e-2);
  useMin_    = lmlist.sublist("Solver").get("Use Smallest Model Iterate",          true);
  useNMSP_   = lmlist.sublist("Solver").get("Use Nonmonotone Search",              false);
  useSimpleSPG_ = !lmlist.sublist("Solver").get("Compute Cauchy Point",            true);
  // Inexactness Information
  ParameterList &glist = list.sublist("General");
  useInexact_.clear();
  useInexact_.push_back(glist.get("Inexact Objective Function",     false));
  useInexact_.push_back(glist.get("Inexact Gradient",               false));
  useInexact_.push_back(glist.get("Inexact Hessian-Times-A-Vector", false));
  // Trust-Region Inexactness Parameters
  ParameterList &ilist = trlist.sublist("Inexact").sublist("Gradient");
  scale0_ = ilist.get("Tolerance Scaling",  static_cast<Real>(0.1));
  scale1_ = ilist.get("Relative Tolerance", static_cast<Real>(2)); 
  // Inexact Function Evaluation Information
  ParameterList &vlist = trlist.sublist("Inexact").sublist("Value");
  scale_       = vlist.get("Tolerance Scaling",                 static_cast<Real>(1.e-1));
  omega_       = vlist.get("Exponent",                          static_cast<Real>(0.9));
  force_       = vlist.get("Forcing Sequence Initial Value",    static_cast<Real>(1.0));
  updateIter_  = vlist.get("Forcing Sequence Update Frequency", static_cast<int>(10));
  forceFactor_ = vlist.get("Forcing Sequence Reduction Factor", static_cast<Real>(0.1));
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
void TrustRegionSPGAlgorithm<Real>::initialize(Vector<Real>          &x,
                                               const Vector<Real>    &g,
                                               Real                   ftol,
                                               Objective<Real>       &obj,
                                               BoundConstraint<Real> &bnd,
                                               std::ostream &outStream) {
  //const Real one(1);
  if (proj_ == nullPtr)
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  nhess_ = 0;
  // Update approximate gradient and approximate objective function.
  proj_->project(x,outStream); state_->nproj++;
  state_->iterateVec->set(x);
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  //obj.gradient(*state_->gradientVec,x,ftol);
  computeGradient(x,*state_->gradientVec,*state_->stepVec,state_->searchSize,obj,true,gtol_,state_->gnorm,outStream);
  state_->ngrad++;
  //state_->stepVec->set(x);
  //state_->stepVec->axpy(-one,state_->gradientVec->dual());
  //proj_->project(*state_->stepVec,outStream); state_->nproj++;
  //state_->stepVec->axpy(-one,x);
  //state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
  // Normalize initial CP step length
  if (normAlpha_) alpha_ /= state_->gradientVec->norm();
  // Compute initial trust region radius if desired.
  if ( state_->searchSize <= static_cast<Real>(0) )
    state_->searchSize = state_->gradientVec->norm();
}

template<typename Real>
Real TrustRegionSPGAlgorithm<Real>::computeValue(Real inTol,
                                                 Real &outTol,
                                                 Real pRed,
                                                 Real &fold,
                                                 int iter,
                                                 const Vector<Real> &x,
                                                 const Vector<Real> &xold,
                                                 Objective<Real> &obj) {
  outTol = std::sqrt(ROL_EPSILON<Real>());
  if ( useInexact_[0] ) {
    if (!(iter%updateIter_) && (iter!=0)) force_ *= forceFactor_;
    const Real one(1);
    Real eta = static_cast<Real>(0.999)*std::min(eta1_,one-eta2_);
    outTol   = scale_*std::pow(eta*std::min(pRed,force_),one/omega_);
    if (inTol > outTol) fold = obj.value(xold,outTol);
  }
  // Evaluate objective function at new iterate
  obj.update(x,UpdateType::Trial);
  Real fval = obj.value(x,outTol);
  return fval;
}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::computeGradient(const Vector<Real> &x,
                                                    Vector<Real> &g,
                                                    Vector<Real> &pwa,
                                                    Real del,
                                                    Objective<Real> &obj,
                                                    bool accept,
                                                    Real &gtol,
                                                    Real &gnorm,
                                                    std::ostream &outStream) const {
  if ( useInexact_[1] ) {
    const Real one(1);
    Real gtol0 = scale0_*del;
    if (accept) gtol  = gtol0 + one;
    else        gtol0 = scale0_*std::min(gnorm,del);
    while ( gtol > gtol0 ) {
      gtol = gtol0;
      obj.gradient(g,x,gtol);
      gnorm = TypeB::Algorithm<Real>::optimalityCriterion(x,g,pwa,outStream);
      gtol0 = scale0_*std::min(gnorm,del);
    }
  }
  else {
    if (accept) {
      gtol = std::sqrt(ROL_EPSILON<Real>());
      obj.gradient(g,x,gtol);
      gnorm = TypeB::Algorithm<Real>::optimalityCriterion(x,g,pwa,outStream);
    }
  }
}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::run(Vector<Real>          &x,
                                        const Vector<Real>    &g, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        std::ostream          &outStream ) {
  const Real zero(0), one(1);
  //Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  Real inTol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(), outTol(inTol);
  Real ftrial(0), pRed(0), rho(1), q(0);
  // Initialize trust-region data
  std::vector<std::string> output;
  initialize(x,g,inTol,obj,bnd,outStream);
  Ptr<Vector<Real>> gmod = g.clone();
  Ptr<Vector<Real>> pwa1 = x.clone(), pwa2 = x.clone();
  Ptr<Vector<Real>> pwa3 = x.clone(), pwa4 = x.clone();
  Ptr<Vector<Real>> pwa5 = x.clone(), pwa6 = x.clone();
  Ptr<Vector<Real>> pwa7 = x.clone();
  Ptr<Vector<Real>> dwa1 = g.clone(), dwa2 = g.clone();
  // Initialize nonmonotone data
  Real rhoNM(0), sigmac(0), sigmar(0);
  Real fr(state_->value), fc(state_->value), fmin(state_->value);
  TRUtils::ETRFlag TRflagNM;
  int L(0);

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Build trust-region model
    model_->setData(obj,*state_->iterateVec,*state_->gradientVec,gtol_);

    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
    q = zero;
    gmod->set(*state_->gradientVec);
    if (useSimpleSPG_)
      dpsg_simple(x,q,*gmod,*state_->iterateVec,state_->searchSize,*model_,
                  *pwa1,*pwa2,*dwa1,outStream);
    else {
      // Compute Cauchy point (TRON notation: x = x[1])
      dcauchy(*state_->stepVec,alpha_,q,*state_->iterateVec,
              state_->gradientVec->dual(),state_->searchSize,
              *model_,*dwa1,*dwa2,outStream); // Solve 1D optimization problem for alpha
      x.plus(*state_->stepVec);               // Set x = x[0] + alpha*g

      // Model gradient at s = x[1] - x[0]
      gmod->plus(*dwa1); // hessVec from Cauchy point computation

      // Apply SPG starting from the Cauchy point
      dpsg(x,q,*gmod,*state_->iterateVec,state_->searchSize,*model_,
           *pwa1,*pwa2,*pwa3,*pwa4,*pwa5,*pwa6,*pwa7,*dwa1,outStream);
    }

    // Update storage and compute predicted reduction
    pRed = -q;
    state_->stepVec->set(x); state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Compute trial objective value
    ftrial = computeValue(inTol,outTol,pRed,state_->value,state_->iter,x,*state_->iterateVec,obj);
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
      computeGradient(x,*state_->gradientVec,*pwa1,state_->searchSize,obj,false,gtol_,state_->gnorm,outStream);
    }
    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
             || (TRflag_ == TRUtils::POSPREDNEG)) { // Step Accepted
      state_->value = ftrial;
      obj.update(x,UpdateType::Accept,state_->iter);
      inTol = outTol;
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
      computeGradient(x,*state_->gradientVec,*pwa1,state_->searchSize,obj,true,gtol_,state_->gnorm,outStream);
      state_->ngrad++;
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
Real TrustRegionSPGAlgorithm<Real>::dgpstep(Vector<Real> &s, const Vector<Real> &w,
                                            const Vector<Real> &x, const Real alpha,
                                            std::ostream &outStream) const {
  s.set(x); s.axpy(alpha,w);
  proj_->project(s,outStream); state_->nproj++;
  s.axpy(static_cast<Real>(-1),x);
  return s.norm();
}

template<typename Real>
Real TrustRegionSPGAlgorithm<Real>::dcauchy(Vector<Real> &s,
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
void TrustRegionSPGAlgorithm<Real>::dpsg_simple(Vector<Real> &y,
                                                Real         &q,
                                                Vector<Real> &gmod,
                                                const Vector<Real> &x,
                                                Real del,
                                                TrustRegionModel_U<Real> &model,
                                                Vector<Real> &pwa,
                                                Vector<Real> &pwa1,
                                                Vector<Real> &dwa,
                                                std::ostream &outStream) {
  // Use SPG to approximately solve TR subproblem:
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)>  subject to y\in C, ||y|| \le del
  //
  //   Inpute:
  //       y = Primal vector
  //       x = Current iterate
  //       g = Current gradient
  const Real half(0.5), one(1), safeguard(1e2*ROL_EPSILON<Real>());
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  Real alpha(1), alphaMax(1), s0s0(0), ss0(0), sHs(0), lambdaTmp(1), snorm(0);
  pwa1.zero();

  // Set y = x
  y.set(x);

  // Compute initial step
  Real coeff  = one/gmod.norm();
  Real lambda = std::max(lambdaMin_,std::min(coeff,lambdaMax_));
  pwa.set(y); pwa.axpy(-lambda,gmod.dual());      // pwa = x - lambda gmod.dual()
  proj_->project(pwa,outStream); state_->nproj++; // pwa = P(x - lambda gmod.dual())
  pwa.axpy(-one,y);                               // pwa = P(x - lambda gmod.dual()) - x = step
  Real gs = gmod.apply(pwa);                      // gs  = <step, model gradient>
  Real ss = pwa.dot(pwa);                         // Norm squared of step
  Real gnorm = std::sqrt(ss);

  // Compute initial projected gradient norm
  const Real gtol = std::min(tol1_,tol2_*gnorm);

  if (verbosity_ > 1)
    outStream << "  Spectral Projected Gradient"          << std::endl;

  SPiter_ = 0;
  while (SPiter_ < maxit_) {
    SPiter_++;

    // Evaluate model Hessian
    model.hessVec(dwa,pwa,x,tol); nhess_++; // dwa = H step
    sHs = dwa.apply(pwa);                   // sHs = <step, H step>

    // Perform line search
    alphaMax = 1;
    if (gnorm >= del-safeguard) { // Trust-region constraint is violated
      ss0      = pwa1.dot(pwa);
      alphaMax = std::min(one, (-ss0 + std::sqrt(ss0*ss0 - ss*(s0s0-del*del)))/ss);
    }
    if (sHs <= safeguard)
      alpha = alphaMax;
    else
      alpha = std::min(alphaMax, -gs/sHs);

    // Update model quantities
    q += alpha * (gs + half * alpha * sHs); // Update model value
    gmod.axpy(alpha,dwa);                   // Update model gradient
    y.axpy(alpha,pwa);                      // New iterate

    // Check trust-region constraint violation
    pwa1.set(y); pwa1.axpy(-one,x);
    s0s0 = pwa1.dot(pwa1);
    snorm = std::sqrt(s0s0);

    if (verbosity_ > 1) {
      outStream << std::endl;
      outStream << "    Iterate:                          " << SPiter_ << std::endl;
      outStream << "    Spectral step length (lambda):    " << lambda  << std::endl;
      outStream << "    Step length (alpha):              " << alpha   << std::endl;
      outStream << "    Model decrease (pRed):            " << -q      << std::endl;
      outStream << "    Optimality criterion:             " << gnorm   << std::endl;
      outStream << "    Step norm:                        " << snorm   << std::endl;
      outStream << std::endl;
    }

    if (snorm >= del - safeguard) { SPflag_ = 2; break; }

    // Compute new spectral step
    lambdaTmp = (sHs <= safeguard ? one/gmod.norm() : ss/sHs);
    lambda = std::max(lambdaMin_,std::min(lambdaTmp,lambdaMax_));
    pwa.set(y); pwa.axpy(-lambda,gmod.dual());
    proj_->project(pwa,outStream); state_->nproj++;
    pwa.axpy(-one,y);
    gs = gmod.apply(pwa);
    ss = pwa.dot(pwa);
    gnorm = std::sqrt(ss);

    if (gnorm <= gtol) { SPflag_ = 0; break; }
  }
  SPflag_ = (SPiter_==maxit_) ? 1 : SPflag_;
}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::dpsg(Vector<Real> &y,
                                         Real         &q,
                                         Vector<Real> &gmod,
                                         const Vector<Real> &x,
                                         Real del,
                                         TrustRegionModel_U<Real> &model,
                                         Vector<Real> &ymin,
                                         Vector<Real> &pwa,
                                         Vector<Real> &pwa1,
                                         Vector<Real> &pwa2,
                                         Vector<Real> &pwa3,
                                         Vector<Real> &pwa4,
                                         Vector<Real> &pwa5,
                                         Vector<Real> &dwa,
                                         std::ostream &outStream) {
  // Use SPG to approximately solve TR subproblem:
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)>  subject to y\in C, ||y|| \le del
  //
  //   Inpute:
  //       y = Cauchy step
  //       x = Current iterate
  //       g = Current gradient
  const Real zero(0), half(0.5), one(1), two(2); //, eps(std::sqrt(ROL_EPSILON<Real>()));
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  Real alpha(1), sHs(0), alphaTmp(1), mmax(0), qmin(0), lambdaTmp(1);
  std::deque<Real> mqueue; mqueue.push_back(q);

  if (useNMSP_ && useMin_) { qmin = q; ymin.set(y); }

  // Compute initial projected gradient norm
  pwa1.set(gmod.dual());
  pwa.set(y); pwa.axpy(-one,pwa1);
  dproj(pwa,x,del,pwa2,pwa3,pwa4,pwa5,outStream);
  pwa.axpy(-one,y);
  Real gnorm = pwa.norm();
  const Real gtol = std::min(tol1_,tol2_*gnorm);

  // Compute initial step
  Real coeff  = one/gmod.norm();
  Real lambda = std::max(lambdaMin_,std::min(coeff,lambdaMax_));
  pwa.set(y); pwa.axpy(-lambda,pwa1);             // pwa = y - lambda gmod.dual()
  dproj(pwa,x,del,pwa2,pwa3,pwa4,pwa5,outStream); // pwa = P(y - lambda gmod.dual())
  pwa.axpy(-one,y);                               // pwa = P(y - lambda gmod.dual()) - y = step
  Real gs = gmod.apply(pwa);                      // gs  = <step, model gradient>
  Real ss = pwa.dot(pwa);                         // Norm squared of step

  if (verbosity_ > 1)
    outStream << "  Spectral Projected Gradient"          << std::endl;

  SPiter_ = 0;
  while (SPiter_ < maxit_) {
    SPiter_++;

    // Evaluate model Hessian
    model.hessVec(dwa,pwa,x,tol); nhess_++; // dwa = H step
    sHs = dwa.apply(pwa);                   // sHs = <step, H step>

    // Perform line search
    if (useNMSP_) { // Nonmonotone
      mmax     = *std::max_element(mqueue.begin(),mqueue.end());
      alphaTmp = (-(one-gamma_)*gs + std::sqrt(std::pow((one-gamma_)*gs,two)-two*sHs*(q-mmax)))/sHs;
    }
    else { // Exact
      alphaTmp = -gs/sHs;
    }
    alpha = (sHs > zero ? std::min(one,std::max(zero,alphaTmp)) : one);

    // Update model quantities
    q += alpha * (gs + half * alpha * sHs); // Update model value
    gmod.axpy(alpha,dwa);                   // Update model gradient
    y.axpy(alpha,pwa);                      // New iterate

    // Update nonmonotone line search information
    if (useNMSP_) {
      if (static_cast<int>(mqueue.size())==maxSize_) mqueue.pop_front();
      mqueue.push_back(q);
      if (useMin_ && q <= qmin) { qmin = q; ymin.set(y); }
    }

    // Compute projected gradient norm
    pwa1.set(gmod.dual());
    pwa.set(y); pwa.axpy(-one,pwa1);
    dproj(pwa,x,del,pwa2,pwa3,pwa4,pwa5,outStream);
    pwa.axpy(-one,y);
    gnorm = pwa.norm();

    if (verbosity_ > 1) {
      outStream << std::endl;
      outStream << "    Iterate:                          " << SPiter_ << std::endl;
      outStream << "    Spectral step length (lambda):    " << lambda  << std::endl;
      outStream << "    Step length (alpha):              " << alpha   << std::endl;
      outStream << "    Model decrease (pRed):            " << -q      << std::endl;
      outStream << "    Optimality criterion:             " << gnorm   << std::endl;
      outStream << std::endl;
    }
    if (gnorm < gtol) break;

    // Compute new spectral step
    //lambda = (sHs<=eps ? lambdaMax_ : std::max(lambdaMin_,std::min(ss/sHs,lambdaMax_)));
    lambdaTmp = (sHs == 0 ? coeff : ss/sHs);
    lambda = std::max(lambdaMin_,std::min(lambdaTmp,lambdaMax_));
    pwa.set(y); pwa.axpy(-lambda,pwa1);
    dproj(pwa,x,del,pwa2,pwa3,pwa4,pwa5,outStream);
    pwa.axpy(-one,y);
    gs = gmod.apply(pwa);
    ss = pwa.dot(pwa);
  }
  if (useNMSP_ && useMin_) { q = qmin; y.set(ymin); }
  SPflag_ = (SPiter_==maxit_) ? 1 : 0;
}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::dproj(Vector<Real> &x,
                                          const Vector<Real> &x0,
                                          Real del,
                                          Vector<Real> &y0,
                                          Vector<Real> &y1,
                                          Vector<Real> &yc,
                                          Vector<Real> &pwa,
                                          std::ostream &outStream) const {
  // Solve ||P(t*x0 + (1-t)*(x-x0))-x0|| = del using Brent's method
  const Real zero(0), half(0.5), one(1), two(2), three(3);
  const Real eps(ROL_EPSILON<Real>()), tol0(1e1*eps), fudge(1.0-1e-2*sqrt(eps));
  Real f0(0), f1(0), fc(0), t0(0), t1(1), tc(0), d1(1), d2(1), tol(1);
  Real p(0), q(0), r(0), s(0), m(0);
  int cnt(state_->nproj);
  y1.set(x);
  proj_->project(y1,outStream); state_->nproj++;
  pwa.set(y1); pwa.axpy(-one,x0);
  f1 = pwa.norm();
  if (f1 <= del) {
    x.set(y1);
    return;
  }
  y0.set(x0);
  tc = t0; fc = f0; yc.set(y0);
  d1 = t1-t0; d2 = d1;
  int code = 0;
  while (true) {
    if (std::abs(fc-del) < std::abs(f1-del)) {
      t0 = t1; t1 = tc; tc = t0;
      f0 = f1; f1 = fc; fc = f0;
      y0.set(y1); y1.set(yc); yc.set(y0);
    }
    tol = two*eps*std::abs(t1) + half*tol0;
    m   = half*(tc - t1);
    if (std::abs(m) <= tol) { code = 1; break; }
    if ((f1 >= fudge*del && f1 <= del)) break;
    if (std::abs(d1) < tol || std::abs(f0-del) <= std::abs(f1-del)) {
      d1 = m; d2 = d1;
    }
    else {
      s = (f1-del)/(f0-del);
      if (t0 == tc) {
        p = two*m*s;
        q = one-s;
      }
      else {
        q = (f0-del)/(fc-del);
        r = (f1-del)/(fc-del);
        p = s*(two*m*q*(q-r)-(t1-t0)*(r-one));
        q = (q-one)*(r-one)*(s-one);
      }
      if (p > zero) q = -q;
      else          p = -p;
      s  = d1;
      d1 = d2;
      if (two*p < three*m*q-std::abs(tol*q) && p < std::abs(half*s*q)) {
        d2 = p/q;
      }
      else {
        d1 = m; d2 = d1;
      }
    }
    t0 = t1; f0 = f1; y0.set(y1);
    if (std::abs(d2) > tol) t1 += d2;
    else if (m > zero)      t1 += tol;
    else                    t1 -= tol;
    y1.set(x); y1.scale(t1); y1.axpy(one-t1,x0);
    proj_->project(y1,outStream); state_->nproj++;
    pwa.set(y1); pwa.axpy(-one,x0);
    f1 = pwa.norm();
    if ((f1 > del && fc > del) || (f1 <= del && fc <= del)) {
      tc = t0; fc = f0; yc.set(y0);
      d1 = t1-t0; d2 = d1;
    }
  }
  if (code==1 && f1>del) x.set(yc);
  else                   x.set(y1);
  if (verbosity_ > 1) {
    outStream << std::endl;
    outStream << "  Trust-Region Subproblem Projection" << std::endl;
    outStream << "    Number of polyhedral projections: " << state_->nproj-cnt << std::endl;
    if (code == 1 && f1 > del) {
      outStream << "    Transformed Multiplier:           " << tc << std::endl;
      outStream << "    Dual Residual:                    " << fc-del << std::endl;
    }
    else {
      outStream << "    Transformed Multiplier:           " << t1 << std::endl;
      outStream << "    Dual Residual:                    " << f1-del << std::endl;
    }
    outStream << "    Exit Code:                        " << code << std::endl;
    outStream << std::endl;
  }
}

// BRACKETING AND BRENTS FOR UNTRANSFORMED MULTIPLIER
//template<typename Real>
//void TrustRegionSPGAlgorithm<Real>::dproj(Vector<Real> &x,
//                                          const Vector<Real> &x0,
//                                          Real del,
//                                          Vector<Real> &y0,
//                                          Vector<Real> &y1,
//                                          Vector<Real> &yc,
//                                          Vector<Real> &pwa,
//                                          std::ostream &outStream) const {
//  // Solve ||P(t*x0 + (1-t)*(x-x0))-x0|| = del using Brent's method
//  const Real zero(0), half(0.5), one(1), two(2), three(3);
//  const Real eps(ROL_EPSILON<Real>()), tol0(1e1*eps), fudge(1.0-1e-2*sqrt(eps));
//  Real f0(0), f1(0), fc(0), u0(0), u1(0), uc(0), t0(1), t1(0), tc(0), d1(1), d2(1), tol(1);
//  Real p(0), q(0), r(0), s(0), m(0);
//  int cnt(state_->nproj);
//  y0.set(x);
//  proj_->project(y0,outStream); state_->nproj++;
//  pwa.set(y0); pwa.axpy(-one,x0);
//  f0 = pwa.norm();
//  if (f0 <= del) {
//    x.set(y0);
//    return;
//  }
//
//  // Bracketing
//  t1 = static_cast<Real>(1e-1);
//  f1 = one+del;
//  while (f1 >= del) {
//    t1 *= static_cast<Real>(5e-2);
//    y1.set(x); y1.scale(t1); y1.axpy(one-t1,x0);
//    proj_->project(y1,outStream); state_->nproj++;
//    pwa.set(y1); pwa.axpy(-one,x0);
//    f1 = pwa.norm();
//  }
//  u1 = (one-t1)/t1;
//
//  // Brents
//  uc = u0; tc = t0; fc = f0; yc.set(y0);
//  d1 = u1-u0; d2 = d1;
//  int code = 0;
//  while (true) {
//    if (std::abs(fc-del) < std::abs(f1-del)) {
//      u0 = u1; u1 = uc; uc = u0;
//      t0 = t1; t1 = tc; tc = t0;
//      f0 = f1; f1 = fc; fc = f0;
//      y0.set(y1); y1.set(yc); yc.set(y0);
//    }
//    tol = two*eps*abs(u1) + half*tol0;
//    m   = half*(uc - u1);
//    if (std::abs(m) <= tol) { code = 1; break; }
//    if ((f1 >= fudge*del && f1 <= del)) break;
//    if (std::abs(d1) < tol || std::abs(f0-del) <= std::abs(f1-del)) {
//      d1 = m; d2 = d1;
//    }
//    else {
//      s = (f1-del)/(f0-del);
//      if (u0 == uc) {
//        p = two*m*s;
//        q = one-s;
//      }
//      else {
//        q = (f0-del)/(fc-del);
//        r = (f1-del)/(fc-del);
//        p = s*(two*m*q*(q-r)-(u1-u0)*(r-one));
//        q = (q-one)*(r-one)*(s-one);
//      }
//      if (p > zero) q = -q;
//      else          p = -p;
//      s  = d1;
//      d1 = d2;
//      if (two*p < three*m*q-std::abs(tol*q) && p < std::abs(half*s*q)) {
//        d2 = p/q;
//      }
//      else {
//        d1 = m; d2 = d1;
//      }
//    }
//    u0 = u1; t0 = t1; f0 = f1; y0.set(y1);
//    if (std::abs(d2) > tol) u1 += d2;
//    else if (m > zero)      u1 += tol;
//    else                    u1 -= tol;
//    t1 = one/(one+u1);
//    y1.set(x); y1.scale(t1); y1.axpy(one-t1,x0);
//    proj_->project(y1,outStream); state_->nproj++;
//    pwa.set(y1); pwa.axpy(-one,x0);
//    f1 = pwa.norm();
//    if ((f1 > del && fc > del) || (f1 <= del && fc <= del)) {
//      uc = u0; tc = t0; fc = f0; yc.set(y0);
//      d1 = u1-u0; d2 = d1;
//    }
//  }
//  if (code==1 && f1>del) x.set(yc);
//  else                   x.set(y1);
//  if (verbosity_ > 1) {
//    outStream << std::endl;
//    outStream << "  Trust-Region Subproblem Projection" << std::endl;
//    outStream << "    Number of polyhedral projections: " << state_->nproj-cnt << std::endl;
//    if (code == 1 && f1 > del) {
//      outStream << "    Multiplier:                       " << uc << std::endl;
//      outStream << "    Transformed Multiplier:           " << tc << std::endl;
//      outStream << "    Dual Residual:                    " << fc-del << std::endl;
//    }
//    else {
//      outStream << "    Multiplier:                       " << u1 << std::endl;
//      outStream << "    Transformed Multiplier:           " << t1 << std::endl;
//      outStream << "    Dual Residual:                    " << f1-del << std::endl;
//    }
//    outStream << "    Exit Code:                        " << code << std::endl;
//    outStream << std::endl;
//  }
//}

// RIDDERS' METHOD FOR TRUST-REGION PROJECTION
//template<typename Real>
//void TrustRegionSPGAlgorithm<Real>::dproj(Vector<Real> &x,
//                                          const Vector<Real> &x0,
//                                          Real del,
//                                          Vector<Real> &y,
//                                          Vector<Real> &y1,
//                                          Vector<Real> &yc,
//                                          Vector<Real> &p,
//                                          std::ostream &outStream) const {
//  // Solve ||P(t*x0 + (1-t)*(x-x0))-x0|| = del using Ridder's method
//  const Real half(0.5), one(1), tol(1e1*ROL_EPSILON<Real>());
//  const Real fudge(1.0-1e-2*std::sqrt(ROL_EPSILON<Real>()));
//  Real e0(0), e1(0), e2(0), e(0), a0(0), a1(0.5), a2(1), a(0);
//  int cnt(state_->nproj);
//  y.set(x);
//  proj_->project(y,outStream); state_->nproj++;
//  p.set(y); p.axpy(-one,x0);
//  e2 = p.norm();
//  if (e2 <= del) {
//    x.set(y);
//    return;
//  }
//  bool code = 1;
//  while (a2-a0 > tol) {
//    a1 = half*(a0+a2);
//    y.set(x); y.scale(a1); y.axpy(one-a1,x0);
//    proj_->project(y,outStream); state_->nproj++;
//    p.set(y); p.axpy(-one,x0);
//    e1 = p.norm();
//    if (e1 >= fudge*del && e1 <= del) break;
//    a = a1-(a1-a0)*(e1-del)/std::sqrt((e1-del)*(e1-del)-(e0-del)*(e2-del));
//    y.set(x); y.scale(a); y.axpy(one-a,x0);
//    proj_->project(y,outStream); state_->nproj++;
//    p.set(y); p.axpy(-one,x0);
//    e = p.norm();
//    if (e < fudge*del) {
//      if (e1 < fudge*del) { e0 = (a < a1 ? e1 : e); a0 = (a < a1 ? a1 : a); }
//      else                { e0 = e; a0 = a; e2 = e1; a2 = a1; };
//    }
//    else if (e > del) {
//      if (e1 < fudge*del) { e0 = e1; a0 = a1; e2 = e; a2 = a; }
//      else                { e2 = (a < a1 ? e : e1); a2 = (a < a1 ? a : a1); }
//    }
//    else {
//      code = 0;
//      break; // Exit if fudge*del <= snorm <= del
//    }
//  }
//  x.set(y);
//  if (verbosity_ > 1) {
//    outStream << std::endl;
//    outStream << "  Trust-Region Subproblem Projection" << std::endl;
//    outStream << "    Number of polyhedral projections: " << state_->nproj-cnt << std::endl;
//    outStream << "    Transformed Multiplier:           " << a1 << std::endl;
//    outStream << "    Dual Residual:                    " << e1-del << std::endl;
//    outStream << "    Exit Code:                        " << code << std::endl;
//    outStream << std::endl;
//  }
//}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(114,'-') << std::endl;
    os << " SPG trust-region method status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)" << std::endl;
    os << "  value   - Objective function value" << std::endl; 
    os << "  gnorm   - Norm of the gradient" << std::endl;
    os << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    os << "  delta   - Trust-Region radius" << std::endl;
    os << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    os << "  #grad   - Number of times the gradient was computed" << std::endl;
    os << "  #hess   - Number of times the Hessian was applied" << std::endl;
    os << "  #proj   - Number of times the projection was computed" << std::endl;
    os << std::endl;
    os << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      os << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    os << std::endl;
    os << "  iterSPG - Number of Spectral Projected Gradient iterations" << std::endl << std::endl;
    os << "  flagSPG - Trust-Region Truncated CG flag" << std::endl;
    os << "    0 - Converged" << std::endl;
    os << "    1 - Iteration Limit Exceeded" << std::endl;
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
  os << std::setw(10) << std::left << "iterSPG";
  os << std::setw(10) << std::left << "flagSPG";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "SPG Trust-Region Method (Type B, Bound Constraints)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void TrustRegionSPGAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
