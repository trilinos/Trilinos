// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_TRUSTREGIONALGORITHM_DEF_HPP
#define ROL_TYPEP_TRUSTREGIONALGORITHM_DEF_HPP

#include <deque>

namespace ROL {
namespace TypeP {


template<typename Real>
TrustRegionAlgorithm<Real>::TrustRegionAlgorithm(ParameterList &list,
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
  verbosity_ = trlist.sublist("General").get("Output Level",      0);
  initProx_  = trlist.get("Apply Prox to Initial Guess",          false);
  t0_        = list.sublist("Status Test").get("Proximal Gradient Parameter", 1.0);
  // Nonmonotone Parameters
  storageNM_ = trlist.get("Nonmonotone Storage Size",             0);
  useNM_     = (storageNM_ <= 0 ? false : true);
  // Algorithm-Specific Parameters
  ROL::ParameterList &lmlist = trlist.sublist("TRN");
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
  // Subsolver (general) parameters
  lambdaMin_ = lmlist.sublist("Solver").get("Minimum Spectral Step Size",    1e-8);
  lambdaMax_ = lmlist.sublist("Solver").get("Maximum Spectral Step Size",    1e8);
  gamma_     = lmlist.sublist("Solver").get("Sufficient Decrease Tolerance", 1e-4);
  maxSize_   = lmlist.sublist("Solver").get("Maximum Storage Size",          10);
  maxit_     = lmlist.sublist("Solver").get("Iteration Limit",               25);
  tol1_      = lmlist.sublist("Solver").get("Absolute Tolerance",            1e-4);
  tol2_      = lmlist.sublist("Solver").get("Relative Tolerance",            1e-2);
  // Subsolver (spectral projected gradient) parameters
  useMin_    = lmlist.sublist("Solver").get("Use Smallest Model Iterate", true);
  useNMSP_   = lmlist.sublist("Solver").get("Use Nonmonotone Search",     false);
  std::string ssname = lmlist.sublist("Solver").get("Subproblem Solver", "SPG");
  algSelect_ = StringToETrustRegionP(ssname);
  // Subsolver (nonlinear conjugate gradient) parameters)
  ncgType_   = lmlist.sublist("Solver").sublist("NCG").get("Nonlinear CG Type",              4);
  etaNCG_    = lmlist.sublist("Solver").sublist("NCG").get("Truncation Parameter for HZ CG", 1e-2);
  desPar_    = lmlist.sublist("Solver").sublist("NCG").get("Descent Parameter",              0.2);
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
void TrustRegionAlgorithm<Real>::initialize(Vector<Real>          &x,
                                            const Vector<Real>    &g,
                                            Real                   ftol,
                                            Objective<Real>       &sobj,
                                            Objective<Real>       &nobj,
                                            Vector<Real>          &px,
                                            Vector<Real>          &dg,
                                            std::ostream          &outStream) {
  // Initialize data
  TypeP::Algorithm<Real>::initialize(x,g);
  nhess_ = 0;
  // Update approximate gradient and approximate objective function.
  if (initProx_){
    nobj.prox(*state_->iterateVec,x,t0_, ftol); state_->nprox++;
    x.set(*state_->iterateVec);
  }
  sobj.update(x,UpdateType::Initial,state_->iter);
  state_->svalue = sobj.value(x,ftol); state_->nsval++;
  nobj.update(x, UpdateType::Initial,state_->iter);
  state_->nvalue = nobj.value(x,ftol); state_->nnval++;
  state_->value = state_->svalue + state_->nvalue;
  computeGradient(x,*state_->gradientVec,px,dg,*state_->stepVec,state_->searchSize,sobj,nobj,true,gtol_,state_->gnorm,outStream);

  state_->snorm = ROL_INF<Real>();
  // Normalize initial CP step length
  if (normAlpha_) alpha_ /= state_->gradientVec->norm();//change here?
  // Compute initial trust region radius if desired.
  if ( state_->searchSize <= static_cast<Real>(0) )
    state_->searchSize = state_->gradientVec->norm();
  SPiter_ = 0;
  SPflag_ = 0;
}

template<typename Real>
Real TrustRegionAlgorithm<Real>::computeValue(Real inTol,
                                              Real &outTol,
                                              Real pRed,
                                              Real &fold,
                                              int iter,
                                              const Vector<Real> &x,
                                              const Vector<Real> &xold,
                                              Objective<Real>    &sobj) {
  outTol = std::sqrt(ROL_EPSILON<Real>());
  if ( useInexact_[0] ) {
    if (!(iter%updateIter_) && (iter!=0)) force_ *= forceFactor_;
    const Real one(1);
    Real eta = static_cast<Real>(0.999)*std::min(eta1_,one-eta2_);
    outTol   = scale_*std::pow(eta*std::min(pRed,force_),one/omega_);
    if (inTol > outTol) {
      fold = sobj.value(xold,outTol); state_->nsval++;
    }
  }
  // Evaluate objective function at new iterate
  sobj.update(x,UpdateType::Trial);
  Real fval = sobj.value(x,outTol); state_->nsval++;
  return fval;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::computeGradient(const Vector<Real> &x,
                                                 Vector<Real> &g,
                                                 Vector<Real> &px,
                                                 Vector<Real> &dg,
                                                 Vector<Real> &step,
                                                 Real del,
                                                 Objective<Real> &sobj,
                                                 Objective<Real> &nobj,
                                                 bool accept,
                                                 Real &gtol,
                                                 Real &gnorm,
                                                 std::ostream &outStream) const {
  if ( useInexact_[1] ) {
    Real gtol0 = scale0_*del;
    if (accept) gtol  = gtol0 + static_cast<Real>(1);
    else        gtol0 = scale0_*std::min(gnorm,del);
    while ( gtol > gtol0 ) {
      gtol = gtol0;
      sobj.gradient(g,x,gtol); state_->ngrad++;
      dg.set(g.dual());
      pgstep(px, step, nobj, x, dg, t0_, gtol0); // change gtol? one or ocScale?
      gnorm = step.norm() / t0_;
      gtol0 = scale0_*std::min(gnorm,del);
    }
  }
  else {
    if (accept) {
      gtol = std::sqrt(ROL_EPSILON<Real>());
      sobj.gradient(g,x,gtol); state_->ngrad++;
      dg.set(g.dual());
      pgstep(px, step, nobj, x, dg, t0_, gtol);
      gnorm = step.norm() / t0_;
    }
  }
}

template<typename Real>
void TrustRegionAlgorithm<Real>::run(Vector<Real>          &x,
                                     const Vector<Real>    &g,
                                     Objective<Real>       &sobj,
                                     Objective<Real>       &nobj,
                                     std::ostream          &outStream ) {
  const Real zero(0), one(1);
  //Real tol0 = std::sqrt(ROL_EPSILON<Real>());
  Real inTol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(), outTol(inTol);
  Real strial(0), ntrial(0), smodel(0), Ftrial(0), pRed(0), rho(1);
  // Initialize trust-region data
  std::vector<std::string> output;
  Ptr<Vector<Real>> gmod = g.clone();
  Ptr<Vector<Real>> px   = x.clone();
  Ptr<Vector<Real>> dg   = x.clone();
  // Initialize Algorithm
  initialize(x,g,inTol,sobj,nobj, *px, *dg, outStream);
  // Initialize storage vectors
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
    model_->setData(sobj,*state_->iterateVec,*state_->gradientVec,gtol_);

    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
    //q = state_->svalue + state_->nvalue;//q is no longer used
    gmod->set(*state_->gradientVec);
    smodel = state_->svalue;
    ntrial = state_->nvalue;
    switch (algSelect_) {
      case TRUSTREGION_P_SPG:
      default:
        // Compute Cauchy point (TRON notation: x = x[1])
        dcauchy(*state_->stepVec,alpha_, smodel, ntrial,
                *state_->iterateVec, *dg, state_->searchSize,
                *model_, nobj, *px, *dwa1, *dwa2, outStream); // Solve 1D optimization problem for alpha
        x.plus(*state_->stepVec);               // Set x = x[0] + alpha*g
        // Model gradient at s = x[1] - x[0]
        gmod->plus(*dwa1); // hessVec from Cauchy point computation

        // Apply SPG starting from the Cauchy point->change input
        dspg(x,smodel,ntrial,*gmod,*state_->iterateVec,state_->searchSize,
             *model_,nobj,*pwa1,*pwa2,*pwa3,*pwa4,*pwa5,*pwa6,*pwa7,*dwa1,
             outStream);
        pRed = state_->value - (smodel+ntrial);
        break;
      case TRUSTREGION_P_SPG2:
        dspg2(x,smodel, ntrial, pRed, *gmod, *state_->iterateVec,
              state_->searchSize, *model_, nobj,
              *pwa1, *pwa2, *px, *dwa1, outStream);
        break;
      case TRUSTREGION_P_NCG:
        dncg(x,smodel,ntrial,*gmod,*state_->iterateVec,state_->searchSize,
             *model_,nobj,*pwa1,*pwa2,*pwa3,*pwa4,*pwa5,*pwa6,*dwa1,
             outStream);
        pRed = state_->value - (smodel+ntrial);
        break;
    }

    // Update storage and compute predicted reduction
    //pRed = -q; // now updated in dcauchy/dspg
    state_->stepVec->set(x); state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Compute trial objective value
    strial = computeValue(inTol,outTol,pRed,state_->svalue,state_->iter,x,*state_->iterateVec,sobj);
    Ftrial = strial + ntrial;

    // Compute ratio of actual and predicted reduction
    TRflag_ = TRUtils::SUCCESS;
    TRUtils::analyzeRatio<Real>(rho,TRflag_,state_->value,Ftrial,pRed,eps_,outStream,verbosity_>1);
    if (useNM_) {
      TRUtils::analyzeRatio<Real>(rhoNM,TRflagNM,fr,Ftrial,pRed+sigmar,eps_,outStream,verbosity_>1);
      TRflag_ = (rho < rhoNM ? TRflagNM : TRflag_);
      rho     = (rho < rhoNM ?    rhoNM :    rho );
    }

    // Update algorithm state
    state_->iter++;
    // Accept/reject step and update trust region radius
    if ((rho < eta0_ && TRflag_ == TRUtils::SUCCESS) || (TRflag_ >= 2)) { // Step Rejected
      x.set(*state_->iterateVec);
      sobj.update(x,UpdateType::Revert,state_->iter);
      nobj.update(x,UpdateType::Revert,state_->iter);
      if (interpRad_ && (rho < zero && TRflag_ != TRUtils::TRNAN)) {
        // Negative reduction, interpolate to find new trust-region radius
        state_->searchSize = TRUtils::interpolateRadius<Real>(*state_->gradientVec,*state_->stepVec,
          state_->snorm,pRed,state_->value,Ftrial,state_->searchSize,gamma0_,gamma1_,eta2_,
          outStream,verbosity_>1);
      }
      else { // Shrink trust-region radius
        state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
      }
      computeGradient(x,*state_->gradientVec,*px,*dg,*pwa1,state_->searchSize,sobj,nobj,false,gtol_,state_->gnorm,outStream);
    }
    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
             || (TRflag_ == TRUtils::POSPREDNEG)) { // Step Accepted
      state_->value  = Ftrial;
      state_->svalue = strial;
      state_->nvalue = ntrial;
      sobj.update(x,UpdateType::Accept,state_->iter);
      nobj.update(x,UpdateType::Accept,state_->iter);
      inTol = outTol;
      if (useNM_) {
        sigmac += pRed; sigmar += pRed;
        if (Ftrial < fmin) { fmin = Ftrial; fc = fmin; sigmac = zero; L = 0; }
        else {
          L++;
          if (Ftrial > fc)     { fc = Ftrial; sigmac = zero;   }
          if (L == storageNM_) { fr = fc;     sigmar = sigmac; }
        }
      }
      // Increase trust-region radius
      if (rho >= eta2_) state_->searchSize = std::min(gamma2_*state_->searchSize, delMax_);
      // Compute gradient at new iterate
      dwa1->set(*state_->gradientVec);
      computeGradient(x,*state_->gradientVec,*px,*dg,*pwa1,state_->searchSize,sobj,nobj,true,gtol_,state_->gnorm,outStream);
      state_->iterateVec->set(x);
      // Update secant information in trust-region model
      model_->update(x,*state_->stepVec,*dwa1,*state_->gradientVec,
                     state_->snorm,state_->iter);
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeP::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
Real TrustRegionAlgorithm<Real>::dcauchy(Vector<Real> &s,
                                         Real &alpha,
                                         Real &sval,
                                         Real &nval,
                                         const Vector<Real> &x,
                                         const Vector<Real> &g,
                                         const Real del,
                                         TrustRegionModel_U<Real> &model,
                                         Objective<Real> &nobj,
                                         Vector<Real> &px,
                                         Vector<Real> &dwa,
                                         Vector<Real> &dwa1,
                                         std::ostream &outStream) {
  const Real half(0.5), sold(sval), nold(nval);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  bool interp = false;
  Real gs(0), snorm(0), Qk(0), pRed(0);
  // Compute s = P(x[0] - alpha g[0]) - x[0]
  pgstep(px, s, nobj, x, g, alpha, tol);
  snorm = s.norm();
  if (snorm > del) {
    interp = true;
  }
  else {
    model.hessVec(dwa,s,x,tol); nhess_++;
    nobj.update(px, UpdateType::Trial);
    nval   = nobj.value(px, tol); state_->nnval++;
    gs     = s.dot(g);
    sval   = sold + gs + half * s.apply(dwa);
    pRed   = (sold + nold) - (sval + nval);
    Qk     = gs + nval - nold;
    interp = (pRed < -mu0_*Qk);
  }
  // Either increase or decrease alpha to find approximate Cauchy point
  int cnt = 0;
  if (interp) {//decrease loop
    bool search = true;
    while (search) {
      alpha *= interpf_;
      pgstep(px, s, nobj, x, g, alpha, tol);
      snorm = s.norm();
      if (snorm <= del) {
        model.hessVec(dwa,s,x,tol); nhess_++;
        nobj.update(px, UpdateType::Trial);
        nval   = nobj.value(px, tol); state_->nnval++;
        gs     = s.dot(g);
        sval   = sold + gs + half * s.apply(dwa);
        pRed   = (sold + nold) - (sval + nval);
        Qk     = gs + nval - nold;
        search = ((pRed < -mu0_*Qk) && (cnt < redlim_)) ;
      }
      cnt++;
    }
  }
  else {
    bool search = true;
    Real alphas = alpha;
    Real mvals  = pRed;
    Real svals  = sval;
    dwa1.set(dwa);
    while (search) {
      alpha *= extrapf_;
      pgstep(px, s, nobj, x, g, alpha, tol);
      snorm = s.norm();
      if (snorm <= del && cnt < explim_){// && mnew < mold + mu0_*Qk) {
        model.hessVec(dwa,s,x,tol); nhess_++;
        nobj.update(px, UpdateType::Trial);
        nval = nobj.value(px, tol); state_->nnval++;
        gs   = s.dot(g);
        sval = sold + gs + half * s.apply(dwa);
        pRed = (sold + nold) - (sval + nval);
        Qk   = gs + nval - nold;
        if (pRed >= -mu0_*Qk && std::abs(pRed-mvals) > qtol_*std::abs(mvals)) {
          dwa1.set(dwa);
          alphas = alpha;
          mvals  = pRed;
          svals  = sval;
          search = true;
        }
        else {
          dwa.set(dwa1);
          pRed   = mvals;
          sval   = svals;
          search = false;
        }
      }
      else {
        search = false;
      }
      cnt++;
    }
    alpha = alphas;
    pgstep(px, s, nobj, x, g, alpha, tol);
    snorm = s.norm();
  }
  if (verbosity_ > 1) {
    outStream << "    Cauchy point"                         << std::endl;
    outStream << "    Step length (alpha):              " << alpha << std::endl;
    outStream << "    Step length (alpha*g):            " << snorm << std::endl;
    outStream << "    Model decrease (pRed):            " << pRed  << std::endl;
    if (!interp)
      outStream << "    Number of extrapolation steps:    " << cnt << std::endl;
  }
  return snorm;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::dspg2(Vector<Real>             &y,
                                       Real                     &sval,
                                       Real                     &nval,
                                       Real                     &pRed,
                                       Vector<Real>             &gmod,
                                       const Vector<Real>       &x,
                                       Real                      del,
                                       TrustRegionModel_U<Real> &model,
                                       Objective<Real>          &nobj,
                                       Vector<Real>             &pwa,
                                       Vector<Real>             &pwa1,
                                       Vector<Real>             &pwa2,
                                       Vector<Real>             &dwa,
                                       std::ostream             &outStream) {
  // Use SPG to approximately solve TR subproblem:
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)>  subject to y\in C, ||y|| \le del
  //
  //   Inpute:
  //       y = Primal vector
  //       x = Current iterate
  //       g = Current gradient
  const Real half(0.5), one(1), safeguard(1e2*ROL_EPSILON<Real>());
  const Real mprev(sval+nval);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  Real coeff(1), alpha(1), alphaMax(1), lambda(1), lambdaTmp(1);
  Real gs(0), ss(0), gnorm(0), s0s0(0), ss0(0), sHs(0), snorm(0), nold(nval);
  pwa1.zero();

  // Set y = x
  y.set(x);

  // Compute initial step
  coeff  = one / gmod.norm();
  lambda = std::max(lambdaMin_,std::min(coeff,lambdaMax_));
  pgstep(pwa2, pwa, nobj, y, gmod.dual(), lambda, tol);
  gs     = gmod.apply(pwa);                      // gs  = <step, model gradient>
  ss     = pwa.dot(pwa);                         // Norm squared of step
  snorm  = std::sqrt(ss);                        // norm(step)
  gnorm  = snorm / lambda;                       // norm(step) / lambda

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

    // Evaluate nonsmooth term
    nobj.update(pwa2,UpdateType::Trial);
    nval = nobj.value(pwa2,tol); state_->nnval++;

    // Perform line search
    alphaMax = one;
    if (snorm >= del-safeguard) { // Trust-region constraint is violated
      ss0      = pwa1.dot(pwa);
      alphaMax = std::min(one, (-ss0 + std::sqrt(ss0*ss0 - ss*(s0s0-del*del)))/ss);
    }
    alpha = (sHs <= safeguard) ? alphaMax : std::min(alphaMax, -(gs + nval - nold)/sHs);

    // Update model quantities
    if (alpha == one) {
      y.set(pwa2);
      sval += gs + half * sHs;
      gmod.plus(dwa);
    }
    else {
      y.axpy(alpha,pwa); // New iterate
      nobj.update(y,UpdateType::Trial);
      nval  = nobj.value(y, tol); state_->nnval++;
      sval += alpha * (gs + half * alpha * sHs);
      gmod.axpy(alpha,dwa);
    }
    nold = nval;
    pRed = mprev - (sval+nval);

    // Check trust-region constraint violation
    pwa1.set(y); pwa1.axpy(-one,x);
    s0s0 = pwa1.dot(pwa1);
    snorm = std::sqrt(s0s0);

    if (verbosity_ > 1) {
      outStream << std::endl;
      outStream << "    Iterate:                          " << SPiter_ << std::endl;
      outStream << "    Spectral step length (lambda):    " << lambda  << std::endl;
      outStream << "    Step length (alpha):              " << alpha   << std::endl;
      outStream << "    Model decrease (pRed):            " << pRed    << std::endl;
      outStream << "    Optimality criterion:             " << gnorm   << std::endl;
      outStream << "    Step norm:                        " << snorm   << std::endl;
      outStream << std::endl;
    }

    if (snorm >= del - safeguard) { SPflag_ = 2; break; }

    // Compute new spectral step
    lambdaTmp = (sHs <= safeguard) ? one/gmod.norm() : ss/sHs;
    lambda    = std::max(lambdaMin_,std::min(lambdaTmp,lambdaMax_));

    pgstep(pwa2, pwa, nobj, y, gmod.dual(), alpha, tol); // pass pwa by reference? *pwa?
    gs    = gmod.apply(pwa);
    ss    = pwa.dot(pwa);
    gnorm = std::sqrt(ss) / lambda;

    if (gnorm <= gtol) { SPflag_ = 0; break; }
  }
  SPflag_ = (SPiter_==maxit_) ? 1 : SPflag_;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::dspg(Vector<Real>             &y,
                                      Real                     &sval,
                                      Real                     &nval,
                                      Vector<Real>             &gmod,
                                      const Vector<Real>       &x,
                                      Real                      del,
                                      TrustRegionModel_U<Real> &model,
                                      Objective<Real>          &nobj,
                                      Vector<Real>             &ymin,
                                      Vector<Real>             &pwa,
                                      Vector<Real>             &pwa1,
                                      Vector<Real>             &pwa2,
                                      Vector<Real>             &pwa3,
                                      Vector<Real>             &pwa4,
                                      Vector<Real>             &pwa5,
                                      Vector<Real>             &dwa,
                                      std::ostream             &outStream) {
  // Use SPG to approximately solve TR subproblem:
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)> + phi(y)  subject to  ||y|| \le del
  //
  //   Inpute:
  //       y = Cauchy step
  //       x = Current iterate
  //       g = Current gradient
  const Real half(0.5), one(1), safeguard(1e2*ROL_EPSILON<Real>());
  const Real mval(sval+nval);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  Real mcomp(0), mval_min(0), sval_min(0), nval_min(0);
  Real alpha(1), coeff(1), lambda(1), lambdaTmp(1);
  Real snew(sval), nnew(nval), mnew(mval);
  Real sold(sval), nold(nval), mold(mval);
  Real sHs(0), ss(0), gs(0), Qk(0), gnorm(0);
  std::deque<Real> mqueue; mqueue.push_back(mold);

  if (useNMSP_ && useMin_) {
    mval_min = mval; sval_min = sval; nval_min = nval; ymin.set(y);
  }

  // Compute initial proximal gradient norm
  pwa1.set(gmod.dual());
  pwa.set(y); pwa.axpy(-t0_,pwa1);
  dprox(pwa,x,t0_,del,nobj,pwa2,pwa3,pwa4,pwa5,outStream);
  pwa.axpy(-one,y);
  gnorm = pwa.norm() / t0_;
  const Real gtol = std::min(tol1_,tol2_*gnorm);

  // Compute initial spectral step size
  coeff  = one / gmod.norm();
  lambda = std::max(lambdaMin_,std::min(coeff,lambdaMax_));

  if (verbosity_ > 1)
    outStream << "  Spectral Projected Gradient"          << std::endl;

  SPiter_ = 0;
  while (SPiter_ < maxit_) {
    SPiter_++;

    // Compuate SPG step
    alpha = one;
    pwa.set(y); pwa.axpy(-lambda,pwa1);                         // pwa = y - lambda gmod.dual()
    dprox(pwa,x,lambda,del,nobj,pwa2,pwa3,pwa4,pwa5,outStream); // pwa = P(y - lambda gmod.dual())
    pwa2.set(pwa);                                              // pwa2 = P(y - lambda gmod.dual())
    pwa.axpy(-one,y);                                           // pwa = P(y - lambda gmod.dual()) - y = step
    ss = pwa.dot(pwa);                                          // Norm squared of step

    // Evaluate model Hessian
    model.hessVec(dwa,pwa,x,tol); nhess_++; // dwa = H step
    nobj.update(pwa2, UpdateType::Trial);
    nnew  = nobj.value(pwa2, tol); state_->nnval++;
    sHs   = dwa.apply(pwa);                 // sHs = <step, H step>
    gs    = gmod.apply(pwa);                // gs  = <step, model gradient>
    snew  = half * sHs + gs + sold;
    mnew  = snew + nnew;
    Qk    = gs + nnew - nold;

    // Perform line search
    //mcomp = useNMSP_ ? *std::max_element(mqueue.begin(),mqueue.end()) : mold;
    //while( mnew > mcomp + mu0_*Qk ) {
    //  alpha *= interpf_;
    //  pwa2.set(y); pwa2.axpy(alpha,pwa);
    //  nobj.update(pwa2, UpdateType::Trial);
    //  nnew  = nobj.value(pwa2, tol); state_->nnval++;
    //  snew  = half * alpha * alpha * sHs + alpha * gs + sold;
    //  mnew  = nnew + snew;
    //  Qk    = alpha * gs + nnew - nold;
    //}
    if (useNMSP_) { // Nonmonotone
      mcomp = *std::max_element(mqueue.begin(),mqueue.end());
      while( mnew > mcomp + mu0_*Qk ) {
        alpha *= interpf_;
        pwa2.set(y); pwa2.axpy(alpha,pwa);
        nobj.update(pwa2, UpdateType::Trial);
        nnew  = nobj.value(pwa2, tol); state_->nnval++;
        snew  = half * alpha * alpha * sHs + alpha * gs + sold;
        mnew  = nnew + snew;
        Qk    = alpha * gs + nnew - nold;
      }
    }
    else {
      alpha = (sHs <= safeguard) ? one : std::min(one,-(gs + nnew - nold)/sHs);
    }

    // Update model quantities
    y.set(pwa2);
    sold = snew;
    nold = nnew;
    mold = mnew;
    gmod.axpy(alpha,dwa);                      // Update model gradient
    nobj.update(y, UpdateType::Accept);

    // Update nonmonotone line search information
    if (useNMSP_) {
      if (static_cast<int>(mqueue.size())==maxSize_) mqueue.pop_front();
      mqueue.push_back(sval+nval);
      if (useMin_ && mval <= mval_min) {
        mval_min = mval; sval_min = sval; nval_min = nval; ymin.set(y);
      }
    }

    // Compute projected gradient norm
    pwa1.set(gmod.dual());
    pwa.set(y); pwa.axpy(-t0_,pwa1);
    dprox(pwa,x,t0_,del,nobj,pwa2,pwa3,pwa4,pwa5,outStream);
    pwa.axpy(-one,y);
    gnorm = pwa.norm() / t0_;

    if (verbosity_ > 1) {
      outStream << std::endl;
      outStream << "    Iterate:                          " << SPiter_   << std::endl;
      outStream << "    Spectral step length (lambda):    " << lambda    << std::endl;
      outStream << "    Step length (alpha):              " << alpha     << std::endl;
      outStream << "    Model decrease (pRed):            " << mval-mold << std::endl;
      outStream << "    Optimality criterion:             " << gnorm     << std::endl;
      outStream << std::endl;
    }
    if (gnorm < gtol) break;

    // Compute new spectral step
    lambdaTmp = (sHs == 0 ? coeff : ss / sHs);
    lambda    = std::max(lambdaMin_, std::min(lambdaTmp, lambdaMax_));
  }
  if (useNMSP_ && useMin_) {
    sval = sval_min; nval = nval_min; y.set(ymin);
  }
  else {
    sval = sold; nval = nold;
  }
  SPflag_ = (SPiter_==maxit_) ? 1 : 0;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::dprox(Vector<Real> &x,
                                       const Vector<Real> &x0,
                                       Real t,
                                       Real del,
                                       Objective<Real> &nobj,
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
  int cnt(state_->nprox);
  nobj.prox(y1, x, t, tol); state_->nprox++;
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
    pwa.set(x); pwa.scale(t1); pwa.axpy(one-t1,x0);
    nobj.prox(y1, pwa, t1*t, tol); state_->nprox++;
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
    outStream << "  Trust-Region Subproblem Proximity Operator" << std::endl;
    outStream << "    Number of proxes:                 " << state_->nprox-cnt << std::endl;
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

template<typename Real>
void TrustRegionAlgorithm<Real>::dbls(Real &alpha, Real &nval, Real &pred,
                                      const Vector<Real> &y,
                                      const Vector<Real> &s,
                                      Real lambda, Real tmax,
                                      Real kappa, Real gs,
                                      Objective<Real> &nobj,
                                      Vector<Real> &pwa) {
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  const Real eps(1e-2*std::sqrt(ROL_EPSILON<Real>())), tol0(1e4*tol);
  const Real eps0(1e2*ROL_EPSILON<Real>());
  const unsigned maxit(50);
  const Real zero(0), half(0.5), one(1), two(2);
  const Real lam(0.5*(3.0-std::sqrt(5.0)));
  const Real nold(nval);
  const Real mu(1e-4);

  // Evaluate model at initial left end point (t = 0)
  Real tL(0), pL(0);
  // Evaluate model at right end point (t = tmax)
  Real tR = tmax;
  pwa.set(y); pwa.axpy(tR, s);
  nobj.update(pwa,UpdateType::Trial);
  Real nR = nobj.value(pwa,tol); state_->nnval++;
  Real pR = tR * (half * tR * kappa + gs) + nR - nold;

  // Compute minimizer of quadratic upper bound
  Real t0 = tR, n0 = nR;
  if (tmax > lambda) {
    t0 = lambda;
    pwa.set(y); pwa.axpy(t0, s);
    nobj.update(pwa,UpdateType::Trial);
    n0 = nobj.value(pwa,tol); state_->nnval++;
  }
  Real ts = (kappa > 0) ? std::min(tR,(((nold - n0) / kappa) / t0) - (gs / kappa)) : tR;
  Real t  = std::min(t0,ts);
  bool useOptT = true;
  if (t <= tL) {
    t = (t0 < tR ? t0 : half*tR);
    useOptT = false;
  }

  // Evaluate model at t (t = minimizer of quadratic upper bound)
  pwa.set(y); pwa.axpy(t, s);
  nobj.update(pwa,UpdateType::Trial);
  Real nt = nobj.value(pwa,tol); state_->nnval++;
  Real Qt = t * gs + nt - nold, Qu = Qt;
  Real pt = half * t * t * kappa + Qt;

  // If phi(x) = phi(x+tmax s) = phi(x+ts), then
  // phi(x+ts) is constant for all t, so the TR
  // model is quadratic---use the minimizer
  if (useOptT && nt == nold && nR == nold) {
    alpha = ts;
    pred  = ts * (half * ts * kappa + gs);
    nval  = nold;
    return;
  }

  // If pt >= max(pL, pR), then the model is concave
  // and the minimum is obtained at the end points
  if (pt >= std::max(pL, pR)) {
    alpha = tR;
    pred  = pR;
    nval  = nR;
    return;
  }

  // Run Brent's method (quadratic interpolation + golden section)
  // to minimize m_k(x+ts) with respect to t
  Real w = t, v = t, pw = pt, pv = pt, d(0), pu(0), nu(0);
  Real u(0), p(0), q(0), r(0), etmp(0), e(0), dL(0), dR(0);
  Real tm   = half * (tL + tR);
  Real tol1 = tol0 * std::abs(t)+eps;
  Real tol2 = two * tol1;
  Real tol3 = eps;
  for (unsigned it = 0u; it < maxit; ++it) {
    dL = tL-t; dR = tR-t;
    if (std::abs(e) > tol1) {
      r = (t-w)*(pt-pv);
      q = (t-v)*(pt-pw);
      p = (t-v)*q-(t-w)*r;
      q = two*(q-r);
      if (q > zero) p = -p;
      q = std::abs(q);
      etmp = e;
      e = d;
      if ( std::abs(p) >= std::abs(half*q*etmp) || p <= q*dL || p >= q*dR ) {
        e = (t > tm ? dL : dR);
        d = lam * e;
      }
      else {
        d = p/q;
        u = t+d;
        if (u-tL < tol2 || tR-u < tol2) d = (tm >= t ? tol1 : -tol1);
      }
    }
    else {
      e = (t > tm ? dL : dR);
      d = lam * e;
    }
    u = t + (std::abs(d) >= tol1 ? d : (d >= zero ? tol1 : -tol1));
    pwa.set(y); pwa.axpy(u, s);
    nobj.update(pwa,UpdateType::Trial);
    nu = nobj.value(pwa,tol); state_->nnval++;
    Qu = u * gs + nu - nold;
    pu = half * u * u * kappa + Qu;
    if (pu <= pt) {
      if (u >= t) tL = t;
      else        tR = t;
      v  = w;  w  = t;  t  = u;
      pv = pw; pw = pt; pt = pu;
      nt = nu; Qt = Qu;
    }
    else {
      if (u < t) tL = u;
      else       tR = u;
      if (pu <= pw || w == t) {
        v  = w;  w  = u;
        pv = pw; pw = pu;
      }
      else if (pu <= pv || v == t || v == w) {
        v = u; pv = pu;
      }
    }
    tm   = half * (tL+tR);
    tol1 = tol0*std::abs(t)+eps;
    tol2 = two*tol1;
    tol3 = eps0 * std::max(std::abs(Qt),one);
    if (pt <= (mu*std::min(zero,Qt)+tol3) && std::abs(t-tm) <= (tol2-half*(tR-tL))) break;
  }
  alpha = t;
  pred  = pt;
  nval  = nt;
}

// NCG Subsolver
template<typename Real>
void TrustRegionAlgorithm<Real>::dncg(Vector<Real>             &y,
                                      Real                     &sval,
                                      Real                     &nval,
                                      Vector<Real>             &gmod,
                                      const Vector<Real>       &x,
                                      Real                      del,
                                      TrustRegionModel_U<Real> &model,
                                      Objective<Real>          &nobj,
                                      Vector<Real>             &s,
                                      Vector<Real>             &pwa1,
                                      Vector<Real>             &pwa2,
                                      Vector<Real>             &pwa3,
                                      Vector<Real>             &pwa4,
                                      Vector<Real>             &pwa5,
                                      Vector<Real>             &dwa,
                                      std::ostream             &outStream) {
  // Use NCG to approximately solve TR subproblem:
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)> + phi(y)  subject to  ||y-x|| \le del
  //
  //   Inpute:
  //       y     = computed iterate
  //       sval  = smooth model value
  //       nval  = nonsmooth value
  //       gmod  = current gradient
  //       x     = current iterate
  //       del   = trust region radius
  //       model = trust region model
  //       nobj  = nonsmooth objective function
  //       s     = the current step
  //       pwa1  = the SPG iterate
  //       pwa2  = the "negative gradient"
  //       pwa3  = y - x
  //       pwa4  = the previous "negative gradient"
  //       pwa5  = temporary storage
  //       dwa   = the Hessian applied to the step
  const Real zero(0), half(0.5), one(1), two(2);
  const Real del2(del*del);
  Real tol(std::sqrt(ROL_EPSILON<Real>())), safeguard(tol);
  Real mold(sval+nval), nold(nval);
  Real snorm(0), snorm0(0), gnorm(0), gnorm0(0), gnorm2(0);
  Real alpha(1), beta(1), lambdaTmp(1), lambda(1), eta(etaNCG_);
  Real alphaMax(1), pred(0), lam1(1);
  Real sy(0), gg(0), sHs(0), gs(0), ds(0), ss(0), ss0(0);
  bool reset(true);

  // Set y = x
  y.set(x);
  pwa3.zero(); // Initially y - x = 0

  // Compute initial spectral step length
  lambdaTmp = t0_ / gmod.norm();
  lambda    = std::max(lambdaMin_,std::min(lambdaTmp,lambdaMax_));
  lam1      = lambda;

  // Compute Cauchy point via SPG
  pgstep(pwa1, pwa2, nobj, y, gmod.dual(), lambda, tol); // pwa1  = prox(x-lambda g), pwa2  = pwa1 - x
  pwa2.scale(one/lambda);                                // pwa2  = (pwa1 - x) / lambda (for smooth: pwa2 = negative gradient)
  s.set(pwa2);                                           // s     = (pwa1 - x) / lambda
  gs     = gmod.apply(s);                                // gs    = <g, prox(x-lambda g)-x> / lambda
  snorm  = s.norm();                                     // snorm = norm(prox(x-lambda g)-x) / lambda
  gnorm  = snorm;
  const Real gtol = std::min(tol1_,tol2_*gnorm);

  if (verbosity_>1) outStream << "  Nonlinear Conjugate Gradient" << std::endl;

  SPiter_ = 0;
  SPflag_ = 1;
  while (SPiter_ < maxit_) {
    SPiter_++;

    // Compute the model curvature
    model.hessVec(dwa,s,x,tol); nhess_++; // dwa = H step
    sHs   = dwa.apply(s);                 // sHs = <step, H step>

    // Compute alpha as the 1D minimize in the s direction
    ss = snorm*snorm;
    ds = s.dot(pwa3);
    alphaMax = (-ds + std::sqrt(ds*ds + ss*(del2 - snorm0*snorm0)))/ss;
    dbls(alpha,nold,pred,y,s,lam1,alphaMax,sHs,gs,nobj,pwa5);
    
    //if (sHs <= safeguard) alpha = alphaMax;
    //else {
    //  pwa5.set(y); pwa5.axpy(alphaMax, s);
    //  nobj.update(pwa5,UpdateType::Trial);
    //  nmax  = nobj.value(pwa5,tol); state_->nnval++;
    //  alpha = alphaMax * std::min(one, -(alphaMax * gs + nmax - nold)/(alphaMax * alphaMax * sHs));
    //  if (alpha <= safeguard) alpha = std::min(one, -(gs + nval - nold)/sHs);
    //}

    // Update quantities to evaluate quadratic model value and gradient
    y.axpy(alpha, s);
    gmod.axpy(alpha, dwa); // need dual here?
    sval += alpha*(gs + half*alpha*sHs);

    // Check step size
    pwa3.set(y); pwa3.axpy(-one,x);
    ss0   += alpha*(alpha*ss + two*ds);
    snorm0 = std::sqrt(ss0); // pwa3.norm();
    
    if (snorm0 >= (one-safeguard)*del) { SPflag_ = 2; break; }

    // Update spectral step length
    lambdaTmp = (sHs <= safeguard) ? t0_/gmod.norm() : ss/sHs;
    lambda    = std::max(lambdaMin_, std::min(lambdaMax_, lambdaTmp));

    // Compute SPG direction
    pwa4.set(pwa2);                                        // store previous "negative gradient"
    pgstep(pwa1, pwa2, nobj, y, gmod.dual(), lambda, tol); // pwa1 = prox(x-lambda g), pwa2 = pwa1 - x
    pwa2.scale(one/lambda);                                // pwa2 = (pwa1 - x) / lambda (for smooth: pwa2 = negative gradient)
    gnorm0 = gnorm;
    gnorm  = pwa2.norm();
    if (gnorm <= gtol) { SPflag_ = 0; break; }

    gnorm2 = gnorm * gnorm;
    switch (ncgType_) {
      case 0: // Fletcher-Reeves
        beta = gnorm2/(gnorm0 * gnorm0);
        break;
      default:
      case 1: // Polyak-Ribiere+
        pwa5.set(pwa4); pwa5.axpy(-one,pwa2);
        beta = std::max(zero, -pwa5.dot(pwa2)/(gnorm0*gnorm0));
        break;
      case 2: // Hager-Zhang
        pwa5.set(pwa4); pwa5.axpy(-one,pwa2);
        sy   = s.dot(pwa5);
        gg   = pwa2.dot(pwa4);
        eta  = -one/(s.norm()*std::min(etaNCG_, gnorm0));
        beta = std::max(eta, (gnorm2-gg-two*pwa2.dot(s)*(gnorm2-two*gg+(gnorm0*gnorm0))/sy)/sy);
        break;
      case 3: // Hestenes-Stiefel+
        pwa5.set(pwa4); pwa5.axpy(-one,pwa2);
        beta = std::max(zero, -pwa2.dot(pwa5)/s.dot(pwa5));
        break;
      case 4: // Dai-Yuan+
        pwa5.set(pwa4); pwa5.axpy(-one,pwa2);
        beta = std::max(zero, gnorm2/s.dot(pwa5));
        break;
      case 5: // Fletcher-Reeves-Polyak-Ribiere
        pwa5.set(pwa4); pwa5.axpy(-one,pwa2);
        beta = std::max(-gnorm2, std::min(gnorm2, -pwa5.dot(pwa2)))/(gnorm0*gnorm0);
        break;
      case 6: //Dai-Yuan-Hestenes-Stiefles
        pwa5.set(pwa4); pwa5.axpy(-one,pwa2);
        beta = std::max(zero, std::min(-pwa2.dot(pwa5), gnorm2)/s.dot(pwa5));
        break;
    }

    reset = true;
    if (beta != zero && beta < ROL_INF<Real>()){
      pwa5.set(pwa2);         // pwa5 = pwa2 (SPG step)
      pwa5.axpy(beta, s);     // pwa5 = pwa2 + beta*s
      pwa4.set(y);            // pwa4 = y
      pwa4.plus(pwa5);        // pwa4 = y + pwa5
      nobj.update(pwa4,UpdateType::Trial);
      nval = nobj.value(pwa4,tol); state_->nnval++;
      gs   = gmod.apply(pwa5);
      if (gs + nval - nold <= -(one - desPar_)*gnorm2) {
        s.set(pwa5);
        lam1  = one;
        reset = false;
      }
    }
    if (reset){ // Reset because either beta=0 or step does not produce descent
      s.set(pwa2);
      gs   = gmod.apply(s);
      lam1 = lambda;
      beta = zero;
    }
    snorm = s.norm();

    if (verbosity_ > 1) {
      outStream << std::endl;
      outStream << "    Iterate:                          " << SPiter_          << std::endl;
      outStream << "    Spectral step length (lambda):    " << lambda           << std::endl;
      outStream << "    Step length (alpha):              " << alpha            << std::endl;
      outStream << "    NCG parameter (beta):             " << beta             << std::endl;
      outStream << "    Model decrease (pRed):            " << mold-(sval+nold) << std::endl;
      outStream << "    Step size:                        " << snorm0           << std::endl;
      outStream << "    Optimality criterion:             " << gnorm            << std::endl;
      outStream << "    Optimality tolerance:             " << gtol             << std::endl;
      outStream << std::endl;
    }
  }
  nval = nold;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if (verbosity_ > 1) {
    os << std::string(114,'-') << std::endl;
    switch (algSelect_) {
      default:
      case TRUSTREGION_P_SPG:  os << " SPG "; break;
      case TRUSTREGION_P_SPG2: os << " Simplified SPG "; break;
      case TRUSTREGION_P_NCG:  os << " NCG "; break;
    }
    os << "trust-region method status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)" << std::endl;
    os << "  value   - Objective function value" << std::endl;
    os << "  gnorm   - Norm of the gradient" << std::endl;
    os << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    os << "  delta   - Trust-Region radius" << std::endl;
    os << "  #sval   - Number of times the smooth objective function was evaluated" << std::endl;
    os << "  #nval   - Number of times the nonsmooth objective function was evaluated" << std::endl;
    os << "  #grad   - Number of times the gradient was computed" << std::endl;
    os << "  #hess   - Number of times the Hessian was applied" << std::endl;
    os << "  #prox   - Number of times the proximal operator was computed" << std::endl;
    os << std::endl;
    os << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      os << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    os << std::endl;
    os << "  iterSP - Number of Spectral Projected Gradient iterations" << std::endl << std::endl;
    os << "  flagSP - Trust-Region Spectral Projected Gradient flag" << std::endl;
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
  os << std::setw(10) << std::left << "#sval";
  os << std::setw(10) << std::left << "#nval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "#hess";
  os << std::setw(10) << std::left << "#prox";
  os << std::setw(10) << std::left << "tr_flag";
  os << std::setw(10) << std::left << "iterSP";
  os << std::setw(10) << std::left << "flagSP";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl;
  switch (algSelect_) {
    default:
    case TRUSTREGION_P_SPG:  os << "SPG "; break;
    case TRUSTREGION_P_SPG2: os << "Simplified SPG "; break;
    case TRUSTREGION_P_NCG:  os << "NCG "; break;
  }
  os << "Trust-Region Method (Type P)" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    os << std::setw(10) << std::left << state_->nsval;
    os << std::setw(10) << std::left << state_->nnval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << nhess_;
    os << std::setw(10) << std::left << state_->nprox;
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
    os << std::setw(10) << std::left << state_->nsval;
    os << std::setw(10) << std::left << state_->nnval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << nhess_;
    os << std::setw(10) << std::left << state_->nprox;
    os << std::setw(10) << std::left << TRflag_;
    os << std::setw(10) << std::left << SPiter_;
    os << std::setw(10) << std::left << SPflag_;
    os << std::endl;
  }
  os.flags(osFlags);
}

} // namespace TypeP
} // namespace ROL

#endif
