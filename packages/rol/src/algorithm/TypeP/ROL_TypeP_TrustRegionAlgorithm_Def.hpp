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
  verbosity_ = trlist.sublist("General").get("Output Level",   0);
  initProx_  = trlist.get("Apply Prox to Initial Guess", false);
  t0_        = list.sublist("Status Test").get("Proximal Gradient Parameter",         1.0);
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
  // Subsolver (spectral projected gradient) parameters
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
  state_->gnorm = computeGradient(x,*state_->gradientVec,px,dg,*state_->stepVec,state_->searchSize,sobj,nobj,outStream);
  
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
Real TrustRegionAlgorithm<Real>::computeGradient(const Vector<Real> &x,
                                                 Vector<Real> &g,
                                                 Vector<Real> &px,
                                                 Vector<Real> &dg,
                                                 Vector<Real> &step,
                                                 Real del,
                                                 Objective<Real> &sobj,
                                                 Objective<Real> &nobj,
                                                 std::ostream &outStream) const {
  Real gnorm(0);
  if ( useInexact_[1] ) {
    const Real one(1);
    Real gtol1 = scale0_*del;
    Real gtol0 = gtol1 + one;
    while ( gtol0 > gtol1 ) {
      sobj.gradient(g,x,gtol1); state_->ngrad++; 
      dg.set(g.dual());
      pgstep(px, step, nobj, x, dg, t0_, gtol1); // change gtol? one or ocScale? 
      gnorm = step.norm() / t0_;
      gtol0 = gtol1;
      gtol1 = scale0_*std::min(gnorm,del);
    }
  }
  else {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    sobj.gradient(g,x,gtol); state_->ngrad++; 
    dg.set(g.dual()); 
    pgstep(px, step, nobj, x, dg, t0_, gtol); 
    gnorm = step.norm() / t0_; 
  }
  return gnorm;
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
    model_->setData(sobj,*state_->iterateVec,*state_->gradientVec);

    /**** SOLVE TRUST-REGION SUBPROBLEM ****/
    //q = state_->svalue + state_->nvalue;//q is no longer used
    gmod->set(*state_->gradientVec);
    if (useSimpleSPG_)
      dpsg_simple(x,state_->svalue, state_->nvalue,pRed,*gmod,*state_->iterateVec,*px, *dg, state_->searchSize,*model_,nobj,
                  *pwa1,*pwa2,*dwa1,outStream);
    else {
      // Compute Cauchy point (TRON notation: x = x[1])
      smodel = state_->svalue;
      ntrial = state_->nvalue;
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
    }

    // Update storage and compute predicted reduction
    //pRed = -q; // now updated in dcauchy/dspg
    state_->stepVec->set(x); state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Compute trial objective value
    strial = computeValue(inTol,outTol,pRed,state_->svalue,state_->iter,x,*state_->iterateVec,sobj);
    nobj.update(x, UpdateType::Trial);
    ntrial = nobj.value(x,outTol); state_->nnval++; 
    Ftrial = strial + ntrial; 

    // Compute ratio of acutal and predicted reduction
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
      state_->gnorm = computeGradient(x,*state_->gradientVec,*px,*dg,*pwa1,state_->searchSize,sobj,nobj,outStream);
      state_->ngrad++;
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
void TrustRegionAlgorithm<Real>::dpsg_simple(Vector<Real> &y,     //x
                                             Real         &sval,
																						 Real         &nval,
																						 Real         &pRed,  //predicted reduction
                                             Vector<Real> &gmod,  // state_->gradientVec
                                             const Vector<Real> &x, //iterateVect
																						 Vector<Real> &px,    // px
																						 Vector<Real> &dg,    // gradient dual? 
                                             Real del,
                                             TrustRegionModel_U<Real> &model,
																						 Objective<Real> &nobj,
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
  Real alpha(1), alphaMax(1), s0s0(0), ss0(0), sHs(0), lambdaTmp(1), snorm(0), mold(sval+nval), mnew(mold);
  pwa1.zero();

  // Set y = x
  y.set(x);

  // Compute initial step
  Real coeff  = one/gmod.norm();
  Real lambda = std::max(lambdaMin_,std::min(coeff,lambdaMax_));
  pgstep(px, pwa, nobj, y, gmod.dual(), alpha, tol); // pass pwa by reference? *pwa?
  //pwa.set(y); pwa.axpy(-lambda,gmod.dual());      // pwa = x - lambda gmod.dual()
  //proj_->project(pwa,outStream); state_->nproj++; // pwa = P(x - lambda gmod.dual())
  //pwa.axpy(-one,y);                               // pwa = P(x - lambda gmod.dual()) - x = step
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
    y.axpy(alpha,pwa); // New iterate
		nval   = nobj.value(y, tol); 
		sval += alpha * (gs + half * alpha * sHs); // Update model value - is this correct? 
    mnew  = sval + nval; 
		gmod.axpy(alpha,dwa);                   // Update model gradient

		pRed = mold - mnew; 

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
    lambdaTmp = (sHs <= safeguard ? one/gmod.norm() : ss/sHs);
    lambda = std::max(lambdaMin_,std::min(lambdaTmp,lambdaMax_));
    
    pgstep(px, pwa, nobj, y, gmod.dual(), alpha, tol); // pass pwa by reference? *pwa?
    //pwa.set(y); pwa.axpy(-lambda,gmod.dual());
    //proj_->project(pwa,outStream); state_->nproj++;
    //pwa.axpy(-one,y);
    gs = gmod.apply(pwa);
    ss = pwa.dot(pwa);
    gnorm = std::sqrt(ss);

    if (gnorm <= gtol) { SPflag_ = 0; break; }
  }
  SPflag_ = (SPiter_==maxit_) ? 1 : SPflag_;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::dspg(Vector<Real> &y,
                                      Real         &sval,
                                      Real         &nval, 
                                      Vector<Real> &gmod,
                                      const Vector<Real> &x,
                                      Real del,
                                      TrustRegionModel_U<Real> &model,
                                      Objective<Real> &nobj,
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
  //   min 1/2 <H(y-x), (y-x)> + <g, (y-x)> + phi(y)  subject to  ||y|| \le del
  //
  //   Inpute:
  //       y = Cauchy step
  //       x = Current iterate
  //       g = Current gradient
  const Real half(0.5), one(1);
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
    pwa.axpy(-one,y);                                           // pwa = P(y - lambda gmod.dual()) - y = step
    pwa2.set(y); pwa2.plus(pwa);                                // pwa2 = P(y - lambda gmod.dual())
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
    mcomp = useNMSP_ ? *std::max_element(mqueue.begin(),mqueue.end()) : mold;
    while( mnew > mcomp + mu0_*Qk) {
      alpha *= interpf_;
      pwa2.set(y); pwa2.axpy(alpha,pwa);
      nobj.update(pwa2, UpdateType::Trial);
      nnew  = nobj.value(pwa2, tol); state_->nnval++; 
      snew  = half * alpha * alpha * sHs + alpha * gs + sold;
      mnew  = nnew + snew;
      Qk    = alpha * gs + nnew - nold;
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

// BRACKETING AND BRENTS FOR UNTRANSFORMED MULTIPLIER
//template<typename Real>
//void TrustRegionSPGAlgorithm<Real>::dprox(Vector<Real> &x,
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
//void TrustRegionSPGAlgorithm<Real>::dprox(Vector<Real> &x,
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
void TrustRegionAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    hist << " SPG trust-region method status output definitions" << std::endl << std::endl;
    hist << "  iter    - Number of iterates (steps taken)" << std::endl;
    hist << "  value   - Objective function value" << std::endl; 
    hist << "  gnorm   - Norm of the gradient" << std::endl;
    hist << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  delta   - Trust-Region radius" << std::endl;
    hist << "  #sval   - Number of times the smooth objective function was evaluated" << std::endl;
    hist << "  #nval   - Number of times the nonsmooth objective function was evaluated" << std::endl;
    hist << "  #grad   - Number of times the gradient was computed" << std::endl;
    hist << "  #hess   - Number of times the Hessian was applied" << std::endl;
    hist << "  #prox   - Number of times the proximal operator was computed" << std::endl;
    hist << std::endl;
    hist << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      hist << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    hist << std::endl;
    hist << "  iterSPG - Number of Spectral Projected Gradient iterations" << std::endl << std::endl;
    hist << "  flagSPG - Trust-Region Spectral Projected Gradient flag" << std::endl;
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
  hist << std::setw(10) << std::left << "#sval";
  hist << std::setw(10) << std::left << "#nval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "#hess";
  hist << std::setw(10) << std::left << "#prox";
  hist << std::setw(10) << std::left << "tr_flag";
  hist << std::setw(10) << std::left << "iterSPG";
  hist << std::setw(10) << std::left << "flagSPG";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "SPG Trust-Region Method (Type P)" << std::endl;
  os << hist.str();
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
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
    hist << std::setw(10) << std::left << state_->nsval;
    hist << std::setw(10) << std::left << state_->nnval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << nhess_;
    hist << std::setw(10) << std::left << state_->nprox;
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
    hist << std::setw(10) << std::left << state_->nsval;
    hist << std::setw(10) << std::left << state_->nnval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << nhess_;
    hist << std::setw(10) << std::left << state_->nprox;
    hist << std::setw(10) << std::left << TRflag_;
    hist << std::setw(10) << std::left << SPiter_;
    hist << std::setw(10) << std::left << SPflag_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeP
} // namespace ROL

#endif
