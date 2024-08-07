// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUSTREGIONALGORITHM_U_DEF_H
#define ROL_TRUSTREGIONALGORITHM_U_DEF_H

#include "ROL_TrustRegion_U_Factory.hpp"

namespace ROL {
namespace TypeU {

template<typename Real>
TrustRegionAlgorithm<Real>::TrustRegionAlgorithm( ParameterList &parlist,
                                                  const Ptr<Secant<Real>> &secant )
    : Algorithm<Real>(), esec_(SECANT_USERDEFINED) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(parlist));

  // Trust-Region Parameters
  ParameterList &slist = parlist.sublist("Step");
  ParameterList &trlist  = slist.sublist("Trust Region");
  state_->searchSize = trlist.get("Initial Radius",            static_cast<Real>(-1));
  delMax_ = trlist.get("Maximum Radius",                       ROL_INF<Real>());
  eta0_   = trlist.get("Step Acceptance Threshold",            static_cast<Real>(0.05));
  eta1_   = trlist.get("Radius Shrinking Threshold",           static_cast<Real>(0.05));
  eta2_   = trlist.get("Radius Growing Threshold",             static_cast<Real>(0.9));
  gamma0_ = trlist.get("Radius Shrinking Rate (Negative rho)", static_cast<Real>(0.0625));
  gamma1_ = trlist.get("Radius Shrinking Rate (Positive rho)", static_cast<Real>(0.25));
  gamma2_ = trlist.get("Radius Growing Rate",                  static_cast<Real>(2.5));
  TRsafe_ = trlist.get("Safeguard Size",                       static_cast<Real>(100.0));
  eps_    = TRsafe_*ROL_EPSILON<Real>();
  // Nonmonotone Information
  NMstorage_ = trlist.get("Nonmonotone Storage Limit", 0);
  useNM_     = (NMstorage_ <= 0 ? false : true);
  // Inexactness Information
  ParameterList &glist = parlist.sublist("General");
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
  // Initialize Trust Region Subproblem Solver Object
  etr_       = StringToETrustRegionU(trlist.get("Subproblem Solver", "Dogleg"));  
  solver_    = TrustRegionUFactory<Real>(parlist);
  verbosity_ = glist.get("Output Level", 0);
  // Secant Information
  useSecantPrecond_ = glist.sublist("Secant").get("Use as Preconditioner", false);
  useSecantHessVec_ = glist.sublist("Secant").get("Use as Hessian",        false);
  if (secant == nullPtr) {
    esec_ = StringToESecant(glist.sublist("Secant").get("Type","Limited-Memory BFGS"));
  }
  // Initialize trust region model
  model_ = makePtr<TrustRegionModel_U<Real>>(parlist,secant);
  printHeader_ = verbosity_ > 2;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::initialize( const Vector<Real> &x,
                                             const Vector<Real> &g,
                                             Vector<Real>       &Bg,
                                             Objective<Real>    &obj,
                                             std::ostream &outStream) {
  const Real zero(0);
  // Initialize data
  Algorithm<Real>::initialize(x,g);
  solver_->initialize(x,g);
  model_->initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(); 
  obj.update(x,UpdateType::Initial,state_->iter);    
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  state_->snorm = ROL_INF<Real>();
  state_->gnorm = ROL_INF<Real>();
  Real Delta = state_->searchSize;
  if (Delta <= zero) state_->searchSize = 1e2*x.norm();
  computeGradient(x,obj,true);
  // Check if inverse Hessian is implemented for dogleg methods
  model_->validate(obj,x,g,etr_);
  // Compute initial trust region radius if desired.
  if ( Delta <= zero ) {
    int nfval = 0;
    state_->searchSize
      = TRUtils::initialRadius<Real>(nfval,x,*state_->gradientVec,Bg,
          state_->value,state_->gnorm,gtol_,obj,*model_,delMax_,
          outStream,(verbosity_>1));
    state_->nfval += nfval;
  }
}

template<typename Real>
Real TrustRegionAlgorithm<Real>::computeValue( const Vector<Real> &x,
                                               Objective<Real>    &obj,
                                               Real               pRed) {
  const Real one(1);
  Real tol(std::sqrt(ROL_EPSILON<Real>())), fval(0);
  if ( useInexact_[0] ) {
    if ( !(state_->iter%updateIter_) && (state_->iter != 0) ) {
      force_ *= forceFactor_;
    }
    Real eta = static_cast<Real>(0.999)*std::min(eta1_,one-eta2_);
    tol      = scale_*std::pow(eta*std::min(pRed,force_),one/omega_);
    state_->value = obj.value(*state_->iterateVec,tol);
    state_->nfval++;
  }
  // Evaluate objective function at new iterate
  obj.update(x,UpdateType::Trial);
  fval = obj.value(x,tol);
  state_->nfval++;
  return fval;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::computeGradient(const Vector<Real> &x,
                                                 Objective<Real>    &obj,
                                                 bool accept) {
  if ( useInexact_[1] ) {
    Real gtol0 = scale0_*state_->searchSize;
    if (accept) gtol_ = gtol0 + static_cast<Real>(1);
    else        gtol0 = scale0_*std::min(state_->gnorm,state_->searchSize);
    while ( gtol_ > gtol0 ) {
      gtol_ = gtol0;
      obj.gradient(*state_->gradientVec,x,gtol_); state_->ngrad++;
      state_->gnorm = state_->gradientVec->norm();
      gtol0 = scale0_*std::min(state_->gnorm,state_->searchSize);
    }
  }
  else {
    if (accept) {
      gtol_ = std::sqrt(ROL_EPSILON<Real>());
      obj.gradient(*state_->gradientVec,x,gtol_); state_->ngrad++;
      state_->gnorm = state_->gradientVec->norm();
    }
  }
}

template<typename Real>
void TrustRegionAlgorithm<Real>::run( Vector<Real>       &x,
                                      const Vector<Real> &g, 
                                      Objective<Real>    &obj,
                                      std::ostream       &outStream ) {
  const Real zero(0);
  // Initialize trust-region data
  Real ftrial(0), pRed(0), rho(0);
  Ptr<Vector<Real>> gvec = g.clone();
  initialize(x,g,*gvec,obj,outStream);
  // Initialize nonmonotone data
  Real rhoNM(0), sigmac(0), sigmar(0);
  Real fr(state_->value), fc(state_->value), fmin(state_->value);
  TRUtils::ETRFlag TRflagNM;
  int L(0);

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Build trust-region model
    model_->setData(obj,x,*state_->gradientVec,gtol_);
    // Minimize trust-region model over trust-region constraint
    pRed = zero;
    SPflag_ = 0; SPiter_ = 0;
    solver_->solve(*state_->stepVec,state_->snorm,pRed,SPflag_,SPiter_,
                   state_->searchSize,*model_);
    // Compute trial objective function value
    x.plus(*state_->stepVec);
    ftrial = computeValue(x,obj,pRed);
    // Compute ratio of actual and predicted reduction
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
    if ((rho < eta0_ && TRflag_ == TRUtils::SUCCESS)
        || (TRflag_ >= 2)) { // Step Rejected
      x.set(*state_->iterateVec);
      obj.update(x,UpdateType::Revert,state_->iter);
      if (rho < zero && TRflag_ != TRUtils::TRNAN) {
        // Negative reduction, interpolate to find new trust-region radius
        state_->searchSize = TRUtils::interpolateRadius<Real>(*state_->gradientVec,*state_->stepVec,
          state_->snorm,pRed,state_->value,ftrial,state_->searchSize,gamma0_,gamma1_,eta2_,
          outStream,verbosity_>1);
      }
      else { // Shrink trust-region radius
        state_->searchSize = gamma1_*std::min(state_->snorm,state_->searchSize);
      }
      computeGradient(x,obj,false);
    }
    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
             || (TRflag_ == TRUtils::POSPREDNEG)) { // Step Accepted
      state_->iterateVec->set(x);
      state_->value = ftrial;
      obj.update(x,UpdateType::Accept,state_->iter);
      if (useNM_) {
        sigmac += pRed; sigmar += pRed;
        if (ftrial < fmin) { fmin = ftrial; fc = fmin; sigmac = zero; L = 0; }
        else {
          L++;
          if (ftrial > fc)     { fc = ftrial; sigmac = zero;   }
          if (L == NMstorage_) { fr = fc;     sigmar = sigmac; }
        }
      }
      // Increase trust-region radius
      if (rho >= eta2_) state_->searchSize = std::min(gamma2_*state_->searchSize, delMax_);
      // Compute gradient at new iterate
      gvec->set(*state_->gradientVec);
      computeGradient(x,obj,true);
      // Update secant information in trust-region model
      model_->update(x,*state_->stepVec,*gvec,*state_->gradientVec,
                     state_->snorm,state_->iter);
    }
    // Update Output
    if (verbosity_ > 0)  writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  if(verbosity_ > 1) {
    os << std::string(114,'-') << std::endl;
    os << "Trust-Region status output definitions" << std::endl << std::endl;
    os << "  iter    - Number of iterates (steps taken)" << std::endl;
    os << "  value   - Objective function value" << std::endl; 
    os << "  gnorm   - Norm of the gradient" << std::endl;
    os << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    os << "  delta   - Trust-Region radius" << std::endl;
    os << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    os << "  #grad   - Number of times the gradient was computed" << std::endl;
    os << std::endl;
    os << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      os << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    if( etr_ == TRUSTREGION_U_TRUNCATEDCG ) {
      os << std::endl;
      os << "  iterCG - Number of Truncated CG iterations" << std::endl << std::endl;
      os << "  flagGC - Trust-Region Truncated CG flag" << std::endl;
      for( int flag = CG_FLAG_SUCCESS; flag != CG_FLAG_UNDEFINED; ++flag ) {
        os << "    " << NumberToString(flag) << " - "
             << ECGFlagToString(static_cast<ECGFlag>(flag)) << std::endl;
      }            
    }
    else if( etr_ == TRUSTREGION_U_SPG ) {
      os << std::endl;
      os << "  iterCG - Number of spectral projected gradient iterations" << std::endl << std::endl;
      os << "  flagGC - Trust-Region spectral projected gradient flag" << std::endl;
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
  os << std::setw(10) << std::left << "tr_flag";
  if ( etr_ == TRUSTREGION_U_TRUNCATEDCG ) {
    os << std::setw(10) << std::left << "iterCG";
    os << std::setw(10) << std::left << "flagCG";
  }
  else if (etr_ == TRUSTREGION_U_SPG) {
    os << std::setw(10) << std::left << "iterSPG";
    os << std::setw(10) << std::left << "flagSPG";
  }
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << ETrustRegionUToString(etr_) << " Trust-Region Solver";
  if ( useSecantPrecond_ || useSecantHessVec_ ) {
    if ( useSecantPrecond_ && !useSecantHessVec_ ) {
      os << " with " << ESecantToString(esec_) << " Preconditioning" << std::endl;
    }
    else if ( !useSecantPrecond_ && useSecantHessVec_ ) {
      os << " with " << ESecantToString(esec_) << " Hessian Approximation" << std::endl;
    }
    else {
      os << " with " << ESecantToString(esec_) << " Preconditioning and Hessian Approximation" << std::endl;
    }
  }
  else {
    os << std::endl;
  }
  os.flags(osFlags);
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeOutput(std::ostream& os, bool print_header) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) {
    writeName(os);
  }
  if ( print_header ) {
    writeHeader(os);
  }
  if ( state_->iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << state_->nfval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(10) << std::left << "---";
    if ( etr_ == TRUSTREGION_U_TRUNCATEDCG || etr_ == TRUSTREGION_U_SPG ) {
      os << std::setw(10) << std::left << "---";
      os << std::setw(10) << std::left << "---";
    }
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
    os << std::setw(10) << std::left << TRflag_;
    if ( etr_ == TRUSTREGION_U_TRUNCATEDCG || etr_ == TRUSTREGION_U_SPG ) {
      os << std::setw(10) << std::left << SPiter_;
      os << std::setw(10) << std::left << SPflag_;
    }
    os << std::endl;
  }
  os.flags(osFlags);
}
} // namespace TypeU
} // namespace ROL

#endif
