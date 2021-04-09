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
  computeGradient(x,obj);
  // Check if inverse Hessian is implemented for dogleg methods
  model_->validate(obj,x,g,etr_);
  // Compute initial trust region radius if desired.
  if ( state_->searchSize <= static_cast<Real>(0) ) {
    int nfval = 0;
    state_->searchSize
      = TRUtils::initialRadius<Real>(nfval,x,*state_->gradientVec,Bg,
          state_->value,state_->gnorm,obj,*model_,delMax_,
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
void TrustRegionAlgorithm<Real>::computeGradient( const Vector<Real> &x,
                                                  Objective<Real>    &obj) {
  if ( useInexact_[1] ) {
    const Real one(1);
    Real gtol1 = scale0_*state_->searchSize;
    Real gtol0 = gtol1 + one;
    while ( gtol0 > gtol1 ) {
      obj.gradient(*state_->gradientVec,x,gtol1);
      state_->gnorm = state_->gradientVec->norm();
      gtol0 = gtol1;
      gtol1 = scale0_*std::min(state_->gnorm,state_->searchSize);
    }
  }
  else {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    obj.gradient(*state_->gradientVec,x,gtol);
    state_->gnorm = state_->gradientVec->norm();
  }
  state_->ngrad++;
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

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Build trust-region model
    model_->setData(obj,x,*state_->gradientVec);
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
      if (useInexact_[1]) computeGradient(x,obj);
    }
    else if ((rho >= eta0_ && TRflag_ != TRUtils::NPOSPREDNEG)
             || (TRflag_ == TRUtils::POSPREDNEG)) { // Step Accepted
      state_->iterateVec->set(x);
      state_->value = ftrial;
      obj.update(x,UpdateType::Accept,state_->iter);
      // Increase trust-region radius
      if (rho >= eta2_) state_->searchSize = std::min(gamma2_*state_->searchSize, delMax_);
      // Compute gradient at new iterate
      gvec->set(*state_->gradientVec);
      computeGradient(x,obj);
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
  std::stringstream hist;
  if(verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    hist << "Trust-Region status output definitions" << std::endl << std::endl;
    hist << "  iter    - Number of iterates (steps taken)" << std::endl;
    hist << "  value   - Objective function value" << std::endl; 
    hist << "  gnorm   - Norm of the gradient" << std::endl;
    hist << "  snorm   - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  delta   - Trust-Region radius" << std::endl;
    hist << "  #fval   - Number of times the objective function was evaluated" << std::endl;
    hist << "  #grad   - Number of times the gradient was computed" << std::endl;
    hist << std::endl;
    hist << "  tr_flag - Trust-Region flag" << std::endl;
    for( int flag = TRUtils::SUCCESS; flag != TRUtils::UNDEFINED; ++flag ) {
      hist << "    " << NumberToString(flag) << " - "
           << TRUtils::ETRFlagToString(static_cast<TRUtils::ETRFlag>(flag)) << std::endl;
    }
    if( etr_ == TRUSTREGION_U_TRUNCATEDCG ) {
      hist << std::endl;
      hist << "  iterCG - Number of Truncated CG iterations" << std::endl << std::endl;
      hist << "  flagGC - Trust-Region Truncated CG flag" << std::endl;
      for( int flag = CG_FLAG_SUCCESS; flag != CG_FLAG_UNDEFINED; ++flag ) {
        hist << "    " << NumberToString(flag) << " - "
             << ECGFlagToString(static_cast<ECGFlag>(flag)) << std::endl;
      }            
    }
    else if( etr_ == TRUSTREGION_U_SPG ) {
      hist << std::endl;
      hist << "  iterCG - Number of spectral projected gradient iterations" << std::endl << std::endl;
      hist << "  flagGC - Trust-Region spectral projected gradient flag" << std::endl;
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
  hist << std::setw(10) << std::left << "tr_flag";
  if ( etr_ == TRUSTREGION_U_TRUNCATEDCG ) {
    hist << std::setw(10) << std::left << "iterCG";
    hist << std::setw(10) << std::left << "flagCG";
  }
  else if (etr_ == TRUSTREGION_U_SPG) {
    hist << std::setw(10) << std::left << "iterSPG";
    hist << std::setw(10) << std::left << "flagSPG";
  }
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << ETrustRegionUToString(etr_) << " Trust-Region Solver";
  if ( useSecantPrecond_ || useSecantHessVec_ ) {
    if ( useSecantPrecond_ && !useSecantHessVec_ ) {
      hist << " with " << ESecantToString(esec_) << " Preconditioning" << std::endl;
    }
    else if ( !useSecantPrecond_ && useSecantHessVec_ ) {
      hist << " with " << ESecantToString(esec_) << " Hessian Approximation" << std::endl;
    }
    else {
      hist << " with " << ESecantToString(esec_) << " Preconditioning and Hessian Approximation" << std::endl;
    }
  }
  else {
    hist << std::endl;
  }
  os << hist.str();
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeOutput(std::ostream& os, bool print_header) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) {
    writeName(os);
  }
  if ( print_header ) {
    writeHeader(os);
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
    hist << std::setw(10) << std::left << "---";
    if ( etr_ == TRUSTREGION_U_TRUNCATEDCG || etr_ == TRUSTREGION_U_SPG ) {
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
    hist << std::setw(10) << std::left << TRflag_;
    if ( etr_ == TRUSTREGION_U_TRUNCATEDCG || etr_ == TRUSTREGION_U_SPG ) {
      hist << std::setw(10) << std::left << SPiter_;
      hist << std::setw(10) << std::left << SPflag_;
    }
    hist << std::endl;
  }
  os << hist.str();
}
} // namespace TypeU
} // namespace ROL

#endif
