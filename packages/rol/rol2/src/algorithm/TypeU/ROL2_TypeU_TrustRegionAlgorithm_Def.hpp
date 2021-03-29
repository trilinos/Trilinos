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

#pragma once
#ifndef ROL2_TYPEU_TRUSTREGIONALGORITHM_DEF_H
#define ROL2_TYPEU_TRUSTREGIONALGORITHM_DEF_H

namespace ROL2 {
namespace TypeU {

template<class Real>
TrustRegionAlgorithm<Real>::TrustRegionAlgorithm(       ParameterList&     parlist,
                                                  const Ptr<Secant<Real>>& secant )
 : Algorithm<Real>(), secantType_( Secant<Real>::Type::UserDefined ),
   solver_(TrustRegion<Real>::create(parlist)) {
  auto& state  = Algorithm<Real>::getState();
  auto& status = Algorithm<Real>::getStatus(); 
  
  // Set status test
  status.reset();
  status.add(makePtr<StatusTest<Real>>(parlist));

  // Trust-Region Parameters
  auto& slist = parlist.sublist("Step");
  auto& trlist  = slist.sublist("Trust Region");
  state.searchSize_ = trlist.get("Initial Radius",             static_cast<Real>(-1));
  delMax_ = trlist.get("Maximum Radius",                       ROL_INF<Real>);
  eta0_   = trlist.get("Step Acceptance Threshold",            static_cast<Real>(0.05));
  eta1_   = trlist.get("Radius Shrinking Threshold",           static_cast<Real>(0.05));
  eta2_   = trlist.get("Radius Growing Threshold",             static_cast<Real>(0.9));
  gamma0_ = trlist.get("Radius Shrinking Rate (Negative rho)", static_cast<Real>(0.0625));
  gamma1_ = trlist.get("Radius Shrinking Rate (Positive rho)", static_cast<Real>(0.25));
  gamma2_ = trlist.get("Radius Growing Rate",                  static_cast<Real>(2.5));
  TRsafe_ = trlist.get("Safeguard Size",                       static_cast<Real>(100.0));
  eps_    = TRsafe_ * ROL_EPSILON<Real>;

  // Inexactness Information
  auto& glist = parlist.sublist("General");
  useInexact_.clear();
  useInexact_.push_back(glist.get("Inexact Objective Function",     false));
  useInexact_.push_back(glist.get("Inexact Gradient",               false));
  useInexact_.push_back(glist.get("Inexact Hessian-Times-A-Vector", false));

  // Trust-Region Inexactness Parameters
  auto& ilist = trlist.sublist("Inexact").sublist("Gradient");
  scale0_ = ilist.get("Tolerance Scaling",  static_cast<Real>(0.1));
  scale1_ = ilist.get("Relative Tolerance", static_cast<Real>(2)); 

  // Inexact Function Evaluation Information
  auto& vlist = trlist.sublist("Inexact").sublist("Value");
  scale_       = vlist.get("Tolerance Scaling",                 static_cast<Real>(1.e-1));
  omega_       = vlist.get("Exponent",                          static_cast<Real>(0.9));
  force_       = vlist.get("Forcing Sequence Initial Value",    static_cast<Real>(1.0));
  updateIter_  = vlist.get("Forcing Sequence Update Frequency", static_cast<int>(10));
  forceFactor_ = vlist.get("Forcing Sequence Reduction Factor", static_cast<Real>(0.1));

  // Initialize Trust Region Subproblem Solver Object
  std::string solvertype = trlist.get("Subproblem Solver", "DogLeg");
  type_   = TrustRegion<Real>::type_dict[solvertype];  
  //std::cout <<  "solver = " << solver_ << std::endl;
  verbosity_ = glist.get("Output Level", 0);

  // Secant Information
  auto& seclist = glist.sublist("Secant");
  useSecantPrecond_ = seclist.get("Use as Preconditioner", false);
  useSecantHessVec_ = seclist.sublist("Secant").get("Use as Hessian", false);

  if (secant == nullPtr) 
    secantType_ = Secant<Real>::type_dict[seclist.get("Type","Limited-Memory BFGS")];

  // Initialize trust region model
  model_ = makePtr<TrustRegionModel<Real>>(parlist,secant);
  printHeader_ = verbosity_ > 2;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::initialize( const Vector<Real>&    x,
                                             const Vector<Real>&    g,
                                                   Vector<Real>&    Bg,
                                                   Objective<Real>& obj,
                                                   std::ostream&    os) {
  // Initialize data
  Algorithm<Real>::initialize(x,g);
  solver_->initialize(x,g);
  model_->initialize(x,g);

  auto& state = Algorithm<Real>::getState();

  // Update approximate gradient and approximate objective function.
  Real ftol = static_cast<Real>(0.1)*ROL_MAX<Real>; 
  obj.update(x,UpdateType::Initial,state.iter_);    
  state.value_ = obj.value(x,ftol); 
  state.nfval_++;
  state.snorm_ = ROL_INF<Real>;
  state.gnorm_ = ROL_INF<Real>;
  computeGradient(x,obj);

  // Check if inverse Hessian is implemented for dogleg methods
  model_->validate(obj,x,g,type_);

  // Compute initial trust region radius if desired.
  if ( state.searchSize_ <= static_cast<Real>(0) ) {
    int nfval = 0;
    state.searchSize_ = TrustRegion<Real>::initialRadius(nfval,x,*state.gradientVec_,Bg,
                       state.value_,state.gnorm_,obj,*model_,delMax_,
                       os,(verbosity_>1));
    state.nfval_ += nfval;
  }
}

template<typename Real>
Real TrustRegionAlgorithm<Real>::computeValue( const Vector<Real>&    x,
                                                     Objective<Real>& obj,
                                                     Real             pRed) {
  const Real one(1);
  Real tol(default_tolerance<Real>()),fval(0);
  auto& state = Algorithm<Real>::getState();

  if ( useInexact_[0] ) {

    if ( !(state.iter_ % updateIter_) && (state.iter_ != 0) ) {
      force_ *= forceFactor_;
    }

    Real eta = static_cast<Real>(0.999)*std::min(eta1_,one-eta2_);
    tol      = scale_*std::pow(eta*std::min(pRed,force_),one/omega_);
    state.value_ = obj.value(*state.iterateVec_,tol);
    state.nfval_++;
  }

  // Evaluate objective function at new iterate
  obj.update(x,UpdateType::Trial);
  fval = obj.value(x,tol);
  state.nfval_++;
  return fval;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::computeGradient(const Vector<Real>&    x,
                                                       Objective<Real>& obj) {
  auto& state = Algorithm<Real>::getState();

  if ( useInexact_[1] ) {
    const Real one(1);
    Real gtol1 = scale0_*state.searchSize_;
    Real gtol0 = gtol1 + one;

    while ( gtol0 > gtol1 ) {
      obj.gradient(*state.gradientVec_,x,gtol1);
      state.gnorm_ = state.gradientVec_->norm();
      gtol0 = gtol1;
      gtol1 = scale0_*std::min(state.gnorm_,state.searchSize_);
    }
  }
  else {
    Real gtol = default_tolerance<Real>();
    obj.gradient(*state.gradientVec_,x,gtol);
    state.gnorm_ = state.gradientVec_->norm();
  }
  state.ngrad_++;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::run(       Vector<Real>&     x,
                                      const Vector<Real>&     g, 
                                             Objective<Real>& obj,
                                             std::ostream&    os ) {
  const Real zero(0);
  auto& state  = Algorithm<Real>::getState();
  auto& status = Algorithm<Real>::getStatus(); 

  // Initialize trust-region data
//  std::vector<std::string> output;
  Real ftrial(0), pRed(0), rho(0);
  auto gvec = g.clone();
  initialize(x,g,*gvec,obj,os);

  // Output
//  output.push_back(print(true));
//  if (verbosity_ > 0) os << print(true);

  while (status.check(state)) {

    // Build trust-region model
    model_->setData(obj,x,*state.gradientVec_);

    // Minimize trust-region model over trust-region constraint
    pRed = zero;
    SPflag_ = 0; SPiter_ = 0;
    solver_->solve(*state.stepVec_,state.snorm_,pRed,SPflag_,SPiter_,
                   state.searchSize_,*model_);

    // Compute trial objective function value
    x.plus(*state.stepVec_);
    ftrial = computeValue(x,obj,pRed);

    // Compute ratio of actual and predicted reduction
    flag_ = TrustRegion<Real>::Flag::Success;
    TrustRegion<Real>::analyzeRatio(rho,flag_,state.value_,ftrial,pRed,eps_,os,verbosity_>1);

    // Update algorithm state
    state.iter_++;

    // Accept/reject step and update trust region radius
    if ((rho < eta0_ && flag_ == TrustRegion<Real>::Flag::Success) || (flag_ >= 2)) { // Step Rejected

      x.set(*state.iterateVec_);
      obj.update(x,UpdateType::Revert,state.iter_);

      if (rho < zero && flag_ != TrustRegion<Real>::Flag::TRNaN) {

        // Negative reduction, interpolate to find new trust-region radius
        state.searchSize_ = TrustRegion<Real>::interpolateRadius(*state.gradientVec_,*state.stepVec_,
        state.snorm_,pRed,state.value_,ftrial,state.searchSize_,gamma0_,gamma1_,eta2_,
        os,verbosity_>1);
      }
      else { // Shrink trust-region radius
        state.searchSize_ = gamma1_*std::min(state.snorm_,state.searchSize_);
      }

      if (useInexact_[1]) computeGradient(x,obj);
    }
    else if ((rho >= eta0_ && flag_ != TrustRegion<Real>::Flag::NPosPredNeg)
             || (flag_ == TrustRegion<Real>::Flag::PosPredNeg)) { // Step Accepted
      state.iterateVec_->set(x);
      state.value_ = ftrial;
      obj.update(x,UpdateType::Accept,state.iter_);

      // Increase trust-region radius
      if (rho >= eta2_) state.searchSize_ = std::min(gamma2_*state.searchSize_, delMax_);

      // Compute gradient at new iterate
      gvec->set(*state.gradientVec_);
      computeGradient(x,obj);

      // Update secant information in trust-region model
      model_->update(x,*state.stepVec_,*gvec,*state.gradientVec_,
                     state.snorm_,state.iter_);
    }

    // Update Output
//    output.push_back(print(printHeader_));

    if (verbosity_ > 0) writeOutput(os,printHeader_);
  }

//  output.push_back(Algorithm<Real>::printExitStatus());
//  if (verbosity_ > 0) os << Algorithm<Real>::printExitStatus();
//  return output;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeHeader( std::ostream& os ) const {

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

    for( auto flag = TrustRegion<Real>::Flag::Success; 
         flag != TrustRegion<Real>::Flag::Undefined; ++flag ) {
      os << "    " << flag << " - "
           << TrustRegion<Real>::flag_dict[flag] << std::endl;
    }

    if( type_ == TrustRegion<Real>::Type::TruncatedCG ) {
      os << std::endl;
      os << "  iterCG - Number of Truncated CG iterations" << std::endl << std::endl;
      os << "  flagGC - Trust-Region Truncated CG flag" << std::endl;
      for( auto flag = ConjugateGradients<Real>::Flag::Success; 
                flag != ConjugateGradients<Real>::Flag::Undefined; ++flag ) {
        os << "    " << flag << " - " << ConjugateGradients<Real>::flag_dict[flag] << std::endl;
      }            
    }

    else if( type_ == TrustRegion<Real>::Type::SPG ) {
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

  if ( type_ == TrustRegion<Real>::Type::TruncatedCG ) {
    os << std::setw(10) << std::left << "iterCG";
    os << std::setw(10) << std::left << "flagCG";
  }
  else if (type_ == TrustRegion<Real>::Type::SPG) {
    os << std::setw(10) << std::left << "iterSPG";
    os << std::setw(10) << std::left << "flagSPG";
  }
  os << std::endl;
}

template<typename Real>
void TrustRegionAlgorithm<Real>::writeName( std::ostream& os ) const {

  os << std::endl << TrustRegion<Real>::type_dict[type_] << " Trust-Region Solver";

  if ( useSecantPrecond_ || useSecantHessVec_ ) {
    if ( useSecantPrecond_ && !useSecantHessVec_ ) {
      os << " with " << Secant<Real>::type_dict[secantType_] << " Preconditioning" << std::endl;
    }
    else if ( !useSecantPrecond_ && useSecantHessVec_ ) {
      os << " with " << Secant<Real>::type_dict[secantType_] << " Hessian Approximation" << std::endl;
    }
    else {
      os << " with " << Secant<Real>::type_dict[secantType_] << " Preconditioning and Hessian Approximation" << std::endl;
    }
  }
  else os << std::endl;
}


template<typename Real>
void TrustRegionAlgorithm<Real>::writeOutput( std::ostream& os, 
                                              bool          print_header) const {
  auto& state = Algorithm<Real>::getState();
  os << std::scientific << std::setprecision(6);
  if( state.iter_ == 0 ) writeName(os);
  if( print_header ) writeHeader(os);
  if( state.iter_ == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state.iter_;
    os << std::setw(15) << std::left << state.value_;
    os << std::setw(15) << std::left << state.gnorm_;
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << state.searchSize_;
    os << std::setw(10) << std::left << state.nfval_;
    os << std::setw(10) << std::left << state.ngrad_;
    os << std::setw(10) << std::left << "---";

    if( type_ == TrustRegion<Real>::Type::TruncatedCG || 
        type_ == TrustRegion<Real>::Type::SPG ) {
      os << std::setw(10) << std::left << "---";
      os << std::setw(10) << std::left << "---";
    }
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state.iter_;
    os << std::setw(15) << std::left << state.value_;
    os << std::setw(15) << std::left << state.gnorm_;
    os << std::setw(15) << std::left << state.snorm_;
    os << std::setw(15) << std::left << state.searchSize_;
    os << std::setw(10) << std::left << state.nfval_;
    os << std::setw(10) << std::left << state.ngrad_;
    os << std::setw(10) << std::left << flag_;

    if( type_ == TrustRegion<Real>::Type::TruncatedCG || type_ == TrustRegion<Real>::Type::SPG ) {
      os << std::setw(10) << std::left << SPiter_;
      os << std::setw(10) << std::left << SPflag_;
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Set/Get component methods

template<class Real>
void TrustRegionAlgorithm<Real>::setModel( const Ptr<TrustRegionModel<Real>>& model ) {
  model_ = model;
}

template<class Real>
void TrustRegionAlgorithm<Real>::setSolver( const Ptr<TrustRegion<Real>>& solver ) {
  solver_ = solver;
}

template<class Real>
TrustRegionModel<Real>& TrustRegionAlgorithm<Real>::getModel() { return *model_; }

template<class Real>
const TrustRegionModel<Real>& TrustRegionAlgorithm<Real>::getModel() const { return *model_; }

template<class Real>
TrustRegion<Real>& TrustRegionAlgorithm<Real>::getSolver() { return *solver_; }

template<class Real>
const TrustRegion<Real>& TrustRegionAlgorithm<Real>::getSolver() const { return *solver_; }

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_TRUSTREGIONALGORITHM_DEF_H

