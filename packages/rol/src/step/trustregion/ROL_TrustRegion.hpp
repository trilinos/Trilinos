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

#ifndef ROL_TRUSTREGION_H
#define ROL_TRUSTREGION_H

/** \class ROL::TrustRegion
    \brief Provides interface for and implements trust-region subproblem solvers.
*/

#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

namespace ROL {

template<class Real>
class TrustRegion {
private:

  Teuchos::RCP<Vector<Real> > xupdate_;
  Teuchos::RCP<Vector<Real> > Hs_;

  Real delmax_;
  Real eta0_;
  Real eta1_;
  Real eta2_;
  Real gamma0_;
  Real gamma1_;
  Real gamma2_;

  Real pRed_;

  Real TRsafe_;
  Real eps_;

  std::vector<bool> useInexact_;

  Real ftol_old_;

  Real scale_;
  Real omega_;
  Real force_;
  int  updateIter_;
  int  forceFactor_;
  int cnt_;

  unsigned verbosity_;

  void updateObj( Vector<Real> &x, int iter, ProjectedObjective<Real> &pObj ) {
    pObj.update(x,true,iter);
  }


public:

  virtual ~TrustRegion() {}

  // Constructor
  TrustRegion( Teuchos::ParameterList & parlist )
    : ftol_old_(ROL_OVERFLOW<Real>()), cnt_(0), verbosity_(0) {
    // Unravel Parameter List
    // Trust-Region Parameters
    Teuchos::ParameterList list = parlist.sublist("Step").sublist("Trust Region");
    delmax_ = list.get("Maximum Radius",5000.0);
    eta0_   = list.get("Step Acceptance Threshold",0.05);
    eta1_   = list.get("Radius Shrinking Threshold",0.05);
    eta2_   = list.get("Radius Growing Threshold",0.9);
    gamma0_ = list.get("Radius Shrinking Rate (Negative rho)",0.0625);
    gamma1_ = list.get("Radius Shrinking Rate (Positive rho)",0.25);
    gamma2_ = list.get("Radius Growing Rate",2.5);
    TRsafe_ = list.get("Safeguard Size",100.0);
    eps_    = TRsafe_*ROL_EPSILON<Real>();

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.sublist("General").get("Inexact Objective Function", false));
    useInexact_.push_back(parlist.sublist("General").get("Inexact Gradient", false));
    useInexact_.push_back(parlist.sublist("General").get("Inexact Hessian-Times-A-Vector", false));
    scale_       = list.sublist("Inexact").sublist("Value").get("Tolerance Scaling",1.e-1);
    omega_       = list.sublist("Inexact").sublist("Value").get("Exponent",0.9);
    force_       = list.sublist("Inexact").sublist("Value").get("Forcing Sequence Initial Value",1.0);
    updateIter_  = list.sublist("Inexact").sublist("Value").get("Forcing Sequence Update Frequency",10);
    forceFactor_ = list.sublist("Inexact").sublist("Value").get("Forcing Sequence Reduction Factor",0.1);

  }

  virtual void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    xupdate_ = x.clone();
    Hs_      = g.clone();
  }

  virtual void update( Vector<Real>      &x,
                       Real              &fnew,
                       Real              &del,
                       int               &nfval,
                       int               &ngrad,
                       ETrustRegionFlag  &flagTR,
                 const Vector<Real>      &s,
                 const Real              snorm,
                 const Real              fold,
                 const Vector<Real>      &g,
                       int               iter,
                       ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>()), one(1), oe4(1.e4), half(0.5), zero(0), p1(0.1);

    // Compute updated iterate vector
    xupdate_->set(x);
    xupdate_->axpy(one,s);
    /***************************************************************************************************/
    // BEGIN OBJECTIVE FUNCTION COMPUTATION
    /***************************************************************************************************/
    // Update inexact objective function
    Real fold1 = fold, ftol = tol, TOL(1.e-2);
    if ( useInexact_[0] ) {
      if ( !(cnt_%updateIter_) && (cnt_ != 0) ) {
        force_ *= forceFactor_;
      }
      Real c = scale_*std::max(TOL,std::min(one,oe4*std::max(pRed_,std::sqrt(ROL_EPSILON<Real>()))));
      ftol   = c*std::pow(std::min(eta1_,one-eta2_)
                *std::min(std::max(pRed_,std::sqrt(ROL_EPSILON<Real>())),force_),one/omega_);
      if ( ftol_old_ > ftol || cnt_ == 0 ) {
        ftol_old_ = ftol;
        fold1 = pObj.value(x,ftol_old_);
      }
      cnt_++;
    }
    // Evaluate objective function at new iterat
    //updateObj(*xupdate_,iter,pObj);
    //pObj.update(*xupdate_,true,iter);
    pObj.update(*xupdate_);
    fnew = pObj.value(*xupdate_,ftol);
    pObj.update(x);

    nfval = 1;
    Real aRed = fold1 - fnew;
    /***************************************************************************************************/
    // FINISH OBJECTIVE FUNCTION COMPUTATION
    /***************************************************************************************************/
    // If constraints are turned on, then compute a different predicted reduction
    if (pObj.isConActivated()) {
      xupdate_->set(x);
      xupdate_->axpy(-one,g.dual());
      pObj.project(*xupdate_);
      xupdate_->axpy(-one,x);
      xupdate_->scale(-one);
 
      pObj.reducedHessVec(*Hs_,s,x,xupdate_->dual(),x,tol);
      pRed_  = -half*s.dot(Hs_->dual());

      Hs_->set(g);
      pObj.pruneActive(*Hs_,xupdate_->dual(),x);
      pRed_ -= s.dot(Hs_->dual());
    }

    if ( verbosity_ > 0 ) {
      std::cout << std::endl;
      std::cout << "  Computation of actual and predicted reduction" << std::endl;
      std::cout << "    Current objective function value:        " << fold1 << std::endl;
      std::cout << "    New objective function value:            " << fnew  << std::endl;
      std::cout << "    Actual reduction:                        " << aRed  << std::endl;
      std::cout << "    Predicted reduction:                     " << pRed_ << std::endl;
    }

    // Compute Ratio of Actual and Predicted Reduction
    aRed  -= eps_*((one > std::abs(fold1)) ? one : std::abs(fold1));
    pRed_ -= eps_*((one > std::abs(fold1)) ? one : std::abs(fold1));
    Real rho(0);
    if ((std::abs(aRed) < eps_) && (std::abs(pRed_) < eps_)) {
      rho = one;
      flagTR = TRUSTREGION_FLAG_SUCCESS;
    }
    else if ( std::isnan(aRed) || std::isnan(pRed_) ) {
      rho = -one;
      flagTR = TRUSTREGION_FLAG_NAN;
    }
    else {
      rho = aRed/pRed_;
      if (pRed_ < zero && aRed > zero) {
        flagTR = TRUSTREGION_FLAG_POSPREDNEG;
      }
      else if (aRed <= zero && pRed_ > zero) {
        flagTR = TRUSTREGION_FLAG_NPOSPREDPOS;
      }
      else if (aRed <= zero && pRed_ < zero) {
        flagTR = TRUSTREGION_FLAG_NPOSPREDNEG;
      }
      else {
        flagTR = TRUSTREGION_FLAG_SUCCESS;
      }
    }

    if ( verbosity_ > 0 ) {
      std::cout << "    Actual reduction with safeguard:         " << aRed   << std::endl;
      std::cout << "    Predicted reduction with safeguard:      " << pRed_  << std::endl;
      std::cout << "    Ratio of actual and predicted reduction: " << rho    << std::endl;
      std::cout << "    Trust-region flag:                       " << flagTR << std::endl;
    }

    // Check Sufficient Decrease in the Reduced Quadratic Model
    bool decr = true;
    if ( pObj.isConActivated() && (std::abs(aRed) > eps_) ) {
      // Compute Criticality Measure || x - P( x - g ) ||
      xupdate_->set(x);
      xupdate_->axpy(-one,g.dual());
      pObj.project(*xupdate_);
      xupdate_->scale(-one);
      xupdate_->plus(x);
      Real pgnorm = xupdate_->norm();
      // Compute Scaled Measure || x - P( x - lam * PI(g) ) ||
      xupdate_->set(g.dual());
      pObj.pruneActive(*xupdate_,g,x);
      Real lam = std::min(one, del/xupdate_->norm());
      xupdate_->scale(-lam);
      xupdate_->plus(x);
      pObj.project(*xupdate_);
      xupdate_->scale(-one);
      xupdate_->plus(x);
      pgnorm *= xupdate_->norm();
      // Sufficient decrease?
      decr = ( aRed >= p1*eta0_*pgnorm );
      flagTR = (!decr ? TRUSTREGION_FLAG_QMINSUFDEC : flagTR);

      if ( verbosity_ > 0 ) {
        std::cout << "    Decrease lower bound (constraints):      " << 0.1*eta0_*pgnorm << std::endl;
        std::cout << "    Trust-region flag (constraints):         " << flagTR << std::endl;
        std::cout << "    Is step feasible:                        " << pObj.isFeasible(x) << std::endl;
      }
    }

    if ( verbosity_ > 0 ) {
      std::cout << std::endl;
    }

    // Accept or Reject Step and Update Trust Region
    if ((rho < eta0_ && flagTR == TRUSTREGION_FLAG_SUCCESS) ||
         flagTR >= 2 || !decr ) { // Step Rejected
      //updateObj(x,iter,pObj);
      //pObj.update(x,true,iter);
      fnew = fold1;
      if (rho < zero) { // Negative reduction, interpolate to find new trust-region radius
        Real gs = s.dot(g.dual());
        pObj.hessVec(*Hs_,s,x,tol);
        Real modelVal = s.dot(Hs_->dual());
        modelVal *= half;
        modelVal += gs + fold1;
        Real theta = (one-eta2_)*gs/((one-eta2_)*(fold1+gs)+eta2_*modelVal-fnew);
        del = std::min(gamma1_*snorm,std::max(gamma0_,theta)*del);
      }
      else { // Shrink trust-region radius
        del = gamma1_*snorm;
      }
    }
    else if ((rho >= eta0_ && flagTR != TRUSTREGION_FLAG_NPOSPREDNEG) ||
      flagTR == TRUSTREGION_FLAG_POSPREDNEG) { // Step Accepted
      x.axpy(one,s);
      pObj.update(x,true,iter);
      if (rho >= eta2_) { // Increase trust-region radius
        del = std::min(gamma2_*del,delmax_);
      }
    }
  }

  void setPredictedReduction(const Real pRed) {
    pRed_ = pRed;
  }
  Real getPredictedReduction(void) const {
    return pRed_;
  }

  virtual void run( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                    const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) = 0;
#if 0
    // If constraints are active, then determine a feasible step
    if ( pObj.isConActivated() && etr_ != TRUSTREGION_CAUCHYPOINT ) {
      // Compute projected step stmp = P( x + alpha*s ) - x and xtmp = x + stmp
      Real alpha = 1.0;
      Teuchos::RCP<Vector<Real> > stmp = x.clone();
      stmp->set(s);
      stmp->scale(alpha);
      pObj.computeProjectedStep(*stmp,x);
      Teuchos::RCP<Vector<Real> > xtmp = x.clone();
      xtmp->set(x);
      xtmp->axpy(1.0,*stmp);
      // Compute model components for alpha = 1.0
      Real tol   = std::sqrt(ROL_EPSILON<Real>());
      Teuchos::RCP<Vector<Real> > Bs = x.clone();
      pObj.hessVec(*Bs,*stmp,x,tol);
      Real sBs   = Bs->dot(*stmp);
      Real gs    = grad.dot(*stmp);
      Real val   = gs + 0.5*sBs;
      Real val0  = val;
      // Backtrack alpha until x+alpha*s is feasible
      int cnt   = 0;
      int maxit = 10;
      while ( val > val0 || !pObj.isFeasible(*xtmp) ) { 
        // Backtrack alpha
        alpha *= 0.5;
        // Computed projected step P( x + alpha*s ) - x and xtmp = x + stmp
        stmp->set(s);
        stmp->scale(alpha);
        pObj.computeProjectedStep(*stmp,x);        
        xtmp->set(x);
        xtmp->axpy(1.0,*stmp);
        // Evaluate Model
        val0 = val;
        pObj.hessVec(*Bs,*stmp,x,tol);
        sBs = Bs->dot(*stmp);
        gs  = grad.dot(*stmp);
        val = gs + 0.5*sBs;
        // Update iteration count
        cnt++;
        if ( cnt >= maxit ) { break; }
      }
      s.set(*stmp);
      pRed_ = -val;
    }
#endif
};

}

#include "ROL_TrustRegionFactory.hpp"

#endif
