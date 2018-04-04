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
#include "ROL_TrustRegionTypes.hpp"
#include "ROL_TrustRegionModel.hpp"
#include "ROL_ColemanLiModel.hpp"
#include "ROL_KelleySachsModel.hpp"

namespace ROL {

template<class Real>
class TrustRegion {
private:

  ROL::Ptr<Vector<Real> > prim_, dual_;

  ETrustRegionModel TRmodel_;

  Real delmax_;
  Real eta0_, eta1_, eta2_;
  Real gamma0_, gamma1_, gamma2_;
  Real pRed_;
  Real TRsafe_, eps_;
  Real mu0_;

  std::vector<bool> useInexact_;

  Real ftol_old_;

  Real scale_, omega_, force_, forceFactor_;
  int updateIter_, cnt_;

  unsigned verbosity_;

public:

  virtual ~TrustRegion() {}

  // Constructor
  TrustRegion( Teuchos::ParameterList &parlist )
    : ftol_old_(ROL_OVERFLOW<Real>()), cnt_(0), verbosity_(0) {
    // Trust-Region Parameters
    Teuchos::ParameterList list = parlist.sublist("Step").sublist("Trust Region");
    TRmodel_ = StringToETrustRegionModel(list.get("Subproblem Model", "Kelley-Sachs"));
    delmax_  = list.get("Maximum Radius",                       static_cast<Real>(5000.0));
    eta0_    = list.get("Step Acceptance Threshold",            static_cast<Real>(0.05));
    eta1_    = list.get("Radius Shrinking Threshold",           static_cast<Real>(0.05));
    eta2_    = list.get("Radius Growing Threshold",             static_cast<Real>(0.9));
    gamma0_  = list.get("Radius Shrinking Rate (Negative rho)", static_cast<Real>(0.0625));
    gamma1_  = list.get("Radius Shrinking Rate (Positive rho)", static_cast<Real>(0.25));
    gamma2_  = list.get("Radius Growing Rate",                  static_cast<Real>(2.5));
    mu0_     = list.get("Sufficient Decrease Parameter",        static_cast<Real>(1.e-4));
    TRsafe_  = list.get("Safeguard Size",                       static_cast<Real>(100.0));
    eps_     = TRsafe_*ROL_EPSILON<Real>();
    // General Inexactness Information
    Teuchos::ParameterList &glist = parlist.sublist("General");
    useInexact_.clear();
    useInexact_.push_back(glist.get("Inexact Objective Function",     false));
    useInexact_.push_back(glist.get("Inexact Gradient",               false));
    useInexact_.push_back(glist.get("Inexact Hessian-Times-A-Vector", false));
    // Inexact Function Evaluation Information
    Teuchos::ParameterList &ilist = list.sublist("Inexact").sublist("Value");
    scale_       = ilist.get("Tolerance Scaling",                 static_cast<Real>(1.e-1));
    omega_       = ilist.get("Exponent",                          static_cast<Real>(0.9));
    force_       = ilist.get("Forcing Sequence Initial Value",    static_cast<Real>(1.0));
    updateIter_  = ilist.get("Forcing Sequence Update Frequency", static_cast<int>(10));
    forceFactor_ = ilist.get("Forcing Sequence Reduction Factor", static_cast<Real>(0.1));
    // Get verbosity level
    verbosity_ = glist.get("Print Verbosity", 0);
  }

  virtual void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    prim_ = x.clone();
    dual_ = g.clone();
  }

  virtual void update( Vector<Real>           &x,
                       Real                   &fnew,
                       Real                   &del,
                       int                    &nfval,
                       int                    &ngrad,
                       ETrustRegionFlag       &flagTR,
                 const Vector<Real>           &s,
                 const Real                   snorm,
                 const Real                   fold,
                 const Vector<Real>           &g,
                       int                    iter,
                       Objective<Real>        &obj,
                       BoundConstraint<Real>  &bnd,
                       TrustRegionModel<Real> &model ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    const Real one(1), zero(0);

    /***************************************************************************************************/
    // BEGIN INEXACT OBJECTIVE FUNCTION COMPUTATION
    /***************************************************************************************************/
    // Update inexact objective function
    Real fold1 = fold, ftol = tol; // TOL(1.e-2);
    if ( useInexact_[0] ) {
      if ( !(cnt_%updateIter_) && (cnt_ != 0) ) {
        force_ *= forceFactor_;
      }
      //const Real oe4(1e4);
      //Real c = scale_*std::max(TOL,std::min(one,oe4*std::max(pRed_,std::sqrt(ROL_EPSILON<Real>()))));
      //ftol   = c*std::pow(std::min(eta1_,one-eta2_)
      //          *std::min(std::max(pRed_,std::sqrt(ROL_EPSILON<Real>())),force_),one/omega_);
      //if ( ftol_old_ > ftol || cnt_ == 0 ) {
      //  ftol_old_ = ftol;
      //  fold1 = obj.value(x,ftol_old_);
      //}
      //cnt_++;
      Real eta = static_cast<Real>(0.999)*std::min(eta1_,one-eta2_);
      ftol     = scale_*std::pow(eta*std::min(pRed_,force_),one/omega_);
      ftol_old_ = ftol;
      fold1 = obj.value(x,ftol_old_);
      cnt_++;
    }
    // Evaluate objective function at new iterate
    prim_->set(x); prim_->plus(s);
    obj.update(*prim_,true);
    fnew = obj.value(*prim_,ftol);

    nfval = 1;
    Real aRed = fold1 - fnew;
    /***************************************************************************************************/
    // FINISH INEXACT OBJECTIVE FUNCTION COMPUTATION
    /***************************************************************************************************/

    /***************************************************************************************************/
    // BEGIN COMPUTE RATIO OF ACTUAL AND PREDICTED REDUCTION
    /***************************************************************************************************/
    // Modify Actual and Predicted Reduction According to Model
    model.updateActualReduction(aRed,s);
    model.updatePredictedReduction(pRed_,s);

    if ( verbosity_ > 0 ) {
      std::cout << std::endl;
      std::cout << "  Computation of actual and predicted reduction" << std::endl;
      std::cout << "    Current objective function value:        "   << fold1 << std::endl;
      std::cout << "    New objective function value:            "   << fnew  << std::endl;
      std::cout << "    Actual reduction:                        "   << aRed  << std::endl;
      std::cout << "    Predicted reduction:                     "   << pRed_ << std::endl;
    }

    // Compute Ratio of Actual and Predicted Reduction
    Real EPS = eps_*((one > std::abs(fold1)) ? one : std::abs(fold1));
    Real aRed_safe = aRed + EPS, pRed_safe = pRed_ + EPS;
    Real rho(0);
    if (((std::abs(aRed_safe) < eps_) && (std::abs(pRed_safe) < eps_)) || aRed == pRed_) {
      rho = one;
      flagTR = TRUSTREGION_FLAG_SUCCESS;
    }
    else if ( std::isnan(aRed_safe) || std::isnan(pRed_safe) ) {
      rho = -one;
      flagTR = TRUSTREGION_FLAG_NAN;
    }
    else {
      rho = aRed_safe/pRed_safe;
      if (pRed_safe < zero && aRed_safe > zero) {
        flagTR = TRUSTREGION_FLAG_POSPREDNEG;
      }
      else if (aRed_safe <= zero && pRed_safe > zero) {
        flagTR = TRUSTREGION_FLAG_NPOSPREDPOS;
      }
      else if (aRed_safe <= zero && pRed_safe < zero) {
        flagTR = TRUSTREGION_FLAG_NPOSPREDNEG;
      }
      else {
        flagTR = TRUSTREGION_FLAG_SUCCESS;
      }
    }

    if ( verbosity_ > 0 ) {
      std::cout << "    Safeguard:                               " << eps_      << std::endl;
      std::cout << "    Actual reduction with safeguard:         " << aRed_safe << std::endl;
      std::cout << "    Predicted reduction with safeguard:      " << pRed_safe << std::endl;
      std::cout << "    Ratio of actual and predicted reduction: " << rho       << std::endl;
      std::cout << "    Trust-region flag:                       " << flagTR    << std::endl;
    }
    /***************************************************************************************************/
    // FINISH COMPUTE RATIO OF ACTUAL AND PREDICTED REDUCTION
    /***************************************************************************************************/


    /***************************************************************************************************/
    // BEGIN CHECK SUFFICIENT DECREASE FOR BOUND CONSTRAINED PROBLEMS
    /***************************************************************************************************/
    bool decr = true;
    if ( bnd.isActivated() && TRmodel_ == TRUSTREGION_MODEL_KELLEYSACHS ) {
      if ( rho >= eta0_ && (std::abs(aRed_safe) > eps_) ) {
        // Compute Criticality Measure || x - P( x - g ) ||
        prim_->set(x);
        prim_->axpy(-one,g.dual());
        bnd.project(*prim_);
        prim_->scale(-one);
        prim_->plus(x);
        Real pgnorm = prim_->norm();
        // Compute Scaled Measure || x - P( x - lam * PI(g) ) ||
        prim_->set(g.dual());
        bnd.pruneActive(*prim_,g,x);
        Real lam = std::min(one, del/prim_->norm());
        prim_->scale(-lam);
        prim_->plus(x);
        bnd.project(*prim_);
        prim_->scale(-one);
        prim_->plus(x);
        pgnorm *= prim_->norm();
        // Sufficient decrease?
        decr = ( aRed_safe >= mu0_*eta0_*pgnorm );
        flagTR = (!decr ? TRUSTREGION_FLAG_QMINSUFDEC : flagTR);

        if ( verbosity_ > 0 ) {
          std::cout << "    Decrease lower bound (constraints):      " << 0.1*eta0_*pgnorm  << std::endl;
          std::cout << "    Trust-region flag (constraints):         " << flagTR            << std::endl;
          std::cout << "    Is step feasible:                        " << bnd.isFeasible(x) << std::endl;
        }
      }
    }
    /***************************************************************************************************/
    // FINISH CHECK SUFFICIENT DECREASE FOR BOUND CONSTRAINED PROBLEMS
    /***************************************************************************************************/

    /***************************************************************************************************/
    // BEGIN STEP ACCEPTANCE AND TRUST REGION RADIUS UPDATE
    /***************************************************************************************************/
    if ( verbosity_ > 0 ) {
      std::cout << "    Norm of step:                            " << snorm << std::endl;
      std::cout << "    Trust-region radius before update:       " << del   << std::endl;
    }
    if ((rho < eta0_ && flagTR == TRUSTREGION_FLAG_SUCCESS) || flagTR >= 2 || !decr ) { // Step Rejected
      fnew = fold1;
      if (rho < zero) { // Negative reduction, interpolate to find new trust-region radius
        Real gs(0);
        if ( bnd.isActivated() ) {
          model.dualTransform(*dual_, *model.getGradient());
          gs = dual_->dot(s.dual());
        }
        else {
          gs = g.dot(s.dual());
        }
        Real modelVal = model.value(s,tol);
        modelVal += fold1;
        Real theta = (one-eta2_)*gs/((one-eta2_)*(fold1+gs)+eta2_*modelVal-fnew);
        del = std::min(gamma1_*std::min(snorm,del),std::max(gamma0_,theta)*del);
        if ( verbosity_ > 0 ) {
          std::cout << "    Interpolation model value:               " << modelVal << std::endl;
          std::cout << "    Interpolation step length:               " << theta    << std::endl;
        }
      }
      else { // Shrink trust-region radius
        del = gamma1_*std::min(snorm,del);
      }
      obj.update(x,true,iter);
    }
    else if ((rho >= eta0_ && flagTR != TRUSTREGION_FLAG_NPOSPREDNEG) ||
             (flagTR == TRUSTREGION_FLAG_POSPREDNEG)) { // Step Accepted
      x.plus(s);
      obj.update(x,true,iter);
      if (rho >= eta2_) { // Increase trust-region radius
        del = std::min(gamma2_*del,delmax_);
      }
    }
    if ( verbosity_ > 0 ) {
      std::cout << "    Trust-region radius after update:        " << del << std::endl;
      std::cout << std::endl;
    }
    /***************************************************************************************************/
    // FINISH STEP ACCEPTANCE AND TRUST REGION RADIUS UPDATE
    /***************************************************************************************************/
  }

  virtual void run( Vector<Real>           &s,           // Step (to be computed)
                    Real                   &snorm,       // Step norm (to be computed)
                    int                    &iflag,       // Exit flag (to be computed)
                    int                    &iter,        // Iteration count (to be computed)
                    const Real              del,         // Trust-region radius
                    TrustRegionModel<Real> &model ) = 0; // Trust-region model

  void setPredictedReduction(const Real pRed) {
    pRed_ = pRed;
  }

  Real getPredictedReduction(void) const {
    return pRed_;
  }
};

}

#include "ROL_TrustRegionFactory.hpp"

#endif
