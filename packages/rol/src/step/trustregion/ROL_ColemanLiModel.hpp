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

#ifndef ROL_COLEMANLIMODEL_HPP
#define ROL_COLEMANLIMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Secant.hpp"

/** @ingroup func_group
    \class ROL::ColemanLiModel
    \brief Provides the interface to evaluate interior trust-region model
    functions from the Coleman-Li bound constrained trust-region algorithm.

    -----
*/

namespace ROL {

template<class Real>
class ColemanLiModel : public TrustRegionModel<Real> {
private:
  ROL::Ptr<BoundConstraint<Real> > bnd_;              // Bound constraint
  ROL::Ptr<Secant<Real> > sec_;                       // Secant storage

  ROL::Ptr<Vector<Real> > prim_, dual_, hv_;          // Auxiliary storage
  ROL::Ptr<Vector<Real> > step_;                      // Step storage
  ROL::Ptr<Vector<Real> > cauchyStep_, cauchyScal_;   // Cauchy point vectors
  ROL::Ptr<Vector<Real> > reflectStep_, reflectScal_; // Reflective step vectors
  ROL::Ptr<Vector<Real> > Dmat_;                      // sqrt(abs(v))
  ROL::Ptr<Vector<Real> > Cmat_;                      // diag(g) * dv/dx

  const bool useSecantPrecond_;                           // Use secant as preconditioner (unused)
  const bool useSecantHessVec_;                           // Use secant as Hessian

  const Real TRradius_, stepBackMax_, stepBackScale_;     // Primal transform parameters
  const bool singleReflect_;                              // Use single reflection
  Real sCs_, pred_;                                       // Actual/predicted reduction

  Elementwise::Multiply<Real> mult_;                      // Elementwise multiply
  Elementwise::Divide<Real>   div_;                       // Elementwise division

  // Apply diagonal D matrix
  void applyD( Vector<Real> &Dv, const Vector<Real> &v ) {
    Dv.set(v);
    Dv.applyBinary(div_,*Dmat_);
  }

  // Apply inverse of diagonal D matrix
  void applyInverseD( Vector<Real> &Dv, const Vector<Real> &v ) {
    Dv.set(v);
    Dv.applyBinary(mult_,*Dmat_);
  }

  // Apply diagonal C matrix
  void applyC( Vector<Real> &Cv, const Vector<Real> &v ) {
    Cv.set(v);
    Cv.applyBinary(mult_, *Cmat_);
  }

  void constructC(void) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > l  = bnd_->getLowerBound();
    const ROL::Ptr<const Vector<Real> > u  = bnd_->getUpperBound();

    // Set Cmat_ to be the sign of the gradient
    Cmat_->set(gc->dual());
    Cmat_->applyUnary(Elementwise::Sign<Real>());
    // If g < 0 and u = INF then set Cmat_ to zero 
    class NegGradInfU : public Elementwise::BinaryFunction<Real> {
    public:
      NegGradInfU(void) {}
      Real apply(const Real &x, const Real &y) const {
        const Real zero(0), one(1), INF(ROL_INF<Real>());
        return (x < zero && y == INF) ? zero : one;
      }
    };
    prim_->set(gc->dual());
    prim_->applyBinary(NegGradInfU(), *u);
    Cmat_->applyBinary(mult_, *prim_);
    // If g >= 0 and l = -INF then set Cmat_ to zero
    class PosGradNinfL : public Elementwise::BinaryFunction<Real> {
    public:
      PosGradNinfL(void) {}
      Real apply(const Real &x, const Real &y) const {
        const Real zero(0), one(1), NINF(ROL_NINF<Real>());
        return (x >= zero && y == NINF) ? zero : one;
      }
    };
    prim_->set(gc->dual());
    prim_->applyBinary(PosGradNinfL(), *l);
    Cmat_->applyBinary(mult_, *prim_);
    // Pointwise multiply Cmat_ with the gradient
    Cmat_->applyBinary(mult_, gc->dual());
  }

  void constructInverseD(void) {
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    const ROL::Ptr<const Vector<Real> > l  = bnd_->getLowerBound();
    const ROL::Ptr<const Vector<Real> > u  = bnd_->getUpperBound();
    const Real zero(0), one(1), INF(ROL_INF<Real>()), NINF(ROL_NINF<Real>());
    const int LESS_THAN    = 0;
    const int EQUAL_TO     = 1;
    const int GREATER_THAN = 2;
    
    Dmat_->zero();
    // CASE (i)
    // Mask for negative gradient (m1 is 1 if g is negative and 0 otherwise)
    reflectStep_->applyBinary(Elementwise::ValueSet<Real>(zero, LESS_THAN),gc->dual());
    // Mask for finite upper bounds (m2 is 1 if u is finite and 0 otherwise)
    reflectScal_->applyBinary(Elementwise::ValueSet<Real>(INF, LESS_THAN),*u);
    // Mask for g_i < 0 and u_i < inf
    reflectScal_->applyBinary(mult_,*reflectStep_);
    // prim_i = { u_i-x_i if g_i < 0 and u_i < inf
    //          { 0       otherwise
    prim_->set(*u); prim_->axpy(-one,*xc);
    prim_->applyBinary(mult_,*reflectScal_);
    // Add to D
    Dmat_->plus(*prim_);

    // CASE (iii)
    // Mask for infinite upper bounds
    reflectScal_->applyBinary(Elementwise::ValueSet<Real>(INF, EQUAL_TO),*u);
    // Mask for g_i < 0 and u_i = inf
    reflectScal_->applyBinary(mult_,*reflectStep_);
    // prim_i = { -1 if g_i < 0 and u_i = inf
    //          { 0  otherwise
    prim_->applyUnary(Elementwise::Fill<Real>(-one)); 
    prim_->applyBinary(mult_,*reflectScal_);
    // Add to D
    Dmat_->plus(*prim_);

    // CASE (ii)
    // m1 = 1-m1
    reflectStep_->scale(-one);
    reflectStep_->applyUnary(Elementwise::Shift<Real>(one));
    // Mask for finite lower bounds
    reflectScal_->applyBinary(Elementwise::ValueSet<Real>(NINF, GREATER_THAN),*l);
    // Zero out elements of Jacobian with l_i=-inf
    reflectScal_->applyBinary(mult_,*reflectStep_);
    // prim_i = { x_i-l_i if g_i >= 0 and l_i > -inf
    //          { 0       otherwise
    prim_->set(*xc); prim_->axpy(-one,*l);
    prim_->applyBinary(mult_,*reflectScal_);
    // Add to D
    Dmat_->plus(*prim_);

    // CASE (iv)
    // Mask for infinite lower bounds
    reflectScal_->applyBinary(Elementwise::ValueSet<Real>(NINF, EQUAL_TO),*l);
    // Mask for g_i>=0 and l_i=-inf
    reflectScal_->applyBinary(mult_,*reflectStep_);
    // prim_i = { 1 if g_i >= 0 and l_i = -inf
    //          { 0 otherwise
    prim_->applyUnary(Elementwise::Fill<Real>(one));
    prim_->applyBinary(mult_,*reflectScal_);
    // Add to D
    Dmat_->plus(*prim_);
  
    // d_i = { u_i-x_i if g_i <  0, u_i<inf
    //       { -1      if g_i <  0, u_i=inf
    //       { x_i-l_i if g_i >= 0, l_i>-inf
    //       { 1       if g_i >= 0, l_i=-inf 
    Dmat_->applyUnary(Elementwise::AbsoluteValue<Real>());
    Dmat_->applyUnary(Elementwise::SquareRoot<Real>());
  }

  // Build diagonal D and C matrices
  void initialize(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g) {
    bnd_ = ROL::makePtrFromRef(bnd);

    prim_ = x.clone();
    dual_ = g.clone();
    hv_   = g.clone();
    step_ = x.clone();
    Dmat_ = x.clone();
    Cmat_ = x.clone();

    cauchyStep_  = x.clone();
    cauchyScal_  = x.clone();
    reflectStep_ = x.clone();
    reflectScal_ = x.clone();

    constructC();
    constructInverseD();
  }

 public:

  ColemanLiModel( Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g)
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      sec_(ROL::nullPtr), useSecantPrecond_(false), useSecantHessVec_(false),
      TRradius_(1), stepBackMax_(0.9999), stepBackScale_(1),
      singleReflect_(true), sCs_(0), pred_(0) {
    initialize(obj,bnd,x,g);
  }

  ColemanLiModel( Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g,
                  const Real TRradius, const Real stepBackMax, const Real stepBackScale,
                  const bool singleReflect = true )
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      sec_(ROL::nullPtr), useSecantPrecond_(false), useSecantHessVec_(false),
      TRradius_(TRradius), stepBackMax_(stepBackMax), stepBackScale_(stepBackScale),
      singleReflect_(singleReflect), sCs_(0), pred_(0) {
    initialize(obj,bnd,x,g);
  }

  ColemanLiModel( Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g,
                  const ROL::Ptr<Secant<Real> > &sec,
                  const bool useSecantPrecond, const bool useSecantHessVec,
                  const Real TRradius, const Real stepBackMax, const Real stepBackScale,
                  const bool singleReflect = true )
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      sec_(sec), useSecantPrecond_(useSecantPrecond), useSecantHessVec_(useSecantHessVec),
      TRradius_(TRradius), stepBackMax_(stepBackMax), stepBackScale_(stepBackScale),
      singleReflect_(singleReflect), sCs_(0), pred_(0) {
    initialize(obj,bnd,x,g);
  }
 
  // Note that s is the \f$\hat{s}\f$ and \f$\psi\f$ is the $\hat\psi$ from the paper
  Real value( const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    // Apply Hessian to s
    hessVec(*hv_, s, s, tol);
    hv_->scale(static_cast<Real>(0.5));
    // Form inv(D) * g
    applyInverseD(*prim_, gc->dual());
    // Add scaled gradient to Hessian in direction s
    hv_->plus(prim_->dual());
    return hv_->dot(s.dual());    
  }

  void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    hessVec(g, s, s, tol);
    applyInverseD(*prim_, gc->dual());
    g.plus(prim_->dual());    
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    // Build B = inv(D) * Hessian * inv(D)
    applyInverseD(*prim_, v);
    if(useSecantHessVec_) {
      sec_->applyB(*dual_, *prim_);
    }
    else {
      const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
      TrustRegionModel<Real>::getObjective()->hessVec(*dual_, *prim_, *xc, tol);   
    }
    applyInverseD(hv, *dual_);
    // Build C = diag(g) J
    applyC(*prim_, v);
    hv.plus(prim_->dual()); 
  }
  
  void dualTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    applyInverseD(tv, v);
  }

  void primalTransform( Vector<Real> &tiv, const Vector<Real> &v ) { 
    Real tol = std::sqrt(ROL_EPSILON<Real>());

    /**************************************************************************/
    /*      PERFORM OPTIMAL SCALING OF TRUST REGION SUBPROBLEM SOLUTION       */
    /**************************************************************************/
    applyInverseD(tiv, v);
    // Get bounds on scalar variable
    Real lowerBoundV(ROL_NINF<Real>()), upperBoundV(ROL_INF<Real>());
    getScalarBounds(lowerBoundV, upperBoundV, tiv);
    // Minimize one dimensional quadratic over bounds
    Real tauV(1);
    Real valueV = minimize1D(tauV, lowerBoundV, upperBoundV, v);

    /**************************************************************************/
    /*      COMPUTE CAUCHY POINT: STORED IN cauchyStep_ AND cauchyScal_       */
    /**************************************************************************/
    Real valueG = computeCauchyPoint();

    /**************************************************************************/
    /*      COMPUTE REFLECTIVE STEP: STORED IN reflectStep_ AND reflectScal_  */
    /**************************************************************************/
    if ( singleReflect_ ) {
      computeReflectiveStep(*reflectStep_, v, tiv);
    }
    else {
      computeFullReflectiveStep(*reflectStep_, v, tiv);
    }
    applyInverseD(*reflectScal_, *reflectStep_);
    // Get bounds on scalar variable
    Real lowerBoundR(ROL_NINF<Real>()), upperBoundR(ROL_INF<Real>());
    getScalarBounds(lowerBoundR, upperBoundR, *reflectScal_);
    // Minimize one dimensional quadratic over bounds
    Real tauR(1);
    Real valueR = minimize1D(tauR, lowerBoundR, upperBoundR, *reflectStep_);

    /**************************************************************************/
    /*      CHOOSE STEP THAT GIVES MOST PREDICTED REDUCTION                   */
    /**************************************************************************/
    Real VALUE(0);
    bool useCauchyPoint = (valueG < valueV);
    if (useCauchyPoint) {
      VALUE = valueG;
      tiv.set(*cauchyScal_);
      // Store unscaled step
      step_->set(*cauchyStep_);
    }
    else {
      VALUE = valueV;
      tiv.scale(tauV);
      // Store unscaled step
      step_->set(v);
      step_->scale(tauV);
    }
    bool useReflectStep = (valueR < VALUE);
    if (useReflectStep) {
      VALUE = valueR;
      tiv.set(*reflectScal_);
      tiv.scale(tauR);
      // Store unscaled step
      step_->set(*reflectStep_);
      step_->scale(tauR);
    }

    /**************************************************************************/
    /*      ENSURE CHOSEN STEP IS STRICTLY FEASIBLE                           */
    /**************************************************************************/
    // Computed predicted reduction based on input step
    if ( !isStrictlyFeasibleStep(tiv) ) {
      Real snorm = step_->norm();
      Real theta = std::max( stepBackMax_, static_cast<Real>(1) - stepBackScale_ * snorm);
      tiv.scale(theta);
      step_->scale(theta);
      // Compute predicted reduction
      pred_ = -value(*step_,tol);
    }
    else { // Theta is one
      // Compute predicted reduction
      pred_ = -VALUE;
    }

    // Compute update for actual reduction
    applyC(*prim_, *step_);
    sCs_ = static_cast<Real>(-0.5) * prim_->dot(*step_);
  }

  void updatePredictedReduction(Real &pred, const Vector<Real> &s) {
    pred = pred_;
  }

  void updateActualReduction(Real &ared, const Vector<Real> &s) {
    ared += sCs_;
  }

  const ROL::Ptr<BoundConstraint<Real> > getBoundConstraint(void) const {
    return bnd_;
  }

private:

  void getScalarBounds( Real &lowerBound, Real &upperBound, const Vector<Real> &p ) {
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    const ROL::Ptr<const Vector<Real> > l  = bnd_->getLowerBound();
    const ROL::Ptr<const Vector<Real> > u  = bnd_->getUpperBound();
    const Real one(1);
    Real pnorm = p.norm();

    // Define elementwise functions
    class PruneNegative : public Elementwise::BinaryFunction<Real> {
    private:
      const Real val_;
    public:
      PruneNegative( const Real val ) : val_(val) {}
      Real apply(const Real &x, const Real &y) const {
        return (y < static_cast<Real>(0)) ? x/y : val_;
      }
    };
    class PrunePositive : public Elementwise::BinaryFunction<Real> {
    private:
      const Real val_;
    public:
      PrunePositive( const Real val ) : val_(val) {}
      Real apply(const Real &x, const Real &y) const {
        return (y > static_cast<Real>(0)) ? x/y : val_;
      }
    };

    // Max of (l-x)/p if p > 0
    prim_->set(*l); prim_->axpy(-one,*xc);
    prim_->applyBinary(PrunePositive(ROL_NINF<Real>()),p);
    Real lowerBound1 = prim_->reduce(Elementwise::ReductionMax<Real>());
    // Max of (u-x)/p if p < 0
    prim_->set(*u); prim_->axpy(-one,*xc);
    prim_->applyBinary(PruneNegative(ROL_NINF<Real>()),p);
    Real lowerBound2 = prim_->reduce(Elementwise::ReductionMax<Real>());
    // Lower bound
    Real lowerBound3 = std::max(lowerBound1, lowerBound2);

    // Min of (u-x)/p if p > 0
    prim_->set(*u); prim_->axpy(-one,*xc);
    prim_->applyBinary(PrunePositive(ROL_INF<Real>()),p);
    Real upperBound1 = prim_->reduce(Elementwise::ReductionMin<Real>());
    // Max of (l-x)/p if p < 0
    prim_->set(*l); prim_->axpy(-one,*xc);
    prim_->applyBinary(PruneNegative(ROL_INF<Real>()),p);
    Real upperBound2 = prim_->reduce(Elementwise::ReductionMin<Real>());
    // Upper bound
    Real upperBound3 = std::min(upperBound1, upperBound2);

    // Adjust for trust-region constraint
    lowerBound = std::max(lowerBound3, -TRradius_/pnorm);
    upperBound = std::min(upperBound3,  TRradius_/pnorm);
  }

  Real minimize1D(Real &tau, const Real lowerBound, const Real upperBound, const Vector<Real> &p) {
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    Real tol = std::sqrt(ROL_EPSILON<Real>());

    // Compute coefficients of one dimensional quadratic
    hessVec(*hv_, p, p, tol);
    Real c2 = static_cast<Real>(0.5) * hv_->dot(p.dual());
    applyInverseD(*prim_, gc->dual());
    Real c1 = prim_->dot(p);

    // Minimize one dimensional quadratic over bounds
    Real lval = (c2 * lowerBound + c1) * lowerBound;
    Real rval = (c2 * upperBound + c1) * upperBound;
    tau  = (lval < rval) ? lowerBound : upperBound;
    if (c2 > static_cast<Real>(0)) {
      Real uncMin = static_cast<Real>(-0.5) * c1/c2;
      tau = (uncMin > lowerBound && uncMin < upperBound) ? uncMin : tau;
    }

    // Return minimal function value
    return (c2 * tau + c1) * tau;
  }

  Real computeCauchyPoint(void) {
    // Set step = -inv(D) g
    const ROL::Ptr<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    applyInverseD(*cauchyStep_, gc->dual());
    cauchyStep_->scale(static_cast<Real>(-1));

    // Scale Cauchy point
    applyInverseD(*cauchyScal_, *cauchyStep_);

    // Scalar bounds
    Real lowerBound(ROL_NINF<Real>()), upperBound(ROL_INF<Real>());
    getScalarBounds(lowerBound, upperBound, *cauchyScal_);

    // Minimize 1D quadratic over bounds
    Real tau(1), value(0);
    value = minimize1D(tau, lowerBound, upperBound, *cauchyStep_);

    // Scale Cauchy point and return minimal function value
    cauchyStep_->scale(tau);
    cauchyScal_->scale(tau);
    return value;
  }

  void computeReflectiveStep(Vector<Real> &Rv, const Vector<Real> &v, const Vector<Real> &Dv) {
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    Real alpha = computeAlpha(Dv);
    Rv.set(v);

    class LowerBound : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return (x == y) ? static_cast<Real>(-1) : static_cast<Real>(1);
      }
    };
    prim_->set(*xc); prim_->axpy(alpha,Dv);
    prim_->applyBinary(LowerBound(),*bnd_->getLowerBound());
    Rv.applyBinary(mult_,*prim_);

    class UpperBound : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return (x == y) ? static_cast<Real>(-1) : static_cast<Real>(1);
      }
    };
    prim_->set(*xc); prim_->axpy(alpha,Dv);
    prim_->applyBinary(UpperBound(),*bnd_->getUpperBound());
    Rv.applyBinary(mult_,*prim_);
  }

  void computeFullReflectiveStep(Vector<Real> &Rv, const Vector<Real> &v, const Vector<Real> &Dv) {
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    Rv.set(v);

    class LowerBound : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return (x < y) ? static_cast<Real>(-1) : static_cast<Real>(1);
      }
    };
    prim_->set(*xc); prim_->plus(Dv);
    prim_->applyBinary(LowerBound(),*bnd_->getLowerBound());
    Rv.applyBinary(mult_,*prim_);

    class UpperBound : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return (x > y) ? static_cast<Real>(-1) : static_cast<Real>(1);
      }
    };
    prim_->set(*xc); prim_->plus(Dv);
    prim_->applyBinary(UpperBound(),*bnd_->getUpperBound());
    Rv.applyBinary(mult_,*prim_);
  }

  Real computeAlpha( const Vector<Real> &d ) {
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
    ROL::Ptr<Vector<Real> > lx = xc->clone(), ux = xc->clone();
    const Real one(1);

    // Define elementwise functions
    class SafeDivide : public Elementwise::BinaryFunction<Real> {
    private:
      const Real val_;
    public:
      SafeDivide( const Real val ) : val_(val) {}
      Real apply(const Real &x, const Real &y) const {
        const Real zero(0);
        return (y == zero) ? val_ : x/y;
      }
    };

    // (l - x) / d
    lx->set(*bnd_->getLowerBound());
    lx->axpy(-one, *xc);
    lx->applyBinary(SafeDivide(ROL_INF<Real>()), d);

    // (u - x) / d
    ux->set(*bnd_->getUpperBound());
    ux->axpy(-one, *xc);
    ux->applyBinary(SafeDivide(ROL_INF<Real>()), d);

    // max{ (l - x) / d, (u - x) / d }
    lx->applyBinary(Elementwise::Max<Real>(),*ux);

    // min{ max{ (l - x) / d, (u - x) / d } }
    return lx->reduce(Elementwise::ReductionMin<Real>());
  }

  bool isStrictlyFeasibleStep( const Vector<Real> &d ) const {
    const ROL::Ptr<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();

    class Greater : public Elementwise::BinaryFunction<Real> {
    public:
      Greater() {}
      Real apply(const Real &x, const Real &y) const {
        return (x > y) ? static_cast<Real>(1) : static_cast<Real>(0);
      }
    };
    prim_->set(*xc); prim_->plus(d);
    prim_->applyBinary(Greater(),*bnd_->getLowerBound());
    Real lowerFeasible = prim_->reduce(Elementwise::ReductionMin<Real>());

    class Lesser : public Elementwise::BinaryFunction<Real> {
    public:
      Lesser() {}
      Real apply(const Real &x, const Real &y) const {
        return (x < y) ? static_cast<Real>(1) : static_cast<Real>(0);
      }
    };
    prim_->set(*xc); prim_->plus(d);
    prim_->applyBinary(Lesser(),*bnd_->getUpperBound());
    Real upperFeasible = prim_->reduce(Elementwise::ReductionMin<Real>());

    return (upperFeasible * lowerFeasible > 0);
  }

}; // class ColemanLiModel

}

#endif // ROL_COLEMANLIMODEL_HPP
