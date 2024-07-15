// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINESEARCH_U_H
#define ROL_LINESEARCH_U_H

/** \class ROL::LineSearch_U
    \brief Provides interface for and implements line searches.
*/

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Types.hpp"
#include "ROL_Objective.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_LineSearch_U_Types.hpp"

namespace ROL {

template<typename Real>
class LineSearch_U {
private:

  ECurvatureConditionU econd_;
  EDescentU            edesc_;

  bool useralpha_;
  bool usePrevAlpha_; // Use the previous step's accepted alpha as an initial guess
  Real alpha0_;
  Real alpha0bnd_;    // Lower bound for initial alpha...if below, set initial alpha to one
  int maxit_;
  Real c1_;
  Real c2_;
  Real c3_;
  Real eps_;
  Real fmin_;         // smallest fval encountered
  Real alphaMin_;     // Alpha that yields the smallest fval encountered
  bool acceptMin_;    // Use smallest fval if sufficient decrease not satisfied
  bool itcond_;       // true if maximum function evaluations reached
  bool FDdirDeriv_;

  Ptr<Vector<Real>> xtst_; 

  Real dirDeriv(const Vector<Real> &x, // current iterate
                const Vector<Real> &s, // current trial step
                const Real alpha,      // current step length
                const Real fnew,       // f(x+alpha*s)
                Objective<Real> &obj) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    Real val(0);
    xtst_->set(x); xtst_->axpy(alpha,s);
    if (FDdirDeriv_) {
      Real snorm = s.norm();
      if (snorm > static_cast<Real>(0)) {
        Real xnorm   = xtst_->norm();
        Real cbrteps = std::cbrt(ROL_EPSILON<Real>());
        Real h       = cbrteps*std::max(xnorm/snorm,static_cast<Real>(1));
        xtst_->axpy(h,s);
        obj.update(*xtst_,UpdateType::Trial);
        Real ftrial = obj.value(*xtst_,tol);
        val = (ftrial - fnew) / h;
      }
    }
    else {
      val = obj.dirDeriv(*xtst_,s,tol);
    }
    return val;
  }

public:

  virtual ~LineSearch_U() {}

  // Constructor
  LineSearch_U( ParameterList &parlist ) {
    Real one(1), p9(0.9), p6(0.6), p4(0.4), oem4(1.e-4), zero(0);
    // Enumerations
    edesc_ = StringToEDescentU(parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").get("Type","Quasi-Newton Method"));
    econd_ = StringToECurvatureConditionU(parlist.sublist("Step").sublist("Line Search").sublist("Curvature Condition").get("Type","Strong Wolfe Conditions"));
    // Linesearch Parameters
    alpha0_       = parlist.sublist("Step").sublist("Line Search").get("Initial Step Size",one);
    alpha0bnd_    = parlist.sublist("Step").sublist("Line Search").get("Lower Bound for Initial Step Size",one);
    useralpha_    = parlist.sublist("Step").sublist("Line Search").get("User Defined Initial Step Size",false);
    usePrevAlpha_ = parlist.sublist("Step").sublist("Line Search").get("Use Previous Step Length as Initial Guess",false);
    acceptMin_    = parlist.sublist("Step").sublist("Line Search").get("Accept Linesearch Minimizer",false);
    maxit_        = parlist.sublist("Step").sublist("Line Search").get("Function Evaluation Limit",20);
    c1_           = parlist.sublist("Step").sublist("Line Search").get("Sufficient Decrease Tolerance",oem4);
    c2_           = parlist.sublist("Step").sublist("Line Search").sublist("Curvature Condition").get("General Parameter",p9);
    c3_           = parlist.sublist("Step").sublist("Line Search").sublist("Curvature Condition").get("Generalized Wolfe Parameter",p6);
 
    fmin_       = std::numeric_limits<Real>::max();
    alphaMin_   = 0; 
    itcond_     = false;
    FDdirDeriv_ = parlist.sublist("Step").sublist("Line Search").get("Finite Difference Directional Derivative",false);

    c1_ = ((c1_ < zero) ? oem4 : c1_);
    c2_ = ((c2_ < zero) ? p9   : c2_);
    c3_ = ((c3_ < zero) ? p9   : c3_);
    if ( c2_ <= c1_ ) {
      c1_ = oem4;
      c2_ = p9;
    }
    if ( edesc_ == DESCENT_U_NONLINEARCG ) {
      c2_ = p4;
      c3_ = std::min(one-c2_,c3_);
    }
  }

  virtual void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    xtst_ = x.clone();
  }

  virtual void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                    const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                    Objective<Real> &obj ) = 0;

  // use this function to modify alpha and fval if the maximum number of iterations
  // are reached
  void setMaxitUpdate(Real &alpha, Real &fnew, const Real &fold) {
    // Use local minimizer
    if( itcond_ && acceptMin_ ) {
      alpha = alphaMin_;
      fnew = fmin_;
    }
    // Take no step
    else if(itcond_ && !acceptMin_) {
      alpha = 0;
      fnew = fold;
    }
    setNextInitialAlpha(alpha);
  }

protected:
  virtual bool status( const ELineSearchU type, int &ls_neval, int &ls_ngrad, const Real alpha, 
                       const Real fold, const Real sgold, const Real fnew, 
                       const Vector<Real> &x, const Vector<Real> &s, 
                       Objective<Real> &obj ) { 
    const Real one(1), two(2);

    // Check Armijo Condition
    bool armijo = false;
    if ( fnew <= fold + c1_*alpha*sgold ) {
      armijo = true;
    }

    // Check Maximum Iteration
    itcond_ = false;
    if ( ls_neval >= maxit_ ) { 
      itcond_ = true;
    }

    // Check Curvature Condition
    bool curvcond = false;
    if ( armijo && ((type != LINESEARCH_U_BACKTRACKING && type != LINESEARCH_U_CUBICINTERP) ||
                    (edesc_ == DESCENT_U_NONLINEARCG)) ) {
      if (econd_ == CURVATURECONDITION_U_GOLDSTEIN) {
        if (fnew >= fold + (one-c1_)*alpha*sgold) {
          curvcond = true;
        }
      }
      else if (econd_ == CURVATURECONDITION_U_NULL) {
        curvcond = true;
      }
      else {
        Real sgnew = dirDeriv(x,s,alpha,fnew,obj); //ls_ngrad++;
        if (    ((econd_ == CURVATURECONDITION_U_WOLFE)       
                     && (sgnew >= c2_*sgold))
             || ((econd_ == CURVATURECONDITION_U_STRONGWOLFE) 
                     && (std::abs(sgnew) <= c2_*std::abs(sgold)))
             || ((econd_ == CURVATURECONDITION_U_GENERALIZEDWOLFE) 
                     && (c2_*sgold <= sgnew && sgnew <= -c3_*sgold))
             || ((econd_ == CURVATURECONDITION_U_APPROXIMATEWOLFE) 
                     && (c2_*sgold <= sgnew && sgnew <= (two*c1_ - one)*sgold)) ) {
          curvcond = true;
        }
      }
    }

    if(fnew<fmin_) {
      fmin_ = fnew;
      alphaMin_ = alpha;
    }

    if (type == LINESEARCH_U_BACKTRACKING || type == LINESEARCH_U_CUBICINTERP) {
      if (edesc_ == DESCENT_U_NONLINEARCG) {
        return ((armijo && curvcond) || itcond_);
      }
      else {
        return (armijo || itcond_);
      }
    }
    else {
      return ((armijo && curvcond) || itcond_);
    }
  }

  virtual Real getInitialAlpha(int &ls_neval, int &ls_ngrad, const Real fval, const Real gs, 
                               const Vector<Real> &x, const Vector<Real> &s, 
                               Objective<Real> &obj) {
    Real val(1);
    if (useralpha_ || usePrevAlpha_ ) {
      val = alpha0_;
    }
    else {
      const Real one(1), half(0.5);
      if (edesc_ == DESCENT_U_STEEPEST || edesc_ == DESCENT_U_NONLINEARCG) {
        Real tol = std::sqrt(ROL_EPSILON<Real>());
        // Evaluate objective at x + s
        xtst_->set(x); xtst_->plus(s);
        obj.update(*xtst_,UpdateType::Trial);
        Real fnew = obj.value(*xtst_,tol);
        ls_neval++;
        // Minimize quadratic interpolate to compute new alpha
        Real denom = (fnew - fval - gs);
        Real alpha = ((denom > ROL_EPSILON<Real>()) ? -half*gs/denom : one);
        val = ((alpha > alpha0bnd_) ? alpha : one);
      }
      else {
        val = one;
      }
    }
    return val;
  }

  void setNextInitialAlpha( Real alpha ) {
    if( usePrevAlpha_ ) {
      alpha0_  = alpha; 
    }
  }

  bool useLocalMinimizer() {
    return itcond_ && acceptMin_;
  }
 
  bool takeNoStep() {
    return itcond_ && !acceptMin_;
  }
}; // class ROL::LineSearch_U

} // namespace ROL

#endif
