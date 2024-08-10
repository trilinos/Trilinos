// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINESEARCH_H
#define ROL_LINESEARCH_H

/** \class ROL::LineSearch
    \brief Provides interface for and implements line searches.
*/

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_ScalarFunction.hpp"

namespace ROL { 

template<class Real>
class LineSearch {
private:

  ECurvatureCondition econd_;
  EDescent            edesc_;

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

  ROL::Ptr<Vector<Real> > xtst_; 
  ROL::Ptr<Vector<Real> > d_;
  ROL::Ptr<Vector<Real> > g_;
  ROL::Ptr<Vector<Real> > grad_;
//  ROL::Ptr<const Vector<Real> > grad_;

public:


  virtual ~LineSearch() {}

  // Constructor
  LineSearch( ROL::ParameterList &parlist ) : eps_(0) {
    Real one(1), p9(0.9), p6(0.6), p4(0.4), oem4(1.e-4), zero(0);
    // Enumerations
    edesc_ = StringToEDescent(parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").get("Type","Quasi-Newton Method"));
    econd_ = StringToECurvatureCondition(parlist.sublist("Step").sublist("Line Search").sublist("Curvature Condition").get("Type","Strong Wolfe Conditions"));
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
 
    fmin_      = std::numeric_limits<Real>::max();
    alphaMin_  = 0; 
    itcond_    = false;

    c1_ = ((c1_ < zero) ? oem4 : c1_);
    c2_ = ((c2_ < zero) ? p9   : c2_);
    c3_ = ((c3_ < zero) ? p9   : c3_);
    if ( c2_ <= c1_ ) {
      c1_ = oem4;
      c2_ = p9;
    }
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      c2_ = p4;
      c3_ = std::min(one-c2_,c3_);
    }
  }


  virtual void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                           Objective<Real> &obj, BoundConstraint<Real> &con ) {
    grad_ = g.clone();
    xtst_ = x.clone();
    d_    = s.clone();
    g_    = g.clone();
  }

  virtual void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                    const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                    Objective<Real> &obj, BoundConstraint<Real> &con ) = 0;

  // Set epsilon for epsilon active sets
  void setData(Real &eps, const Vector<Real> &g) {
    eps_ = eps;
    grad_->set(g);
  }

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
  virtual bool status( const ELineSearch type, int &ls_neval, int &ls_ngrad, const Real alpha, 
                       const Real fold, const Real sgold, const Real fnew, 
                       const Vector<Real> &x, const Vector<Real> &s, 
                       Objective<Real> &obj, BoundConstraint<Real> &con ) { 
    Real tol = std::sqrt(ROL_EPSILON<Real>()), one(1), two(2);

    // Check Armijo Condition
    bool armijo = false;
    if ( con.isActivated() ) {
      Real gs(0);
      if ( edesc_ == DESCENT_STEEPEST ) {
        updateIterate(*d_,x,s,alpha,con);
        d_->scale(-one);
        d_->plus(x);
        gs = -s.dot(*d_);
      }
      else {
        d_->set(s);
        d_->scale(-one);
        con.pruneActive(*d_,grad_->dual(),x,eps_);
        gs = alpha*(grad_)->dot(d_->dual());
        d_->zero();
        updateIterate(*d_,x,s,alpha,con);
        d_->scale(-one);
        d_->plus(x);
        con.pruneInactive(*d_,grad_->dual(),x,eps_);
        gs += d_->dot(grad_->dual());
      }
      if ( fnew <= fold - c1_*gs ) {
        armijo = true;
      }
    }
    else {
      if ( fnew <= fold + c1_*alpha*sgold ) {
        armijo = true;
      }
    }

    // Check Maximum Iteration
    itcond_ = false;
    if ( ls_neval >= maxit_ ) { 
      itcond_ = true;
    }

    // Check Curvature Condition
    bool curvcond = false;
    if ( armijo && ((type != LINESEARCH_BACKTRACKING && type != LINESEARCH_CUBICINTERP) ||
                    (edesc_ == DESCENT_NONLINEARCG)) ) {
      if (econd_ == CURVATURECONDITION_GOLDSTEIN) {
        if (fnew >= fold + (one-c1_)*alpha*sgold) {
          curvcond = true;
        }
      }
      else if (econd_ == CURVATURECONDITION_NULL) {
        curvcond = true;
      }
      else {
        updateIterate(*xtst_,x,s,alpha,con);
        obj.update(*xtst_);
        obj.gradient(*g_,*xtst_,tol);
        Real sgnew(0);
        if ( con.isActivated() ) {
          d_->set(s);
          d_->scale(-alpha);
          con.pruneActive(*d_,s,x);
          sgnew = -d_->dot(g_->dual());
        }
        else {
          sgnew = s.dot(g_->dual());
        }
        ls_ngrad++;
   
        if (    ((econd_ == CURVATURECONDITION_WOLFE)       
                     && (sgnew >= c2_*sgold))
             || ((econd_ == CURVATURECONDITION_STRONGWOLFE) 
                     && (std::abs(sgnew) <= c2_*std::abs(sgold)))
             || ((econd_ == CURVATURECONDITION_GENERALIZEDWOLFE) 
                     && (c2_*sgold <= sgnew && sgnew <= -c3_*sgold))
             || ((econd_ == CURVATURECONDITION_APPROXIMATEWOLFE) 
                     && (c2_*sgold <= sgnew && sgnew <= (two*c1_ - one)*sgold)) ) {
          curvcond = true;
        }
      }
    }

    if(fnew<fmin_) {
      fmin_ = fnew;
      alphaMin_ = alpha;
    }

    if (type == LINESEARCH_BACKTRACKING || type == LINESEARCH_CUBICINTERP) {
      if (edesc_ == DESCENT_NONLINEARCG) {
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
                               Objective<Real> &obj, BoundConstraint<Real> &con) {
    Real val(1), one(1), half(0.5);
    if (useralpha_ || usePrevAlpha_ ) {
      val = alpha0_;
    }
    else {
      if (edesc_ == DESCENT_STEEPEST || edesc_ == DESCENT_NONLINEARCG) {
        Real tol = std::sqrt(ROL_EPSILON<Real>());
        // Evaluate objective at x + s
        updateIterate(*d_,x,s,one,con);
        obj.update(*d_);
        Real fnew = obj.value(*d_,tol);
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

  void updateIterate(Vector<Real> &xnew, const Vector<Real> &x, const Vector<Real> &s, Real alpha,
                     BoundConstraint<Real> &con ) {

    xnew.set(x);
    xnew.axpy(alpha,s);

    if ( con.isActivated() ) {
      con.project(xnew);
    }

  }

  bool useLocalMinimizer() {
    return itcond_ && acceptMin_;
  }
 
  bool takeNoStep() {
    return itcond_ && !acceptMin_;
  }
};

}

#include "ROL_LineSearchFactory.hpp"

#endif
