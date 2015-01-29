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

#ifndef ROL_LINESEARCH_H
#define ROL_LINESEARCH_H

/** \class ROL::LineSearch
    \brief Provides interface for and implements line searches.
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL { 

template<class Real>
class LineSearch {
private:

  ECurvatureCondition econd_;
  EDescent            edesc_;

  bool useralpha_;
  Real alpha0_;
  int maxit_;
  Real c1_;
  Real c2_;
  Real c3_;
  Real eps_;

  Teuchos::RCP<Vector<Real> > xtst_; 
  Teuchos::RCP<Vector<Real> > d_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<const Vector<Real> > grad_;

public:

  virtual ~LineSearch() {}

  // Constructor
  LineSearch( Teuchos::ParameterList &parlist ) : eps_(0.0) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    econd_ = StringToECurvatureCondition( parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    // Linesearc Parameters
    maxit_     = parlist.get("Maximum Number of Function Evaluations",            20);
    c1_        = parlist.get("Sufficient Decrease Parameter",                     1.e-4);
    c2_        = parlist.get("Curvature Conditions Parameter",                    0.9);
    c3_        = parlist.get("Curvature Conditions Parameter: Generalized Wolfe", 0.6);
    alpha0_    = parlist.get("Initial Linesearch Parameter",1.0);
    useralpha_ = parlist.get("User Defined Linesearch Parameter",false);

    if ( c1_ < 0.0 ) {
      c1_ = 1.e-4;
    }
    if ( c2_ < 0.0 ) {
      c2_ = 0.9;
    }
    if ( c3_ < 0.0 ) {
      c3_ = 0.9;
    }
    if ( c2_ <= c1_ ) {
      c1_ = 1.e-4;
      c2_ = 0.9;
    }
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      c2_ = 0.4;
      c3_ = std::min(1.0-c2_,c3_);
    }
  }

  virtual void initialize( const Vector<Real> &x, const Vector<Real> &g,
                           Objective<Real> &obj, BoundConstraint<Real> &con ) {
    grad_ = Teuchos::rcp(&g, false);
    xtst_ = x.clone();
    d_    = x.clone();
    g_    = g.clone();
  }

  void setData(Real &eps) {
    eps_ = eps;
  }

  virtual bool status( const ELineSearch type, int &ls_neval, int &ls_ngrad, const Real alpha, 
                       const Real fold, const Real sgold, const Real fnew, 
                       const Vector<Real> &x, const Vector<Real> &s, 
                       Objective<Real> &obj, BoundConstraint<Real> &con ) { 
    Real tol = std::sqrt(ROL_EPSILON);

    // Check Armijo Condition
    bool armijo = false;
    if ( con.isActivated() ) {
      Real gs = 0.0;
      if ( edesc_ == DESCENT_STEEPEST ) {
        updateIterate(*d_,x,s,alpha,con);
        d_->scale(-1.0);
        d_->plus(x);
        gs = -s.dot(*d_);
      }
      else {
        d_->set(s);
        d_->scale(-1.0);
        con.pruneActive(*d_,*(grad_),x,eps_);
        gs = alpha*(grad_)->dot(*d_);
        d_->zero();
        updateIterate(*d_,x,s,alpha,con);
        d_->scale(-1.0);
        d_->plus(x);
        con.pruneInactive(*d_,*(grad_),x,eps_);
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
    bool itcond = false;
    if ( ls_neval >= maxit_ ) { 
      itcond = true;
    }

    // Check Curvature Condition
    bool curvcond = false;
    if ( armijo && ((type != LINESEARCH_BACKTRACKING && type != LINESEARCH_CUBICINTERP) ||
                    (edesc_ == DESCENT_NONLINEARCG)) ) {
      if (econd_ == CURVATURECONDITION_GOLDSTEIN) {
        if (fnew >= fold + (1.0-c1_)*alpha*sgold) {
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
        Real sgnew = 0.0;
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
                     && (c2_*sgold <= sgnew && sgnew <= (2.0*c1_ - 1.0)*sgold)) ) {
          curvcond = true;
        }
      }
    }

    if (type == LINESEARCH_BACKTRACKING || type == LINESEARCH_CUBICINTERP) {
      if (edesc_ == DESCENT_NONLINEARCG) {
        return ((armijo && curvcond) || itcond);
      }
      else {
        return (armijo || itcond);
      }
    }
    else {
      return ((armijo && curvcond) || itcond);
    }
  }

  virtual void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                    const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                    Objective<Real> &obj, BoundConstraint<Real> &con ) = 0;

  virtual Real getInitialAlpha(int &ls_neval, int &ls_ngrad, const Real fval, const Real gs, 
                               const Vector<Real> &x, const Vector<Real> &s, 
                               Objective<Real> &obj, BoundConstraint<Real> &con) {
    Real val = 1.0;
    if (useralpha_) {
      val = alpha0_;
    }
    else {
      if (edesc_ == DESCENT_STEEPEST || edesc_ == DESCENT_NONLINEARCG) {
        Real tol = std::sqrt(ROL_EPSILON);
        Real alpha = 1.0;
        // Evaluate objective at x + s
        updateIterate(*d_,x,s,alpha,con);
        obj.update(*d_);
        Real fnew = obj.value(*d_,tol);
        ls_neval++;
        // Minimize quadratic interpolate to compute new alpha
        alpha = -gs/(2.0*(fnew-fval-gs));
        // Evaluate objective at x + alpha s 
        updateIterate(*d_,x,s,alpha,con);
        obj.update(*d_);
        fnew = obj.value(*d_,tol);
        ls_neval++;
        // Ensure that sufficient decrease and curvature conditions are satisfied
        bool stat = status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,alpha,fval,gs,fnew,x,s,obj,con);
        if ( !stat ) {
          alpha = 1.0;
        }
        val = alpha;
      }
      else {
        val = 1.0;
      }
    }
    return val;
  }

  void updateIterate(Vector<Real> &xnew, const Vector<Real> &x, const Vector<Real> &s, Real alpha,
                     BoundConstraint<Real> &con ) {
    xnew.set(x);
    xnew.axpy(alpha,s);
    if ( con.isActivated() ) {
      con.project(xnew);
    }
  }
};

}

#include "ROL_IterationScaling.hpp"
#include "ROL_PathBasedTargetLevel.hpp"
#include "ROL_BackTracking.hpp"
#include "ROL_CubicInterp.hpp"
#include "ROL_Bisection.hpp"
#include "ROL_GoldenSection.hpp"
#include "ROL_Brents.hpp"

#endif
