// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ScalarMinimizationLineSearch_H
#define ROL_ScalarMinimizationLineSearch_H

/** \class ROL::ScalarMinimizationLineSearch
    \brief Implements line search methods that attempt to minimize the
           scalar function \f$\phi(t) := f(x+ts)\f$.
*/

#include "ROL_LineSearch.hpp"
#include "ROL_BrentsScalarMinimization.hpp"
#include "ROL_BisectionScalarMinimization.hpp"
#include "ROL_GoldenSectionScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_Bracketing.hpp"

namespace ROL { 

template<class Real>
class ScalarMinimizationLineSearch : public LineSearch<Real> {
private:
  ROL::Ptr<Vector<Real> >             xnew_; 
  ROL::Ptr<Vector<Real> >             g_;
  ROL::Ptr<ScalarMinimization<Real> > sm_;
  ROL::Ptr<Bracketing<Real> >         br_;
  ROL::Ptr<ScalarFunction<Real> >     sf_;

  ECurvatureCondition econd_;
  Real c1_;
  Real c2_;
  Real c3_;
  int max_nfval_;

  class Phi : public ScalarFunction<Real> {
  private:
    const ROL::Ptr<Vector<Real> > xnew_;
    const ROL::Ptr<Vector<Real> > g_;
    const ROL::Ptr<const Vector<Real> > x_;
    const ROL::Ptr<const Vector<Real> > s_;
    const ROL::Ptr<Objective<Real> > obj_;
    const ROL::Ptr<BoundConstraint<Real> > con_;
    Real ftol_;
    void updateIterate(Real alpha) {
      xnew_->set(*x_);
      xnew_->axpy(alpha,*s_);
      if ( con_->isActivated() ) {
        con_->project(*xnew_);
      }
    }
  public:
    Phi(const ROL::Ptr<Vector<Real> > &xnew,
        const ROL::Ptr<Vector<Real> > &g,
        const ROL::Ptr<const Vector<Real> > &x,
        const ROL::Ptr<const Vector<Real> > &s,
        const ROL::Ptr<Objective<Real> > &obj,
        const ROL::Ptr<BoundConstraint<Real> > &con)
     : xnew_(xnew), g_(g), x_(x), s_(s), obj_(obj), con_(con),
       ftol_(std::sqrt(ROL_EPSILON<Real>())) {}
    Real value(const Real alpha) {
      updateIterate(alpha);
      obj_->update(*xnew_);
      return obj_->value(*xnew_,ftol_);
    }
    Real deriv(const Real alpha) {
      updateIterate(alpha);
      obj_->update(*xnew_);
      obj_->gradient(*g_,*xnew_,ftol_);
      return s_->dot(g_->dual()); 
    }
  };

  class LineSearchStatusTest : public ScalarMinimizationStatusTest<Real> {
  private:
    ROL::Ptr<ScalarFunction<Real> > phi_;

    const Real f0_;
    const Real g0_;

    const Real c1_;
    const Real c2_;
    const Real c3_;
    const int max_nfval_;
    const ECurvatureCondition econd_;


  public:
    LineSearchStatusTest(const Real f0, const Real g0,
                         const Real c1, const Real c2, const Real c3,
                         const int max_nfval, ECurvatureCondition econd,
                         const ROL::Ptr<ScalarFunction<Real> > &phi)
      : phi_(phi), f0_(f0), g0_(g0), c1_(c1), c2_(c2), c3_(c3),
        max_nfval_(max_nfval), econd_(econd) {}

    bool check(Real &x, Real &fx, Real &gx,
               int &nfval, int &ngval, const bool deriv = false) {
      Real one(1), two(2);
      bool armijo = (fx <= f0_ + c1_*x*g0_);
//      bool itcond = (nfval >= max_nfval_);
      bool curvcond = false;
//      if (armijo && !itcond) {
      if (armijo) {
        if (econd_ == CURVATURECONDITION_GOLDSTEIN) {
          curvcond = (fx >= f0_ + (one-c1_)*x*g0_);
        }
        else if (econd_ == CURVATURECONDITION_NULL) {
          curvcond = true;
        }
        else {
          if (!deriv) {
            gx = phi_->deriv(x); ngval++;
          }
          if (econd_ == CURVATURECONDITION_WOLFE) {
            curvcond = (gx >= c2_*g0_);
          }
          else if (econd_ == CURVATURECONDITION_STRONGWOLFE) {
            curvcond = (std::abs(gx) <= c2_*std::abs(g0_));
          }
          else if (econd_ == CURVATURECONDITION_GENERALIZEDWOLFE) {
            curvcond = (c2_*g0_ <= gx && gx <= -c3_*g0_);
          }
          else if (econd_ == CURVATURECONDITION_APPROXIMATEWOLFE) {
            curvcond = (c2_*g0_ <= gx && gx <= (two*c1_ - one)*g0_);
          }
        }
      }
      //return (armijo && curvcond) || itcond;
      return (armijo && curvcond);
    }
  };

public:
  // Constructor
  ScalarMinimizationLineSearch( ROL::ParameterList &parlist, 
    const ROL::Ptr<ScalarMinimization<Real> > &sm = ROL::nullPtr,
    const ROL::Ptr<Bracketing<Real> > &br = ROL::nullPtr,
    const ROL::Ptr<ScalarFunction<Real> > &sf  = ROL::nullPtr )
    : LineSearch<Real>(parlist) {
    Real zero(0), p4(0.4), p6(0.6), p9(0.9), oem4(1.e-4), oem10(1.e-10), one(1);
    ROL::ParameterList &list0 = parlist.sublist("Step").sublist("Line Search");
    ROL::ParameterList &list  = list0.sublist("Line-Search Method");
    // Get Bracketing Method
    if( br == ROL::nullPtr ) {
      br_ = ROL::makePtr<Bracketing<Real>>();
    }
    else {
      br_ = br;
    }
    // Get ScalarMinimization Method
    std::string type = list.get("Type","Brent's");
    Real tol         = list.sublist(type).get("Tolerance",oem10);
    int niter        = list.sublist(type).get("Iteration Limit",1000);
    ROL::ParameterList plist;
    plist.sublist("Scalar Minimization").set("Type",type);
    plist.sublist("Scalar Minimization").sublist(type).set("Tolerance",tol);
    plist.sublist("Scalar Minimization").sublist(type).set("Iteration Limit",niter);

    if( sm == ROL::nullPtr ) { // No user-provided ScalarMinimization object

      if ( type == "Brent's" ) {
        sm_ = ROL::makePtr<BrentsScalarMinimization<Real>>(plist);
      }
      else if ( type == "Bisection" ) {
        sm_ = ROL::makePtr<BisectionScalarMinimization<Real>>(plist);
      }
      else if ( type == "Golden Section" ) {
        sm_ = ROL::makePtr<GoldenSectionScalarMinimization<Real>>(plist);
      }
      else {
        ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          ">>> (ROL::ScalarMinimizationLineSearch): Undefined ScalarMinimization type!");
      }
    }
    else {
      sm_ = sm;
    }

    sf_ = sf;


    // Status test for line search
    econd_ = StringToECurvatureCondition(list0.sublist("Curvature Condition").get("Type","Strong Wolfe Conditions"));
    max_nfval_ = list0.get("Function Evaluation Limit",20);
    c1_        = list0.get("Sufficient Decrease Tolerance",oem4);
    c2_        = list0.sublist("Curvature Condition").get("General Parameter",p9);
    c3_        = list0.sublist("Curvature Condition").get("Generalized Wolfe Parameter",p6);
    // Check status test inputs
    c1_ = ((c1_ < zero) ? oem4 : c1_);
    c2_ = ((c2_ < zero) ? p9   : c2_);
    c3_ = ((c3_ < zero) ? p9   : c3_);
    if ( c2_ <= c1_ ) {
      c1_ = oem4;
      c2_ = p9;
    }
    EDescent edesc = StringToEDescent(list0.sublist("Descent Method").get("Type","Quasi-Newton Method"));
    if ( edesc == DESCENT_NONLINEARCG ) {
      c2_ = p4;
      c3_ = std::min(one-c2_,c3_);
    }
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                   Objective<Real> &obj, BoundConstraint<Real> &con ) {
    LineSearch<Real>::initialize(x,s,g,obj,con);
    xnew_ = x.clone();
    g_    = g.clone();
  }

  // Find the minimum of phi(alpha) = f(x + alpha*s) using Brent's method
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    ls_neval = 0; ls_ngrad = 0;

    // Get initial line search parameter
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con);

    // Build ScalarFunction and ScalarMinimizationStatusTest
    ROL::Ptr<const Vector<Real> > x_ptr = ROL::makePtrFromRef(x);
    ROL::Ptr<const Vector<Real> > s_ptr = ROL::makePtrFromRef(s);
    ROL::Ptr<Objective<Real> > obj_ptr = ROL::makePtrFromRef(obj);
    ROL::Ptr<BoundConstraint<Real> > bnd_ptr = ROL::makePtrFromRef(con);


    ROL::Ptr<ScalarFunction<Real> > phi;

    if( sf_ == ROL::nullPtr ) {
      phi = ROL::makePtr<Phi>(xnew_,g_,x_ptr,s_ptr,obj_ptr,bnd_ptr);
    }
    else {
      phi = sf_;
    }

    ROL::Ptr<ScalarMinimizationStatusTest<Real> > test
      = ROL::makePtr<LineSearchStatusTest>(fval,gs,c1_,c2_,c3_,max_nfval_,econd_,phi);

    // Run Bracketing
    int nfval = 0, ngrad = 0;
    Real A(0),      fA = fval;
    Real B = alpha, fB = phi->value(B);
    br_->run(alpha,fval,A,fA,B,fB,nfval,ngrad,*phi,*test); 
    B = alpha;
    ls_neval += nfval; ls_ngrad += ngrad;

    // Run ScalarMinimization
    nfval = 0, ngrad = 0;
    sm_->run(fval, alpha, nfval, ngrad, *phi, A, B, *test);
    ls_neval += nfval; ls_ngrad += ngrad;

    LineSearch<Real>::setNextInitialAlpha(alpha); 
  }
};

}

#endif
