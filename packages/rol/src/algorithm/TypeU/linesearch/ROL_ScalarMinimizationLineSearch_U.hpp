// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARMINIMIZATIONLINESEARCH_U_H
#define ROL_SCALARMINIMIZATIONLINESEARCH_U_H

/** \class ROL::ScalarMinimizationLineSearch_U
    \brief Implements line search methods that attempt to minimize the
           scalar function \f$\phi(t) := f(x+ts)\f$.
*/

#include "ROL_LineSearch_U.hpp"
#include "ROL_BrentsScalarMinimization.hpp"
#include "ROL_BisectionScalarMinimization.hpp"
#include "ROL_GoldenSectionScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_Bracketing.hpp"

namespace ROL {

template<typename Real>
class ScalarMinimizationLineSearch_U : public LineSearch_U<Real> {
private:
  Ptr<Vector<Real>>             xnew_; 
  Ptr<Vector<Real>>             g_;
  Ptr<ScalarMinimization<Real>> sm_;
  Ptr<Bracketing<Real>>         br_;
  Ptr<ScalarFunction<Real>>     sf_;

  ECurvatureConditionU econd_;
  Real c1_;
  Real c2_;
  Real c3_;
  int max_nfval_;

  bool FDdirDeriv_;

  class Phi : public ScalarFunction<Real> {
  private:
    const Ptr<Vector<Real>> xnew_;
    const Ptr<const Vector<Real>> x_, s_;
    const Ptr<Objective<Real>> obj_;
    Real ftol_, alpha_, val_;
    bool FDdirDeriv_;
  public:
    Phi(const Ptr<Vector<Real>> &xnew,
        const Ptr<const Vector<Real>> &x,
        const Ptr<const Vector<Real>> &s,
        const Ptr<Objective<Real>> &obj,
        const bool FDdirDeriv = false)
     : xnew_(xnew), x_(x), s_(s), obj_(obj),
       ftol_(std::sqrt(ROL_EPSILON<Real>())),
       alpha_(ROL_INF<Real>()), val_(ROL_INF<Real>()),
       FDdirDeriv_(FDdirDeriv) {}
    Real value(const Real alpha) {
      if (alpha_ != alpha) {
        alpha_ = alpha;
        xnew_->set(*x_); xnew_->axpy(alpha,*s_);
        obj_->update(*xnew_,UpdateType::Trial);
        val_ = obj_->value(*xnew_,ftol_);
      }
      return val_;
    }
    Real deriv(const Real alpha) {
      Real tol = std::sqrt(ROL_EPSILON<Real>());
      Real val(0);
      xnew_->set(*x_); xnew_->axpy(alpha,*s_);
      if (FDdirDeriv_) {
        Real snorm = s_->norm();
        if (snorm > static_cast<Real>(0)) {
          Real xnorm   = xnew_->norm();
          Real cbrteps = std::cbrt(ROL_EPSILON<Real>());
          Real h       = cbrteps*std::max(xnorm/snorm,static_cast<Real>(1));
          Real fnew    = value(alpha);
          xnew_->axpy(h,*s_);
          obj_->update(*xnew_,UpdateType::Trial);
          Real ftrial = obj_->value(*xnew_,tol);
          val = (ftrial - fnew) / h;
        }
      }
      else {
        val = obj_->dirDeriv(*xnew_,*s_,tol);
      }
      return val;
    }
  };

  class StatusTest : public ScalarMinimizationStatusTest<Real> {
  private:
    Ptr<ScalarFunction<Real>> phi_;

    const Real f0_;
    const Real g0_;

    const Real c1_;
    const Real c2_;
    const Real c3_;
    const int max_nfval_;
    const ECurvatureConditionU econd_;

  public:
    StatusTest(const Real f0, const Real g0,
               const Real c1, const Real c2, const Real c3,
               const int max_nfval, ECurvatureConditionU econd,
               const Ptr<ScalarFunction<Real>> &phi)
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
        if (econd_ == CURVATURECONDITION_U_GOLDSTEIN) {
          curvcond = (fx >= f0_ + (one-c1_)*x*g0_);
        }
        else if (econd_ == CURVATURECONDITION_U_NULL) {
          curvcond = true;
        }
        else {
          if (!deriv) {
            gx = phi_->deriv(x); ngval++;
          }
          if (econd_ == CURVATURECONDITION_U_WOLFE) {
            curvcond = (gx >= c2_*g0_);
          }
          else if (econd_ == CURVATURECONDITION_U_STRONGWOLFE) {
            curvcond = (std::abs(gx) <= c2_*std::abs(g0_));
          }
          else if (econd_ == CURVATURECONDITION_U_GENERALIZEDWOLFE) {
            curvcond = (c2_*g0_ <= gx && gx <= -c3_*g0_);
          }
          else if (econd_ == CURVATURECONDITION_U_APPROXIMATEWOLFE) {
            curvcond = (c2_*g0_ <= gx && gx <= (two*c1_ - one)*g0_);
          }
        }
      }
      //return (armijo && curvcond) || itcond;
      return (armijo && curvcond);
    }
  };

  using LineSearch_U<Real>::getInitialAlpha;
  using LineSearch_U<Real>::setNextInitialAlpha;

public:
  // Constructor
  ScalarMinimizationLineSearch_U( ParameterList &parlist, 
    const Ptr<ScalarMinimization<Real>> &sm = nullPtr,
    const Ptr<Bracketing<Real>>         &br = nullPtr,
    const Ptr<ScalarFunction<Real>>     &sf = nullPtr )
    : LineSearch_U<Real>(parlist) {
    const Real zero(0), p4(0.4), p6(0.6), p9(0.9), oem4(1.e-4), oem10(1.e-10), one(1);
    ParameterList &list0 = parlist.sublist("Step").sublist("Line Search");
    FDdirDeriv_ = list0.get("Finite Difference Directional Derivative",false);
    ParameterList &list  = list0.sublist("Line-Search Method");
    // Get Bracketing Method
    if( br == nullPtr ) {
      br_ = makePtr<Bracketing<Real>>();
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

    if( sm == nullPtr ) { // No user-provided ScalarMinimization object
      if ( type == "Brent's" ) {
        sm_ = makePtr<BrentsScalarMinimization<Real>>(plist);
      }
      else if ( type == "Bisection" ) {
        sm_ = makePtr<BisectionScalarMinimization<Real>>(plist);
      }
      else if ( type == "Golden Section" ) {
        sm_ = makePtr<GoldenSectionScalarMinimization<Real>>(plist);
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
    econd_ = StringToECurvatureConditionU(list0.sublist("Curvature Condition").get("Type","Strong Wolfe Conditions"));
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
    EDescentU edesc = StringToEDescentU(list0.sublist("Descent Method").get("Type","Quasi-Newton Method"));
    if ( edesc == DESCENT_U_NONLINEARCG ) {
      c2_ = p4;
      c3_ = std::min(one-c2_,c3_);
    }
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) override {
    LineSearch_U<Real>::initialize(x,g);
    xnew_ = x.clone();
    g_    = g.clone();
  }

  // Find the minimum of phi(alpha) = f(x + alpha*s) using Brent's method
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj ) override {
    ls_neval = 0; ls_ngrad = 0;

    // Get initial line search parameter
    alpha = getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj);

    // Build ScalarFunction and ScalarMinimizationStatusTest
    Ptr<const Vector<Real>> x_ptr = ROL::makePtrFromRef(x);
    Ptr<const Vector<Real>> s_ptr = ROL::makePtrFromRef(s);
    Ptr<Objective<Real>> obj_ptr = ROL::makePtrFromRef(obj);

    Ptr<ScalarFunction<Real>> phi;
    if( sf_ == ROL::nullPtr ) {
      phi = ROL::makePtr<Phi>(xnew_,x_ptr,s_ptr,obj_ptr,FDdirDeriv_);
    }
    else {
      phi = sf_;
    }

    Ptr<ScalarMinimizationStatusTest<Real>> test
      = makePtr<StatusTest>(fval,gs,c1_,c2_,c3_,max_nfval_,econd_,phi);

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

    setNextInitialAlpha(alpha); 
  }
}; // class ROL::ScalarMinimization_U

} // namespace ROL

#endif
