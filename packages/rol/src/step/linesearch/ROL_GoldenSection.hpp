// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GOLDENSECTION_H
#define ROL_GOLDENSECTION_H

/** \class ROL::GoldenSection
    \brief Implements a golden section line search.
*/

#include "ROL_LineSearch.hpp"
#include "ROL_BackTracking.hpp"

namespace ROL { 

template<class Real>
class GoldenSection : public LineSearch<Real> {
private:
  Real tol_;
  ROL::Ptr<Vector<Real> > xnew_; 
  ROL::Ptr<LineSearch<Real> > btls_;

public:

  virtual ~GoldenSection() {}

  // Constructor
  GoldenSection( ROL::ParameterList &parlist ) : LineSearch<Real>(parlist) {
    Real oem8(1.e-8);
    tol_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Bracketing Tolerance",oem8);
    btls_ = ROL::makePtr<BackTracking<Real>>(parlist);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                   Objective<Real> &obj, BoundConstraint<Real> &con ) {
    LineSearch<Real>::initialize(x,s,g,obj,con);
    xnew_ = x.clone();
    btls_->initialize(x,s,g,obj,con);
  }

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>()), zero(0), one(1), two(2), five(5);
    ls_neval = 0; 
    ls_ngrad = 0;
    // Get initial line search parameter
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con);
    
    // Reciprocal of golden ratio
    Real c = two/(one+sqrt(five));

    // Compute value phi(0)
    Real tl  = zero;
    Real val_tl = fval;

    // Compute value phi(alpha)
    Real tr  = alpha;
    LineSearch<Real>::updateIterate(*xnew_,x,s,tr,con);
    obj.update(*xnew_);
    Real val_tr = obj.value(*xnew_,tol);
    ls_neval++;

    // Check if phi(alpha) < phi(0)
    if ( val_tr < val_tl ) {
      if ( LineSearch<Real>::status(LINESEARCH_GOLDENSECTION,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj,con) ) {
        alpha = tr;
        fval  = val_tr;
        return;
      }
    }

    // Compute min( phi(0), phi(alpha) )
    Real t     = zero;
    Real val_t = zero;
    if ( val_tl < val_tr ) {
      t     = tl;
      val_t = val_tl;
    }
    else {
      t     = tr;
      val_t = val_tr;
    }

    // Compute value phi(t1)
    Real tc1 = c*tl + (one-c)*tr;
    LineSearch<Real>::updateIterate(*xnew_,x,s,tc1,con);
    obj.update(*xnew_);
    Real val_tc1 = obj.value(*xnew_,tol);
    ls_neval++;

    // Compute value phi(t2)
    Real tc2 = (one-c)*tl + c*tr;
    LineSearch<Real>::updateIterate(*xnew_,x,s,tc2,con);
    obj.update(*xnew_);
    Real val_tc2 = obj.value(*xnew_,tol);
    ls_neval++;

    // Compute min( phi(0), phi(t1), phi(t2), phi(alpha) )
    if ( val_tl <= val_tc1 && val_tl <= val_tc2 && val_tl <= val_tr ) {
      val_t = val_tl;
      t     = tl;
    }
    else if ( val_tc1 <= val_tl && val_tc1 <= val_tc2 && val_tc1 <= val_tr ) {
      val_t = val_tc1;
      t     = tc1;
    }
    else if ( val_tc2 <= val_tl && val_tc2 <= val_tc1 && val_tc2 <= val_tr ) {
      val_t = val_tc2;
      t     = tc2;
    }
    else {
      val_t = val_tr;
      t     = tr;
    }

    while (    !LineSearch<Real>::status(LINESEARCH_GOLDENSECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con) 
            && (std::abs(tl-tr) >= tol_) ) {
      if ( val_tc1 > val_tc2 ) {
        tl      = tc1;
        val_tl  = val_tc1;
        tc1     = tc2;
        val_tc1 = val_tc2;
 
        tc2     = (one-c)*tl + c*tr;     
        LineSearch<Real>::updateIterate(*xnew_,x,s,tc2,con);
        obj.update(*xnew_);
        val_tc2 = obj.value(*xnew_,tol);
        ls_neval++;
      }
      else {
        tr      = tc2;
        val_tr  = val_tc2;
        tc2     = tc1;
        val_tc2 = val_tc1;

        tc1     = c*tl + (one-c)*tr;
        LineSearch<Real>::updateIterate(*xnew_,x,s,tc1,con);
        obj.update(*xnew_);
        val_tc1 = obj.value(*xnew_,tol);
        ls_neval++;
      }

      if ( val_tl <= val_tc1 && val_tl <= val_tc2 && val_tl <= val_tr ) {
        val_t = val_tl;
        t     = tl;
      }
      else if ( val_tc1 <= val_tl && val_tc1 <= val_tc2 && val_tc1 <= val_tr ) {
        val_t = val_tc1;
        t     = tc1;
      }
      else if ( val_tc2 <= val_tl && val_tc2 <= val_tc1 && val_tc2 <= val_tr ) {
        val_t = val_tc2;
        t     = tc2;
      }
      else {
        val_t = val_tr;
        t     = tr;
      }
    }
    alpha = t;
    fval  = val_t;  

    if ( alpha < ROL_EPSILON<Real>() ) {
      btls_->run(alpha,fval,ls_neval,ls_ngrad,gs,s,x,obj,con);
    }
  }
};

}

#endif
