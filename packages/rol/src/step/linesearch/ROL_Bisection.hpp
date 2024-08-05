// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BISECTION_H
#define ROL_BISECTION_H

/** \class ROL::Bisection
    \brief Implements a bisection line search.
*/

#include "ROL_LineSearch.hpp"
#include "ROL_BackTracking.hpp"

namespace ROL { 

template<class Real>
class Bisection : public LineSearch<Real> {
private:
  Real tol_;
  ROL::Ptr<Vector<Real> > xnew_; 
  ROL::Ptr<LineSearch<Real> > btls_;

public:

  virtual ~Bisection() {}

  // Constructor
  Bisection( ROL::ParameterList &parlist ) : LineSearch<Real>(parlist) {
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
    Real tol = std::sqrt(ROL_EPSILON<Real>()), half(0.5);
    ls_neval = 0;
    ls_ngrad = 0;
    // Get initial line search parameter
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con);

    // Compute value phi(0)
    Real tl(0);
    Real val_tl = fval;

    // Compute value phi(alpha)
    Real tr = alpha;
    LineSearch<Real>::updateIterate(*xnew_,x,s,tr,con);
    obj.update(*xnew_);
    Real val_tr = obj.value(*xnew_,tol); 
    ls_neval++;

    // Check if phi(alpha) < phi(0)
    if ( val_tr < val_tl ) {
      if ( LineSearch<Real>::status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj,con) ) {
        alpha = tr;
        fval  = val_tr;
        return;
      }
    }

    // Get min( phi(0), phi(alpha) )
    Real t(0);
    Real val_t(0);
    if ( val_tl < val_tr ) {
      t     = tl;
      val_t = val_tl;
    }
    else {
      t     = tr;
      val_t = val_tr;
    }

    // Compute value phi(midpoint)
    Real tc = half*(tl+tr);
    LineSearch<Real>::updateIterate(*xnew_,x,s,tc,con);
    Real val_tc = obj.value(*xnew_,tol);
    ls_neval++;

    // Get min( phi(0), phi(alpha), phi(midpoint) )
    if ( val_tc < val_t ) {
      t     = tc;
      val_t = val_tc;
    }

    Real t1(0), val_t1(0);
    Real t2(0), val_t2(0);

    while (    !LineSearch<Real>::status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con)  
            && std::abs(tr - tl) > tol_ ) {
      t1 = half*(tl+tc);
      LineSearch<Real>::updateIterate(*xnew_,x,s,t1,con);
      obj.update(*xnew_);
      val_t1 = obj.value(*xnew_,tol);
      ls_neval++;

      t2 = half*(tr+tc);
      LineSearch<Real>::updateIterate(*xnew_,x,s,t2,con);
      obj.update(*xnew_);
      val_t2 = obj.value(*xnew_,tol);
      ls_neval++;

      if (    ( (val_tl <= val_tr) && (val_tl <= val_t1) && (val_tl <= val_t2) && (val_tl <= val_tc) ) 
           || ( (val_t1 <= val_tr) && (val_t1 <= val_tl) && (val_t1 <= val_t2) && (val_t1 <= val_tc) ) ) {
        if ( val_tl < val_t1 ) {
          t     = tl;
          val_t = val_tl;
        }
        else {
          t     = t1;
          val_t = val_t1;
        }

        tr     = tc;
        val_tr = val_tc;
        tc     = t1;
        val_tc = val_t1;
      }
      else if ( ( (val_tc <= val_tr) && (val_tc <= val_tl) && (val_tc <= val_t1) && (val_tc <= val_t2) ) ) { 
        t     = tc;
        val_t = val_tc;

        tl     = t1;
        val_tl = val_t1;
        tr     = t2;
        val_tr = val_t2;
      }
      else if (    ( (val_t2 <= val_tr) && (val_t2 <= val_tl) && (val_t2 <= val_t1) && (val_t2 <= val_tc) ) 
                || ( (val_tr <= val_tl) && (val_tr <= val_t1) && (val_tr <= val_t2) && (val_tr <= val_tc) ) ) {
        if ( val_tr < val_t2 ) {
          t     = tr;
          val_t = val_tr;
        }
        else {
          t     = t2;
          val_t = val_t2;
        }

        tl     = tc;
        val_tl = val_tc;
        tc     = t2;
        val_tc = val_t2;
      }
    }

    fval  = val_t;
    alpha = t;
   
    if ( alpha < ROL_EPSILON<Real>() ) {
      btls_->run(alpha,fval,ls_neval,ls_ngrad,gs,s,x,obj,con);
    }
  }
};

}

#endif
