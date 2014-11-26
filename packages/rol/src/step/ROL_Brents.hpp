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

#ifndef ROL_BRENTS_H
#define ROL_BRENTS_H

/** \class ROL::Brents
    \brief Implements a Brent's method line search.
*/

#include "ROL_LineSearch.hpp"
#include "ROL_BackTracking.hpp"

namespace ROL { 

template<class Real>
class Brents : public LineSearch<Real> {
private:
  Real tol_;
  Teuchos::RCP<Vector<Real> > xnew_; 
  Teuchos::RCP<LineSearch<Real> > btls_;

public:

  virtual ~Brents() {}

  // Constructor
  Brents( Teuchos::ParameterList &parlist ) : LineSearch<Real>(parlist) {
    tol_ = parlist.get("Bracketing Tolerance",1.e-8);
    btls_ = Teuchos::rcp(new BackTracking<Real>(parlist));
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &g,
                   Objective<Real> &obj, BoundConstraint<Real> &con ) {
    LineSearch<Real>::initialize(x,g,obj,con);
    xnew_ = x.clone();
    btls_->initialize(x,g,obj,con);
  }

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);
    ls_neval = 0;
    ls_ngrad = 0;
    // Get initial line search parameter
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con);

    // Compute value phi(0)
    Real tl = 0.0;         // Left interval point
    Real val_tl = fval;

    // Initialize value phi(t)
    Real tc = 0.0;        // Center interval point
    Real val_tc = 0.0;

    // Compute value phi(alpha)
    Real tr = alpha;      // Right interval point
    LineSearch<Real>::updateIterate(*xnew_,x,s,tr,con);
    obj.update(*xnew_);
    Real val_tr = obj.value(*xnew_,tol);
    ls_neval++;

    // Check if phi(alpha) < phi(0)
    if ( val_tr < val_tl ) {
      if ( LineSearch<Real>::status(LINESEARCH_BRENTS,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj,con) ) {
        alpha = tr;
        fval  = val_tr;
        return;
      }
    }

    // Compute min( phi(0), phi(alpha) )
    Real t     = 0.0;
    Real val_t = 0.0;
    if ( val_tl < val_tr ) {
      t     = tl;
      val_t = val_tl;
    }
    else {
      t     = tr;
      val_t = val_tr;
    }

    // Determine bracketing triple
    const Real gr                = (1.0+sqrt(5.0))/2.0;
    const Real inv_gr2           = 1.0/(gr*gr);
    const Real goldinv           = 1.0/(1.0+gr);
    const Real tiny              = sqrt(ROL_EPSILON);
    const Real max_extrap_factor = 100.0;
    Real tmp    = 0.0;
    Real q      = 0.0;
    Real r      = 0.0; 
    Real tm     = 0.0;
    Real tlim   = 0.0; 
    Real val_tm = 0.0;

    int itbt = 0;
    while ( val_tr > val_tl && itbt < 8 ) {
      tc     = tr;
      val_tc = val_tr;

      tr     = goldinv * (tc + gr*tl);
      LineSearch<Real>::updateIterate(*xnew_,x,s,tr,con);
      obj.update(*xnew_);
      val_tr = obj.value(*xnew_,tol);
      ls_neval++;

      itbt++;
    }
    if ( val_tr > val_tl ) {
      tmp    = tl;
      tl     = tr;
      tr     = tmp;
      tmp    = val_tr;
      val_tr = val_tl;
      val_tl = tmp;
      tc     = 0.0;
    }

    if ( std::abs(tc) < ROL_EPSILON ) {
      tc = tl + (gr-1.0)*(tr-tl);
      LineSearch<Real>::updateIterate(*xnew_,x,s,tc,con);
      obj.update(*xnew_);
      val_tc = obj.value(*xnew_,tol);
      ls_neval++;
    }

    if ( val_tl <= val_tr && val_tl <= val_tc ) {
      t     = tl;
      val_t = val_tl;
    }
    else if ( val_tc <= val_tr && val_tc <= val_tl ) {
      t     = tc;
      val_t = val_tc;
    }
    else {
      t     = tr;
      val_t = val_tr;
    }

    if ( LineSearch<Real>::status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con) ) {
      alpha = t;
      fval  = val_t;
      return;
    }
    
    int itb = 0;
    while ( val_tr >= val_tc && itb < 8 ) {
      q = ( val_tr-val_tl ) * (tr - tc);
      r = ( val_tr-val_tc ) * (tr - tl);
      tmp = fabs(q-r);
      tmp = (tmp > tiny ? tmp : -tmp);
      tm  = tr - (q*(tr-tc) - r*(tr-tl))/(2.0*tmp);

      tlim = tl + max_extrap_factor * (tc-tr);

      if ( (tr-tm)*(tm-tc) > 0.0 ) {
        LineSearch<Real>::updateIterate(*xnew_,x,s,tm,con);
        obj.update(*xnew_);
        val_tm = obj.value(*xnew_,tol);
        ls_neval++;
        if ( val_tm < val_tc ) {
          tl     = tr;
          val_tl = val_tr;
          tr     = tm;
          val_tr = val_tm;
        }
        else if ( val_tm > val_tr) {
          tc     = tm;
          val_tc = val_tm;
        }
        tm = tc + gr*(tc-tr);
        LineSearch<Real>::updateIterate(*xnew_,x,s,tm,con);
        obj.update(*xnew_);
        val_tm = obj.value(*xnew_,tol);
        ls_neval++;
      }
      else if ( (tc - tm)*(tm -tlim) > 0.0 ) {
        LineSearch<Real>::updateIterate(*xnew_,x,s,tm,con);
        obj.update(*xnew_);
        val_tm = obj.value(*xnew_,tol);
        ls_neval++;
        if ( val_tm < val_tc ) {
          tr     = tc;
          val_tr = val_tc;

          tc     = tm;
          val_tc = val_tm;

          tm     = tc + gr*(tc-tr);
          LineSearch<Real>::updateIterate(*xnew_,x,s,tm,con);
          obj.update(*xnew_);
          val_tm = obj.value(*xnew_,tol);
          ls_neval++;
        }
      }
      else if ( (tm-tlim)*(tlim-tc) >= 0.0 ) {
        tm = tlim;
        LineSearch<Real>::updateIterate(*xnew_,x,s,tm,con);
        obj.update(*xnew_);
        val_tm = obj.value(*xnew_,tol);
        ls_neval++;
      }
      else {
        tm = tc + gr*(tc-tr);
        LineSearch<Real>::updateIterate(*xnew_,x,s,tm,con);
        obj.update(*xnew_);
        val_tm = obj.value(*xnew_,tol);
        ls_neval++;
      }
      tl     = tr;
      val_tl = val_tr;
      tr     = tc;
      val_tr = val_tc;
      tc     = tm;
      val_tc = val_tm;
      itb++;
    }
     
    if ( val_tl <= val_tr && val_tl <= val_tc ) {
      t     = tl;
      val_t = val_tl;
    }
    else if ( val_tc <= val_tr && val_tc <= val_tl ) {
      t     = tc;
      val_t = val_tc;
    }
    else {
      t     = tr;
      val_t = val_tr;
    }

    if ( LineSearch<Real>::status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con) ) {
      alpha = t;
      fval  = val_t;
      return;
    }
 
    // Run Brent's using the triple (tl,tr,tc)
    Real a     = 0.0;
    Real b     = 0.0;
    Real d     = 0.0;
    Real e     = 0.0;
    Real etemp = 0.0; 
    Real fu    = 0.0; 
    Real fv    = 0.0;
    Real fw    = 0.0;
    Real ft    = 0.0;
    Real p     = 0.0;
    Real u     = 0.0;
    Real v     = 0.0;
    Real w     = 0.0;
    int it     = 0;
 
    fw = (val_tl<val_tc ? val_tl : val_tc);
    if ( fw == val_tl ) {
      w  = tl;
      v  = tc;
      fv = val_tc;
    }
    else {
      w  = tc;
      v  = tl;
      fv = val_tl; 
    }
    t  = tr;
    ft = val_tr;
    a  = (tr < tc ? tr : tc);
    b  = (tr > tc ? tr : tc);

    while (    !LineSearch<Real>::status(LINESEARCH_BRENTS,ls_neval,ls_ngrad,t,fval,gs,ft,x,s,obj,con)
            && std::abs(t - tm) > tol_*(b-a) ) {
      if ( it < 2 ) {
        e = 2.0*(b-a);
      }
      tm = (a+b)/2.0;

      Real tol1 = tol_*std::abs(t) + tiny;
      Real tol2 = 2.0*tol1;

      if ( std::abs(e) > tol1 || it < 2 ) {
        r     = (t-w)*(ft-fv);
        q     = (t-v)*(ft-fw);
        p     = (t-v)*q-(t-w)*r;
        q     = 2.0*(q-r);
        if ( q > 0.0 ) {
          p = -p;
        }
        q     = std::abs(q);
        etemp = e;
        e     = d;
        if ( std::abs(p) > std::abs(0.5*q*etemp) || p <= q*(a-t) || p >= q*(b-t) ) {
          d = inv_gr2*(e=(t>=tm ? a-t : b-t));  
        }
        else {
          d = p/q;
          u = t+d;
          if ( u-a < tol2 || b-u < tol2 ) {
            d = ( tm-t > 0.0 ? std::abs(tol1) : -std::abs(tol1) );
          }
        }
      }
      else  {
        d = inv_gr2*(e = (t>=tm ? a-t : b-t) );
      }
      u = (std::abs(d)>=tol1 ? t+d : t+(d>=0.0 ? std::abs(tol1) : -std::abs(tol1)));
      LineSearch<Real>::updateIterate(*xnew_,x,s,u,con);
      obj.update(*xnew_);
      fu = obj.value(*xnew_,tol);
      ls_neval++;

      if ( fu <= ft ) {
        if ( u >= t ) {
          a = t;
        }
        else {
          b = t;
        }
        v  = w;
        fv = fw;
        w  = t;
        fw = ft;
        t  = u;
        ft = fu;
      }
      else {
        if ( u < t ) {
          a = u;
        }
        else {
          b = u;
        }
        if ( fu <= fw || w == t ) {
          v  = w;
          fv = fw;
          w  = u;
          fw = fu;
        }
        else if ( fu <= fv || v == t || v == w ) {
          v  = u;
          fv = fu;
        }
      }
      it++;
    }
    alpha = t;
    fval  = ft;

    if ( alpha < ROL_EPSILON ) {
      btls_->run(alpha,fval,ls_neval,ls_ngrad,gs,s,x,obj,con);
    }
  }
};

}

#endif
