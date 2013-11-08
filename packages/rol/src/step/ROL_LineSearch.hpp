//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_LINESEARCH_H
#define ROL_LINESEARCH_H

/** \class ROL::LineSearch
    \brief Provides interfrace for different line-searches.
*/

#include <Teuchos_ScalarTraits.hpp>

namespace ROL { 

enum LineSearchStepType {
  LineSearchStep_Gradient = 0,
  LineSearchStep_NonlinearCG,
  LineSearchStep_Secant,
  LineSearchStep_Newton,
  LineSearchStep_NewtonKrylov,
  LineSearchStep_NewtonKrylovSecantPreconditioning
};

enum LineSearchType {
  LineSearchType_Backtracking = 0,   // Backtracking
  LineSearchType_SimpleBacktracking, // Simple Backtracking
  LineSearchType_Brents,             // Brents
  LineSearchType_Bisection,          // Bisection
  LineSearchType_GoldenSection       // GoldenSection
};

enum LineSearchCondition {
  LineSearchCondition_Wolfe = 0,   // Wolfe
  LineSearchCondition_StrongWolfe, // Strong Wolfe
  LineSearchCondition_Goldstein    // Goldstein
};

template<class Real>
class LineSearch {
private:

  LineSearchType      type_;
  LineSearchCondition cond_;
  LineSearchStepType  step_;

  int maxit_;
  Real c1_;
  Real c2_;
  Real tol_;
  Real rho_;
  Real alpha0_;
  Real gs_old_;
  Real alpha_old_;
  bool firstIt_;

public:

  virtual ~LineSearch() {}

  // Constructor
  LineSearch( LineSearchType type, LineSearchCondition cond, LineSearchStepType step,  
              int maxit = 20, Real c1 = 1.e-4, Real c2 = 0.9, Real tol = 1.e-8, Real rho = 0.5 ) : 
    type_(type), cond_(cond), step_(step), maxit_(maxit), c1_(c1), c2_(c2), tol_(tol), rho_(rho) {
    if ( c2_ <= c1_ ) {
      c1_ = 1.e-4;
      c2_ = 0.9;
      if ( step_ == LineSearchStep_NonlinearCG ) {
        c2_ = 0.1;
      }
    }

    gs_old_    = 0.0;
    alpha_old_ = 0.0;
    firstIt_   = true;
  }

  bool status( const LineSearchType type, int &ls_neval, int &ls_ngrad, const Real alpha, 
               const Real fold, const Real sgold, const Real fnew, 
               const Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj ) { 
    // Check Armijo Condition
    bool armijo = false;
    if ( fnew <= fold + c1_*alpha*sgold ) {
      armijo = true;
    }

    // Check Maximum Iteration
    bool itcond = false;
    if ( ls_neval >= maxit_ ) { 
      itcond = true;
    }

    // Check Curvature Condition
    bool curvcond = false;
    if ( armijo && (type != LineSearchType_Backtracking && type != LineSearchType_SimpleBacktracking) ) {
      if ( (cond_ == LineSearchCondition_Goldstein) && (fnew >= fold + (1.0-c1_)*alpha*sgold) ) {
        curvcond = true;
      }
      else { 
        Teuchos::RCP<Vector<Real> > xnew = x.clone();
        xnew->set(x);
        xnew->axpy(alpha,s);   
        Teuchos::RCP<Vector<Real> > grad = x.clone();
        obj.gradient(*grad,*xnew);
        Real sgnew = grad->dot(s);
        ls_ngrad++;
   
        if (    ((cond_ == LineSearchCondition_Wolfe)       && (sgnew >= c2_*sgold))
             || ((cond_ == LineSearchCondition_StrongWolfe) && (std::abs(sgnew) <= c2_*std::abs(sgold))) ) {
          curvcond = true;
        }
      }
    }

    if (type == LineSearchType_Backtracking || type == LineSearchType_SimpleBacktracking) {
      return armijo;
    }
    else {
      return ((armijo && curvcond) || itcond);
    }
  }

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    // Determine Initial Step Length
    alpha0_ = 1.0;
    if ( firstIt_ == true ) {
      firstIt_ = false;
    }
    else {
      if ( step_ == LineSearchStep_Gradient || step_ == LineSearchStep_NonlinearCG ) {
        alpha0_ = alpha_old_*gs_old_/gs;      
      }
    }

    // Run Linesearch
    ls_neval = 0;
    ls_ngrad = 0;
    if ( type_ == LineSearchType_SimpleBacktracking ) {
      simplebacktracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Backtracking ) {
      backtracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Brents ) {
      brents( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Bisection ) {
      bisection( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_GoldenSection ) {
      goldensection( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }

    // Update Storage
    gs_old_    = gs;
    alpha_old_ = alpha;
  }

  void simplebacktracking( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                           const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    alpha = alpha0_;

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(alpha, s);

    Real fold = fval;
    fval = obj.value(*xnew);
    ls_neval++;

    while ( !status(LineSearchType_SimpleBacktracking,ls_neval,ls_ngrad,alpha,fold,gs,fval,*xnew,s,obj) ) {
      alpha *= rho_;
      xnew->set(x);
      xnew->axpy(alpha, s);
      fval = obj.value(*xnew);
      ls_neval++;
    }
  }

  void backtracking( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                     const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    alpha = alpha0_;

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(alpha, s);

    Real fold = fval;
    fval = obj.value(*xnew);
    ls_neval++;
    Real fvalp = 0.0;

    Real alpha1 = 0.0;
    Real alpha2 = 0.0;
    Real a      = 0.0;
    Real b      = 0.0;
    Real x1     = 0.0;
    Real x2     = 0.0;

    while ( !status(LineSearchType_Backtracking,ls_neval,ls_ngrad,alpha,fold,gs,fval,x,s,obj) ) {
      if ( alpha == alpha0_ ) {
        alpha1 = -gs*alpha*alpha/(2.0*(fval-fold-gs*alpha));
      }
      else {
        x1 = fval-fold-alpha*gs;
        x2 = fvalp-fval-alpha2*gs;
        a = (1.0/(alpha - alpha2))*( x1/(alpha*alpha) - x2/(alpha2*alpha2));
        b = (1.0/(alpha - alpha2))*(-x1*alpha2/(alpha*alpha) + x2*alpha/(alpha2*alpha2));
        if ( std::abs(a) < Teuchos::ScalarTraits<Real>::eps() ) {
          alpha1 = -gs/(2.0*b);
        }
        else {
          alpha1 = (-b+sqrt(b*b-3.0*a*gs))/(3.0*a);
        }
        if ( alpha1 > 0.5*alpha ) {
          alpha1 = 0.5*alpha;
        }
      }

      alpha2 = alpha;
      fvalp  = fval;

      if ( alpha1 <= 0.1*alpha ) {
        alpha = 0.1*alpha;
      }
      else {
        alpha = alpha1;
      }

      xnew->set(x);
      xnew->axpy(alpha, s);
      fval = obj.value(*xnew);
      ls_neval++;
    }
  }

  void bisection( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
               const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    Real tl = 0.0;
    Real tr = alpha0_;
 
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(tr,s);

    Real val_tr = obj.value(*xnew); 
    ls_neval++;
    Real val_tl = fval;

    if ( status(LineSearchType_Bisection,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj) ) {
      alpha = tr;
      fval  = val_tr;
      return;
    }

    Real tc = (tl+tr)/2.0;
    xnew->set(x);
    xnew->axpy(tc,s);
    Real val_tc = obj.value(*xnew);
    ls_neval++;

    Real t1     = 0.0;
    Real val_t1 = 0.0;
    Real t2     = 0.0;
    Real val_t2 = 0.0;

    while (    !status(LineSearchType_Bisection,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj)  
            && std::abs(tr - tl) > tol_ ) {
      t1 = (tl+tc)/2.0;
      xnew->set(x);
      xnew->axpy(t1,s);
      val_t1 = obj.value(*xnew);
      ls_neval++;

      t2 = (tr+tc)/2.0;
      xnew->set(x);
      xnew->axpy(t2,s);
      val_t2 = obj.value(*xnew);
      ls_neval++;

      if (    ( (val_tl <= val_tr) && (val_tl <= val_t1) && (val_tl <= val_t2) && (val_tl <= val_tc) ) 
           || ( (val_t1 <= val_tr) && (val_t1 <= val_tl) && (val_t1 <= val_t2) && (val_t1 <= val_tc) ) ) {
        tr     = tc;
        val_tr = val_tc;
        tc     = t1;
        val_tc = val_t1;
      }
      else if ( ( (val_tc <= val_tr) && (val_tc <= val_tl) && (val_tc <= val_t1) && (val_tc <= val_t2) ) ) { 
        tl     = t1;
        val_tl = val_t1;
        tr     = t2;
        val_tr = val_t2;
      }
      else if (    ( (val_t2 <= val_tr) && (val_t2 <= val_tl) && (val_t2 <= val_t1) && (val_t2 <= val_tc) ) 
                || ( (val_tr <= val_tl) && (val_tr <= val_t1) && (val_tr <= val_t2) && (val_tr <= val_tc) ) ) {
        tl     = tc;
        val_tl = val_tc;
        tc     = t2;
        val_tc = val_t2;
      }
    }

    fval  = val_tr;
    alpha = tr;
  }

  void goldensection( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                      const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    Teuchos::RCP<Vector<Real> > grad = x.clone();
    Real c   = (-1.0+sqrt(5.0))/2.0;
    Real tl  = 0.0;
    Real tr  = alpha0_;

    xnew->set(x);
    xnew->axpy(tr,s);
    Real val_tr = obj.value(*xnew);
    ls_neval++;
    Real val_tl = fval;

    if ( status(LineSearchType_GoldenSection,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj) ) {
      alpha = tr;
      fval  = val_tr;
      return;
    }

    Real tc1 = c*tl + (1.0-c)*tr;
    Real tc2 = (1.0-c)*tl + c*tr;
   
    xnew->set(x);
    xnew->axpy(tc1,s);
    Real val_tc1 = obj.value(*xnew);
    ls_neval++;

    xnew->set(x);
    xnew->axpy(tc2,s);
    Real val_tc2 = obj.value(*xnew);
    ls_neval++;

    while (    !status(LineSearchType_GoldenSection,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj) 
            && (std::abs(tl-tr) >= tol_) ) {
      if ( val_tc1 > val_tc2 ) {
        tl      = tc1;
        val_tl  = val_tc1;
        tc1     = tc2;
        val_tc1 = val_tc2;
 
        tc2     = (1.0-c)*tl + c*tr;     
        xnew->set(x);
        xnew->axpy(tc2,s);
        val_tc2 = obj.value(*xnew);
        ls_neval++;
      }
      else {
        tr      = tc2;
        val_tr  = val_tc2;
        tc2     = tc1;
        val_tc2 = val_tc1;

        tc1     = c*tl + (1.0-c)*tr;
        xnew->set(x);
        xnew->axpy(tc1,s);
        val_tc1 = obj.value(*xnew);
        ls_neval++;
      }
    }
    alpha = tr;
    fval  = val_tr;  
  }


  void brents( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
               const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    Real tl = 0.0;         // Left interval point
    Real tr = alpha0_;     // Right interval point
    Real tc = 0.0;         // Center interval point

    xnew->set(x);
    xnew->axpy(tr, s);
 
    Real val_tl = fval;
    Real val_tr = obj.value(*xnew);
    Real val_tc = 0.0;
    ls_neval++;

    if ( status(LineSearchType_Brents,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj) ) {
      alpha = tr;
      fval  = val_tr;
      return;
    }

    // Determine bracketing triple
    const Real gr                = (1.0+sqrt(5.0))/2.0;
    const Real inv_gr2           = 1.0/(gr*gr);
    const Real goldinv           = 1.0/(1.0+gr);
    const Real tiny              = sqrt(Teuchos::ScalarTraits<Real>::eps());
    const Real max_extrap_factor = 100.0;
    const int max_backtrack      = 8;
    Real tmp    = 0.0;
    Real q      = 0.0;
    Real r      = 0.0; 
    Real tm     = 0.0;
    Real tlim   = 0.0; 
    Real val_tm = 0.0;

    int itbt = 0;
    while ( val_tr > val_tl ) {
      if ( itbt <= max_backtrack ) {
        tc     = tr;
        val_tc = val_tr;

        tr     = goldinv * (tc + gr*tl);
        xnew->set(x);
        xnew->axpy(tr,s);
        val_tr = obj.value(*xnew);
        ls_neval++;
      }
      else {
        tmp    = tl;
        tl     = tr;
        tr     = tmp;
        tmp    = val_tr;
        val_tr = val_tl;
        val_tl = tmp;
        tc     = 0.0;
      }
      itbt++;
    }
    if ( std::abs(tc) < Teuchos::ScalarTraits<Real>::eps() ) {
      tc = tr + (gr-1.0)*(tr-tl);
      xnew->set(x);
      xnew->axpy(tc,s);
      val_tc = obj.value(*xnew);
      ls_neval++;
    }

    int itb = 0;
    while ( val_tr >= val_tc ) {
      q = ( val_tr-val_tl ) * (tr - tc);
      r = ( val_tr-val_tc ) * (tr - tl);
      tmp = fabs(q-r);
      tmp = (tmp > tiny ? tmp : -tmp);
      tm  = tr - (q*(tr-tc) - r*(tr-tl))/(2.0*tmp);

      tlim = tl + max_extrap_factor * (tc-tr);

      if ( (tr-tm)*(tm-tc) > 0.0 ) {
        xnew->set(x);
        xnew->axpy(tm,s);
        val_tm = obj.value(*xnew);
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
        xnew->set(x);
        xnew->axpy(tm,s);
        val_tm = obj.value(*xnew);
        ls_neval++;
      }
      else if ( (tc - tm)*(tm -tlim) > 0.0 ) {
        xnew->set(x);
        xnew->axpy(tm,s);
        val_tm = obj.value(*xnew);
        ls_neval++;
        if ( val_tm < val_tc ) {
          tr     = tc;
          val_tr = val_tc;

          tc     = tm;
          val_tc = val_tm;

          tm     = tc + gr*(tc-tr);
          xnew->set(x);
          xnew->axpy(tm,s);
          val_tm = obj.value(*xnew);
          ls_neval++;
        }
      }
      else if ( (tm-tlim)*(tlim-tc) >= 0.0 ) {
        tm = tlim;
        xnew->set(x);
        xnew->axpy(tm,s);
        val_tm = obj.value(*xnew);
        ls_neval++;
      }
      else {
        tm = tc + gr*(tc-tr);
        xnew->set(x);
        xnew->axpy(tm,s);
        val_tm = obj.value(*xnew);
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
    Real t     = 0.0;
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

    while (    !status(LineSearchType_Brents,ls_neval,ls_ngrad,t,fval,gs,ft,x,s,obj)
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
      xnew->set(x);
      xnew->axpy(u,s);
      fu = obj.value(*xnew);
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

    // Safeguard if Brent's does not find a minimizer
    if ( std::abs(alpha) < Teuchos::ScalarTraits<Real>::eps() ) {
      simplebacktracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
  }

};

}

#endif
