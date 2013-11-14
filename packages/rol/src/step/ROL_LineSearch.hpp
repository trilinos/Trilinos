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

namespace ROL { 

template<class Real>
class LineSearch {
private:

  ELineSearch         els_;
  ECurvatureCondition econd_;
  EDescent            edesc_;

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
  LineSearch( ELineSearch els, ECurvatureCondition econd, EDescent edesc,  
              int maxit = 20, Real c1 = 1.e-4, Real c2 = 0.9, Real tol = 1.e-8, Real rho = 0.5 ) : 
    els_(els), econd_(econd), edesc_(edesc), maxit_(maxit), c1_(c1), c2_(c2), tol_(tol), rho_(rho) {
    if ( c2_ <= c1_ ) {
      c1_ = 1.e-4;
      c2_ = 0.9;
      if ( edesc_ == DESCENT_NONLINEARCG ) {
        c2_ = 0.1;
      }
    }

    gs_old_    = 0.0;
    alpha_old_ = 0.0;
    firstIt_   = true;
  }

  bool status( const ELineSearch type, int &ls_neval, int &ls_ngrad, const Real alpha, 
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
    if ( armijo && (type != LINESEARCH_BACKTRACKING && type != LINESEARCH_CUBICINTERP) ) {
      if ( (econd_ == CURVATURECONDITION_GOLDSTEIN) && (fnew >= fold + (1.0-c1_)*alpha*sgold) ) {
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
   
        if (    ((econd_ == CURVATURECONDITION_WOLFE)       && (sgnew >= c2_*sgold))
             || ((econd_ == CURVATURECONDITION_STRONGWOLFE) && (std::abs(sgnew) <= c2_*std::abs(sgold))) ) {
          curvcond = true;
        }
      }
    }

    if (type == LINESEARCH_BACKTRACKING || type == LINESEARCH_CUBICINTERP) {
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
      if ( edesc_ == DESCENT_STEEPEST || edesc_ == DESCENT_NONLINEARCG ) {
        alpha0_ = alpha_old_*gs_old_/gs;      
      }
    }

    // Run Linesearch
    ls_neval = 0;
    ls_ngrad = 0;
    if ( els_ == LINESEARCH_BACKTRACKING ) {
      simplebacktracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( els_ == LINESEARCH_CUBICINTERP ) {
      backtracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( els_ == LINESEARCH_BRENTS ) {
      brents( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( els_ == LINESEARCH_BISECTION ) {
      bisection( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
    else if ( els_ == LINESEARCH_GOLDENSECTION ) {
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

    while ( !status(LINESEARCH_BACKTRACKING,ls_neval,ls_ngrad,alpha,fold,gs,fval,*xnew,s,obj) ) {
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

    while ( !status(LINESEARCH_CUBICINTERP,ls_neval,ls_ngrad,alpha,fold,gs,fval,x,s,obj) ) {
      if ( alpha == alpha0_ ) {
        alpha1 = -gs*alpha*alpha/(2.0*(fval-fold-gs*alpha));
      }
      else {
        x1 = fval-fold-alpha*gs;
        x2 = fvalp-fval-alpha2*gs;
        a = (1.0/(alpha - alpha2))*( x1/(alpha*alpha) - x2/(alpha2*alpha2));
        b = (1.0/(alpha - alpha2))*(-x1*alpha2/(alpha*alpha) + x2*alpha/(alpha2*alpha2));
        if ( std::abs(a) < ROL_EPSILON ) {
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

    if ( status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj) ) {
      alpha = t;
      fval  = val_t;
      return;
    }

    Real tc = (tl+tr)/2.0;
    xnew->set(x);
    xnew->axpy(tc,s);
    Real val_tc = obj.value(*xnew);
    ls_neval++;

    if ( val_tc < val_t ) {
      t     = tc;
      val_t = val_tc;
    }

    Real t1     = 0.0;
    Real val_t1 = 0.0;
    Real t2     = 0.0;
    Real val_t2 = 0.0;

    while (    !status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj)  
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

    if ( status(LINESEARCH_GOLDENSECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj) ) {
      alpha = t;
      fval  = val_t;
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

    while (    !status(LINESEARCH_GOLDENSECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj) 
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

    if ( status(LINESEARCH_BRENTS,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj) ) {
      alpha = t;
      fval  = val_t;
      return;
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
      xnew->set(x);
      xnew->axpy(tr,s);
      val_tr = obj.value(*xnew);
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
      xnew->set(x);
      xnew->axpy(tc,s);
      val_tc = obj.value(*xnew);
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

    if ( status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj) ) {
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

    if ( status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj) ) {
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

    while (    !status(LINESEARCH_BRENTS,ls_neval,ls_ngrad,t,fval,gs,ft,x,s,obj)
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
    if ( std::abs(alpha) < ROL_EPSILON ) {
      simplebacktracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj ); 
    }
  }

};

}

#endif
