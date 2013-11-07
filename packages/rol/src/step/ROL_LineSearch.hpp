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
  LineSearchStep_NewtonKrylov
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

  bool status( const LineSearchType type, const int it, const Real alpha, 
               const Real fold, const Real sgold, const Real fnew, 
               const Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj ) { 
    if ( it >= maxit_ ) { 
      return true;
    }
    if ( type == LineSearchType_SimpleBacktracking || type == LineSearchType_Backtracking ) {
      if ( fnew <= fold + c1_*alpha*sgold ) {
        return true;
      }
    }
    else {
      if ( cond_ == LineSearchCondition_Wolfe ) {
        Teuchos::RCP<Vector<Real> > xnew = x.clone();
        xnew->set(x);
        xnew->axpy(alpha,s);   
        Teuchos::RCP<Vector<Real> > grad = x.clone();
        obj.gradient(*grad,*xnew);
        Real sgnew = grad->dot(s);

        if ( (fnew <= fold + c1_*alpha*sgold) && (sgnew >= c2_*sgold) ) {
          return true;
        }
      }
      else if ( cond_ == LineSearchCondition_StrongWolfe ) {
        Teuchos::RCP<Vector<Real> > xnew = x.clone();
        xnew->set(x);
        xnew->axpy(alpha,s);   
        Teuchos::RCP<Vector<Real> > grad = x.clone();
        obj.gradient(*grad,*xnew);
        Real sgnew = grad->dot(s);

        if ( (fnew <= fold + c1_*alpha*sgold) && (std::abs(sgnew) <= c2_*std::abs(sgold)) ) {
          return true;
        }
      }
      else if ( cond_ == LineSearchCondition_Goldstein ) {
        if ( (fnew <= fold + c1_*alpha*sgold) && (fnew >= fold + (1.0-c1_)*alpha*sgold) ) {
          return true;
        }
      }
    }
    return false;
  }

  void run( Real &alpha, Real &fval, int &it, 
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
    if ( type_ == LineSearchType_SimpleBacktracking ) {
      simplebacktracking( alpha, fval, it, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Backtracking ) {
      backtracking( alpha, fval, it, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Brents ) {
      brents( alpha, fval, it, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Bisection ) {
      bisection( alpha, fval, it, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_GoldenSection ) {
      goldensection( alpha, fval, it, gs, s, x, obj ); 
    }

    // Update Storage
    gs_old_    = gs;
    alpha_old_ = alpha;
  }

  void simplebacktracking( Real &alpha, Real &fval, int &it, 
                           const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    it    = 0;
    alpha = alpha0_;

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(alpha, s);

    Real fold = fval;
    fval = obj.value(*xnew);

    while ( !status(LineSearchType_SimpleBacktracking,it,alpha,fold,gs,fval,*xnew,s,obj) ) {
      alpha *= rho_;
      xnew->set(x);
      xnew->axpy(alpha, s);
      fval = obj.value(*xnew);
      it++;
    }
  }

  void backtracking( Real &alpha, Real &fval, int &it, 
                     const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    it    = 0;
    alpha = alpha0_;

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(alpha, s);

    Real fold = fval;
    fval = obj.value(*xnew);
    Real fvalp = 0.0;

    Real alpha1 = 0.0;
    Real alpha2 = 0.0;
    Real a      = 0.0;
    Real b      = 0.0;
    Real x1     = 0.0;
    Real x2     = 0.0;

    while ( !status(LineSearchType_Backtracking,it,alpha,fold,gs,fval,x,s,obj) ) {
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

      it++;
    }
  }

  void bisection( Real &alpha, Real &fval, int &it,
               const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    it = 0;

    Real tl = 0.0;
    Real tr = alpha0_;
 
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(tr,s);

    Real val_tr = obj.value(*xnew); 
    Real val_tl = fval;

    if ( status(LineSearchType_Bisection,it,tr,fval,gs,val_tr,x,s,obj) ) {
      alpha = tr;
      fval  = val_tr;
      return;
    }

    Real tc = (tl+tr)/2.0;
    xnew->set(x);
    xnew->axpy(tc,s);
    Real val_tc = obj.value(*xnew);

    Real t1     = 0.0;
    Real val_t1 = 0.0;
    Real t2     = 0.0;
    Real val_t2 = 0.0;

    while (    !status(LineSearchType_Bisection,it,tr,fval,gs,val_tr,x,s,obj)  
            && std::abs(tr - tl) > tol_ ) {
      t1 = (tl+tc)/2.0;
      xnew->set(x);
      xnew->axpy(t1,s);
      val_t1 = obj.value(*xnew);

      t2 = (tr+tc)/2.0;
      xnew->set(x);
      xnew->axpy(t2,s);
      val_t2 = obj.value(*xnew);

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
      it++;
    }

    fval  = val_tr;
    alpha = tr;
  }

  void goldensection( Real &alpha, Real &fval, int &it,
                      const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    Teuchos::RCP<Vector<Real> > grad = x.clone();
    Real c   = (-1.0+sqrt(5.0))/2.0;
    Real tl  = 0.0;
    Real tr  = alpha0_;

    it = 0;
    xnew->set(x);
    xnew->axpy(tr,s);
    Real val_tr = obj.value(*xnew);
    Real val_tl = fval;

    if ( status(LineSearchType_GoldenSection,it,tr,fval,gs,val_tr,x,s,obj) ) {
      alpha = tr;
      fval  = val_tr;
      return;
    }

    Real tc1 = c*tl + (1.0-c)*tr;
    Real tc2 = (1.0-c)*tl + c*tr;
   
    xnew->set(x);
    xnew->axpy(tc1,s);
    Real val_tc1 = obj.value(*xnew);

    xnew->set(x);
    xnew->axpy(tc2,s);
    Real val_tc2 = obj.value(*xnew);

    while (    !status(LineSearchType_GoldenSection,it,tr,fval,gs,val_tr,x,s,obj) 
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
      }
      it++;
    }
    alpha = tr;
    fval  = val_tr;  
  }


  void brents( Real &alpha, Real &fval, int &it, 
               const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
/*
    Real tl = 0.0;         // Left interval point
    Real tr = alpha0_;     // Right interval point

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(tr, s);
 
    Real val_tl = fval;
    Real val_tr = obj.value(*xnew);
 
    it = 0;
  
    if ( status(LineSearchType_Brents,it,tr,fval,gs,val_tr,x,s,obj) ) {
      alpha = tr;
      fval  = fnew;
      return;
    }

    Real tc = (tl+tr)/2.0; // Center invertval point
    xnew->set(x);
    xnew->axpy(tc, s);
    Real val_tc = obj.value(*xnew);

    Real t      = tr;
    Real val_t  = val_tr;
 
    while (    !status(LineSearchType_Brents,it,t,fval,gs,fnew,x,s,obj)
            && std::abs(tl - tr) > tol_ ) {

      

      if (    std::abs(val_tr - val_tc) > Teuchos::ScalarTraits<Real>::eps() 
           && std::abs(val_tl - val_tc) > Teuchos::ScalarTraits<Real>::eps() ) {
        t =   tl*val_tr*val_tc/((val_tl-val_tr)*(val_tl-val_tc)) 
            + tr*val_tl*val_tc/((val_tr-val_tl)*(val_tr-val_tc)) 
            + tc*val_tl*val_tr/((val_tc-val_tr)*(val_tc-val_tl));
      }
      else {
        t = tr - val_tr * (tr - tl) / (val_tr - val_tl);
      }
   
      Real tmp2 = (3.0 * tr + tl)/4.0;
      if ( (!(((t > tmp2) && (t < tr)) || ((t < tmp2) && (t > tr)))) 
          || (mflag && (std::abs(t-tr) >= (std::abs(tr-tc)/2.0)))
          || (!mflag && (std::abs(t-tr) >= (std::abs(tc-t0)/2.0)))) {
        t = (tl+tr)/2.0;
        mflag = true;
      }
      else {
        if ( (mflag && (std::abs(tr-tc) < tol_)) 
            || (!mflag && (std::abs(tc-t0) < tol_)) ) {
          t = (tl+tr)/2;
          mflag = true;
        }
        else {
          mflag = false;
        }
      }

      xnew->set(x);
      xnew->axpy(t,s);
      obj.gradient(*grad,*xnew);
      val_t = grad->dot(s);
      fnew = obj.value(*xnew);

      t0     = tc;
      tc     = tr;
      val_tc = val_tr;
      if ( val_tl * val_t < 0 ) {
        tr     = t; 
        val_tr = val_t;
      }
      else {
        tl     = t;
        val_tl = val_t;
      }

      if ( std::abs(val_tl) < std::abs(val_tr) ) {
        Real tmp = tl;
        tl       = tr;
        tr       = tmp;
        tmp      = val_tl;
        val_tl   = val_tr;
        val_tr   = tmp;
      }
      it++;
    } 
    alpha = t;
    fval  = obj.value(*xnew);
*/
  }

};

}

#endif
