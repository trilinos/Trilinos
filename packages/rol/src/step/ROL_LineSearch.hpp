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
  LineSearchType_Backtracking = 0, // Backtracking
  LineSearchType_Brents,           // Brents
  LineSearchType_Bisection         // Bisection
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

  bool status( const LineSearchType type, const int it, const Real alpha, const Real fold, const Real sgold, const Real fnew, const Real sgnew = 0.0 ) { 
    if ( it >= maxit_ ) { 
      return true;
    }
    if ( type == LineSearchType_Backtracking ) {
      if ( fnew <= fold + c1_*alpha*sgold ) {
        return true;
      }
    }
    else {
      if ( cond_ == LineSearchCondition_Wolfe ) {
        if ( (fnew <= fold + c1_*alpha*sgold) && (sgnew >= c2_*sgold) ) {
          return true;
        }
      }
      else if ( cond_ == LineSearchCondition_StrongWolfe ) {
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
      if ( step_ == LineSearchStep_Gradient ) {
        alpha0_ = alpha_old_*gs_old_/gs;      
      }
    }

    // Run Linesearch
    if ( type_ == LineSearchType_Backtracking ) {
      backtracking( alpha, fval, it, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Brents ) {
      brents( alpha, fval, it, gs, s, x, obj ); 
    }
    else if ( type_ == LineSearchType_Bisection ) {
      bisection( alpha, fval, it, gs, s, x, obj ); 
    }

    // Update Storage
    gs_old_    = gs;
    alpha_old_ = alpha;
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

    while ( !status(LineSearchType_Backtracking,it,alpha,fold,gs,fval) ) {
      alpha *= rho_;
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

    Teuchos::RCP<Vector<Real> > grad = x.clone();
    obj.gradient(*grad,*xnew);
    Real fnew = obj.value(*xnew); 

    Real val_tl = gs;
    Real val_tr = grad->dot(s);

    if ( val_tl*val_tr >= Teuchos::ScalarTraits<Real>::eps() ) {
      if ( status(LineSearchType_Bisection,it,tr,fval,gs,fnew,val_tr) ) {
        alpha = tr;
        fval  = fnew;
        return;
      }
      else {
        backtracking(alpha,fval,it,gs,s,x,obj);
        //alpha = tl; 
        return;
      }
    }

    Real tc = (tl+tr)/2.0;
    xnew->set(x);
    xnew->axpy(tc,s);
    obj.gradient(*grad,*xnew);
    Real val_tc = grad->dot(s);
    fnew        = obj.value(*xnew);
    it++;

    while (    !status(LineSearchType_Bisection,it,tc,fval,gs,fnew,val_tc)  
            && std::abs(val_tc) > Teuchos::ScalarTraits<Real>::eps()
            && std::abs(tr - tl)/2.0 > tol_ ) {
      if (    (val_tc < 0.0  && val_tl < 0.0) 
           || (val_tc > 0.0  && val_tl > 0.0) 
           || (val_tc == 0.0 && val_tl == 0.0) ) {
        tl = tc;
      }
      else {
        tr = tc;
      }

      tc = (tl + tr)/2.0;
     
      xnew->set(x);
      xnew->axpy(tc,s);
      obj.gradient(*grad,*xnew);
      val_tc = grad->dot(s);
      fnew   = obj.value(*xnew);
      it++;
    }

    fval  = obj.value(*xnew);
    alpha = tc;
  }

  void goldensection( Real &alpha, Real &fval, int &it,
                      const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    Real c   = (-1.0+sqrt(5.0))/2.0;
    Real tl  = 0.0;
    Real tr  = alpha0_;

    Real tc1 = c*tl + (1.0-c)*tr;
    Real tc2 = (1.0-c)*tl + c*tr;
    
  }


  void brents( Real &alpha, Real &fval, int &it, 
               const Real &gs, const Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj ) {
    Real tl = 0.0;     // Left interval point
    Real tr = alpha0_; // Right interval point
    Real tc = 0.0;     // Center invertval point
    Real t0 = Teuchos::ScalarTraits<Real>::rmax();

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(tr, s);
 
    Real val_tl = gs;
    Teuchos::RCP<Vector<Real> > grad = x.clone();
    obj.gradient(*grad,*xnew);
    Real val_tr = grad->dot(s);
    Real fnew = obj.value(*xnew);
 
    Real val_tc = 0.0;
    Real t      = tr;
    Real val_t  = val_tr;
 
    it = 0;
  
    if ( val_tl*val_tr >= Teuchos::ScalarTraits<Real>::eps() ) {
      if ( status(LineSearchType_Brents,it,tr,fval,gs,fnew,val_tr) ) {
        alpha = tr;
        fval  = fnew;
        return;
      }
      else {
        backtracking(alpha,fval,it,gs,s,x,obj);
        //alpha = tl; 
        return;
      }
    }

    if ( std::abs(val_tl) < std::abs(val_tr) ) {
      Real tmp = tl;
      tl       = tr;
      tr       = tmp;
      tmp      = val_tl;
      val_tl   = val_tr;
      val_tr   = tmp;
    }

    tc     = tl;
    val_tc = val_tl;
    bool mflag = true;
 
    while (    !status(LineSearchType_Brents,it,t,fval,gs,fnew,val_t)
            && std::abs(val_tr)  > Teuchos::ScalarTraits<Real>::eps() 
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
  }

};

}

#endif
