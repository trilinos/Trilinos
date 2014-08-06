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
  Real c3_;
  Real tol_;
  Real rho_;
  Real alpha0_;
  bool useralpha_;
  int algo_iter_;

  Teuchos::RCP<Vector<Real> > grad_;
  Real eps_;

  void updateIterate(Vector<Real> &xnew, const Vector<Real> &x, const Vector<Real> &s, Real alpha, 
                     BoundConstraint<Real> &con ) {
    xnew.set(x); 
    xnew.axpy(alpha,s);
    if ( con.isActivated() ) {
      con.project(xnew);
    }
  }

public:

  virtual ~LineSearch() {}

  // Constructor
  LineSearch( Teuchos::ParameterList &parlist ) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    els_   = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    econd_ = StringToECurvatureCondition( parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    // Linesearc Parameters
    maxit_     = parlist.get("Maximum Number of Function Evaluations",            20);
    c1_        = parlist.get("Sufficient Decrease Parameter",                     1.e-4);
    c2_        = parlist.get("Curvature Conditions Parameter",                    0.9);
    c3_        = parlist.get("Curvature Conditions Parameter: Generalized Wolfe", 0.6);
    tol_       = parlist.get("Bracketing Tolerance",                              1.e-8);
    rho_       = parlist.get("Backtracking Rate",                                 0.5);
    alpha0_    = parlist.get("Initial Linesearch Parameter",                      1.0);
    useralpha_ = parlist.get("User Defined Linesearch Parameter",                 false);
    algo_iter_ = 0;

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

  void setData( const Teuchos::RCP<Vector<Real> > &grad, const Real eps = 0.0 ) { 
    this->grad_ = grad->clone();
    this->grad_->set(*grad);
    this->eps_ = eps;
  }

  bool status( const ELineSearch type, int &ls_neval, int &ls_ngrad, const Real alpha, 
               const Real fold, const Real sgold, const Real fnew, 
               const Vector<Real> &x, const Vector<Real> &s, 
               Objective<Real> &obj, BoundConstraint<Real> &con ) { 
    Real tol = std::sqrt(ROL_EPSILON);

    // Check Armijo Condition
    bool armijo = false;
    if ( con.isActivated() ) {
      Real gs = 0.0;
      Teuchos::RCP<Vector<Real> > d = x.clone();
      if ( this->edesc_ == DESCENT_STEEPEST ) {
        this->updateIterate(*d,x,s,alpha,con);
        d->scale(-1.0);
        d->plus(x);
        gs = -s.dot(*d);
      }
      else {
        d->set(s);
        d->scale(-1.0);
        con.pruneActive(*d,*(this->grad_),x,this->eps_);
        gs = alpha*(this->grad_)->dot(*d);
        d->zero();
        this->updateIterate(*d,x,s,alpha,con);
        d->scale(-1.0);
        d->plus(x);
        con.pruneInactive(*d,*(this->grad_),x,this->eps_);
        gs += (this->grad_)->dot(*d);
      }
      if ( fnew <= fold - this->c1_*gs ) {
        armijo = true;
      }
    }
    else {
      if ( fnew <= fold + this->c1_*alpha*sgold ) {
        armijo = true;
      }
    }

    // Check Maximum Iteration
    bool itcond = false;
    if ( ls_neval >= this->maxit_ ) { 
      itcond = true;
    }

    // Check Curvature Condition
    bool curvcond = false;
    if ( armijo && ((type != LINESEARCH_BACKTRACKING && type != LINESEARCH_CUBICINTERP) ||
                    (this->edesc_ == DESCENT_NONLINEARCG)) ) {
      if (this->econd_ == CURVATURECONDITION_GOLDSTEIN) {
        if (fnew >= fold + (1.0-this->c1_)*alpha*sgold) {
          curvcond = true;
        }
      }
      else if (this->econd_ == CURVATURECONDITION_NULL) {
        curvcond = true;
      }
      else { 
        Teuchos::RCP<Vector<Real> > xnew = x.clone();
        this->updateIterate(*xnew,x,s,alpha,con);
        Teuchos::RCP<Vector<Real> > grad = x.clone();
        obj.update(*xnew);
        obj.gradient(*grad,*xnew,tol);
        Real sgnew = 0.0;
        if ( con.isActivated() ) {
          Teuchos::RCP<Vector<Real> > d = x.clone();
          d->set(s);
          d->scale(-alpha);
          con.pruneActive(*d,s,x);
          sgnew = -grad->dot(*d);
        }
        else {
          sgnew = grad->dot(s);
        }
        ls_ngrad++;
   
        if (    ((this->econd_ == CURVATURECONDITION_WOLFE)       
                     && (sgnew >= this->c2_*sgold))
             || ((this->econd_ == CURVATURECONDITION_STRONGWOLFE) 
                     && (std::abs(sgnew) <= this->c2_*std::abs(sgold)))
             || ((this->econd_ == CURVATURECONDITION_GENERALIZEDWOLFE) 
                     && (this->c2_*sgold <= sgnew && sgnew <= -this->c3_*sgold))
             || ((this->econd_ == CURVATURECONDITION_APPROXIMATEWOLFE) 
                     && (this->c2_*sgold <= sgnew && sgnew <= (2.0*this->c1_ - 1.0)*sgold)) ) {
          curvcond = true;
        }
      }
    }

    if (type == LINESEARCH_BACKTRACKING || type == LINESEARCH_CUBICINTERP) {
      if (this->edesc_ == DESCENT_NONLINEARCG) {
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

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    // Determine Initial Step Length
    if (this->useralpha_ || this->els_ == LINESEARCH_ITERATIONSCALING) {
      alpha = this->alpha0_;
    }
    else if ( this->edesc_ == DESCENT_STEEPEST || this->edesc_ == DESCENT_NONLINEARCG ) {
      this->updateIterate(*xnew,x,s,1.0,con);
      Real ftol = 0.0;
      // TODO: Think about reusing for efficiency!
      obj.update(*xnew);
      Real fnew = obj.value(*xnew, ftol);
      alpha = -gs/(2.0*(fnew-fval-gs));
      this->updateIterate(*xnew,x,s,alpha,con);
      obj.update(*xnew);
      fnew = obj.value(*xnew, ftol);
      bool stat = status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,alpha,fval,gs,fnew,x,s,obj,con);
      if (!stat) {
        alpha = 1.0;
      }
      ls_neval++;
    }
    else { // Newton-type methods
      alpha = 1.0;
    }

    // Run Linesearch
    ls_neval = 0;
    ls_ngrad = 0;
    if ( this->els_ == LINESEARCH_ITERATIONSCALING ) {
      this->iterationscaling( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con );
    }
    else if ( this->els_ == LINESEARCH_BACKTRACKING ) {
      this->simplebacktracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con ); 
    }
    else if ( this->els_ == LINESEARCH_CUBICINTERP ) {
      this->backtracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con ); 
    }
    else if ( this->els_ == LINESEARCH_BRENTS ) {
      this->brents( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con ); 
    }
    else if ( this->els_ == LINESEARCH_BISECTION ) {
      this->bisection( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con ); 
    }
    else if ( this->els_ == LINESEARCH_GOLDENSECTION ) {
      this->goldensection( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con ); 
    }
  }

  void iterationscaling( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                         const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                         Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);

    this->algo_iter_++;
    alpha /= this->algo_iter_;

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    this->updateIterate(*xnew,x,s,alpha,con);

    obj.update(*xnew);
    fval = obj.value(*xnew,tol);
    ls_neval++;
  }

  void simplebacktracking( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                           const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                           Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    this->updateIterate(*xnew,x,s,alpha,con);

    Real fold = fval;
    obj.update(*xnew);
    fval = obj.value(*xnew,tol);
    ls_neval++;

    while ( !status(LINESEARCH_BACKTRACKING,ls_neval,ls_ngrad,alpha,fold,gs,fval,*xnew,s,obj,con) ) {
      alpha *= this->rho_;
      this->updateIterate(*xnew,x,s,alpha,con);
      obj.update(*xnew);
      fval = obj.value(*xnew,tol);
      ls_neval++;
    }
  }

  void backtracking( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                     const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                     Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);

    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    this->updateIterate(*xnew,x,s,alpha,con);

    Real fold = fval;
    obj.update(*xnew);
    fval = obj.value(*xnew,tol);
    ls_neval++;
    Real fvalp = 0.0;

    Real alpha1 = 0.0;
    Real alpha2 = 0.0;
    Real a      = 0.0;
    Real b      = 0.0;
    Real x1     = 0.0;
    Real x2     = 0.0;

    bool first_iter = true;

    while ( !status(LINESEARCH_CUBICINTERP,ls_neval,ls_ngrad,alpha,fold,gs,fval,x,s,obj,con) ) {
      if ( first_iter ) {
        alpha1 = -gs*alpha*alpha/(2.0*(fval-fold-gs*alpha));
        first_iter = false;
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
        alpha *= 0.1;
      }
      else if ( alpha1 >= 0.5*alpha ) {
        alpha *= 0.5;
      }
      else {
        alpha = alpha1;
      }

      this->updateIterate(*xnew,x,s,alpha,con);
      obj.update(*xnew);
      fval = obj.value(*xnew,tol);
      ls_neval++;
    }
  }

  void bisection( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
                  const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                  Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);

    // Compute value phi(0)
    Real tl = 0.0;
    Real val_tl = fval;

    // Compute value phi(alpha)
    Real tr = alpha;
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    this->updateIterate(*xnew,x,s,tr,con);
    obj.update(*xnew);
    Real val_tr = obj.value(*xnew,tol); 
    ls_neval++;

    // Check if phi(alpha) < phi(0)
    if ( val_tr < val_tl ) {
      if ( status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj,con) ) {
        alpha = tr;
        fval  = val_tr;
        return;
      }
    }

    // Get min( phi(0), phi(alpha) )
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

    // Compute value phi(midpoint)
    Real tc = (tl+tr)/2.0;
    this->updateIterate(*xnew,x,s,tc,con);
    Real val_tc = obj.value(*xnew,tol);
    ls_neval++;

    // Get min( phi(0), phi(alpha), phi(midpoint) )
    if ( val_tc < val_t ) {
      t     = tc;
      val_t = val_tc;
    }

    Real t1     = 0.0;
    Real val_t1 = 0.0;
    Real t2     = 0.0;
    Real val_t2 = 0.0;

    while (    !status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con)  
            && std::abs(tr - tl) > this->tol_ ) {
      t1 = (tl+tc)/2.0;
      this->updateIterate(*xnew,x,s,t1,con);
      obj.update(*xnew);
      val_t1 = obj.value(*xnew,tol);
      ls_neval++;

      t2 = (tr+tc)/2.0;
      this->updateIterate(*xnew,x,s,t2,con);
      obj.update(*xnew);
      val_t2 = obj.value(*xnew,tol);
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
                      const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
                      Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);

    Teuchos::RCP<Vector<Real> > grad = x.clone();
    Real c   = (-1.0+sqrt(5.0))/2.0;

    // Compute value phi(0)
    Real tl  = 0.0;
    Real val_tl = fval;

    // Compute value phi(alpha)
    Real tr  = alpha;
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    this->updateIterate(*xnew,x,s,tr,con);
    obj.update(*xnew);
    Real val_tr = obj.value(*xnew,tol);
    ls_neval++;

    // Check if phi(alpha) < phi(0)
    if ( val_tr < val_tl ) {
      if ( status(LINESEARCH_GOLDENSECTION,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj,con) ) {
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

    // Compute value phi(t1)
    Real tc1 = c*tl + (1.0-c)*tr;
    this->updateIterate(*xnew,x,s,tc1,con);
    obj.update(*xnew);
    Real val_tc1 = obj.value(*xnew,tol);
    ls_neval++;

    // Compute value phi(t2)
    Real tc2 = (1.0-c)*tl + c*tr;
    this->updateIterate(*xnew,x,s,tc2,con);
    obj.update(*xnew);
    Real val_tc2 = obj.value(*xnew,tol);
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

    while (    !status(LINESEARCH_GOLDENSECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con) 
            && (std::abs(tl-tr) >= this->tol_) ) {
      if ( val_tc1 > val_tc2 ) {
        tl      = tc1;
        val_tl  = val_tc1;
        tc1     = tc2;
        val_tc1 = val_tc2;
 
        tc2     = (1.0-c)*tl + c*tr;     
        this->updateIterate(*xnew,x,s,tc2,con);
        obj.update(*xnew);
        val_tc2 = obj.value(*xnew,tol);
        ls_neval++;
      }
      else {
        tr      = tc2;
        val_tr  = val_tc2;
        tc2     = tc1;
        val_tc2 = val_tc1;

        tc1     = c*tl + (1.0-c)*tr;
        this->updateIterate(*xnew,x,s,tc1,con);
        obj.update(*xnew);
        val_tc1 = obj.value(*xnew,tol);
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
               const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
               Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON);

    // Compute value phi(0)
    Real tl = 0.0;         // Left interval point
    Real val_tl = fval;

    // Initialize value phi(t)
    Real tc = 0.0;        // Center interval point
    Real val_tc = 0.0;

    // Compute value phi(alpha)
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    Real tr = alpha;      // Right interval point
    this->updateIterate(*xnew,x,s,tr,con);
    obj.update(*xnew);
    Real val_tr = obj.value(*xnew,tol);
    ls_neval++;

    // Check if phi(alpha) < phi(0)
    if ( val_tr < val_tl ) {
      if ( status(LINESEARCH_BRENTS,ls_neval,ls_ngrad,tr,fval,gs,val_tr,x,s,obj,con) ) {
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
      this->updateIterate(*xnew,x,s,tr,con);
      obj.update(*xnew);
      val_tr = obj.value(*xnew,tol);
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
      this->updateIterate(*xnew,x,s,tc,con);
      obj.update(*xnew);
      val_tc = obj.value(*xnew,tol);
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

    if ( status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con) ) {
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
        this->updateIterate(*xnew,x,s,tm,con);
        obj.update(*xnew);
        val_tm = obj.value(*xnew,tol);
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
        this->updateIterate(*xnew,x,s,tm,con);
        obj.update(*xnew);
        val_tm = obj.value(*xnew,tol);
        ls_neval++;
      }
      else if ( (tc - tm)*(tm -tlim) > 0.0 ) {
        this->updateIterate(*xnew,x,s,tm,con);
        obj.update(*xnew);
        val_tm = obj.value(*xnew,tol);
        ls_neval++;
        if ( val_tm < val_tc ) {
          tr     = tc;
          val_tr = val_tc;

          tc     = tm;
          val_tc = val_tm;

          tm     = tc + gr*(tc-tr);
          this->updateIterate(*xnew,x,s,tm,con);
          obj.update(*xnew);
          val_tm = obj.value(*xnew,tol);
          ls_neval++;
        }
      }
      else if ( (tm-tlim)*(tlim-tc) >= 0.0 ) {
        tm = tlim;
        this->updateIterate(*xnew,x,s,tm,con);
        obj.update(*xnew);
        val_tm = obj.value(*xnew,tol);
        ls_neval++;
      }
      else {
        tm = tc + gr*(tc-tr);
        this->updateIterate(*xnew,x,s,tm,con);
        obj.update(*xnew);
        val_tm = obj.value(*xnew,tol);
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

    if ( status(LINESEARCH_BISECTION,ls_neval,ls_ngrad,t,fval,gs,val_t,x,s,obj,con) ) {
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

    while (    !status(LINESEARCH_BRENTS,ls_neval,ls_ngrad,t,fval,gs,ft,x,s,obj,con)
            && std::abs(t - tm) > this->tol_*(b-a) ) {
      if ( it < 2 ) {
        e = 2.0*(b-a);
      }
      tm = (a+b)/2.0;

      Real tol1 = this->tol_*std::abs(t) + tiny;
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
      this->updateIterate(*xnew,x,s,u,con);
      obj.update(*xnew);
      fu = obj.value(*xnew,tol);
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
      simplebacktracking( alpha, fval, ls_neval, ls_ngrad, gs, s, x, obj, con ); 
    }
  }

};

}

#endif
