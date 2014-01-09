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

/** \file
    \brief  Contains definitions for helper functions in ROL.
        \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_HELPERFUNCTIONS_HPP
#define ROL_HELPERFUNCTIONS_HPP

#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraints.hpp"
#include "ROL_Secant.hpp"

namespace ROL {

  template<class Real>
  void applyHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, 
                     Objective<Real> &obj, Real tol = 0.0,
                     Teuchos::RCP<Secant<Real> > &secant = Teuchos::null, bool useSecant = false ) {
    if ( secant != Teuchos::null && useSecant ) {
      secant->applyB( Hv, v, x );
    }
    else {
      obj.hessVec( Hv, v, x, tol );
    }
  }

  template<class Real>
  void applyReducedHessVec( Vector<Real> &Hp, const Vector<Real> &p, const Vector<Real> &g, const Vector<Real> &x,
                            Constraints<Real> &con, Objective<Real> &obj, Real tol = 0.0, 
                            Teuchos::RCP<Secant<Real> > &secant = Teuchos::null, bool useSecant = false ) {
    if ( con.isActivated() ) {
      Teuchos::RCP<Vector<Real> > pnew = x.clone();
      pnew->set(p);
      con.pruneActive(*pnew,g,x);
      applyHessVec(Hp,*pnew,x,obj,tol,secant,useSecant);
      con.pruneActive(Hp,g,x);
      pnew->set(p);
      con.pruneInactive(*pnew,g,x);
      Hp.plus(*pnew);
    }
    else {
      applyHessVec(Hp,p,x,obj,tol,secant,useSecant);
    }
  }

  template<class Real>
  void applyInvHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, 
                        Objective<Real> &obj, Real tol = 0.0, 
                        Teuchos::RCP<Secant<Real> > &secant = Teuchos::null, bool useSecant = false ) {
    if ( secant != Teuchos::null && useSecant ) {
      secant->applyH(Hv,v,x);
    }
    else {
      obj.invHessVec(Hv,v,x,tol);
    }
  }

  template<class Real>
  void applyReducedInvHessVec( Vector<Real> &Hp, const Vector<Real> &p, const Vector<Real> &g, const Vector<Real> &x,
                               Constraints<Real> &con, Objective<Real> &obj, Real tol = 0.0, 
                               Teuchos::RCP<Secant<Real> > &secant = Teuchos::null, bool useSecant = false ) {
    if ( con.isActivated() ) {
      Teuchos::RCP<Vector<Real> > pnew = x.clone();
      pnew->set(p);
      con.pruneActive(*pnew,g,x);
      applyInvHessVec(Hp,*pnew,x,obj,tol,secant,useSecant);
      con.pruneActive(Hp,g,x);
      pnew->set(p);
      con.pruneInactive(*pnew,g,x);
      Hp.plus(*pnew);
    }
    else {
      applyInvHessVec(Hp,p,x,obj,tol,secant,useSecant);
    }
  }

  template<class Real>
  void applyPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &x, 
                     Objective<Real> &obj, Real tol = 0.0, 
                     Teuchos::RCP<Secant<Real> > &secant = Teuchos::null, bool useSecant = false ) {
    if ( secant != Teuchos::null && useSecant ) {
      secant->applyH( Mv, v, x );
    }
    else {
      obj.precond( Mv, v, x );
    }
  }

  template<class Real>
  void applyReducedPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x,
                            Constraints<Real> &con, Objective<Real> &obj, Real tol = 0.0, 
                            Teuchos::RCP<Secant<Real> > &secant = Teuchos::null, bool useSecant = false ) {
    if ( con.isActivated() ) {
      Teuchos::RCP<Vector<Real> > vnew = x.clone();
      vnew->set(v);
      con.pruneActive(*vnew,g,x);
      applyPrecond(Mv,*vnew,x,obj,tol,secant,useSecant);
      con.pruneActive(Mv,g,x);
      vnew->set(v);
      con.pruneInactive(*vnew,g,x);
      Mv.plus(*vnew);
    }
    else {
      applyPrecond(Mv,v,x,obj,tol,secant,useSecant);
    }
  }


  template<class Real> 
  class projectedObjective {
  private:
    Teuchos::RCP<Objective<Real> >   obj_;
    Teuchos::RCP<Constraints<Real> > con_;
    Teuchos::RCP<Secant<Real> >      secant_;
    bool useSecantPrecond_;
    bool useSecantHessVec_;

  public:
    projectedObjective( Objective<Real> &obj, Constraints<Real> &con, Secant<Real> &secant, 
                        bool useSecantPrecond = false, bool useSecantHessVec = false ) {
      obj_              = Teuchos::rcp(&obj,    false);
      con_              = Teuchos::rcp(&con,    false);
      secant_           = Teuchos::rcp(&secant, false);
      useSecantPrecond_ = useSecantPrecond;
      useSecantHessVec_ = useSecantHessVec;
    }

    void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
      this->obj_->update(x,flag,iter);
      this->con_->update(x,flag,iter);
    }

    Real value( const Vector<Real> &x, Real &tol ) {
      return this->obj_->value(x,tol);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      this->obj_->gradient(g,x,tol);
    }

    Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
      return this->obj_->dirDeriv(x,d,tol);
    }

    void hessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      if ( this->useSecantHessVec_ ) {
        this->secant_->applyB( Hv, v, x );
      }
      else {
        this->obj->hessVec( Hv, v, x, tol );
      }
    }

    void reducedHessVec( Vector<Real> &Hp, const Vector<Real> &p, const Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > pnew = x.clone();
        pnew->set(p);
        this->con_->pruneActive(*pnew,g,x);
        this->hessVec(Hp,*pnew,x,tol);
        this->con_->pruneActive(Hp,g,x);
        pnew->set(p);
        this->con_->pruneInactive(*pnew,g,x);
        Hp.plus(*pnew);
      }
      else {
        this->hessVec(Hp,p,x,tol);
      }
    }
 
    void invHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) { 
      if ( this->useSecantHessVec_ ) {
        this->secant_->applyH(Hv,v,x);
      }
      else {
        this->obj_->invHessVec(Hv,v,x,tol);
      }
    }

    void reducedInvHessVec( Vector<Real> &Hp, const Vector<Real> &p, const Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > pnew = x.clone();
        pnew->set(p);
        this->con_->pruneActive(*pnew,g,x);
        this->invHessVec(Hp,*pnew,x,tol);
        this->con_->pruneActive(Hp,g,x);
        pnew->set(p);
        this->con_->pruneInactive(*pnew,g,x);
        Hp.plus(*pnew);
      }
      else {
        this->invHessVec(Hp,p,x,tol);
      }
    }

    void precond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      if ( this->useSecantPrecond_ ) {
        this->secant_->applyH( Mv, v, x );
      }
      else {
        this->obj_->precond( Mv, v, x );
      }
    }

    void reducedPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        vnew->set(v);
        this->con_->pruneActive(*vnew,g,x);
        this->precond(Mv,*vnew,x,tol);
        this->con_->pruneActive(Mv,g,x);
        vnew->set(v);
        this->con_->pruneInactive(*vnew,g,x);
        Mv.plus(*vnew);
      }
      else {
        this->precond(Mv,v,x,tol);
      }
    }

    void project( Vector<Real> &x ) {
      this->con_->project(x);
    } 

    void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x ) {
      this->con_->pruneActive(v,g,x);
    }

    void pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x ) {
      this->con_->pruneInactive(v,g,x);
    }

    bool isFeasible( const Vector<Real> &v ) {
      return this->con_->isFeasible();
    }

    bool isActivated(void) {
      return this->con_->isActivated();
    }
  }; 

} // namespace ROL

#endif 
