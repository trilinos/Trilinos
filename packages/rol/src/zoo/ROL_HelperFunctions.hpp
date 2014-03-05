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
  class ProjectedObjective {
  private:
    Teuchos::RCP<Objective<Real> >   obj_;
    Teuchos::RCP<Constraints<Real> > con_;
    Teuchos::RCP<Secant<Real> >      secant_;
    bool useSecantPrecond_;
    bool useSecantHessVec_;
    Real eps_;

  public:
    ProjectedObjective( Objective<Real> &obj, Constraints<Real> &con, Teuchos::RCP<Secant<Real> > &secant, 
                        bool useSecantPrecond = false, bool useSecantHessVec = false, Real eps = 0.0 ) {
      obj_              = Teuchos::rcp(&obj,    false);
      con_              = Teuchos::rcp(&con,    false);
      secant_           = secant; //Teuchos::rcp(&secant, false);
      useSecantPrecond_ = useSecantPrecond;
      useSecantHessVec_ = useSecantHessVec;
      eps_              = eps;
    }

    void update( const Vector<Real> &x, bool flag = true, int iter = -1, Real eps = 0.0 ) {
      this->obj_->update(x,flag,iter);
      this->con_->update(x,flag,iter);
      if ( std::abs(eps) >= ROL_EPSILON ) { this->eps_ = eps; }
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
        this->obj_->hessVec( Hv, v, x, tol );
      }
    }

    /** \brief Apply the reduced Hessian to a vector, v.  
               The reduced Hessian first removes elements of v 
               corresponding to the feasible indices from 
               the point p in the direction -d.
                   Hv   the Hessian times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   d    negative of search direction
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p, 
                         const Vector<Real> &d, const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        // Set vnew to v
        vnew->set(v);
        // Remove elements of vnew corresponding to binding set
        this->con_->pruneActive(*vnew,d,p,this->eps_);
        // Apply full Hessian to reduced vector
        this->hessVec(Hv,*vnew,x,tol);
        // Remove elements of Hv corresponding to binding set
        this->con_->pruneActive(Hv,d,p,this->eps_);
        // Set vnew to v
        vnew->set(v);                             
        // Remove Elements of vnew corresponding to complement of binding set
        this->con_->pruneInactive(*vnew,d,p,this->eps_); 
        // Fill complement of binding set elements in Hp with v
        Hv.plus(*vnew);                           
      }
      else {
        this->hessVec(Hv,v,x,tol);
      }
    }
 
    /** \brief Apply the reduced Hessian to a vector, v.  
               The reduced Hessian first removes elements of v 
               corresponding to the feasible indices from 
               the point p.
                   Hv   the Hessian times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p, 
                         const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        // Set vnew to v
        vnew->set(v);
        // Remove elements of vnew corresponding to binding set
        this->con_->pruneActive(*vnew,p,this->eps_);
        // Apply full Hessian to reduced vector
        this->hessVec(Hv,*vnew,x,tol);
        // Remove elements of Hv corresponding to binding set
        this->con_->pruneActive(Hv,p,this->eps_);
        // Set vnew to v
        vnew->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        this->con_->pruneInactive(*vnew,p,this->eps_);
        // Fill complement of binding set elements in Hp with v
        Hv.plus(*vnew);
      }
      else {
        this->hessVec(Hv,v,x,tol);
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

    /** \brief Apply the reduced inverse Hessian to a vector, v.  
               The reduced inverse Hessian first removes elements 
               of v corresponding to the feasible indices from 
               the point p in the direction -d.
                   Hv   the inverse Hessian times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   d    negative of search direction
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedInvHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p, 
                            const Vector<Real> &d, const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        // Set vnew to v
        vnew->set(v);
        // Remove elements of vnew corresponding to binding set
        this->con_->pruneActive(*vnew,d,p,this->eps_);
        // Apply full Hessian to reduced vector
        this->invHessVec(Hv,*vnew,x,tol);
        // Remove elements of Hv corresponding to binding set
        this->con_->pruneActive(Hv,d,p,this->eps_);
        // Set vnew to v
        vnew->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        this->con_->pruneInactive(*vnew,d,p,this->eps_);
        // Fill complement of binding set elements in Hv with v
        Hv.plus(*vnew);
      }
      else {
        this->invHessVec(Hv,v,x,tol);
      }
    }

    /** \brief Apply the reduced inverse Hessian to a vector, v.  
               The reduced inverse Hessian first removes elements 
               of v corresponding to the feasible indices from 
               the point p.
                   Hv   the inverse Hessian times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedInvHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p, 
                            const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        // Set vnew to v
        vnew->set(v);
        // Remove elements of vnew corresponding to binding set
        this->con_->pruneActive(*vnew,p,this->eps_);
        // Apply full Hessian to reduced vector
        this->invHessVec(Hv,*vnew,x,tol);
        // Remove elements of Hv corresponding to binding set
        this->con_->pruneActive(Hv,p,this->eps_);
        // Set vnew to v
        vnew->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        this->con_->pruneInactive(*vnew,p,this->eps_);
        // Fill complement of binding set elements in Hv with v
        Hv.plus(*vnew);
      }
      else {
        this->invHessVec(Hv,v,x,tol);
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

    /** \brief Apply the reduced preconditioner to a vector, v.  
               The reduced preconditioner first removes elements 
               of v corresponding to the feasible indices from 
               the point p in the direction -d.
                   Hv   the preconditioner times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   d    negative of search direction
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &p, 
                         const Vector<Real> &d, const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        // Set vnew to v
        vnew->set(v);
        // Remove elements of vnew corresponding to binding set
        this->con_->pruneActive(*vnew,d,p,this->eps_);
        // Apply full Hessian to reduced vector
        this->precond(Mv,*vnew,x,tol);
        // Remove elements of Mv corresponding to binding set
        this->con_->pruneActive(Mv,d,p,this->eps_);
        // Set vnew to v
        vnew->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        this->con_->pruneInactive(*vnew,d,p,this->eps_);
        // Fill complement of binding set elements in Mv with v
        Mv.plus(*vnew);
      }
      else {
        this->precond(Mv,v,x,tol);
      }
    }

    /** \brief Apply the reduced preconditioner to a vector, v.  
               The reduced preconditioner first removes elements 
               of v corresponding to the feasible indices from 
               the point p.
                   Hv   the preconditioner times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &p, 
                         const Vector<Real> &x, Real &tol ) {
      if ( this->con_->isActivated() ) {
        Teuchos::RCP<Vector<Real> > vnew = x.clone();
        // Set vnew to v
        vnew->set(v);
        // Remove elements of vnew corresponding to binding set
        this->con_->pruneActive(*vnew,p,this->eps_);
        // Apply full Hessian to reduced vector
        this->precond(Mv,*vnew,x,tol);
        // Remove elements of Mv corresponding to binding set
        this->con_->pruneActive(Mv,p,this->eps_);
        // Set vnew to v
        vnew->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        this->con_->pruneInactive(*vnew,p,this->eps_);
        // Fill complement of binding set elements in Mv with v
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
      this->con_->pruneActive(v,g,x,this->eps_);
    }

    void pruneActive( Vector<Real> &v, const Vector<Real> &x ) {
      this->con_->pruneActive(v,x,this->eps_);
    }

    void pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x ) {
      this->con_->pruneInactive(v,g,x,this->eps_);
    }

    void pruneInactive( Vector<Real> &v, const Vector<Real> &x ) {
      this->con_->pruneInactive(v,x,this->eps_);
    }

    bool isFeasible( const Vector<Real> &v ) {
      return this->con_->isFeasible(v);
    }

    bool isConActivated(void) {
      return this->con_->isActivated();
    }

    void computeProjectedStep( Vector<Real> &v, const Vector<Real> &x ) {
      this->con_->computeProjectedStep(v,x);
    } 
  }; 

} // namespace ROL

#endif 
