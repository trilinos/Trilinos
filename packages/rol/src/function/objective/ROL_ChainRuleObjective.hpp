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

#ifndef ROL_CHAIN_RULE_OBJECTIVE_HPP
#define ROL_CHAIN_RULE_OBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::ChainRuleObjective
    \brief Defines an objective of the form f(g(x)) where

    \f$g:\mathcal{X}\to\mathcal{Y}\f$ and \f$f:\mathcal{Y}\to\mathbb{R}\f$

    It is assumed that both $f$ and $g$ are twice differentiable and that
    the mapping performed by \f$g\f$ is implemented by a ROL::Constraint,
    while a ROL::Objective implements the mapping performed by \f$f\f$.

---
*/


namespace ROL {

template<typename Real>
class ChainRuleObjective : public ROL::Objective<Real> {
public:

  /** \brief Constructor
      @param[in]         obj    is the objective function that performs the mapping \f$f:\mathcal{Y}\to\mathbb{R}\f$
      @param[in]         con    is the constraint function that performs the mapping \f$g:\mathcal{X}\to\mathbb{Y}\f$
      @param[in]         x      is an optimization space vector (\f$x\in\mathcal{X}\f$) provided for allocating memory for intermediate computations
      @param[in]               is a constraint space dual vector (\f$l\in\mathcal{Y}^\ast\f$) provided for allocating memory for intermediate computations
  */
  ChainRuleObjective( const Ptr<Objective<Real>>&  obj,
                      const Ptr<Constraint<Real>>& con,
                      const Vector<Real>& x,
                      const Vector<Real>& l ) 
  : obj_(obj), con_(con), g_(l.clone()), y_(l.dual().clone()), Jv_(l.dual().clone()), 
    HJv_(l.clone()), JtHJv_(x.dual().clone()), tol_(0) {}

  virtual ~ChainRuleObjective() = default;

    /** \brief Update objective function. 

      This function updates the objective function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          type   is the type of update requested.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    con_->update(x,type,iter);
    con_->value(*y_,x,tol_);
    obj_->update(*y_,type,iter);
  }

  /** \brief Update objective function. 

      This function updates the objective function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
    con_->value(*y_,x,tol_);
    obj_->update(*y_,flag,iter);
  }

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual Real value( const Vector<Real> &x, Real &tol ) {
    con_->value(*y_,x,tol);
    return obj_->value(*y_,tol);
  }

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.

    */
  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    con_->value(*y_,x,tol); 
    obj_->gradient(*g_,*y_,tol);
    con_->applyAdjointJacobian(g,*g_,x,tol);
  }

  /** \brief Apply Hessian approximation to vector.

      This function applies the Hessian of the objective function to the vector \f$v\f$.
      @param[out]         hv  is the the action of the Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    con_->value(*y_,x,tol);
    obj_->gradient(*g_,*y_,tol);
    con_->applyJacobian(*Jv_,v,x,tol);
    obj_->hessVec(*HJv_,*Jv_,*y_,tol);
    con_->applyAdjointJacobian(*JtHJv_,*HJv_,x,tol);
    con_->applyAdjointHessian(hv,*g_,v,x,tol);
    hv.plus(*JtHJv_);
  }

private:

  const Ptr<Objective<Real>>  obj_;
  const Ptr<Constraint<Real>> con_;
  Ptr<Vector<Real>>           g_, y_, Jv_, HJv_, JtHJv_;    
  Real                        tol_;

}; // class ChainRuleObjective

} // namespace ROL

#endif // ROL_CHAIN_RULE_OBJECTIVE_HPP
