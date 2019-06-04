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

#ifndef ROL_OBJECTIVE_FROM_CONSTRAINT_H
#define ROL_OBJECTIVE_FROM_CONSTRAINT_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::ObjectiveFromConstraint
    \brief Form an objective function from a ROL::Constraint and a
           vector in the dual constraint space \f$\lambda\in \mathcal{C}^\ast\f$

    \f[ f(x;\lambda) = \langle \lambda, c(x)\rangle_{\mathcal{C}^*,\mathcal{C}} \f]
*/


namespace ROL {

template <class Real>
class ObjectiveFromConstraint : public Objective<Real> {

private:

  ROL::Ptr<Constraint<Real> > con_;
  ROL::Ptr<Vector<Real> >     l_;      // Lagrange multiplier 
  ROL::Ptr<Vector<Real> >     c_;      // Constraint vector


public:

  ObjectiveFromConstraint( const ROL::Ptr<Constraint<Real> > &con, 
                           const Vector<Real> &l ) :
    con_(con), l_(l.clone()), c_(l.dual().clone()) {
    l_->set(l);
  }

  virtual ~ObjectiveFromConstraint() {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    con_->update(x,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    con_->value(*c_,x,tol);
    return l_->dot(c_->dual());  
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    con_->applyAdjointJacobian(g,*l_,x,tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    con_->applyAdjointHessian(hv,*l_,v,x,tol);
  }

  void updateMultiplier( const Vector<Real> &l ) {
    l_->set(l);
  }

}; // class ObjectiveFromConstraint

} // namespace ROL

#endif // ROL_OBJECTIVE_FROM_CONSTRAINT_H
