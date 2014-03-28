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


#ifndef ROL_INEQUALITYCONSTRAINT_DEF_H
#define ROL_INEQUALITYCONSTRAINT_DEF_H

/** \class ROL::InequalityConstraint
    \brief Provides the definition of the inequality constraint function interface.
*/

namespace ROL {

template <class Real>
void InequalityConstraint<Real>::applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, const bool &adj, Real &tol ) {
  Real ctol = std::sqrt(ROL_EPSILON);

  if (!adj) {
    // Get step length.
    Real h = std::max(1.0,x.norm()/v.norm())*tol;
    //Real h = 2.0/(v.norm()*v.norm())*tol;

    // Compute constraint at x.
    Teuchos::RCP<Vector<Real> > c = jv.clone();
    this->value(*c,x,ctol);

    // Compute new step x + h*v.
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(h,v);
    this->update(*xnew);

    // Compute constraint at x + h*v.
    jv.zero();
    this->value(jv,*xnew,ctol);

    // Compute Newton quotient.
    jv.axpy(-1.0,*c);
    jv.scale(1.0/h);
  }
  else {
    Real h = 0.0;
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    Teuchos::RCP<Vector<Real> > e    = x.clone();
    Teuchos::RCP<Vector<Real> > J    = v.clone();
    Teuchos::RCP<Vector<Real> > c    = v.clone();
    this->value(*c,x,ctol);
    jv.zero();
    for ( unsigned i = 0; i < (unsigned)x.dimension(); i++ ) {
      e = x.basis(i);
      h = std::max(1.0,x.norm()/e->norm())*tol;
      xnew->set(x);
      xnew->axpy(h,*e);
      this->update(*xnew);
      this->value(*J,*xnew,ctol);
      J->axpy(-1.0,*c);
      J->scale(1.0/h);
      jv.axpy(J->dot(v),*e);
    }
  }
}

template <class Real>
void InequalityConstraint<Real>::applyHessian( Vector<Real> &huv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  Real jtol = std::sqrt(ROL_EPSILON);

  // Get step length.
  Real h = std::max(1.0,x.norm()/v.norm())*tol;
  //Real h = 2.0/(v.norm()*v.norm())*tol;

  // Compute constraint Jacobian at x.
  Teuchos::RCP<Vector<Real> > ju = huv.clone();
  bool adj = false;
  this->applyJacobian(*ju,u,x,adj,jtol);

  // Compute new step x + h*v.
  Teuchos::RCP<Vector<Real> > xnew = x.clone();
  xnew->set(x);
  xnew->axpy(h,v);
  this->update(*xnew);

  // Compute constraint Jacobian at x + h*v.
  huv.zero();
  this->applyJacobian(huv,u,*xnew,adj,jtol);

  // Compute Newton quotient.
  huv.axpy(-1.0,*ju);
  huv.scale(1.0/h);
}

} // namespace ROL

#endif
