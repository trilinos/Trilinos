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

#include "ROL_Constraint.hpp"
#include "ROL_LinearOperator.hpp"


#ifndef ROL_LINEAR_CONSTRAINT_H
#define ROL_LINEAR_CONSTRAINT_H


/** @ingroup func_group
    \class ROL::LinearConstraint
    \brief Provides the interface to evaluate linear constraints.

    This class implements the linear constraint
    \f[
       c(x) = Ax-b
    \f]

    Where A is a linear operator

    ---
*/

namespace ROL {

template <class Real>
class LinearConstraint : public Constraint<Real> {
private:
  const ROL::Ptr<const LinearOperator<Real> > A_;
  const ROL::Ptr<const LinearOperator<Real> > Atrans_;
  const ROL::Ptr<const Vector<Real> > b_;
  bool  isSymmetric_;
public:
  // Nonsymmetric case
  LinearConstraint( const ROL::Ptr<const LinearOperator<Real> > &A,
                            const ROL::Ptr<const LinearOperator<Real> > &Atrans,
                            const ROL::Ptr<const Vector<Real> &b ) :
      A_(A), Atrans_(Atrans), b_(b), isSymmetric_(false) {
  }
  // Symmetric case
  LinearConstraint( const ROL::Ptr<const LinearOperator<Real> > &A,
                            const ROL::Ptr<const Vector<Real> &b ) : 
      A_(A), Atrans_(A), b_(b), isSymmetric_(true) {
  }

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    A_->apply(c,x,tol);
    c_->axpy(-1.0,*b_);
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, 
                     const Vector<Real> &x, Real &tol) {
    A_->apply(jv,v,tol);
  }
 
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                            const Vector<Real> &x, Real &tol) {
    Atrans_->apply(ajv,v,tol);
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                           const Vector<Real> &x, Real &tol) {
    ahuv.zero();
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    A_->update(x,flag,iter);
    if( !isSymmetric_ ) {
      A_->update(x,flag,iter);
    }
  }
}; // class LinearConstraint


} // namespace ROL

#endif //ROL_LINEAR_EQUALITY_CONSTRAINT_H
