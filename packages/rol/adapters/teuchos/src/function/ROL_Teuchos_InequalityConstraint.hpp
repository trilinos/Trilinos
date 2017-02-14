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

#ifndef ROL_TEUCHOS_INEQUALITYCONSTRAINT_HPP 
#define ROL_TEUCHOS_INEQUALITYCONSTRAINT_HPP

#include "ROL_TeuchosEqualityConstraint.hpp"

/**  @ingroup func_group
     \class ROL::TeuchosInequalityConstraint 
     \brief Provides a unique argument for inequality constraints using
            Teuchos::SerialDenseVector types, which otherwise behave exactly as equality constraints
*/

namespace ROL {

template<class Ordinal, class Real> 
class TeuchosInequalityConstraint : public virtual TeuchosEqualityConstraint<Ordinal, Real>,
                                    public virtual InequalityConstraint<Real> {

  typedef TeuchosEqualityConstraint<Real>  TEC;
  typedef Vector<Real>                     V;  

public:

  using EqualityConstraint<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    TEC::update(x,flag,iter);
  }

  using EqualityConstraint<Real>::value;
  void value(V &c, const V &x, Real &tol ) {
    TEC::value(c,x,tol);
  }

  using EqualityConstraint<Real>::applyJacobian;
  void applyJacobian(V &jv, const V &v, const V &x, Real &tol) {
    TEC::applyJacobian(jv, v, x, tol);
  }

  using EqualityConstraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(V &aju, const V &u, const V &x, Real &tol) {
    TEC::applyAdjointJacobian(aju, u, x, tol);
  }

  using EqualityConstraint<Real>::applyAdjointHessian;
  void applyAdjointHessian(V &ahuv, const V &u, const V &v, const V &x, Real &tol) {
    TEC::applyAdjointHessian(ahuv, u, v, x, tol);  

  }

  using EqualityConstraint<Real>::solveAugmentedSystem;
  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                         const Vector<Real> &b1, const Vector<Real> &b2,
                                         const Vector<Real> &x, Real &tol) {
    return TEC::solveAugmentedSystem(v1,v2,b1,b2,x,tol);
  }

  using EqualityConstraint<Real>::applyPreconditioner;
  void applyPreconditioner(Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &x,
                           const Vector<Real> &g, Real &tol) {
    TEC::applyPreconditioner(pv,v,x,g,tol);
  }
}; // class TeuchosInequalityConstraint

} // namespace ROL

#endif // ROL_TEUCHOS_INEQUALITYCONSTRAINT_HPP
