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

#ifndef ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_H
#define ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::ScalarLinearConstraint
    \brief This equality constraint defines an affine hyperplane.

    ROL's scalar linear equality constraint interface implements
    \f[
       c(x) := \langle a, x\rangle_{\mathcal{X}^*,\mathcal{X}} - b = 0
    \f]
    where \f$a\in\mathcal{X}^*\f$ and \f$b\in\mathbb{R}\f$.  The range space of
    \f$c\f$ is an ROL::StdVector with dimension 1.

    Note: If \f$a\neq 0\f$ then there exists an explicit solution of the
    augmented system.  Namely,
    \f[
       v_1 = I^{-1}(b_1-av_2)
         \quad\text{and}\quad
       v_2 = \frac{(\langle a,I^{-1}b_1\rangle_{\mathcal{X}^*,\mathcal{X}}
               - b_2)}{\|a\|_{\mathcal{X}^*}^2}\,.
    \f]
    Moreover, note that \f$I^{-1}v\f$ for any \f$v\in\mathcal{X}^*\f$ is
    implemented in ROL as v.dual(). 

    ---
*/


namespace ROL {

template <class Real>
class ScalarLinearConstraint : public Constraint<Real> {
private:
  const Teuchos::RCP<Vector<Real> > a_; ///< Dual vector defining hyperplane
  const Real b_;                        ///< Affine shift

public:
  ScalarLinearConstraint(const Teuchos::RCP<Vector<Real> > &a,
                                     const Real b)
    : a_(a), b_(b) {}

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    StdVector<Real>    cc = Teuchos::dyn_cast<StdVector<Real> >(c);
    std::vector<Real> &cp = *(cc.getVector());
    cp[0] = a_->dot(x.dual()) - b_;
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
               const Vector<Real> &x,  Real &tol) {
    StdVector<Real>    jc = Teuchos::dyn_cast<StdVector<Real> >(jv);
    std::vector<Real> &jp = *(jc.getVector());
    jp[0] = a_->dot(v.dual());
  }

  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                      const Vector<Real> &x,   Real &tol) {
    const StdVector<Real>    vc = Teuchos::dyn_cast<const StdVector<Real> >(v);
    const std::vector<Real> &vp = *(vc.getVector());
    ajv.set(*a_);
    ajv.scale(vp[0]);
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                     const Vector<Real> &v,    const Vector<Real> &x,
                           Real &tol) {
    ahuv.zero();
  }

  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                   const Vector<Real> &b1, const Vector<Real> &b2,
                                   const Vector<Real> &x,  Real &tol) {
    StdVector<Real>    v2c = Teuchos::dyn_cast<StdVector<Real> >(v2);
    std::vector<Real> &v2p = *(v2c.getVector());
    const StdVector<Real>    b2c = Teuchos::dyn_cast<const StdVector<Real> >(b2);
    const std::vector<Real> &b2p = *(b2c.getVector());

    v2p[0] = (a_->dot(b1.dual()) - b2p[0])/a_->dot(*a_);
    v1.set(b1.dual());
    v1.axpy(-v2p[0],a_->dual());

    std::vector<Real> out;
    return out;
  }

}; // class ScalarLinearConstraint

} // namespace ROL

#endif
