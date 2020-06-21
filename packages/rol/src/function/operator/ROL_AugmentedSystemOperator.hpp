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

#ifndef ROL_AUGMENTED_SYSTEM_OPERATOR_H
#define ROL_AUGMENTED_SYSTEM_OPERATOR_H

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::AugmentedSystemOperator
    \brief Apply the augmented system operator.

    ---
*/

namespace ROL {

template <class Real>
class AugmentedSystemOperator : public LinearOperator<Real> {
private:
  const Ptr<Constraint<Real>> con_;
  const Ptr<Vector<Real>>     x_;
  const Real                  delta_;

public:
  virtual ~AugmentedSystemOperator() {}
  AugmentedSystemOperator(const Ptr<Constraint<Real>> &con,
                          const Ptr<Vector<Real>>     &x,
                          const Real                   delta = -1.0)
    : con_(con), x_(x), delta_(delta) {}

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    PartitionedVector<Real>      &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
    const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

    con_->applyAdjointJacobian(*(Hvp.get(0)), *(vp.get(1)), *x_, tol);
    Hvp.get(0)->plus(*(vp.get(0)));

    con_->applyJacobian(*(Hvp.get(1)), *(vp.get(0)), *x_, tol);
    if ( delta_ > static_cast<Real>(0) ) {
      Hvp.get(1)->axpy(-delta_*delta_, *(vp.get(1)));
    }
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    apply(Hv,v,tol);
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemOperator::applyInverse : Not implemented!");
  }

  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemOperator::applyAdjointInverse : Not implemented!");
  }

}; // class AugmentedSystemOperator

} // namespace ROL

#endif
