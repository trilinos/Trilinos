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

#ifndef ROL_AUGMENTED_SYSTEM_PREC_OPERATOR_H
#define ROL_AUGMENTED_SYSTEM_PREC_OPERATOR_H

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::AugmentedSystemPrecOperator
    \brief Implements a preconditioner for the augmented system.

    ---
*/

namespace ROL {

template <class Real>
class AugmentedSystemPrecOperator : public LinearOperator<Real> {
private:
  const Ptr<Constraint<Real>>   con_;
  const Ptr<const Vector<Real>> x_;

public:
  virtual ~AugmentedSystemPrecOperator() {}
  AugmentedSystemPrecOperator(const Ptr<Constraint<Real>>   &con,
                              const Ptr<const Vector<Real>> &x)
    : con_(con), x_(x) {}

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemPrecOperator::apply : Not Implemented!");
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemPrecOperator::applyAdjoint : Not Implemented!");
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Real zero(0);
    PartitionedVector<Real>      &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
    const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

    Hvp.set(0, *(vp.get(0)));
    // Second x should be dual, but unused?
    con_->applyPreconditioner(*(Hvp.get(1)),*(vp.get(1)),*x_,*x_, zero);
  }

  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    applyInverse(Hv,v,tol);
  }

}; // class AugmentedSystemPrecOperator

} // namespace ROL

#endif
