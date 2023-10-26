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

#ifndef ROL_L1OBJECTIVE_H
#define ROL_L1OBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::l1Objective
    \brief Provides the interface to evaluate the weighted/shifted l1 objective function.

        ---
*/


namespace ROL {

template<typename Real>
class l1Objective : public Objective<Real> {
private:
  const Ptr<Vector<Real>> weights_, shift_;
  Ptr<Vector<Real>> tmp_;

  struct ProjSymBnd : public Elementwise::BinaryFunction<Real> {
       Real apply(const Real &xc, const Real &yc) const { return std::min(yc, std::max(-yc, xc)); }
  } psb_;
 
public:

  l1Objective(const Ptr<Vector<Real>> &weights)
    : weights_(weights), shift_(weights->dual().clone()) {
    shift_->zero();
    tmp_ = shift_->clone();
  }

  l1Objective(const Ptr<Vector<Real>> &weights, const Ptr<Vector<Real>> &shift)
    : weights_(weights),  shift_(shift) {
    tmp_ = shift_->clone();
  }
  
  Real value( const Vector<Real> &x, Real &tol ) {
    tmp_->set(x);
    tmp_->axpy(static_cast<Real>(-1),*shift_);
    tmp_->applyUnary(Elementwise::AbsoluteValue<Real>());
    return weights_->apply(*tmp_);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.set(x);
    g.axpy(static_cast<Real>(-1),*shift_);
    g.applyUnary(Elementwise::Sign<Real>());
    g.applyBinary(Elementwise::Multiply<Real>(), *weights_);
  }

  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
    gradient(*tmp_, x, tol);
    return tmp_->apply(d);
  }

  void prox( Vector<Real> &Pv, const Vector<Real> &v, Real t, Real &tol){
    Pv.set(*shift_);
    Pv.axpy(static_cast<Real>(-1), v);
    Pv.scale(static_cast<Real>(1) / t);
    Pv.applyBinary(psb_, *weights_);
    Pv.scale(t);
    Pv.plus(v);
  } 
}; // class l1Objective

} // namespace ROL

#endif
