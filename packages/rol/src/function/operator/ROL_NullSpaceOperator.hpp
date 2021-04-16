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

#ifndef ROL_NULL_SPACE_OPERATOR_H
#define ROL_NULL_SPACE_OPERATOR_H

#include "ROL_LinearOperator.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::NullSpaceOperator
    \brief Projects on to the null space of a linear constraint.

    ---
*/

namespace ROL {

template <class Real>
class NullSpaceOperator : public LinearOperator<Real> {
private:
  const Ptr<Constraint<Real>> con_;
  const Ptr<Vector<Real>> x_;
  const bool useAugSys_;

  mutable Ptr<Vector<Real>> b1_;
  mutable Ptr<Vector<Real>> b1dual_;
  mutable Ptr<Vector<Real>> b2_;
  mutable Ptr<Vector<Real>> mul_;

  int dim_;
  Real b1sqr_;

public:
  virtual ~NullSpaceOperator() {}
  NullSpaceOperator(const Ptr<Constraint<Real>> &con,
                    const Vector<Real> &dom,
                    const Vector<Real> &ran,
                    const bool useAugSys = false)
    : con_(con), x_(dom.clone()), useAugSys_(useAugSys) {
    x_->set(dom);
    dim_ = ran.dimension();
    if (dim_==1 && !useAugSys_) {
      Real tol = std::sqrt(ROL_EPSILON<Real>());
      b1_     = dom.dual().clone();
      b1dual_ = dom.clone();
      b2_     = ran.dual().clone(); b2_->setScalar(1.0);
      con_->applyAdjointJacobian(*b1_,*b2_,dom,tol);
      b1dual_->set(b1_->dual());
      b1sqr_ = b1_->dot(*b1_);
    }
    else {
      b1_  = dom.dual().clone();
      b2_  = ran.clone();
      mul_ = ran.dual().clone();
    }
  }

  NullSpaceOperator(const Ptr<Constraint<Real>>   &con,
                    const Ptr<const Vector<Real>> &dom,
                    const Ptr<const Vector<Real>> &ran)
    : NullSpaceOperator(con,*dom,*ran) {}

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    x_->set(x);
    if (dim_==1 && !useAugSys_) {
      Real tol = std::sqrt(ROL_EPSILON<Real>());
      con_->applyAdjointJacobian(*b1_,*b2_,x,tol);
      b1dual_->set(b1_->dual());
      b1sqr_ = b1_->dot(*b1_);
    }
  }

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    if (dim_==1 && !useAugSys_) {
      Real dot = v.dot(*b1dual_);
      Hv.set(v);
      Hv.axpy(-dot/b1sqr_,*b1dual_);
    }
    else {
      b1_->set(v.dual()); b2_->zero();
      con_->solveAugmentedSystem(Hv,*mul_,*b1_,*b2_,*x_,tol);
    }
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    apply(Hv,v,tol);
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> NullSpaceOperator::applyInverse : Not Implemented!");
  }

  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> NullSpaceOperator::applyAdjointInverse : Not Implemented!");
  }

}; // class NullSpaceOperator

} // namespace ROL

#endif
