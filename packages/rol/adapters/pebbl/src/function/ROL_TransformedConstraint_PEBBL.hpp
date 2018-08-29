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

#ifndef ROL_TRANSFORMEDCONSTRAINT_PEBBL_H
#define ROL_TRANSFORMEDCONSTRAINT_PEBBL_H

#include "ROL_Transform_PEBBL.hpp"

/** @ingroup func_group
    \class ROL::TransformedConstraint_PEBBL
    \brief Defines the pebbl transformed constraint operator interface.

    ROL's pebbl transformed constraint interface is designed to set individual
    components of a vector to a fixed value.

    ---
*/


namespace ROL {

template <class Real>
class TransformedConstraint_PEBBL : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>>      con_;
  const Ptr<Transform_PEBBL<Real>> trans_;
  Ptr<Vector<Real>>                Tx_;
  Ptr<Vector<Real>>                Tv_;
  Ptr<Vector<Real>>                dx_;
  bool                             isInit_;
  bool                             isTransformed_;

  void initialize(const Vector<Real> &x) {
    if (!isInit_) {
      Tx_ = x.clone();
      Tv_ = x.clone();
      dx_ = x.dual().clone();
      isInit_ = true;
    }
  }

  void transform(const Vector<Real> &x, Real &tol) {
    if (!isTransformed_) {
      trans_->value(*Tx_,x,tol);
      isTransformed_ = true;
    }
  }

public:
  TransformedConstraint_PEBBL(const Ptr<Constraint<Real>>      &con,
                              const Ptr<Transform_PEBBL<Real>> &trans)
    : con_(con), trans_(trans), isInit_(false), isTransformed_(false) {}

  TransformedConstraint_PEBBL(const TransformedConstraint_PEBBL &con)
    : con_(con.con_), trans_(con.trans_), isInit_(false), isTransformed_(false) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    isTransformed_ = false;
    initialize(x);
    Tx_->zero(); Tv_->zero(); dx_->zero();
    transform(x,tol);
    con_->update(*Tx_,flag,iter);
  }

  void value(Vector<Real> &c,
       const Vector<Real> &x,
             Real &tol) {
    initialize(x);
    transform(x,tol);
    con_->value(c,*Tx_,tol);    
  }

  void applyJacobian(Vector<Real> &jv,
               const Vector<Real> &v,
               const Vector<Real> &x,
                     Real &tol) {
    initialize(x);
    transform(x,tol);
    trans_->applyJacobian(*Tv_,v,x,tol);
    con_->applyJacobian(jv,*Tv_,*Tx_,tol);
  }

  void applyAdjointJacobian(Vector<Real> &ajv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                            Real &tol) {
    initialize(x);
    transform(x,tol);
    con_->applyAdjointJacobian(*dx_,v,*Tx_,tol);
    trans_->applyAdjointJacobian(ajv,*dx_,x,tol);
  }

  void applyAdjointHessian(Vector<Real> &ahuv,
                     const Vector<Real> &u,
                     const Vector<Real> &v,
                     const Vector<Real> &x,
                           Real &tol) {
    initialize(x);
    transform(x,tol);
    trans_->applyJacobian(*Tv_,v,x,tol);
    con_->applyAdjointHessian(*dx_,u,*Tv_,*Tx_,tol);
    trans_->applyAdjointJacobian(ahuv,*dx_,x,tol);
  }

  bool isEmpty(void) const {
    return trans_->isEmpty();
  }

  void reset(void) {
    trans_->clear();
  }

  void add(const std::map<int,Real> &in) {
    trans_->add(in);
  }

  void add(const std::pair<int,Real> &in) {
    trans_->add(in);
  }

}; // class TransformedConstraint_PEBBL

} // namespace ROL

#endif
