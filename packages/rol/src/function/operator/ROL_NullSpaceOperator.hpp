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

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_Krylov.hpp"
#include "ROL_AugmentedSystemOperator.hpp"
#include "ROL_AugmentedSystemPrecOperator.hpp"

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
  const bool useInexact_;

  Ptr<LinearOperator<Real>> augsys_, augsysprec_;
  Ptr<Krylov<Real>> krylov_;
  mutable int iterKrylov_;
  mutable int flagKrylov_;

  mutable Ptr<Vector<Real>> v1_;
  mutable Ptr<Vector<Real>> v2_;
  mutable Ptr<PartitionedVector<Real>> vv_;
  mutable Ptr<Vector<Real>> b1_;
  mutable Ptr<Vector<Real>> b2_;
  mutable Ptr<PartitionedVector<Real>> bb_;
  mutable Ptr<Vector<Real>> w1_;
  mutable Ptr<Vector<Real>> w2_;
  mutable Ptr<PartitionedVector<Real>> ww_;
  mutable Ptr<Vector<Real>> mul_;

  int dim_;
  Real b1sqr_;

  void solveAugmentedSystem(Vector<Real> &v,
                            Vector<Real> &b,
                            Real &tol,
                            bool refine = false) const {
    if( refine ) {
      // TODO: Make sure this tol is actually ok...
      Real origTol = tol;
      ww_->set(v);
      augsys_->apply(*vv_, *ww_, tol);
      tol = origTol;
      b.axpy( static_cast<Real>(-1), *vv_ );
    }
    vv_->zero();
    // If inexact, change tolerance
    if( useInexact_ ) {
      krylov_->resetAbsoluteTolerance(tol);
    }

    flagKrylov_ = 0;
    tol = krylov_->run(*vv_,*augsys_,b,*augsysprec_,iterKrylov_,flagKrylov_);

    if( refine ) {
      v.plus(*vv_);
    }
    else {
      v.set(*vv_);
    }
  }

public:
  virtual ~NullSpaceOperator() {}
  NullSpaceOperator(const Ptr<Constraint<Real>> &con,
                    const Ptr<Vector<Real>>     &dom,
                    const Ptr<Vector<Real>>     &ran)
    : con_(con), useInexact_(false) {
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    dim_        = ran->dimension();
    if (dim_==1) {
      Real tol = std::sqrt(ROL_EPSILON<Real>());
      b1_ = dom->dual().clone();
      b2_ = ran->clone(); b2_->setScalar(1.0);
      con_->applyAdjointJacobian(*b1_,*b2_,*dom,tol);
      b1sqr_ = b1_->dot(*b1_);
    }
    else {
      ParameterList list;
      Real atol = static_cast<Real>(1e-12);
      Real rtol = static_cast<Real>(1e-2);
      list.sublist("General").sublist("Krylov").set("Type", "GMRES");
      list.sublist("General").sublist("Krylov").set("Absolute Tolerance", atol);
      list.sublist("General").sublist("Krylov").set("Relative Tolerance", rtol);
      list.sublist("General").sublist("Krylov").set("Iteration Limit", 200);
      krylov_ = KrylovFactory<Real>(list);

      augsys_ = makePtr<AugmentedSystemOperator<Real>>(con,dom);
      augsysprec_ = makePtr<AugmentedSystemPrecOperator<Real>>(con,dom);

      v1_ = dom->dual().clone();
      v2_ = ran->dual().clone();
      vv_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({v1_, v2_}));

      w1_ = dom->dual().clone();
      w2_ = ran->dual().clone();
      ww_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({w1_, w2_}));

      b1_ = dom->dual().clone();
      b2_ = ran->clone();
      bb_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({b1_, b2_}));

      mul_ = ran->dual().clone();
    }
  }

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    if (dim_==1) {
      Real dot = v.dot(*b1_);
      Hv.set(v);
      Hv.axpy(-dot/b1sqr_,*b1_);
    }
    else {
      b1_->set(v); b2_->zero();
      Ptr<PartitionedVector<Real>> sol = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({makePtrFromRef(Hv),mul_}));
      solveAugmentedSystem(*sol,*bb_,tol);
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
