// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUNDS_DEF_H
#define ROL_BOUNDS_DEF_H

#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class Bounds
    \brief Provides the elementwise interface to apply upper and lower bound
           constraints.

*/

namespace ROL {

template<typename Real>
Bounds<Real>::Bounds(const Vector<Real> &x, bool isLower, Real scale, Real feasTol)
  : scale_(scale), feasTol_(feasTol), mask_(x.clone()), min_diff_(ROL_INF<Real>()) {
  lower_ = x.clone();
  upper_ = x.clone();
  if (isLower) {
    lower_->set(x);
    upper_->applyUnary(Elementwise::Fill<Real>( BoundConstraint<Real>::computeInf(x)));
    BoundConstraint<Real>::activateLower();
  }
  else {
    lower_->applyUnary(Elementwise::Fill<Real>(-BoundConstraint<Real>::computeInf(x)));
    upper_->set(x);
    BoundConstraint<Real>::activateUpper();
  }
}

template<typename Real>
Bounds<Real>::Bounds(const Ptr<Vector<Real>> &x_lo, const Ptr<Vector<Real>> &x_up,
                     const Real scale, const Real feasTol)
  : scale_(scale), feasTol_(feasTol), mask_(x_lo->clone()) {
  lower_ = x_lo;
  upper_ = x_up;
  const Real half(0.5), one(1);
  // Compute difference between upper and lower bounds
  mask_->set(*upper_);
  mask_->axpy(-one,*lower_);
  // Compute minimum difference
  min_diff_ = mask_->reduce(minimum_);
  min_diff_ *= half;
}

template<typename Real>
void Bounds<Real>::project( Vector<Real> &x ) {
  struct Lesser : public Elementwise::BinaryFunction<Real> {
    Real apply(const Real &xc, const Real &yc) const { return xc<yc ? xc : yc; }
  } lesser;

  struct Greater : public Elementwise::BinaryFunction<Real> {
    Real apply(const Real &xc, const Real &yc) const { return xc>yc ? xc : yc; }
  } greater;

  if (BoundConstraint<Real>::isUpperActivated()) {
    x.applyBinary(lesser, *upper_); // Set x to the elementwise minimum of x and upper_
  }
  if (BoundConstraint<Real>::isLowerActivated()) {
    x.applyBinary(greater,*lower_); // Set x to the elementwise maximum of x and lower_
  }
}

template<typename Real>
void Bounds<Real>::projectInterior( Vector<Real> &x ) {
  // Make vector strictly feasible
  // Lower feasibility
  if (BoundConstraint<Real>::isLowerActivated()) {
    class LowerFeasible : public Elementwise::BinaryFunction<Real> {
    private:
      const Real eps_;
      const Real diff_;
    public:
      LowerFeasible(const Real eps, const Real diff)
        : eps_(eps), diff_(diff) {}
      Real apply( const Real &xc, const Real &yc ) const {
        const Real tol = static_cast<Real>(100)*ROL_EPSILON<Real>();
        const Real one(1);
        Real val = ((yc <-tol) ? yc*(one-eps_)
                 : ((yc > tol) ? yc*(one+eps_)
                 : yc+eps_));
        val = std::min(yc+eps_*diff_, val);
        return xc < val ? val : xc;
      }
    };
    x.applyBinary(LowerFeasible(feasTol_,min_diff_), *lower_);
  }
  // Upper feasibility
  if (BoundConstraint<Real>::isUpperActivated()) {
    class UpperFeasible : public Elementwise::BinaryFunction<Real> {
    private:
      const Real eps_;
      const Real diff_;
    public:
      UpperFeasible(const Real eps, const Real diff)
        : eps_(eps), diff_(diff) {}
      Real apply( const Real &xc, const Real &yc ) const {
        const Real tol = static_cast<Real>(100)*ROL_EPSILON<Real>();
        const Real one(1);
        Real val = ((yc <-tol) ? yc*(one+eps_)
                 : ((yc > tol) ? yc*(one-eps_)
                 : yc-eps_));
        val = std::max(yc-eps_*diff_, val);
        return xc > val ? val : xc;
      }
    };
    x.applyBinary(UpperFeasible(feasTol_,min_diff_), *upper_);
  }
}

template<typename Real>
void Bounds<Real>::pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
  if (BoundConstraint<Real>::isUpperActivated()) {
    Real one(1), epsn(std::min(scale_*eps,static_cast<Real>(0.1)*min_diff_));

    mask_->set(*upper_);
    mask_->axpy(-one,x);

    Active op(epsn);
    v.applyBinary(op,*mask_);
  }
}

template<typename Real>
void Bounds<Real>::pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) {
  if (BoundConstraint<Real>::isUpperActivated()) {
    Real one(1), epsn(std::min(scale_*xeps,static_cast<Real>(0.1)*min_diff_));

    mask_->set(*upper_);
    mask_->axpy(-one,x);

    UpperBinding op(epsn,geps);
    mask_->applyBinary(op,g);

    v.applyBinary(prune_,*mask_);
  }
}

template<typename Real>
void Bounds<Real>::pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
  if (BoundConstraint<Real>::isLowerActivated()) {
    Real one(1), epsn(std::min(scale_*eps,static_cast<Real>(0.1)*min_diff_));

    mask_->set(x);
    mask_->axpy(-one,*lower_);

    Active op(epsn);
    v.applyBinary(op,*mask_);
  }
}

template<typename Real>
void Bounds<Real>::pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps ) {
  if (BoundConstraint<Real>::isLowerActivated()) {
    Real one(1), epsn(std::min(scale_*xeps,static_cast<Real>(0.1)*min_diff_));

    mask_->set(x);
    mask_->axpy(-one,*lower_);

    LowerBinding op(epsn,geps);
    mask_->applyBinary(op,g);

    v.applyBinary(prune_,*mask_);
  }
}

template<typename Real>
bool Bounds<Real>::isFeasible( const Vector<Real> &v ) {
  const Real half(0.5);
  bool flagU = false, flagL = false;
  if (BoundConstraint<Real>::isUpperActivated()) {
    mask_->set(v);
    mask_->applyBinary(isGreater_,*upper_);
    flagU = mask_->reduce(maximum_) > half ? true : false;
  }
  if (BoundConstraint<Real>::isLowerActivated()) {
    mask_->set(*lower_);
    mask_->applyBinary(isGreater_,v);
    flagL = mask_->reduce(maximum_) > half ? true : false;
  }
  return ((flagU || flagL) ? false : true);
}

template<typename Real>
void Bounds<Real>::buildScalingFunction(Vector<Real> &d, const Vector<Real> &x, const Vector<Real> &g) const {
  // TODO: Cache values?

  const Real zero(0), one(1);

  // This implementation handles -l and/or u = infinity properly.

  /*
  The first if statement defining d_II (Equation 4.4 of [Ulbrich et al. 1999])
  is equivalent to x - l <= min(u - x, -g) = - max(g, x - u).
    * Our x, l, u represent u, a, b in the paper.
    * Indeed, if x - l = -g < u - x, then min{|g|,c} = min{x - l, u - x, c}.
  */

  d.set(x);
  d.axpy(-one,*upper_);
  d.applyBinary(Elementwise::Max<Real>(),g);
  d.plus(x);
  d.axpy(-one,*lower_);  // = x - l + max(g, x - u)

  mask_->set(x);
  mask_->axpy(-one,*lower_);
  mask_->applyBinary(Elementwise::Min<Real>(),g);
  mask_->plus(x);
  mask_->axpy(-one,*upper_);
  mask_->scale(-one);   // = u - x - min(g, x - l)

  mask_->applyBinary(Elementwise::Min<Real>(),d);

  d.setScalar(-one);
  d.applyBinary(Active(zero),*mask_);
  mask_->setScalar(one);
  d.plus(*mask_);
  // d[i] =   1    when one of the if conditions in (4.4) are true else 0

  mask_->set(g);
  mask_->applyUnary(Elementwise::AbsoluteValue<Real>());
  d.applyBinary(Elementwise::Multiply<Real>(),*mask_);
  // d[i] = |g[i]| when one of the if conditions in (4.4) are true else 0

  // Compute min(x - l, u - x).
  // * We use the identity min(p, q) = min(p + r, q + r) - r to handle the case
  //   where -l or u = infinity.
  mask_->set(x);
  mask_->axpy(-one,*lower_);
  mask_->plus(x);
  mask_->applyBinary(Elementwise::Min<Real>(),*upper_);
  mask_->axpy(-one,x);  // = min(x - l, u - x)

  // When one of the if conditions in (4.4) are true |g|[i] >= (*mask_)[i].
  // Thus by taking
  d.applyBinary(Elementwise::Max<Real>(),*mask_);
  // d_II follows as min(d, c), i.e.,
  mask_->set(*upper_);
  mask_->axpy(-one,*lower_);
  mask_->applyUnary(buildC_);
  d.applyBinary(Elementwise::Min<Real>(),*mask_);
}

template<typename Real>
void Bounds<Real>::applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
  Bounds<Real>::buildScalingFunction(dv, x, g);
  dv.applyBinary(Elementwise::DivideAndInvert<Real>(),v);
}

template<typename Real>
void Bounds<Real>::applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
  const Real one(1), two(2), three(3);

  // This implementation builds Equation (5.1) of [Ulbrich et al. 1999].

  Bounds<Real>::buildScalingFunction(dv, x, g);  // dv = d_II

  mask_->set(*upper_);
  mask_->axpy(-one,*lower_);
  mask_->applyUnary(buildC_);        // = c
  mask_->axpy(-one,dv);              // = c - d_II
  dv.setScalar(one);
  dv.applyBinary(Active(0),*mask_);  // = \chi

  mask_->setScalar(three);
  dv.applyBinary(Elementwise::Multiply<Real>(),*mask_);
  dv.axpy(-one,*mask_);
  // dv[i] = 0 if \chi[i] = 1 else -3

  mask_->set(g);
  mask_->applyUnary(Elementwise::Sign<Real>());
  dv.plus(*mask_);  // dv[i] = sgn(g[i]) if \chi[i] = 1 else dv[i] <= -2

  // Set the dv elements that = 0 to sgn(u + l - 2x).
  mask_->set(*upper_);
  mask_->plus(*lower_);
  mask_->axpy(-two,x);
  mask_->applyUnary(Elementwise::Sign<Real>());
  dv.applyBinary(setZeroEntry_,*mask_);

  // Set the dv elements that = 0 to 1.
  mask_->setScalar(one);
  dv.applyBinary(setZeroEntry_,*mask_);

  // Set the dv elements that are <= -2 to 0.
  mask_->set(dv);
  dv.applyBinary(Active(-two),*mask_);

  dv.applyBinary(Elementwise::Multiply<Real>(), g);
  dv.applyBinary(Elementwise::Multiply<Real>(), v);
}

} // namespace ROL

#endif
