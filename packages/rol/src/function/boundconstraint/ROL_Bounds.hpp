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

#ifndef ROL_BOUNDS_H
#define ROL_BOUNDS_H

#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::Bounds
    \brief Provides the elementwise interface to apply upper and lower bound
           constraints.

*/

namespace ROL {

template <class Real>
class Bounds : public BoundConstraint<Real> {
private:
  const ROL::Ptr<Vector<Real> > x_lo_;
  const ROL::Ptr<Vector<Real> > x_up_;
  const Real scale_;

  ROL::Ptr<Vector<Real> > mask_;

  Real min_diff_;

  Elementwise::ReductionMin<Real> minimum_;

  class Active : public Elementwise::BinaryFunction<Real> {
    public:
    Active(Real offset) : offset_(offset) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y <= offset_) ? 0 : x);
    }
    private:
    Real offset_;
  };

  class UpperBinding : public Elementwise::BinaryFunction<Real> {
    public:
    UpperBinding(Real offset) : offset_(offset) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y < 0 && x <= offset_) ? 0 : 1);
    }
    private:
    Real offset_;
  };

  class LowerBinding : public Elementwise::BinaryFunction<Real> {
    public:
    LowerBinding(Real offset) : offset_(offset) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y > 0 && x <= offset_) ? 0 : 1);
    }
    private:
    Real offset_;
  };

  class PruneBinding : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return ((y == 1) ? x : 0);
      }
  } prune_;

public:

  Bounds(const Vector<Real> &x,
         bool isLower = true,
         Real scale = 1)
    : x_lo_(x.clone()),  x_up_(x.clone()),
      scale_(scale), mask_(x.clone()),
      min_diff_(ROL_INF<Real>()) {
    if (isLower) {
      x_lo_->set(x);
      x_up_->applyUnary(Elementwise::Fill<Real>(ROL_INF<Real>()));
      BoundConstraint<Real>::activateLower();
    }
    else {
      x_lo_->applyUnary(Elementwise::Fill<Real>(ROL_NINF<Real>()));
      x_up_->set(x);
      BoundConstraint<Real>::activateUpper();
    }
  }

  Bounds(const ROL::Ptr<Vector<Real> > &x_lo,
         const ROL::Ptr<Vector<Real> > &x_up,
         const Real scale = 1)
    : x_lo_(x_lo), x_up_(x_up), scale_(scale), mask_(x_lo->clone()) {
    const Real half(0.5), one(1);
    // Compute difference between upper and lower bounds
    mask_->set(*x_up_);
    mask_->axpy(-one,*x_lo_);
    // Compute minimum difference
    min_diff_ = mask_->reduce(minimum_);
    min_diff_ *= half;
  }

  virtual void project( Vector<Real> &x ) {
    struct Lesser : public Elementwise::BinaryFunction<Real> {
      Real apply(const Real &xc, const Real &yc) const { return xc<yc ? xc : yc; }
    } lesser;

    struct Greater : public Elementwise::BinaryFunction<Real> {
      Real apply(const Real &xc, const Real &yc) const { return xc>yc ? xc : yc; }
    } greater;

    if (BoundConstraint<Real>::isUpperActivated()) {
      x.applyBinary(lesser, *x_up_); // Set x to the elementwise minimum of x and x_up_
    }
    if (BoundConstraint<Real>::isLowerActivated()) {
      x.applyBinary(greater,*x_lo_); // Set x to the elementwise maximum of x and x_lo_
    }
  }

  virtual void projectInterior( Vector<Real> &x ) {
    // Make vector strictly feasible
    const Real eps(1e-1);
    // Lower feasibility
    if (BoundConstraint<Real>::isLowerActivated()) {
      class LowerFeasible : public Elementwise::BinaryFunction<Real> {
      private:
        const Real eps_;
        const Real diff_;
      public:
        LowerFeasible(const Real eps, const Real diff)
          : eps_(eps), diff_(diff) {}
        Real apply( const Real &x, const Real &y ) const {
          const Real tol = static_cast<Real>(100)*ROL_EPSILON<Real>();
          const Real one(1);
          Real val = ((y <-tol) ? y*(one-eps_)
                   : ((y > tol) ? y*(one+eps_)
                   : y+eps_));
          val = std::min(y+eps_*diff_, val);
          return (x < y+tol) ? val : x;
        }
      };
      x.applyBinary(LowerFeasible(eps,min_diff_), *x_lo_);
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
        Real apply( const Real &x, const Real &y ) const {
          const Real tol = static_cast<Real>(100)*ROL_EPSILON<Real>();
          const Real one(1);
          Real val = ((y <-tol) ? y*(one+eps_)
                   : ((y > tol) ? y*(one-eps_)
                   : y-eps_));
          val = std::max(y-eps_*diff_, val);
          return (x > y-tol) ? val : x;
        }
      };
      x.applyBinary(UpperFeasible(eps,min_diff_), *x_up_);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if (BoundConstraint<Real>::isUpperActivated()) {
      Real one(1), epsn(std::min(scale_*eps,min_diff_));

      mask_->set(*x_up_);
      mask_->axpy(-one,x);

      Active op(epsn);
      v.applyBinary(op,*mask_);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if (BoundConstraint<Real>::isUpperActivated()) {
      Real one(1), epsn(std::min(scale_*eps,min_diff_));

      mask_->set(*x_up_);
      mask_->axpy(-one,x);

      UpperBinding op(epsn);
      mask_->applyBinary(op,g);

      v.applyBinary(prune_,*mask_);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if (BoundConstraint<Real>::isLowerActivated()) {
      Real one(1), epsn(std::min(scale_*eps,min_diff_));

      mask_->set(x);
      mask_->axpy(-one,*x_lo_);

      Active op(epsn);
      v.applyBinary(op,*mask_);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if (BoundConstraint<Real>::isLowerActivated()) {
      Real one(1), epsn(std::min(scale_*eps,min_diff_));

      mask_->set(x);
      mask_->axpy(-one,*x_lo_);

      LowerBinding op(epsn);
      mask_->applyBinary(op,g);

      v.applyBinary(prune_,*mask_);
    }
  }

  const ROL::Ptr<const Vector<Real> > getLowerBound( void ) const {
    return x_lo_;
  }

  const ROL::Ptr<const Vector<Real> > getUpperBound( void ) const {
    return x_up_;
  }

  virtual bool isFeasible( const Vector<Real> &v ) { 
    const Real one(1);
    bool flagU = false, flagL = false;
    if (BoundConstraint<Real>::isUpperActivated()) {
      mask_->set(*x_up_);
      mask_->axpy(-one,v);
      Real uminusv = mask_->reduce(minimum_);
      flagU = ((uminusv<0) ? true : false);
    }
    if (BoundConstraint<Real>::isLowerActivated()) {
      mask_->set(v);
      mask_->axpy(-one,*x_lo_);
      Real vminusl = mask_->reduce(minimum_);

      flagL = ((vminusl<0) ? true : false);
    }
    return ((flagU || flagL) ? false : true);
  }

}; // class Bounds

} // namespace ROL

#endif
