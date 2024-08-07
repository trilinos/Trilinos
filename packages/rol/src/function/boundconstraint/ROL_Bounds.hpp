// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

template<typename Real>
class Bounds : public BoundConstraint<Real> {
private:
  const Real scale_;
  const Real feasTol_;

  using BoundConstraint<Real>::lower_;
  using BoundConstraint<Real>::upper_;

  Ptr<Vector<Real>> mask_;

  Real min_diff_;

  Elementwise::ReductionMin<Real> minimum_;
  Elementwise::ReductionMax<Real> maximum_;

  class isGreater : public Elementwise::BinaryFunction<Real> {
  public:
    isGreater() {}
    Real apply(const Real &x, const Real &y) const {
      return (x > y) ? static_cast<Real>(1) : static_cast<Real>(0);
    }
  } isGreater_;

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
    UpperBinding(Real xeps, Real geps) : xeps_(xeps), geps_(geps) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y < -geps_ && x <= xeps_) ? 0 : 1);
    }
    private:
    Real xeps_, geps_;
  };

  class LowerBinding : public Elementwise::BinaryFunction<Real> {
    public:
    LowerBinding(Real xeps, Real geps) : xeps_(xeps), geps_(geps) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y > geps_ && x <= xeps_) ? 0 : 1);
    }
    private:
    Real xeps_, geps_;
  };

  class PruneBinding : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return ((y == 1) ? x : 0);
      }
  } prune_;

  class BuildC : public Elementwise::UnaryFunction<Real> {
    public:
      Real apply( const Real &x ) const {
        const Real zeta(0.5), kappa(1);
        return std::min(zeta * x, kappa);
      }
  } buildC_;

  class SetZeroEntry : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply(const Real &x, const Real &y) const {
        const Real zero(0);
        return (x==zero ? y : x);
      }
  } setZeroEntry_;

  void buildScalingFunction(Vector<Real> &d, const Vector<Real> &x, const Vector<Real> &g) const;

public:

  Bounds(const Vector<Real> &x,
         bool isLower = true,
         Real scale = 1,
         Real feasTol = std::sqrt(ROL_EPSILON<Real>()));

  Bounds(const Ptr<Vector<Real>> &x_lo,
         const Ptr<Vector<Real>> &x_up,
         const Real scale = 1,
         const Real feasTol = std::sqrt(ROL_EPSILON<Real>()));

  void project( Vector<Real> &x ) override;

  void projectInterior( Vector<Real> &x ) override;

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) override;

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) override;

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) override;

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) override;

  bool isFeasible( const Vector<Real> &v ) override;

  void applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const override;

  void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const override;
}; // class Bounds

} // namespace ROL

#include "ROL_Bounds_Def.hpp"

#endif
