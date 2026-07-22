// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DOUGLASRACHFORDPROJECTION_H
#define ROL_DOUGLASRACHFORDPROJECTION_H

#include "ROL_PolyhedralProjection.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<typename Real>
class DouglasRachfordProjection : public PolyhedralProjection<Real> {
private:
  int dim_;
  Ptr<Vector<Real>> tmp_, p_, q_, y_, z_;
  Real b_, cdot_;

  Real DEFAULT_atol_, DEFAULT_rtol_;
  int DEFAULT_maxit_, DEFAULT_verbosity_;
  Real DEFAULT_alpha1_, DEFAULT_gamma_, DEFAULT_t0_;

  Real atol_, rtol_;
  int maxit_, verbosity_;
  Real alpha1_, alpha2_, gamma_, t0_;

  using PolyhedralProjection<Real>::bnd_;
  using PolyhedralProjection<Real>::con_;
  using PolyhedralProjection<Real>::xprim_;
  using PolyhedralProjection<Real>::xdual_;
  using PolyhedralProjection<Real>::mul_;
  using PolyhedralProjection<Real>::res_;

public:

  DouglasRachfordProjection(const Vector<Real>               &xprim,
                    const Vector<Real>               &xdual,
                    const Ptr<BoundConstraint<Real>> &bnd,
                    const Ptr<Constraint<Real>>      &con,
                    const Vector<Real>               &mul,
                    const Vector<Real>               &res);

  DouglasRachfordProjection(const Vector<Real>               &xprim,
                    const Vector<Real>               &xdual,
                    const Ptr<BoundConstraint<Real>> &bnd,
                    const Ptr<Constraint<Real>>      &con,
                    const Vector<Real>               &mul,
                    const Vector<Real>               &res,
                    ParameterList                    &list);

  void project(Vector<Real> &x, std::ostream &stream = std::cout) override;

private:

  Real residual_1d(const Vector<Real> &x) const;

  void residual_nd(Vector<Real> &r, const Vector<Real> &y) const;

  void project_bnd(Vector<Real> &x, const Vector<Real> &y) const;

  void project_con(Vector<Real> &x, const Vector<Real> &y) const;

  void project_DouglasRachford(Vector<Real> &x, std::ostream &stream = std::cout) const;

}; // class DouglasRachfordProjection

} // namespace ROL

#include "ROL_DouglasRachfordProjection_Def.hpp"

#endif
