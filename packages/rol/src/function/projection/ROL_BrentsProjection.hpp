// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BRENTSPROJECTION_H
#define ROL_BRENTSPROJECTION_H

#include "ROL_PolyhedralProjection.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<typename Real>
class BrentsProjection : public PolyhedralProjection<Real> {
private:
  int dim_;
  Ptr<Vector<Real>> xnew_, Px_;
  Real b_, mul1_, dlam1_, cdot_;

  Real DEFAULT_atol_, DEFAULT_rtol_, DEFAULT_ltol_;
  int DEFAULT_maxit_, DEFAULT_verbosity_;

  Real atol_, rtol_, ltol_;
  int maxit_, verbosity_;

  Real ctol_;

  using PolyhedralProjection<Real>::bnd_;
  using PolyhedralProjection<Real>::con_;
  using PolyhedralProjection<Real>::xprim_;
  using PolyhedralProjection<Real>::xdual_;
  using PolyhedralProjection<Real>::mul_;
  using PolyhedralProjection<Real>::res_;

  void initialize(const Vector<Real>               &xprim,
                  const Vector<Real>               &xdual,
                  const Ptr<BoundConstraint<Real>> &bnd,
                  const Ptr<Constraint<Real>>      &con,
                  const Vector<Real>               &mul,
                  const Vector<Real>               &res);

public:

  BrentsProjection(const Vector<Real>               &xprim,
                   const Vector<Real>               &xdual,
                   const Ptr<BoundConstraint<Real>> &bnd,
                   const Ptr<Constraint<Real>>      &con,
                   const Vector<Real>               &mul,
                   const Vector<Real>               &res);

  BrentsProjection(const Vector<Real>               &xprim,
                   const Vector<Real>               &xdual,
                   const Ptr<BoundConstraint<Real>> &bnd,
                   const Ptr<Constraint<Real>>      &con,
                   const Vector<Real>               &mul,
                   const Vector<Real>               &res,
                   ParameterList                    &list);

  void project(Vector<Real> &x, std::ostream &stream = std::cout) override;

private:

  Real residual(const Vector<Real> &x) const;

  void update_primal(Vector<Real> &y, const Vector<Real> &x, const Real lam) const;

  void project_df(Vector<Real> &x, Real &lam, Real &dlam, std::ostream &stream = std::cout) const;

}; // class RiddersProjection

} // namespace ROL

#include "ROL_BrentsProjection_Def.hpp"

#endif
