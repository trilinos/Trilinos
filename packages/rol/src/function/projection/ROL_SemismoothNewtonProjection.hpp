// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SEMISMOOTHNEWTONPROJECTION_H
#define ROL_SEMISMOOTHNEWTONPROJECTION_H

#include "ROL_PolyhedralProjection.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_KrylovFactory.hpp"

namespace ROL {

template<typename Real>
class SemismoothNewtonProjection : public PolyhedralProjection<Real> {
private:
  int dim_;
  Ptr<Krylov<Real>> krylov_;
  Ptr<Vector<Real>> xnew_, lnew_, dlam_;

  Real DEFAULT_atol_, DEFAULT_rtol_, DEFAULT_stol_;
  Real DEFAULT_decr_, DEFAULT_factor_, DEFAULT_regscale_, DEFAULT_errscale_;
  int DEFAULT_maxit_, DEFAULT_lstype_, DEFAULT_verbosity_;
  bool DEFAULT_useproj_;

  Real atol_, rtol_, stol_, decr_, factor_, regscale_, errscale_, ctol_;
  int maxit_, lstype_, verbosity_;
  bool useproj_;

  using PolyhedralProjection<Real>::bnd_;
  using PolyhedralProjection<Real>::con_;
  using PolyhedralProjection<Real>::xprim_;
  using PolyhedralProjection<Real>::xdual_;
  using PolyhedralProjection<Real>::mul_;
  using PolyhedralProjection<Real>::res_;

public:

  SemismoothNewtonProjection(const Vector<Real>               &xprim,
                             const Vector<Real>               &xdual,
                             const Ptr<BoundConstraint<Real>> &bnd,
                             const Ptr<Constraint<Real>>      &con,
                             const Vector<Real>               &mul,
                             const Vector<Real>               &res);

  SemismoothNewtonProjection(const Vector<Real>               &xprim,
                             const Vector<Real>               &xdual,
                             const Ptr<BoundConstraint<Real>> &bnd,
                             const Ptr<Constraint<Real>>      &con,
                             const Vector<Real>               &mul,
                             const Vector<Real>               &res,
                             ParameterList                    &list);

  void project(Vector<Real> &x, std::ostream &stream = std::cout) override;

private:

  // Jv = inv(A DP(y) A^T) v
  // y  = P(x + A^T lam)
  class Jacobian : public LinearOperator<Real> {
    private:
      const Ptr<Constraint<Real>>      con_;
      const Ptr<BoundConstraint<Real>> bnd_;
      const Ptr<const Vector<Real>>    y_;
      const Ptr<Vector<Real>> xdual_, xprim_;
      const Real alpha_;
    public:
      Jacobian(const Ptr<Constraint<Real>>      &con,
               const Ptr<BoundConstraint<Real>> &bnd,
               const Ptr<const Vector<Real>>    &y,
               const Ptr<Vector<Real>>          &xdual,
               const Ptr<Vector<Real>>          &xprim,
               const Real                        alpha = 1e-4)
        : con_(con), bnd_(bnd), y_(y), xdual_(xdual), xprim_(xprim), alpha_(alpha) {}
      void apply(Vector<Real> &Jx, const Vector<Real> &x, Real &tol) const {
        con_->applyAdjointJacobian(*xdual_,x.dual(),*y_,tol);
        xprim_->set(xdual_->dual());
        bnd_->pruneActive(*xprim_,*y_);
        con_->applyJacobian(Jx,*xprim_,*y_,tol);
        // This is a hack to make the Jacobian invertible
        Jx.axpy(alpha_,x);
      }
  };

  class Precond : public LinearOperator<Real> {
  private:
    const Real alpha_;
  public:
    Precond(Real alpha = 1e-4) : alpha_(alpha) {}
    void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      const Real one(1);
      Hv.set(v.dual());
      Hv.scale(one+alpha_);
    }
    void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      const Real one(1);
      Hv.set(v.dual());
      Hv.scale(one/(one+alpha_));
    }
  };

  Real residual(Vector<Real> &r, const Vector<Real> &y) const;

  void solve_newton_system(Vector<Real>       &s,
                           const Vector<Real> &r,
                           const Vector<Real> &y,
                           const Real          mu,
                           const Real          rho,
                           int                &iter,
                           int                &flag) const;

  void update_primal(Vector<Real>       &y,
                     const Vector<Real> &x,
                     const Vector<Real> &lam) const;

  void project_ssn(Vector<Real> &x,
                   Vector<Real> &lam,
                   Vector<Real> &dlam,
                   std::ostream &stream = std::cout) const;

   Real compute_tolerance() const;

}; // class SemismoothNewtonProjection

} // namespace ROL

#include "ROL_SemismoothNewtonProjection_Def.hpp"

#endif
