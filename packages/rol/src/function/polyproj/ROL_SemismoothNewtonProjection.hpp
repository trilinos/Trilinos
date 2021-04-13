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
  public:
    void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      Hv.set(v.dual());
    }
    void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      Hv.set(v.dual());
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

}; // class SemismoothNewtonProjection

} // namespace ROL

#include "ROL_SemismoothNewtonProjection_Def.hpp"

#endif
