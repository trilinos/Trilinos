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


#ifndef ROL_FLETCHEROBJECTIVEE_H
#define ROL_FLETCHEROBJECTIVEE_H

#include "ROL_FletcherObjectiveBase.hpp"

namespace ROL {

template<typename Real>
class FletcherObjectiveE : public FletcherObjectiveBase<Real> {
private:
  // Temporaries
  Ptr<Vector<Real>> Tv_;       // Temporary for matvecs
  Ptr<Vector<Real>> w_;        // first component of augmented system solve solution
  Ptr<Vector<Real>> v_;        // second component of augmented system solve solution
  Ptr<Vector<Real>> wdual_;    // first component of augmented system solve solution in dual space
  Ptr<Vector<Real>> wg_;       // first component of augmented system solve solution for gradient
  Ptr<Vector<Real>> vg_;       // second component of augmented system solve solution for gradient
  Ptr<Vector<Real>> xzeros_;   // zero vector
  Ptr<Vector<Real>> czeros_;   // zero vector

  // Required for Fletcher penalty function definition
  using FletcherObjectiveBase<Real>::obj_;
  using FletcherObjectiveBase<Real>::con_;
  using FletcherObjectiveBase<Real>::sigma_;
  using FletcherObjectiveBase<Real>::delta_;    // regularization parameter
  using FletcherObjectiveBase<Real>::quadPenaltyParameter_;
  using FletcherObjectiveBase<Real>::useInexact_;
  using FletcherObjectiveBase<Real>::HessianApprox_;
  using FletcherObjectiveBase<Real>::nfval_;
  using FletcherObjectiveBase<Real>::ngval_;
  using FletcherObjectiveBase<Real>::ncval_;
  using FletcherObjectiveBase<Real>::fPhi_;     // value of penalty function
  using FletcherObjectiveBase<Real>::gPhi_;     // gradient of penalty function
  using FletcherObjectiveBase<Real>::y_;        // multiplier estimate
  using FletcherObjectiveBase<Real>::fval_;     // value of objective function
  using FletcherObjectiveBase<Real>::g_;        // gradient of objective value
  using FletcherObjectiveBase<Real>::c_;        // constraint value
  using FletcherObjectiveBase<Real>::scaledc_;  // sigma_ * c_
  using FletcherObjectiveBase<Real>::gL_;       // gradient of Lagrangian (g - A*y)
  using FletcherObjectiveBase<Real>::gLdual_;   // dual gradient of Lagrangian (g - A*y)
  using FletcherObjectiveBase<Real>::cnorm_;    // norm of constraint violation
  using FletcherObjectiveBase<Real>::xprim_;
  using FletcherObjectiveBase<Real>::xdual_;
  using FletcherObjectiveBase<Real>::cprim_;
  using FletcherObjectiveBase<Real>::cdual_;
  using FletcherObjectiveBase<Real>::multSolverError_; // Error from augmented system solve in value()
  using FletcherObjectiveBase<Real>::gradSolveError_;  // Error from augmented system solve in gradient()
  using FletcherObjectiveBase<Real>::krylov_;
  using FletcherObjectiveBase<Real>::iterKrylov_;
  using FletcherObjectiveBase<Real>::flagKrylov_;
  using FletcherObjectiveBase<Real>::v1_;
  using FletcherObjectiveBase<Real>::v2_;
  using FletcherObjectiveBase<Real>::vv_;
  using FletcherObjectiveBase<Real>::w1_;
  using FletcherObjectiveBase<Real>::w2_;
  using FletcherObjectiveBase<Real>::ww_;
  using FletcherObjectiveBase<Real>::b1_;
  using FletcherObjectiveBase<Real>::b2_;
  using FletcherObjectiveBase<Real>::bb_;

  class AugSystem : public LinearOperator<Real> {
  private:
    const Ptr<Constraint<Real>> con_;
    const Ptr<const Vector<Real>> x_;
    const Real delta_;
  public:
    AugSystem(const Ptr<Constraint<Real>> &con,
              const Ptr<const Vector<Real>> &x,
              const Real delta) : con_(con), x_(x), delta_(delta) {}

    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

      con_->applyAdjointJacobian(*Hvp.get(0), *vp.get(1), *x_, tol);
      Hvp.get(0)->plus(vp.get(0)->dual());

      con_->applyJacobian(*Hvp.get(1), *vp.get(0), *x_, tol);
      Hvp.get(1)->axpy(-delta_*delta_, vp.get(1)->dual());
    }
  };

  class AugSystemPrecond : public LinearOperator<Real> {
  private:
    const Ptr<Constraint<Real>> con_;
    const Ptr<const Vector<Real>> x_;
    const Ptr<const Vector<Real>> g_;
  public:
    AugSystemPrecond(const Ptr<Constraint<Real>> con,
                     const Ptr<const Vector<Real>> x,
                     const Ptr<const Vector<Real>> g) : con_(con), x_(x), g_(g) {}

    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual());
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

      Hvp.set(0, vp.get(0)->dual());
      con_->applyPreconditioner(*(Hvp.get(1)),*(vp.get(1)),*x_,*g_, tol); 
    }
  };

public:
  FletcherObjectiveE(const ROL::Ptr<Objective<Real>> &obj,
                     const ROL::Ptr<Constraint<Real>> &con,
                     const Vector<Real> &xprim,
                     const Vector<Real> &xdual,
                     const Vector<Real> &cprim,
                     const Vector<Real> &cdual,
                     ROL::ParameterList &parlist);

  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

protected:

  void solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2, const Vector<Real> &b1, const Vector<Real> &b2, const Vector<Real> &x, Real &tol, bool refine = false) override;

}; // class FletcherObjectiveE

} // namespace ROL

#include "ROL_FletcherObjectiveE_Def.hpp"

#endif
