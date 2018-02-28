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


#ifndef ROL_FLETCHER_H
#define ROL_FLETCHER_H

#include "ROL_FletcherBase.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Krylov.hpp"
#include "ROL_PartitionedVector.hpp"
#include <iostream>

namespace ROL {

template <class Real>
class Fletcher : public FletcherBase<Real> {
private:
  // Required for Fletcher penalty function definition
  using FletcherBase<Real>::obj_;
  using FletcherBase<Real>::con_;

  using FletcherBase<Real>::penaltyParameter_;

  // Evaluation counters
  using FletcherBase<Real>::nfval_;
  using FletcherBase<Real>::ngval_;
  using FletcherBase<Real>::ncval_;

  using FletcherBase<Real>::fPhi_;                   // value of penalty function
  using FletcherBase<Real>::gPhi_;     // gradient of penalty function

  using FletcherBase<Real>::y_;        // multiplier estimate

  using FletcherBase<Real>::fval_;                   // value of objective function
  using FletcherBase<Real>::g_;        // gradient of objective value
  using FletcherBase<Real>::c_;        // constraint value
  using FletcherBase<Real>::scaledc_;  // penaltyParameter_ * c_
  using FletcherBase<Real>::gL_;       // gradient of Lagrangian (g - A*y)

  using FletcherBase<Real>::isValueComputed_;
  using FletcherBase<Real>::isGradientComputed_;
  using FletcherBase<Real>::isMultiplierComputed_;
  using FletcherBase<Real>::isObjValueComputed_;
  using FletcherBase<Real>::isObjGradComputed_;
  using FletcherBase<Real>::isConValueComputed_;

  using FletcherBase<Real>::delta_;                  // regularization parameter

  int HessianApprox_;

  // Temporaries
  Ptr<Vector<Real> > Tv_;       // Temporary for matvecs
  Ptr<Vector<Real> > w_;        // first component of augmented system solve solution
  Ptr<Vector<Real> > v_;        // second component of augmented system solve solution

  Ptr<Vector<Real> > xzeros_;   // zero vector
  Ptr<Vector<Real> > czeros_;   // zero vector

  using FletcherBase<Real>::useInexact_;

  using FletcherBase<Real>::multSolverError_;         // Error from augmented system solve in value()
  using FletcherBase<Real>::gradSolveError_;          // Error from augmented system solve in gradient()

  // For Augmented system solves
  using FletcherBase<Real>::krylov_;
  using FletcherBase<Real>::iterKrylov_;
  using FletcherBase<Real>::flagKrylov_;
  using FletcherBase<Real>::v1_;
  using FletcherBase<Real>::v2_;
  using FletcherBase<Real>::vv_;
  using FletcherBase<Real>::b1_;
  using FletcherBase<Real>::b2_;
  using FletcherBase<Real>::bb_;

  class AugSystem : public LinearOperator<Real> {
  private:
    const Ptr<Constraint<Real> > con_;
    const Ptr<const Vector<Real> > x_;
    const Real delta_;
  public:
    AugSystem(const Ptr<Constraint<Real> > &con,
              const Ptr<const Vector<Real> > &x,
              const Real delta) : con_(con), x_(x), delta_(delta) {}

    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

      con_->applyAdjointJacobian(*(Hvp.get(0)), *(vp.get(1)), *x_, tol);
      Hvp.get(0)->plus(*(vp.get(0)));

      con_->applyJacobian(*(Hvp.get(1)), *(vp.get(0)), *x_, tol);
      Hvp.get(1)->axpy(-delta_*delta_, *(vp.get(1)));
    }
  };

  class AugSystemPrecond : public LinearOperator<Real> {
  private:
    const Ptr<Constraint<Real> > con_;
    const Ptr<const Vector<Real> > x_;
  public:
    AugSystemPrecond(const Ptr<Constraint<Real> > con,
                     const Ptr<const Vector<Real> > x) : con_(con), x_(x) {}

    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual());
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Real zero(0);
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

      Hvp.set(0, *(vp.get(0)));
      // Second x should be dual, but unused?
      con_->applyPreconditioner(*(Hvp.get(1)),*(vp.get(1)),*x_,*x_, zero); 
    }
  };

public:
  Fletcher(const ROL::Ptr<Objective<Real> > &obj,
           const ROL::Ptr<Constraint<Real> > &con,
           const Vector<Real> &optVec,
           const Vector<Real> &conVec,
           Teuchos::ParameterList &parlist)
      : FletcherBase<Real>(obj, con) {

      gPhi_    = optVec.dual().clone();
      y_       = conVec.dual().clone();
      g_       = optVec.dual().clone();
      gL_      = optVec.dual().clone();
      c_       = conVec.clone();
      scaledc_ = conVec.clone();

      Tv_ = optVec.dual().clone();
      w_ = optVec.dual().clone();
      v_ = conVec.dual().clone();

      xzeros_ = optVec.dual().clone();
      xzeros_->zero();
      czeros_ = conVec.clone();
      czeros_->zero();

      v1_ = optVec.dual().clone();
      v2_ = conVec.dual().clone();
      vv_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>> >({v1_, v2_}));

      b1_ = optVec.dual().clone();
      b2_ = conVec.clone();
      bb_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>> >({b1_, b2_}));

      Teuchos::ParameterList& sublist = parlist.sublist("Step").sublist("Fletcher");
      HessianApprox_ = sublist.get("Level of Hessian Approximation",  0);
      penaltyParameter_ = sublist.get("Penalty Parameter", 1.0);

      delta_ = sublist.get("Regularization Parameter", 0.0);

      useInexact_ = sublist.get("Inexact Solves", false);
      
      Teuchos::ParameterList krylovList;
      Real atol = static_cast<Real>(1e-12);
      Real rtol = static_cast<Real>(1e-2);
      krylovList.sublist("General").sublist("Krylov").set("Type", "GMRES");
      krylovList.sublist("General").sublist("Krylov").set("Absolute Tolerance", atol);
      krylovList.sublist("General").sublist("Krylov").set("Relative Tolerance", rtol);
      krylovList.sublist("General").sublist("Krylov").set("Iteration Limit", 200);
      krylov_ = KrylovFactory<Real>(krylovList);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    con_->update(x,flag,iter);
    isValueComputed_ = (flag ? false : isValueComputed_);
    isGradientComputed_ = (flag ? false : isGradientComputed_);
    isMultiplierComputed_ = (flag ? false : isMultiplierComputed_);
    isObjValueComputed_ = (flag ? false : isObjValueComputed_);
    isObjGradComputed_ = (flag ? false : isObjGradComputed_);
    isConValueComputed_ = (flag ? false : isConValueComputed_);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    if( isValueComputed_ )
      return fPhi_;

    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;

    FletcherBase<Real>::objValue(x, tol2); tol2 = origTol;
    multSolverError_ = origTol / static_cast<Real>(2);
    computeMultipliers(x, multSolverError_);
    tol = multSolverError_;

    fPhi_ = fval_ - c_->dot(y_->dual());
    isValueComputed_ = true;

    return fPhi_;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    if( isGradientComputed_ && gradSolveError_ <= tol) {
      tol = gradSolveError_;
      g.set(*gPhi_);
      return;
    }

    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;

    gradSolveError_ = origTol / static_cast<Real>(2);
    computeMultipliers(x, gradSolveError_);

    // gPhi = sum y_i H_i w + sigma w + sum v_i H_i gL - H w + gL
    solveAugmentedSystem( *w_, *v_, *xzeros_, *c_, x, gradSolveError_ );
    gradSolveError_ += multSolverError_;
    tol = gradSolveError_;

    con_->applyAdjointHessian( *gPhi_, *y_, *w_, x, tol2 ); tol2 = origTol;
    gPhi_->axpy( penaltyParameter_, *w_ );

    obj_->hessVec( *Tv_, *w_, x, tol2 ); tol2 = origTol;
    gPhi_->axpy( static_cast<Real>(-1), *Tv_ );

    con_->applyAdjointHessian( *Tv_, *v_, *gL_, x, tol2 ); tol2 = origTol;
    gPhi_->plus( *Tv_ );

    gPhi_->plus( *gL_ );

    g.set(*gPhi_);

    isGradientComputed_ = true;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;

    computeMultipliers(x, tol);

    obj_->hessVec( hv, v, x, tol2 ); tol2 = origTol;
    con_->applyAdjointHessian( *Tv_, *y_, v, x, tol2 ); tol2 = origTol;
    hv.axpy(static_cast<Real>(-1), *Tv_ );

    tol2 = tol;
    solveAugmentedSystem( *w_, *v_, hv, *czeros_, x, tol2 ); tol2 = origTol;
    hv.scale( static_cast<Real>(-1) );
    hv.plus( *w_ );

    Tv_->set(v);
    tol2 = tol;
    solveAugmentedSystem( *w_, *v_, *Tv_, *czeros_, x, tol2 ); tol2 = origTol;
    hv.axpy(static_cast<Real>(-2)*penaltyParameter_, *w_);

    obj_->hessVec( *Tv_, *w_, x, tol2 ); tol2 = origTol;
    hv.plus( *Tv_ );
    con_->applyAdjointHessian( *Tv_, *y_, *w_, x, tol2 ); tol2 = origTol;
    hv.axpy( static_cast<Real>(-1), *Tv_ );

    hv.axpy(static_cast<Real>(2)*penaltyParameter_, v);
  }

  void solveAugmentedSystem(Vector<Real> &v1,
                            Vector<Real> &v2,
                            const Vector<Real> &b1,
                            const Vector<Real> &b2,
                            const Vector<Real> &x,
                            Real &tol) {
    // Ignore tol for now
    ROL::Ptr<LinearOperator<Real> > K
      = ROL::makePtr<AugSystem>(con_, makePtrFromRef(x), delta_);
    ROL::Ptr<LinearOperator<Real> > P
      = ROL::makePtr<AugSystemPrecond>(con_, makePtrFromRef(x));

    v1_->set(v1);
    v2_->set(v2);

    b1_->set(b1);
    b2_->set(b2);

    // If inexact, change tolerance
    if( useInexact_ ) {
      krylov_->resetAbsoluteTolerance(tol);
    }

    flagKrylov_ = 0;
    tol = krylov_->run(*vv_,*K,*bb_,*P,iterKrylov_,flagKrylov_);
    v1.set(*v1_);
    v2.set(*v2_);
  }

  void computeMultipliers(const Vector<Real>& x, const Real tol) {
    if( isMultiplierComputed_ && multSolverError_ <= tol) {
      return;
    }
    Real tol2 = tol;
    FletcherBase<Real>::objGrad(x, tol2); tol2 = tol;
    FletcherBase<Real>::conValue(x, tol2);

    multSolverError_ = tol;
    solveAugmentedSystem(*gL_, *y_, *g_, *scaledc_, x, multSolverError_);

    isMultiplierComputed_ = true;
  }

}; // class Fletcher

} // namespace ROL

#endif
