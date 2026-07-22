// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUNDFLETCHER_H
#define ROL_BOUNDFLETCHER_H

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
class BoundFletcher : public FletcherBase<Real> {
private:
  // Required for Fletcher penalty function definition
  using FletcherBase<Real>::obj_;
  using FletcherBase<Real>::con_;
  Ptr<const Vector<Real> > low_;
  Ptr<const Vector<Real> > upp_;

  using FletcherBase<Real>::penaltyParameter_;
  using FletcherBase<Real>::quadPenaltyParameter_;

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

  using FletcherBase<Real>::cnorm_;    // norm of constraint violation

  using FletcherBase<Real>::isValueComputed_;
  using FletcherBase<Real>::isGradientComputed_;
  using FletcherBase<Real>::isMultiplierComputed_;
  using FletcherBase<Real>::isObjValueComputed_;
  using FletcherBase<Real>::isObjGradComputed_;
  using FletcherBase<Real>::isConValueComputed_;

  using FletcherBase<Real>::delta_;                  // regularization parameter

  Ptr<Vector<Real> > Q_;
  Ptr<Vector<Real> > umx_;
  Ptr<Vector<Real> > DQ_;
  Ptr<Vector<Real> > Qsqrt_;

  int HessianApprox_;
  // Flag to determine type of augmented system to solve
  // AugSolve_ = 0 : Symmetric system
  //   [   I      Q^{1/2} A][w] = [b1]
  //   [A'Q^{1/2}          ][v]   [b2]
  // AugSolve_ = 1 : Nonsymmetric system
  //   [ I  A][w] = [b1]
  //   [A'Q  ][v]   [b2]
  int AugSolve_;

  Ptr<Vector<Real> > QsgL_;     // scaled gradient of Lagrangian Q^{1/2}*(g-A*y)
  Ptr<Vector<Real> > QgL_;      // scaled gradient of Lagrangian Q*(g-A*y)
  Ptr<Vector<Real> > Qsg_;      // Scaled gradient of objective Q^{1/2}*g
  Ptr<Vector<Real> > DQgL_;     // gradient of Lagrangian scaled by DQ, DQ*(g-A*y)
  Ptr<Vector<Real> > Qv_;       // used in augmented system solve

  // Temporaries
  Ptr<Vector<Real> > Tv_;       // Temporary for matvecs
  Ptr<Vector<Real> > w_;        // first component of augmented system solve solution
  Ptr<Vector<Real> > v_;        // second component of augmented system solve solution
  Ptr<Vector<Real> > htmp1_;    // Temporary for rhs
  Ptr<Vector<Real> > htmp2_;    // Temporary for rhs

  Ptr<Vector<Real> > xzeros_;   // zero vector
  Ptr<Vector<Real> > czeros_;   // zero vector

  using FletcherBase<Real>::useInexact_;

  bool isQComputed_;
  bool isDQComputed_;

  using FletcherBase<Real>::multSolverError_;         // Error from augmented system solve in value()
  using FletcherBase<Real>::gradSolveError_;          // Error from augmented system solve in gradient()

  // For Augmented system solves
  using FletcherBase<Real>::krylov_;
  using FletcherBase<Real>::iterKrylov_;
  using FletcherBase<Real>::flagKrylov_;
  using FletcherBase<Real>::v1_;
  using FletcherBase<Real>::v2_;
  using FletcherBase<Real>::vv_;
  using FletcherBase<Real>::w1_;
  using FletcherBase<Real>::w2_;
  using FletcherBase<Real>::ww_;
  using FletcherBase<Real>::b1_;
  using FletcherBase<Real>::b2_;
  using FletcherBase<Real>::bb_;

  class DiffLower : public Elementwise::BinaryFunction<Real> {
  public:
    DiffLower(void) {}
    Real apply(const Real& x, const Real& y) const {
      const Real NINF(ROL_NINF<Real>());
      return (y <= NINF ? static_cast<Real>(-1.) : x - y);
    }
  };

  class DiffUpper : public Elementwise::BinaryFunction<Real> {
  public:
    DiffUpper(void) {}
    Real apply(const Real& x, const Real& y) const {
      const Real INF(ROL_INF<Real>());
      return (y >= INF ? static_cast<Real>(-1.) : y - x);
    }
  };

  class FormQ : public Elementwise::BinaryFunction<Real> {
  public:
    FormQ(void) {}
    Real apply(const Real& x, const Real& y) const {
      Real zero(0.);
      if( x < zero  && y < zero) {
        return static_cast<Real>(1);
      }
      if( x < zero  && y >= zero ) {
        return y;
      }
      if( x >= zero && y < zero ) {
        return x;
      }
      return std::min(x, y);
    }
  };

  class FormDQ : public Elementwise::BinaryFunction<Real> {
  public:
    FormDQ(void) {}
    Real apply(const Real& x, const Real& y) const {
      Real zero(0.), one(1.), mone(-1.);
      if( x < zero  && y < zero) {
        return zero;
      }
      if( x < zero  && y >= zero ) {
        return mone;
      }
      if( x >= zero && y < zero ) {
        return one;
      }
      if( x < y ) {
        return one;
      } else if( y < x) {
        return mone;
      } else {
        return zero;
      }
    }
  };

  class AugSystemSym : public LinearOperator<Real> {
  private:
    const Ptr<Constraint<Real> > con_;
    const Ptr<const Vector<Real> > x_;
    const Ptr<Vector<Real> > Qsqrt_;
    const Ptr<Vector<Real> > Qv_;
    const Real delta_;
  public:
    AugSystemSym(const Ptr<Constraint<Real> > &con,
                 const Ptr<const Vector<Real> > &x,
                 const Ptr<Vector<Real> > &Qsqrt,
                 const Ptr<Vector<Real> > &Qv,
                 const Real delta) : con_(con), x_(x), Qsqrt_(Qsqrt), Qv_(Qv), delta_(delta) {}

    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

      con_->applyAdjointJacobian(*(Hvp.get(0)), *(vp.get(1)), *x_, tol);
      Hvp.get(0)->applyBinary(Elementwise::Multiply<Real>(), *Qsqrt_);
      Hvp.get(0)->plus(*(vp.get(0)));

      Qv_->set(*(vp.get(0)));
      Qv_->applyBinary(Elementwise::Multiply<Real>(), *Qsqrt_);
      con_->applyJacobian(*(Hvp.get(1)), *(Qv_), *x_, tol);
      Hvp.get(1)->axpy(-delta_*delta_, *(vp.get(1)));
    }
  };

  class AugSystemNonSym : public LinearOperator<Real> {
  private:
    const Ptr<Constraint<Real> > con_;
    const Ptr<const Vector<Real> > x_;
    const Ptr<Vector<Real> > Q_;
    const Ptr<Vector<Real> > Qv_;
    const Real delta_;
  public:
    AugSystemNonSym(const Ptr<Constraint<Real> > &con,
                    const Ptr<const Vector<Real> > &x,
                    const Ptr<Vector<Real> > &Q,
                    const Ptr<Vector<Real> > &Qv,
                    const Real delta) : con_(con), x_(x), Q_(Q), Qv_(Qv), delta_(delta) {}

    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

      con_->applyAdjointJacobian(*(Hvp.get(0)), *(vp.get(1)), *x_, tol);
      Hvp.get(0)->plus(*(vp.get(0)));

      Qv_->set(*(vp.get(0)));
      Qv_->applyBinary(Elementwise::Multiply<Real>(), *Q_);
      con_->applyJacobian(*(Hvp.get(1)), *(Qv_), *x_, tol);
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
  BoundFletcher(const ROL::Ptr<Objective<Real> > &obj,
                const ROL::Ptr<Constraint<Real> > &con,
                const ROL::Ptr<BoundConstraint<Real> > &bnd,
                const Vector<Real> &optVec,
                const Vector<Real> &conVec,
                ROL::ParameterList &parlist)
  : FletcherBase<Real>(obj, con), isQComputed_(false), isDQComputed_(false) {
      
      low_ = bnd->getLowerBound();
      upp_ = bnd->getUpperBound();

      gPhi_    = optVec.dual().clone();
      y_       = conVec.dual().clone();
      g_       = optVec.dual().clone();
      gL_      = optVec.dual().clone();
      c_       = conVec.clone();
      scaledc_ = conVec.clone();

      Q_ = optVec.clone();
      DQ_ = optVec.clone();
      umx_ = optVec.clone();
      Qsqrt_ = optVec.clone();
      Qv_ = optVec.dual().clone();
      QsgL_ = optVec.dual().clone();
      QgL_ = optVec.dual().clone();
      Qsg_ = optVec.dual().clone();
      DQgL_ = optVec.dual().clone();

      Tv_ = optVec.dual().clone();
      w_ = optVec.dual().clone();
      v_ = conVec.dual().clone();
      htmp1_ = optVec.dual().clone();
      htmp2_ = conVec.clone();

      xzeros_ = optVec.dual().clone();
      xzeros_->zero();
      czeros_ = conVec.clone();
      czeros_->zero();

      v1_ = optVec.dual().clone();
      v2_ = conVec.dual().clone();
      vv_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>> >({v1_, v2_}));

      w1_ = optVec.dual().clone();
      w2_ = conVec.dual().clone();
      ww_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>> >({w1_, w2_}));

      b1_ = optVec.dual().clone();
      b2_ = conVec.clone();
      bb_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>> >({b1_, b2_}));

      ROL::ParameterList& sublist = parlist.sublist("Step").sublist("Fletcher");
      HessianApprox_ = sublist.get("Level of Hessian Approximation",  0);

      AugSolve_ = sublist.get("Type of Augmented System Solve",  0);
      AugSolve_ = (0 < AugSolve_ && AugSolve_ < 2) ? AugSolve_ : 0; 

      penaltyParameter_ = sublist.get("Penalty Parameter", 1.0);
      quadPenaltyParameter_ = sublist.get("Quadratic Penalty Parameter", 0.0);

      delta_ = sublist.get("Regularization Parameter", 0.0);

      useInexact_ = sublist.get("Inexact Solves", false);
      
      ROL::ParameterList krylovList;
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
    isQComputed_ = (flag ? false : isQComputed_);
    isDQComputed_ = (flag ? false : isDQComputed_);
    multSolverError_ = (flag ? ROL_INF<Real>() : multSolverError_);
    gradSolveError_ = (flag ? ROL_INF<Real>() : gradSolveError_);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    if( isValueComputed_ && multSolverError_*cnorm_ <= tol) {
      tol = multSolverError_*cnorm_;
      return fPhi_;
    }

    Real zero(0);

    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;

    FletcherBase<Real>::objValue(x, tol2); tol2 = origTol;
    multSolverError_ = origTol / (static_cast<Real>(2) * std::max(static_cast<Real>(1), cnorm_));
    computeMultipliers(x, multSolverError_);
    tol = multSolverError_;

    fPhi_ = fval_ - c_->dot(y_->dual());

    if( quadPenaltyParameter_ > zero ) {
      fPhi_ = fPhi_ + Real(0.5)*quadPenaltyParameter_*(c_->dot(c_->dual()));
    }

    isValueComputed_ = true;

    return fPhi_;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    if( isGradientComputed_ && gradSolveError_ <= tol) {
      tol = gradSolveError_;
      g.set(*gPhi_);
      return;
    }

    Real zero(0);

    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;

    gradSolveError_ = origTol / static_cast<Real>(2);
    computeMultipliers(x, gradSolveError_);

    bool refine = isGradientComputed_;

    switch( AugSolve_ ) {
      case 0: {
        solveAugmentedSystem( *w_, *v_, *xzeros_, *c_, x, gradSolveError_, refine );
        gradSolveError_ += multSolverError_;
        tol = gradSolveError_;

        w_->applyBinary(Elementwise::Multiply<Real>(), *Qsqrt_);
        con_->applyAdjointHessian( *gPhi_, *y_, *w_, x, tol2 ); tol2 = origTol;
        obj_->hessVec( *Tv_, *w_, x, tol2 ); tol2 = origTol;
        gPhi_->axpy( static_cast<Real>(-1), *Tv_ );

        con_->applyAdjointJacobian( *Tv_, *v_, x, tol2); tol2 = origTol;
        gPhi_->axpy( -penaltyParameter_, *Tv_);
        
        Tv_->applyBinary(Elementwise::Multiply<Real>(), *DQgL_);
        gPhi_->plus( *Tv_ );

        con_->applyAdjointHessian( *Tv_, *v_, *QgL_, x, tol2 ); tol2 = origTol;
        gPhi_->plus( *Tv_ );

        gPhi_->plus( *gL_ );
        break;
      }
      case 1: {
        solveAugmentedSystem( *w_, *v_, *xzeros_, *c_, x, gradSolveError_, refine );
        gradSolveError_ += multSolverError_;
        tol = gradSolveError_;

        gPhi_->set( *w_ );
        gPhi_->scale( penaltyParameter_ );
        Tv_->set( *w_ );
        Tv_->applyBinary( Elementwise::Multiply<Real>(), *DQgL_ );
        gPhi_->axpy(static_cast<Real>(-1), *Tv_);
        
        w_->applyBinary( Elementwise::Multiply<Real>(), *Q_ );
        obj_->hessVec( *Tv_, *w_, x, tol2); tol2 = origTol;
        gPhi_->axpy( static_cast<Real>(-1), *Tv_ );
        con_->applyAdjointHessian( *Tv_, *y_, *w_, x, tol2 ); tol2 = origTol;
        gPhi_->plus( *Tv_ );

        con_->applyAdjointHessian( *Tv_, *v_, *QgL_, x, tol2 ); tol2 = origTol;
        gPhi_->plus( *Tv_ );

        gPhi_->plus( *gL_ );

        break;
      }
    }

    if( quadPenaltyParameter_ > zero ) {
      con_->applyAdjointJacobian( *Tv_, *c_, x, tol2 ); tol2 = origTol;
      gPhi_->axpy( quadPenaltyParameter_, *Tv_ );
    }

    g.set(*gPhi_);
    isGradientComputed_ = true;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Real zero(0);

    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;

    // Make sure everything is already computed
    value(x, tol2); tol2 = origTol;
    computeMultipliers(x, tol2); tol2 = origTol;
    gradient(*Tv_, x, tol2); tol2 = origTol;

    switch( AugSolve_ ) {
      case 0: {
        // hv <- HL*v
        obj_->hessVec( hv, v, x, tol2 ); tol2 = origTol;
        con_->applyAdjointHessian( *Tv_, *y_, v, x, tol2 ); tol2 = origTol;
        hv.axpy(static_cast<Real>(-1), *Tv_ );

        // htmp1_ <- Q^{1/2}*hv
        htmp1_->set(hv);
        htmp1_->applyBinary(Elementwise::Multiply<Real>(), *Qsqrt_);
        htmp1_->scale(static_cast<Real>(-1));
        // htmp2_ <- A'*(R(gL) - sigma*I)*v
        Tv_->set( *DQgL_ );
        Tv_->applyBinary( Elementwise::Multiply<Real>(), v );
        Tv_->axpy(-penaltyParameter_, v);
        con_->applyJacobian( *htmp2_, *Tv_, x, tol2); tol2 = origTol;
        // v_ <- - (A'QA)^-1 [A' (R(gL)-sigma I) u + A' Q HL u]
        solveAugmentedSystem( *w_, *v_, *htmp1_, *htmp2_, x, tol2 ); tol2 = origTol;
        con_->applyAdjointJacobian( *Tv_, *v_, x, tol2 ); tol2 = origTol;
        hv.plus(*Tv_);

        con_->applyJacobian( *htmp2_, v, x, tol2 ); tol2 = origTol;
        solveAugmentedSystem( *w_, *v_, *xzeros_, *htmp2_, x, tol2 ); tol2 = origTol;
        con_->applyAdjointJacobian( *Tv_, *v_, x, tol2 ); tol2 = origTol;
        hv.axpy( -penaltyParameter_, *Tv_);
        Tv_->applyBinary( Elementwise::Multiply<Real>(), *DQgL_ );
        hv.plus( *Tv_ );

        w_->applyBinary( Elementwise::Multiply<Real>(), *Qsqrt_ );
        obj_->hessVec( *Tv_, *w_, x, tol2 ); tol2 = origTol;
        hv.axpy( static_cast<Real>(-1), *Tv_);
        con_->applyAdjointHessian( *Tv_, *y_, *w_, x, tol2 ); tol2 = origTol;
        hv.plus( *Tv_ );
        break;
      }
      case 1: {
        // hv <- HL*v
        obj_->hessVec( hv, v, x, tol2 ); tol2 = origTol;
        con_->applyAdjointHessian( *Tv_, *y_, v, x, tol2 ); tol2 = origTol;
        hv.axpy(static_cast<Real>(-1), *Tv_ );

        htmp1_->set(hv);
        Tv_->set(v);
        Tv_->applyBinary( Elementwise::Multiply<Real>(), *DQgL_ );
        Tv_->axpy( -penaltyParameter_, v );
        Tv_->scale( static_cast<Real>(-1) );
        con_->applyJacobian( *htmp2_, *Tv_, x, tol2 ); tol2 = origTol;
        solveAugmentedSystem( *w_, *v_, *htmp1_, *htmp2_, x, tol2 ); tol2 = origTol;
        hv.set( *w_ );

        con_->applyJacobian( *htmp2_, v, x, tol2 ); tol2 = origTol;
        solveAugmentedSystem( *w_, *v_, *xzeros_, *htmp2_, x, tol2 ); tol2 = origTol;
        hv.axpy( penaltyParameter_, *w_ );
        Tv_->set( *w_ );
        Tv_->applyBinary( Elementwise::Multiply<Real>(), *DQgL_ );
        hv.axpy( static_cast<Real>(-1), *Tv_ );

        w_->applyBinary( Elementwise::Multiply<Real>(), *Q_ );
        obj_->hessVec( *Tv_, *w_, x, tol2 ); tol2 = tol;
        hv.axpy( static_cast<Real>(-1), *Tv_);
        con_->applyAdjointHessian( *Tv_, *y_, *w_, x, tol2 ); tol2 = origTol;
        hv.plus( *Tv_ );
        break;
      }
    }

    if( quadPenaltyParameter_ > zero ) {
      con_->applyJacobian( *b2_, v, x, tol2 ); tol2 = origTol;
      con_->applyAdjointJacobian( *Tv_, *b2_, x, tol2 ); tol2 = origTol;
      hv.axpy( quadPenaltyParameter_, *Tv_ );
      con_->applyAdjointHessian( *Tv_, *c_, v, x, tol2); tol2 = origTol;
      hv.axpy( -quadPenaltyParameter_, *Tv_ );
    }

  }

  void solveAugmentedSystem(Vector<Real> &v1,
                            Vector<Real> &v2,
                            const Vector<Real> &b1,
                            const Vector<Real> &b2,
                            const Vector<Real> &x,
                            Real &tol,
                            bool refine = false) {
    // Ignore tol for now
    ROL::Ptr<LinearOperator<Real> > K;
    switch( AugSolve_ ) {
      case 0: {
        K = ROL::makePtr<AugSystemSym>(con_, makePtrFromRef(x), Qsqrt_, Qv_, delta_);
        break;
      }
      case 1: {
        K = ROL::makePtr<AugSystemNonSym>(con_, makePtrFromRef(x), Q_, Qv_, delta_);
        break;        
      }
    }
    ROL::Ptr<LinearOperator<Real> > P
      = ROL::makePtr<AugSystemPrecond>(con_, makePtrFromRef(x));

    b1_->set(b1);
    b2_->set(b2);

    if( refine ) {
      // TODO: Make sure this tol is actually ok...
      Real origTol = tol;
      w1_->set(v1);
      w2_->set(v2);
      K->apply(*vv_, *ww_, tol); tol = origTol;

      b1_->axpy( static_cast<Real>(-1), *v1_ );
      b2_->axpy( static_cast<Real>(-1), *v2_ );
    }

    v1_->zero();
    v2_->zero();

    // If inexact, change tolerance
    if( useInexact_ ) {
      krylov_->resetAbsoluteTolerance(tol);
    }

    flagKrylov_ = 0;
    tol = krylov_->run(*vv_,*K,*bb_,*P,iterKrylov_,flagKrylov_);

    if( refine ) {
      v1.plus(*v1_);
      v2.plus(*v2_);
    } else {
      v1.set(*v1_);
      v2.set(*v2_);
    }
  }

  void computeMultipliers(const Vector<Real>& x, const Real tol) {
    if( isMultiplierComputed_ && multSolverError_ <= tol) {
      return;
    }

    if( !isMultiplierComputed_ ) {
      Real tol2 = tol;
      FletcherBase<Real>::objGrad(x, tol2); tol2 = tol;
      FletcherBase<Real>::conValue(x, tol2); tol2 = tol;
      cnorm_ = c_->norm();
      computeQ(x);
      computeDQ(x);
    }

    bool refine = isMultiplierComputed_;

    switch( AugSolve_ ) {
      case 0: {
        Qsg_->set(*g_);
        Qsg_->applyBinary(Elementwise::Multiply<Real>(), *Qsqrt_);

        multSolverError_ = tol;
        solveAugmentedSystem(*QsgL_, *y_, *Qsg_, *scaledc_, x, multSolverError_, refine);

        gL_->set(*QsgL_);
        gL_->applyBinary(Elementwise::Divide<Real>(), *Qsqrt_);
        QgL_->set(*QsgL_);
        QgL_->applyBinary(Elementwise::Multiply<Real>(), *Qsqrt_);
        break;
      }
      case 1: {
        multSolverError_ = tol;
        solveAugmentedSystem(*gL_, *y_, *g_, *scaledc_, x, multSolverError_, refine);
        QgL_->set(*gL_);
        QgL_->applyBinary(Elementwise::Multiply<Real>(), *Q_);
        break;
      }
    }
    
    DQgL_->set(*gL_);
    DQgL_->applyBinary(Elementwise::Multiply<Real>(), *DQ_);

    isMultiplierComputed_ = true;
  }

  void computeQ(const Vector<Real>& x) {
    if( isQComputed_ ) {
      return;
    }

    Q_->set(x);
    Q_->applyBinary(DiffLower(), *low_);
    umx_->set(x);
    umx_->applyBinary(DiffUpper(), *upp_);
    Q_->applyBinary(FormQ(), *umx_);
    Qsqrt_->set(*Q_);
    Qsqrt_->applyUnary(Elementwise::SquareRoot<Real>());

    isQComputed_ = true;
  }

  void computeDQ(const Vector<Real> &x) {
    if( isDQComputed_ ) {
      return;
    }

    DQ_->set(x);
    DQ_->applyBinary(DiffLower(), *low_);
    umx_->set(x);
    umx_->applyBinary(DiffUpper(), *upp_);   
    DQ_->applyBinary(FormDQ(), *umx_);

    isDQComputed_ = true;
  }

}; // class Fletcher

} // namespace ROL

#endif
