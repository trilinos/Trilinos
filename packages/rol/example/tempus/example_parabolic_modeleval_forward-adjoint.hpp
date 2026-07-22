// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_parabolic.cpp
    \brief Shows how to solve a control problem based on the 1D semilinear parabolic equation,
           \f[
               u_t(t,x) - u_{xx}(t,x) + f(u(t,x),x) = z(t,x) \quad t\in (0,T], \; x\in (0,1)
           \f]
           with boundary conditions
           \f[
               u_x(t,0) = 0, \; u_x(t,1) = 0
           \f]
           and initial condition
           \f[
               u(0,x) = 0.
           \f]
           using Thyra::ModelEvaluator for operators and derivatives and Tempus for time stepping.
           The objective function is given by
           \f[
              J(u,z) = \frac{1}{2} \int_0^T \int_0^1 (u(t,x)-\bar{u}(x))^2\,\mathrm{d}x\,\mathrm{d}t
                       \frac{\alpha}{2} \int_0^T z(t)^2\,\mathrm{d}t.
           \f]
           We first define the Thyra::ModelEvaluator Parabolic model, followed by the constraint equation,
           followed by the objective function.
*/

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_DynamicConstraint.hpp"
#include "ROL_DynamicObjective.hpp"
#include "ROL_ThyraVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Thyra_ModelEvaluator.hpp"              // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterList.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"


/***********************************************************************
***********************  PARABOLIC MODEL AND ***************************
*********************  MODEL EVALUATOR WRAPPER  ************************
************************************************************************/

/// 1. ParabolicModel
/// -----------------
/// Let's start with the definition of the ParabolicModel.
/// The purpose of this class is to define linear and nonlinear operations
/// that will be used to implement a Thyra::ModelEvaluator wrapper object,
/// which is needed for the ROL::TempusDynamicConstraint interface.
/// Note that the public member functions operate on Thyra::VectorBase objects.
/// -----------------
template<class Real>
class ParabolicModel {

  typedef typename std::vector<Real>::size_type uint;

private: // ParabolicModel data members

  Real eps1_;
  Real eps2_;
  uint nx_;
  Real dx_;
  int  type_;
  std::vector<Real> pts_;
  std::vector<Real> wts_;

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> u_space_;
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> c_space_;
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> z_space_;

public: // ParabolicModel methods

  ParabolicModel(ROL::ParameterList &pl) {
    const Real one(1);
    type_     = pl.get("Reaction Type",             1);
    eps1_     = pl.get("Diffusion Scale",         1.0);
    eps2_     = pl.get("Reaction Scale",          1.0);
    nx_       = pl.get("Spatial Discretization",  128);
    dx_       = one/(static_cast<Real>(nx_)-one);

    pts_.resize(4,0.0);
    pts_[0] = -0.8611363115940526;
    pts_[1] = -0.3399810435848563;
    pts_[2] =  0.3399810435848563;
    pts_[3] =  0.8611363115940526;

    wts_.resize(4,0.0);
    wts_[0] = 0.3478548451374538;
    wts_[1] = 0.6521451548625461;
    wts_[2] = 0.6521451548625461;
    wts_[3] = 0.3478548451374538;

    Real sum = wts_[0]+wts_[1]+wts_[2]+wts_[3];
    wts_[0] *= 2.0/sum;
    wts_[1] *= 2.0/sum;
    wts_[2] *= 2.0/sum;
    wts_[3] *= 2.0/sum;

    // Create state space, constraint space and control space.
    u_space_ = Thyra::defaultSpmdVectorSpace<Real>(nx_);
    c_space_ = Thyra::defaultSpmdVectorSpace<Real>(nx_);
    z_space_ = Thyra::defaultSpmdVectorSpace<Real>(nx_+2);
  }

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> getStateSpace() const {
    return u_space_;
  }

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> getConstraintSpace() const {
    return c_space_;
  }

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> getControlSpace() const {
    return z_space_;
  }

  void apply_c_udot(const ROL::Ptr<Thyra::VectorBase<Real>> &Mu_ptr,
                    const ROL::Ptr<const Thyra::VectorBase<Real>> &u_ptr ) const {
    // This is the application of df/dxdot in ModelEvaluator notation.
    const Real zero(0), two(2), six(6);
    const Thyra::ConstDetachedVectorView<Real> u(u_ptr);
    Thyra::put_scalar(zero, Mu_ptr.ptr());
    { // scope for Mu eval
      Thyra::DetachedVectorView<Real> Mu(Mu_ptr);
      for (uint n = 0; n < nx_; n++) {
        if ( n < nx_-1 ) {
          Mu[n] += dx_/six*(two*u[n]+u[n+1]);
        }
        if ( n > 0 ) {
          Mu[n] += dx_/six*(u[n-1]+two*u[n]);
        }
      }
    } // end scope for Mu eval
  }

  void apply_c_u(const ROL::Ptr<Thyra::VectorBase<Real>> &jv_ptr,
                 const ROL::Ptr<const Thyra::VectorBase<Real>> &v_ptr,
                 const ROL::Ptr<const Thyra::VectorBase<Real>> &u_ptr) const {
    // This is the application of df/dx in ModelEvaluator notation.
    const Real zero(0), half(0.5);
    const Thyra::ConstDetachedVectorView<Real> v(v_ptr);
    const Thyra::ConstDetachedVectorView<Real> u(u_ptr);
    Thyra::put_scalar(zero, jv_ptr.ptr());
    Real phi1(0), phi2(0), f(0), x(0), w(0);
    { // scope for jv eval
      Thyra::DetachedVectorView<Real> jv(jv_ptr);
      for (uint n = 0; n < nx_; n++) {
        if ( n < nx_-1 ) {
          jv[n] += eps1_/dx_*(v[n]-v[n+1]); // Stiffness
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n+1);
            w = half*dx_*wts_[j];
            f = evaluate_nonlinearity(x,u,1);
            // Diagonal contribution
            phi1 = (static_cast<Real>(n+1)*dx_-x)/dx_;
            jv[n]+= w*f*phi1*phi1*v[n];
            // Off diagonal contribution
            phi2 = (x-static_cast<Real>(n)*dx_)/dx_;
            jv[n]+= w*f*phi1*phi2*v[n+1];
          }
        }
        if ( n > 0 ) {
          jv[n] += eps1_/dx_*(v[n]-v[n-1]); // Stiffness
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n-1);
            w = half*dx_*wts_[j];
            f = evaluate_nonlinearity(x,u,1);
            // Diagonal contribution
            phi1 = (x-static_cast<Real>(n-1)*dx_)/dx_;
            jv[n]+= w*f*phi1*phi1*v[n];
            // Off diagonal contribution
            phi2 = (static_cast<Real>(n)*dx_-x)/dx_;
            jv[n]+= w*f*phi1*phi2*v[n-1];
          }
        }
      }
    } // end scope for jv eval
  }

  void compute_alpha_c_udot_beta_c_u(std::vector<Real> &d,
                                     std::vector<Real> &o,
                                     const Thyra::VectorBase<Real> &u_vb,
                                     const Real &alpha,
                                     const Real &beta) const {
    const Real half(0.5), two(2), three(3), four(4), six(6);
    const Thyra::ConstDetachedVectorView<Real> u(u_vb);
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(nx_,alpha*four*dx_/six + beta*eps1_*two/dx_);
    d[0]     = alpha*dx_/three + beta*eps1_/dx_;
    d[nx_-1] = alpha*dx_/three + beta*eps1_/dx_;
    o.clear();
    o.resize(nx_-1,alpha*dx_/six - beta*eps1_/dx_);
    // Contribution from nonlinearity
    Real phi1(0), phi2(0), f(0), x(0), w(0);
    for (uint i=0; i<nx_; i++) {
      if (i<nx_-1) {
        for (uint j=0; j<4; j++) {
          x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*i+1);
          w = half*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = (static_cast<Real>(i+1)*dx_-x)/dx_;
          d[i]+= beta*w*f*phi1*phi1;
          // Off diagonal contribution
          phi2 = (x-static_cast<Real>(i)*dx_)/dx_;
          o[i]+= beta*w*f*phi1*phi2;
        }
      }
      if (i>0) {
        for (uint j=0; j<4; j++) {
          x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*i-1);
          w = half*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = (x-static_cast<Real>(i-1)*dx_)/dx_;
          d[i]+= beta*w*f*phi1*phi1;
        }
      }
    }
  }

  void apply_c_z(const ROL::Ptr<Thyra::VectorBase<Real>> &jv_ptr,
                 const ROL::Ptr<const Thyra::VectorBase<Real>> &v_ptr,
                 bool adjoint = false) const {
    const Real zero(0), four(4), six(6);
    const Thyra::ConstDetachedVectorView<Real> v(v_ptr);
    uint dim = ((adjoint == true) ? nx_+2 : nx_);
    Thyra::put_scalar(zero, jv_ptr.ptr());
    { // scope for jv eval
      Thyra::DetachedVectorView<Real> jv(jv_ptr);
      for (uint n = 0; n < dim; n++) {
        if ( adjoint ) {
          if ( n == 0 ) {
            jv[n] = -dx_/six*v[n];
          }
          else if ( n == 1 ) {
            jv[n] = -dx_/six*(four*v[n-1]+v[n]);
          }
          else if ( n == nx_ ) {
            jv[n] = -dx_/six*(four*v[n-1]+v[n-2]);
          }
          else if ( n == nx_+1 ) {
            jv[n] = -dx_/six*v[n-2];
          }
          else {
            jv[n] = -dx_/six*(v[n-2]+four*v[n-1]+v[n]);
          }
        }
        else {
          jv[n] -= dx_/six*(v[n]+four*v[n+1]+v[n+2]);
        }
      }
    } // end scope for jv eval
  }

  Real evaluate_solution(const Real x,
                         const Thyra::ConstDetachedVectorView<Real> &u) const {
    // Determine u(x)
    Real pt(0), val(0);
    for (uint n = 0; n < nx_; n++) {
      if (x <= static_cast<Real>(n+1)*dx_ && x >= static_cast<Real>(n)*dx_) {
        pt  = static_cast<Real>(n+1)*dx_;
        val = u[n]*(pt-x)/dx_;
        pt  = static_cast<Real>(n)*dx_;
        val+= u[n+1]*(x-pt)/dx_;
        break;
      }
      else if (x <= static_cast<Real>(n)*dx_ && x >= static_cast<Real>(n-1)*dx_) {
        pt  = static_cast<Real>(n)*dx_;
        val = u[n-1]*(pt-x)/dx_;
        pt  = static_cast<Real>(n-1)*dx_;
        val+= u[n]*(x-pt)/dx_;
        break;
      }
    }
    return val;
  }

  Real evaluate_nonlinearity(const Real x,
                             const Thyra::ConstDetachedVectorView<Real> &u,
                             const int deriv = 0) const {
    Real nl(0);
    Real val = evaluate_solution(x,u);
    if (type_==0) {       // Linear reaction
      const Real zero(0), one(1);
      nl = (deriv==0 ? val : (deriv==1 ? one : zero));
    }
    else if (type_==1) { // Cahn-Hillard reaction
      const Real one(1), two(2), three(3), six(6);
      nl = (deriv==0 ? (std::pow(val,three)-val)
         : (deriv==1 ? (three*std::pow(val,two)-one)
         : (six*val)));
    }
    else if (type_==2) { // Poisson-Boltzmann reaction
      const Real half(0.5);
      Real epu = std::exp(val), enu = std::exp(-val);
      nl = (deriv==0 ? half*(epu-enu)
         : (deriv==1 ? half*(epu+enu)
         : half*(epu-enu)));
    }
    else {               // Default linear reaction
      const Real zero(0), one(1);
      nl = (deriv==0 ? val : (deriv==1 ? one : zero));
    }
    return eps2_*nl;
  }

  void compute_residual(const ROL::Ptr<Thyra::VectorBase<Real>> &r_ptr,
                        const Thyra::VectorBase<Real> &udot_vb,
                        const Thyra::VectorBase<Real> &u_vb,
                        const Thyra::VectorBase<Real> &z_vb) const {
    const Real zero(0), half(0.5), two(2), four(4), six(6);
    Thyra::put_scalar(zero, r_ptr.ptr());
    Real x(0), w(0);
    const Thyra::ConstDetachedVectorView<Real> udot(udot_vb);
    const Thyra::ConstDetachedVectorView<Real> u(u_vb);
    const Thyra::ConstDetachedVectorView<Real> z(z_vb);
    { // scope for r eval
      Thyra::DetachedVectorView<Real> r(r_ptr);
      for (uint n = 0; n < nx_; n++) {
        if ( n < nx_-1 ) {
          r[n] += dx_/six*(two*udot[n]+udot[n+1]) + eps1_/dx_*(u[n]-u[n+1]); // Mass & stiffness
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n+1);
            w = half*dx_*wts_[j];
            r[n] += w*evaluate_nonlinearity(x,u,0)*(static_cast<Real>(n+1)*dx_-x)/dx_;
          }
        }
        if ( n > 0 ) {
          r[n] += dx_/six*(two*udot[n]+udot[n-1]) + eps1_/dx_*(u[n]-u[n-1]); // Mass & stiffness
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n-1);
            w = half*dx_*wts_[j];
            r[n] += w*evaluate_nonlinearity(x,u,0)*(x-static_cast<Real>(n-1)*dx_)/dx_;
          }
        }
        // Control
        r[n] -= dx_/six*(z[n]+four*z[n+1]+z[n+2]);
      }
    } // end scope for r eval

  }

  Real compute_norm(const Thyra::VectorBase<Real> &v_vb) const {
    return Thyra::norm_2(v_vb);
  }

  void linear_solve(const ROL::Ptr<Thyra::VectorBase<Real>> &u_ptr,
                    std::vector<Real> &d,
                    std::vector<Real> &o,
                    const Thyra::VectorBase<Real> &r_vb) const {
    Thyra::copy(r_vb, u_ptr.ptr());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int nx = static_cast<int>(nx_);
    int info;
    int ldb  = nx;
    int nhrs = 1;
    lp.PTTRF(nx,&d[0],&o[0],&info);
    { // scope for u eval
      Thyra::DetachedVectorView<Real> u(u_ptr);
      lp.PTTRS(nx,nhrs,&d[0],&o[0],u.values(),ldb,&info);
    } // end scope for u eval
  }

  void run_newton(const ROL::Ptr<Thyra::VectorBase<Real>> &u,
                  const Thyra::VectorBase<Real> &up,
                  const Thyra::VectorBase<Real> &z) const {
    const Real half(0.5), one(1), lstol(1e-4);
    // Set initial guess
    Thyra::copy(up, u.ptr());
    // Compute residual and residual norm
    ROL::Ptr<Thyra::VectorBase<Real>> r = u->clone_v();
    Thyra::put_scalar(0.0, r.ptr());
    compute_residual(r,up,*u,z);
    Real rnorm = compute_norm(*r);
    // Define tolerances
    Real tol   = static_cast<Real>(1e2)*ROL::ROL_EPSILON<Real>();
    Real maxit = 100;
    // Initialize Jacobian storage
    std::vector<Real> d(nx_);
    std::vector<Real> o(nx_-1);
    // Iterate Newton's method
    Real alpha(1), tmp(0);
    ROL::Ptr<Thyra::VectorBase<Real>> s = u->clone_v();
    Thyra::put_scalar(0.0, s.ptr());
    ROL::Ptr<Thyra::VectorBase<Real>> utmp = u->clone_v();
    Thyra::put_scalar(0.0, utmp.ptr());
    for (uint i = 0; i < maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      compute_pde_jacobian(d,o,*u);
      // Solve Newton system
      linear_solve(s,d,o,*r);
      // Perform line search
      tmp = rnorm;
      alpha = one;
      Thyra::copy(*u, utmp.ptr());
      Thyra::Vp_StV(utmp.ptr(), -alpha, *s);
      compute_residual(r,up,*utmp,z);
      rnorm = compute_norm(*r);
      while ( rnorm > (one-lstol*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON<Real>()) ) {
        alpha *= half;
        Thyra::copy(*u, utmp.ptr());
        Thyra::Vp_StV(utmp.ptr(), -alpha, *s);
        compute_residual(r,up,*utmp,z);
        rnorm = compute_norm(*r);
      }
      // Update iterate
      Thyra::copy(*utmp, u.ptr());
      if ( rnorm < tol ) {
        break;
      }
    }
  }

};



/// 2. Linear Operators
/// --------------------------
/// These classes provide implementations of LinearOp and LinearOpWithSolve
/// objects required by the ModelEvaluator.  The goal is to reuse or wrap
/// the operations already provided by the ParabolicModel class.
/// --------------------------
template<class Real>
class ParabolicModelWopForward : virtual public Thyra::LinearOpWithSolveBase<Real> {
private:
  const ROL::Ptr<ParabolicModel<Real>> model_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> space_;
  ROL::Ptr<const Thyra::VectorBase<Real>> u_;
  Real alpha_;
  Real beta_;

public:

  ParabolicModelWopForward(const ROL::Ptr<ParabolicModel<Real>> &model) :
    model_(model),
    space_(model->getStateSpace()) {}

  void set_u(const ROL::Ptr<const Thyra::VectorBase<Real>> &u) {
    u_ = u;
  }

  void set_alpha(const Real &alpha) {
    alpha_ = alpha;
  }

  void set_beta(const Real &beta) {
    beta_ = beta;
  }

  /** @name Overridden from LinearOpBase */
  //@{

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> range() const { return space_; }

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> domain() const { return space_; }

  ROL::Ptr<const Thyra::LinearOpBase<Real>> clone() const {
    ROL::Ptr<ParabolicModelWopForward<Real>> op =
      ROL::makePtr<ParabolicModelWopForward<Real>>(model_);
    return op;
  }
  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const {
    if (M_trans == Thyra::NOTRANS) {
      return true;
    }
    else {
      return false;
    }
  }

  void applyImpl(const Thyra::EOpTransp M_trans,
                 const Thyra::MultiVectorBase<Real>& X,
                 const Teuchos::Ptr<Thyra::MultiVectorBase<Real> >& Y,
                 const Real a,
                 const Real b) const {
    // assert(opSupportedImpl(M_trans));
    ROL::Ptr<const Thyra::VectorBase<Real>> x = X.col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> y = Y->col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> Mx = y->clone_v();
    ROL::Ptr<Thyra::VectorBase<Real>> cudotx = y->clone_v();
    Thyra::put_scalar(0.0, Mx.ptr());
    Thyra::put_scalar(0.0, cudotx.ptr());

    model_->apply_c_u(Mx, x, u_);
    Thyra::scale(beta_, Mx.ptr());
    model_->apply_c_udot(cudotx, x);
    Thyra::Vp_StV(Mx.ptr(), alpha_, *cudotx);  // Mx = alpha_*c_udot*x _ beta_*c_u*x
    Thyra::scale(b, y.ptr());
    Thyra::Vp_StV(y.ptr(), a, *Mx);  // y = a*Mx + b*y
  }
  //@}

  /** @name Overridden from LinearOpWithSolveBase */
  //@{
  bool solveSupportsImpl(Thyra::EOpTransp M_trans) const {
    if (M_trans == Thyra::NOTRANS) {
      return true;
    }
    else {
      return false;
    }
  }

  bool solveSupportsNewImpl(
    Thyra::EOpTransp M_trans,
    const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria) const
  { return false; }

  bool solveSupportsSolveMeasureTypeImpl(
    Thyra::EOpTransp M_trans,
    const Thyra::SolveMeasureType &solveMeasureType) const
  { return false; }

  Thyra::SolveStatus<Real> solveImpl(const Thyra::EOpTransp M_trans,
                                     const Thyra::MultiVectorBase<Real>& B,
                                     const Teuchos::Ptr<Thyra::MultiVectorBase<Real> >& X,
                                     const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria) const
  {
    assert(solveSupportsImpl(M_trans));
    ROL::Ptr<const Thyra::VectorBase<Real>> b = B.col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> x = X->col(0);

    // Initialize Jacobian storage.
    int nx = space_->dim();
    std::vector<Real> d(nx);
    std::vector<Real> o(nx-1);
    // Get Jacobian.
    model_->compute_alpha_c_udot_beta_c_u(d,o,*u_, alpha_,beta_);
    // Solve system.
    model_->linear_solve(x,d,o,*b);

    // Return status.
    Thyra::SolveStatus<Real> solveStatus;
    solveStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
    return solveStatus;
  }
  //@}

};


template<class Real>
class ParabolicModelWopAdjoint : virtual public Thyra::LinearOpWithSolveBase<Real> {
private:
  const ROL::Ptr<ParabolicModel<Real>> model_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> space_;
  ROL::Ptr<const Thyra::VectorBase<Real>> u_;
  Real alpha_;
  Real beta_;

public:

  ParabolicModelWopAdjoint(const ROL::Ptr<ParabolicModel<Real>> &model) :
    model_(model),
    space_(model->getStateSpace()) {}

  void set_u(const ROL::Ptr<const Thyra::VectorBase<Real>> &u) {
    u_ = u;
  }

  void set_alpha(const Real &alpha) {
    alpha_ = alpha;
  }

  void set_beta(const Real &beta) {
    beta_ = beta;
  }

  /** @name Overridden from LinearOpBase */
  //@{

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> range() const { return space_; }

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> domain() const { return space_; }

  ROL::Ptr<const Thyra::LinearOpBase<Real>> clone() const {
    ROL::Ptr<ParabolicModelWopAdjoint<Real>> op =
      ROL::makePtr<ParabolicModelWopAdjoint<Real>>(model_);
    return op;
  }
  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const {
    if (M_trans == Thyra::NOTRANS) {
      return true;
    }
    else {
      return false;
    }
  }

  void applyImpl(const Thyra::EOpTransp M_trans,
                 const Thyra::MultiVectorBase<Real>& X,
                 const Teuchos::Ptr<Thyra::MultiVectorBase<Real> >& Y,
                 const Real a,
                 const Real b) const {
    assert(opSupportedImpl(M_trans));
    ROL::Ptr<const Thyra::VectorBase<Real>> x = X.col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> y = Y->col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> Mx = y->clone_v();
    ROL::Ptr<Thyra::VectorBase<Real>> cudotx = y->clone_v();
    Thyra::put_scalar(0.0, Mx.ptr());
    Thyra::put_scalar(0.0, cudotx.ptr());

    model_->apply_c_u(Mx, x, u_);
    Thyra::scale(beta_, Mx.ptr());
    model_->apply_c_udot(cudotx, x);
    Thyra::Vp_StV(Mx.ptr(), alpha_, *cudotx);  // Mx = alpha_*c_udot*x _ beta_*c_u*x
    Thyra::scale(b, y.ptr());
    Thyra::Vp_StV(y.ptr(), a, *Mx);  // y = a*Mx + b*y
  }
  //@}

  /** @name Overridden from LinearOpWithSolveBase */
  //@{
  bool solveSupportsImpl(Thyra::EOpTransp M_trans) const {
    if (M_trans == Thyra::NOTRANS) {
      return true;
    }
    else {
      return false;
    }
  }

  bool solveSupportsNewImpl(
    Thyra::EOpTransp M_trans,
    const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria) const
  { return false; }

  bool solveSupportsSolveMeasureTypeImpl(
    Thyra::EOpTransp M_trans,
    const Thyra::SolveMeasureType &solveMeasureType) const
  { return false; }

  Thyra::SolveStatus<Real> solveImpl(const Thyra::EOpTransp M_trans,
                                     const Thyra::MultiVectorBase<Real>& B,
                                     const Teuchos::Ptr<Thyra::MultiVectorBase<Real> >& X,
                                     const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria) const
  {
    assert(solveSupportsImpl(M_trans));
    ROL::Ptr<const Thyra::VectorBase<Real>> b = B.col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> x = X->col(0);

    // Initialize Jacobian storage.
    int nx = space_->dim();
    std::vector<Real> d(nx);
    std::vector<Real> o(nx-1);
    // Get Jacobian.
    model_->compute_alpha_c_udot_beta_c_u(d,o,*u_, alpha_,beta_);
    // Solve system.
    model_->linear_solve(x,d,o,*b);

    // Return status.
    Thyra::SolveStatus<Real> solveStatus;
    solveStatus.solveStatus = Thyra::SOLVE_STATUS_CONVERGED;
    return solveStatus;
  }
  //@}

};


template<class Real>
class ParabolicModelDfDpop : virtual public Thyra::LinearOpBase<Real> {
private:
  const ROL::Ptr<ParabolicModel<Real>> model_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> controlspace_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> constraintspace_;

public:

  ParabolicModelDfDpop(const ROL::Ptr<ParabolicModel<Real>> &model) :
    model_(model),
    controlspace_(model->getControlSpace()),
    constraintspace_(model->getConstraintSpace()) {}

  /** @name Overridden from LinearOpBase */
  //@{

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> range() const { return constraintspace_; }

  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> domain() const { return controlspace_; }

  ROL::Ptr<const Thyra::LinearOpBase<Real>> clone() const {
    ROL::Ptr<ParabolicModelDfDpop<Real>> op =
      ROL::makePtr<ParabolicModelDfDpop<Real>>(model_);
    return op;
  }
  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const { return true; }

  void applyImpl(const Thyra::EOpTransp M_trans,
                 const Thyra::MultiVectorBase<Real>& X,
                 const Teuchos::Ptr<Thyra::MultiVectorBase<Real> >& Y,
                 const Real a,
                 const Real b) const {
    assert(opSupportedImpl(M_trans));
    // Note:  The operator is non-square, so we have to distinguish
    // between two cases.
    ROL::Ptr<const Thyra::VectorBase<Real>> x = X.col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> y = Y->col(0);
    ROL::Ptr<Thyra::VectorBase<Real>> Mx = y->clone_v();

    if (M_trans == Thyra::NOTRANS) {
      model_->apply_c_z(Mx, x, false);
    } else if (M_trans == Thyra::TRANS) {
      model_->apply_c_z(Mx, x, true);
    }
    // Mx = c_z*x

    Thyra::scale(b, y.ptr());
    Thyra::Vp_StV(y.ptr(), a, *Mx); // y = a*Mx + b*y
  }
  //@}

};



/// 3. ParabolicModelMEWrapper
/// --------------------------
/// This is the Thyra::ModelEvaluator wrapper for the ParabolicModel class.
/// The purpose of this class is to provide the implementation of an
/// interface that is compatible with the Tempus::Integrator.
/// --------------------------

// Forward wrapper.
template<class Real>
class ParabolicModelMEWrapperForward : public Thyra::StateFuncModelEvaluatorBase<Real> {

public: // ParabolicModelMEWrapperForward methods.

  ParabolicModelMEWrapperForward(const ROL::Ptr<ParabolicModel<Real>> &model);
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> get_x_space() const;
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> get_f_space() const;
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> get_p_space(int l) const;
  Thyra::ModelEvaluatorBase ::InArgs<Real>     getNominalValues() const;
  Thyra::ModelEvaluatorBase::InArgs<Real>      createInArgs() const;
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> create_W() const;
  ROL::Ptr<Thyra::LinearOpBase<Real>>          create_W_op() const;
  ROL::Ptr<Thyra::LinearOpBase<Real>>          create_DfDp_op_impl(int l) const;


private: // ParabolicModelMEWrapperForward methods.

  void setupInOutArgs_() const;

  // methods overridden from the base class
  Thyra::ModelEvaluatorBase::OutArgs<Real> createOutArgsImpl() const;
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Real>  &inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<Real> &outArgs) const;

private: // ParabolicModelMEWrapperForward data members.

  const ROL::Ptr<ParabolicModel<Real>> model_;

  const int Np_;

  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> x_space_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> f_space_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> p_space_;

  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Real>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Real>  nominalValues_;

};

/// ------------------------------------------
/// Implementation of ParabolicModelMEWrapperForward.
/// ------------------------------------------

template<class Real>
ParabolicModelMEWrapperForward<Real>::ParabolicModelMEWrapperForward(const ROL::Ptr<ParabolicModel<Real>> &model) :
  model_(model),
  Np_(1),
  x_space_(model->getStateSpace()),
  f_space_(model->getConstraintSpace()),
  p_space_(model->getControlSpace()) {
  isInitialized_ = false;
}


template<class Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> ParabolicModelMEWrapperForward<Real>::get_x_space() const {
  return x_space_;
}


template<class Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> ParabolicModelMEWrapperForward<Real>::get_f_space() const {
  return f_space_;
}


template<class Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> ParabolicModelMEWrapperForward<Real>::get_p_space(int l) const {
  return p_space_;
}

template<class Real>
Thyra::ModelEvaluatorBase::InArgs<Real> ParabolicModelMEWrapperForward<Real>::getNominalValues() const {
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template<class Real>
Thyra::ModelEvaluatorBase::InArgs<Real> ParabolicModelMEWrapperForward<Real>::createInArgs() const {
  setupInOutArgs_();
  return inArgs_;
}

template<class Real>
Thyra::ModelEvaluatorBase::OutArgs<Real> ParabolicModelMEWrapperForward<Real>::createOutArgsImpl() const {
  setupInOutArgs_();
  return outArgs_;
}

template<class Real>
ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> ParabolicModelMEWrapperForward<Real>::create_W() const {
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> W =
    ROL::makePtr<ParabolicModelWopForward<Real>>(model_);
  return W;
}

template<class Real>
ROL::Ptr<Thyra::LinearOpBase<Real>> ParabolicModelMEWrapperForward<Real>::create_W_op() const {
  ROL::Ptr<Thyra::LinearOpBase<Real>> W_op =
    ROL::makePtr<ParabolicModelWopForward<Real>>(model_);
  return W_op;
}

template<class Real>
ROL::Ptr<Thyra::LinearOpBase<Real>> ParabolicModelMEWrapperForward<Real>::create_DfDp_op_impl(int l) const {
  ROL::Ptr<Thyra::LinearOpBase<Real>> DfDp_op =
    ROL::makePtr<ParabolicModelDfDpop<Real>>(model_);
  return DfDp_op;
}

template<class Real>
void ParabolicModelMEWrapperForward<Real>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Real> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Real> &outArgs) const {

  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");

  ROL::Ptr<Thyra::VectorBase<Real> >         r = outArgs.get_f();
  ROL::Ptr<const Thyra::VectorBase<Real>> udot = inArgs.get_x_dot();
  ROL::Ptr<const Thyra::VectorBase<Real>>    u = inArgs.get_x();
  ROL::Ptr<const Thyra::VectorBase<Real>>    z = inArgs.get_p(0);

  if (!is_null(r)) {
    model_->compute_residual(r, *udot, *u, *z);
  }

  Real alpha = inArgs.get_alpha();
  Real beta  = inArgs.get_beta();
  const ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> W = outArgs.get_W();
  const ROL::Ptr<Thyra::LinearOpBase<Real>> W_op = outArgs.get_W_op();
  //const ROL::Ptr<Thyra::LinearOpBase<Real>> DfDp_op = outArgs.get_DfDp_op(0).getLinearOp(); // not needed in this evaluator

  if (!is_null(W)) {
    (ROL::dynamicPtrCast<ParabolicModelWopForward<Real>>(W))->set_u(u);
    (ROL::dynamicPtrCast<ParabolicModelWopForward<Real>>(W))->set_alpha(alpha);
    (ROL::dynamicPtrCast<ParabolicModelWopForward<Real>>(W))->set_beta(beta);
  }

  if (!is_null(W_op)) {
    (ROL::dynamicPtrCast<ParabolicModelWopForward<Real>>(W_op))->set_u(u);
    (ROL::dynamicPtrCast<ParabolicModelWopForward<Real>>(W_op))->set_alpha(alpha);
    (ROL::dynamicPtrCast<ParabolicModelWopForward<Real>>(W_op))->set_beta(beta);
  }

  //if (!is_null(DfDp_op)) { // not needed in this evaluator
  //  (ROL::dynamicPtrCast<ParabolicModelWop<Real>>(DfDp_op))->set_u(u);
  //}

}

template<class Real>
void ParabolicModelMEWrapperForward<Real>::setupInOutArgs_() const {
  if (isInitialized_) {
    return;
  }

  { // Set up prototypical InArgs.
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_t );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_beta );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x_dot );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_alpha );
    inArgs.set_Np(Np_);
    inArgs_ = inArgs;
  }

  { // Set up prototypical OutArgs,
    Thyra::ModelEvaluatorBase::OutArgsSetup<Real> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.set_Np_Ng(Np_,0);
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f );
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_W_op );
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_W );
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,0,
                         Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP );
    outArgs_ = outArgs;
  }

  // Set up nominal values.
  nominalValues_ = inArgs_;
  const Teuchos::RCP<Thyra::VectorBase<Real>> x_ic = createMember(x_space_);
  Thyra::put_scalar(0.0, x_ic.ptr());
  nominalValues_.set_x(x_ic);
  nominalValues_.set_x_dot(x_ic);

  isInitialized_ = true;
}

// Adjoint wrapper.
template<class Real>
class ParabolicModelMEWrapperAdjoint : public Thyra::StateFuncModelEvaluatorBase<Real> {

public: // ParabolicModelMEWrapperAdjoint methods.

  ParabolicModelMEWrapperAdjoint(const ROL::Ptr<ParabolicModel<Real>> &model);
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> get_x_space() const;
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> get_f_space() const;
  ROL::Ptr<const Thyra::VectorSpaceBase<Real>> get_p_space(int l) const;
  Thyra::ModelEvaluatorBase ::InArgs<Real>     getNominalValues() const;
  Thyra::ModelEvaluatorBase::InArgs<Real>      createInArgs() const;
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> create_W() const;
  ROL::Ptr<Thyra::LinearOpBase<Real>>          create_W_op() const;
  //ROL::Ptr<Thyra::LinearOpBase<Real>>          create_DfDp_op_impl(int l) const;


private: // ParabolicModelMEWrapperAdjoint methods.

  void setupInOutArgs_() const;

  // methods overridden from the base class
  Thyra::ModelEvaluatorBase::OutArgs<Real> createOutArgsImpl() const;
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Real>  &inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<Real> &outArgs) const;

private: // ParabolicModelMEWrapperAdjoint data members.

  const ROL::Ptr<ParabolicModel<Real>> model_;

  const int Np_;

  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> x_space_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> f_space_;
  const ROL::Ptr<const Thyra::VectorSpaceBase<Real>> p_space_;

  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Real>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Real>  nominalValues_;

};

/// ------------------------------------------
/// Implementation of ParabolicModelMEWrapperAdjoint.
/// ------------------------------------------

template<class Real>
ParabolicModelMEWrapperAdjoint<Real>::ParabolicModelMEWrapperAdjoint(const ROL::Ptr<ParabolicModel<Real>> &model) :
  model_(model),
  Np_(1),
  x_space_(model->getStateSpace()),
  f_space_(model->getConstraintSpace()),
  p_space_(model->getControlSpace()) {
  isInitialized_ = false;
}


template<class Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> ParabolicModelMEWrapperAdjoint<Real>::get_x_space() const {
  return x_space_;
}


template<class Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> ParabolicModelMEWrapperAdjoint<Real>::get_f_space() const {
  return f_space_;
}


template<class Real>
Teuchos::RCP<const Thyra::VectorSpaceBase<Real>> ParabolicModelMEWrapperAdjoint<Real>::get_p_space(int l) const {
  return p_space_;
}

template<class Real>
Thyra::ModelEvaluatorBase::InArgs<Real> ParabolicModelMEWrapperAdjoint<Real>::getNominalValues() const {
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template<class Real>
Thyra::ModelEvaluatorBase::InArgs<Real> ParabolicModelMEWrapperAdjoint<Real>::createInArgs() const {
  setupInOutArgs_();
  return inArgs_;
}

template<class Real>
Thyra::ModelEvaluatorBase::OutArgs<Real> ParabolicModelMEWrapperAdjoint<Real>::createOutArgsImpl() const {
  setupInOutArgs_();
  return outArgs_;
}

template<class Real>
ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> ParabolicModelMEWrapperAdjoint<Real>::create_W() const {
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> W =
    ROL::makePtr<ParabolicModelWopAdjoint<Real>>(model_);
  return W;
}

template<class Real>
ROL::Ptr<Thyra::LinearOpBase<Real>> ParabolicModelMEWrapperAdjoint<Real>::create_W_op() const {
  ROL::Ptr<Thyra::LinearOpBase<Real>> W_op =
    ROL::makePtr<ParabolicModelWopAdjoint<Real>>(model_);
  return W_op;
}

//template<class Real>
//ROL::Ptr<Thyra::LinearOpBase<Real>> ParabolicModelMEWrapperAdjoint<Real>::create_DfDp_op_impl(int l) const {
//  ROL::Ptr<Thyra::LinearOpBase<Real>> DfDp_op =
//    ROL::makePtr<ParabolicModelDfDpop<Real>>(model_);
//  return DfDp_op;
//}

template<class Real>
void ParabolicModelMEWrapperAdjoint<Real>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Real> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Real> &outArgs) const {

  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");

  //ROL::Ptr<Thyra::VectorBase<Real> >         r = outArgs.get_f();
  ROL::Ptr<const Thyra::VectorBase<Real>> udot = inArgs.get_x_dot();
  ROL::Ptr<const Thyra::VectorBase<Real>>    u = inArgs.get_x();
  ROL::Ptr<const Thyra::VectorBase<Real>>    z = inArgs.get_p(0);

  //if (!is_null(r)) {
  //  model_->compute_residual(r, *udot, *u, *z);
  //}

  Real alpha = inArgs.get_alpha();
  Real beta  = inArgs.get_beta();
  const ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> W = outArgs.get_W();
  const ROL::Ptr<Thyra::LinearOpBase<Real>> W_op = outArgs.get_W_op();
  //const ROL::Ptr<Thyra::LinearOpBase<Real>> DfDp_op = outArgs.get_DfDp_op(0).getLinearOp(); // not needed in this evaluator

  if (!is_null(W)) {
    (ROL::dynamicPtrCast<ParabolicModelWopAdjoint<Real>>(W))->set_u(u);
    (ROL::dynamicPtrCast<ParabolicModelWopAdjoint<Real>>(W))->set_alpha(alpha);
    (ROL::dynamicPtrCast<ParabolicModelWopAdjoint<Real>>(W))->set_beta(beta);
  }

  if (!is_null(W_op)) {
    (ROL::dynamicPtrCast<ParabolicModelWopAdjoint<Real>>(W_op))->set_u(u);
    (ROL::dynamicPtrCast<ParabolicModelWopAdjoint<Real>>(W_op))->set_alpha(alpha);
    (ROL::dynamicPtrCast<ParabolicModelWopAdjoint<Real>>(W_op))->set_beta(beta);
  }

  //if (!is_null(DfDp_op)) { // not needed in this evaluator
  //  (ROL::dynamicPtrCast<ParabolicModelWop<Real>>(DfDp_op))->set_u(u);
  //}

}

template<class Real>
void ParabolicModelMEWrapperAdjoint<Real>::setupInOutArgs_() const {
  if (isInitialized_) {
    return;
  }

  { // Set up prototypical InArgs.
    Thyra::ModelEvaluatorBase::InArgsSetup<Real> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_t );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_beta );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x_dot );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_alpha );
    inArgs.set_Np(Np_);
    inArgs_ = inArgs;
  }

  { // Set up prototypical OutArgs,
    Thyra::ModelEvaluatorBase::OutArgsSetup<Real> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.set_Np_Ng(Np_,0);
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f );
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_W_op );
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_W );
    //outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,0,
    //                     Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP );
    outArgs_ = outArgs;
  }

  // Set up nominal values.
  nominalValues_ = inArgs_;
  const Teuchos::RCP<Thyra::VectorBase<Real>> x_ic = createMember(x_space_);
  Thyra::put_scalar(0.0, x_ic.ptr());
  nominalValues_.set_x(x_ic);
  nominalValues_.set_x_dot(x_ic);

  isInitialized_ = true;
}




/***********************************************************************
****************************  OBJECTIVE  *******************************
************************************************************************/

template<class Real>
class Objective_ParabolicControl : public ROL::DynamicObjective<Real> {

  typedef typename std::vector<Real>::size_type uint;

private:
  Real alpha_;
  Real theta_;
  uint nx_;
  Real dx_;
  Real dt_;
  int  type_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  Real dot(const Thyra::VectorBase<Real> &x_vb,
           const Thyra::VectorBase<Real> &y_vb) const {
    const Real two(2), four(4), six(6);
    const Thyra::ConstDetachedVectorView<Real> x(x_vb);
    const Thyra::ConstDetachedVectorView<Real> y(y_vb);
    Real ip(0);
    Real c = ((static_cast<uint>(x.subDim())==nx_) ? four : two);
    for (uint i=0; i<static_cast<uint>(x.subDim()); i++) {
      if ( i == 0 ) {
        ip += dx_/six*(c*x[i] + x[i+1])*y[i];
      }
      else if ( i == static_cast<uint>(x.subDim())-1 ) {
        ip += dx_/six*(x[i-1] + c*x[i])*y[i];
      }
      else {
        ip += dx_/six*(x[i-1] + four*x[i] + x[i+1])*y[i];
      }
    }
    return ip;
  }

  void apply_mass(ROL::Ptr<Thyra::VectorBase<Real>> &Mu_ptr,
                  const Thyra::VectorBase<Real> &u_vb) const {
    const Real two(2), four(4), six(6);
    Thyra::put_scalar(0.0, Mu_ptr.ptr());
    const Thyra::ConstDetachedVectorView<Real> u(u_vb);
    Real c = ((static_cast<uint>(u.subDim())==nx_) ? four : two);
    { // scope for Mu eval
      Thyra::DetachedVectorView<Real> Mu(Mu_ptr);
      for (uint i=0; i<static_cast<uint>(u.subDim()); i++) {
        if ( i == 0 ) {
          Mu[i] = dx_/six*(c*u[i] + u[i+1]);
        }
        else if ( i == static_cast<uint>(u.subDim())-1 ) {
          Mu[i] = dx_/six*(u[i-1] + c*u[i]);
        }
        else {
          Mu[i] = dx_/six*(u[i-1] + four*u[i] + u[i+1]);
        }
      }
    } // end scope for Mu eval
  }

  Real evaluate_target(Real x) const {
    const Real zero(0), half(0.5), eight(8), pi(M_PI);
    Real val(0);
    switch (type_) {
      case 1:  val = (x<half ? half : zero);                 break;
      case 2:  val = half;                                   break;
      case 3:  val = half*std::abs(std::sin(eight*pi*x));    break;
      case 4:  val = half*std::exp(-half*(x-half)*(x-half)); break;
      default: val = half;                                   break;
    }
    return val;
  }

  ROL::Ptr<const Thyra::VectorBase<Real>> getVector( const ROL::Vector<Real>& x ) const {
    return dynamic_cast<const ROL::ThyraVector<Real>&>(x).getVector();
  }

  ROL::Ptr<Thyra::VectorBase<Real>> getVector( ROL::Vector<Real>& x ) const {
    return dynamic_cast<ROL::ThyraVector<Real>&>(x).getVector();
  }
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_ParabolicControl(ROL::ParameterList &pl) {
    const Real one(1);
    uint nt   = pl.get("Temporal Discretization",  100);
    Real  T   = pl.get("End Time",                 1.0);
    type_     = pl.get("Target Type",                1);
    alpha_    = pl.get("Control Penalty",         1e-4);
    nx_       = pl.get("Spatial Discretization",   128);
    theta_    = pl.get("Integration Factor",       0.5);
    dx_       = one/(static_cast<Real>(nx_)-one);
    dt_       = T/(static_cast<Real>(nt)-one);
  }

  Real value(const ROL::Vector<Real> &uold,
             const ROL::Vector<Real> &unew,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts ) const {
    const Real half(0.5), one(1);
    ROL::Ptr<const Thyra::VectorBase<Real>>  z_ptr = getVector(z);
    ROL::Ptr<const Thyra::VectorBase<Real>> uo_ptr = getVector(uold);
    ROL::Ptr<const Thyra::VectorBase<Real>> un_ptr = getVector(unew);
    const Thyra::ConstDetachedVectorView<Real> uo(uo_ptr);
    const Thyra::ConstDetachedVectorView<Real> un(un_ptr);
    // Compute Norm Squared of State
    Real target(0);
    ROL::Ptr<Thyra::VectorBase<Real>> uOld_ptr = uo_ptr->clone_v();
    ROL::Ptr<Thyra::VectorBase<Real>> uNew_ptr = un_ptr->clone_v();
    Thyra::put_scalar(0.0, uOld_ptr.ptr());
    Thyra::put_scalar(0.0, uNew_ptr.ptr());
    { // scope for uOld, uNew eval
      Thyra::DetachedVectorView<Real> uOld(uOld_ptr);
      Thyra::DetachedVectorView<Real> uNew(uNew_ptr);
      for (uint n = 0; n < nx_; n++) {
        target = evaluate_target(static_cast<Real>(n)*dx_);
        uOld[n] = uo[n] - target;
        uNew[n] = un[n] - target;
      }
    } // end scope for uOld, uNew eval
    Real uoval = dot(*uOld_ptr, *uOld_ptr);
    Real unval = dot(*uNew_ptr, *uNew_ptr);
    // Add Norm Squared of Control
    Real zval = dot(*z_ptr, *z_ptr);
    return half * dt_ * (theta_*uoval + (one-theta_)*unval + alpha_ * zval);
  }

  void gradient_uo(ROL::Vector<Real> &g,
                   const ROL::Vector<Real> &uold,
                   const ROL::Vector<Real> &unew,
                   const ROL::Vector<Real> &z,
                   const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<Thyra::VectorBase<Real>>        g_ptr = getVector(g);
    ROL::Ptr<const Thyra::VectorBase<Real>> uo_ptr = getVector(uold);
    const Thyra::ConstDetachedVectorView<Real> uo(uo_ptr);
    ROL::Ptr<Thyra::VectorBase<Real>> uDiff_ptr = uo_ptr->clone_v();
    Thyra::put_scalar(0.0, uDiff_ptr.ptr());
    Real target(0);
    { // scope for uDiff
      Thyra::DetachedVectorView<Real> uDiff(uDiff_ptr);
      for (uint n = 0; n < nx_; n++) {
        target = evaluate_target(static_cast<Real>(n)*dx_);
        uDiff[n] = uo[n] - target;
      }
    } // end scope for uDiff
    apply_mass(g_ptr, *uDiff_ptr);
    g.scale(theta_*dt_);
  }

  void gradient_un(ROL::Vector<Real> &g,
                   const ROL::Vector<Real> &uold,
                   const ROL::Vector<Real> &unew,
                   const ROL::Vector<Real> &z,
                   const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    ROL::Ptr<Thyra::VectorBase<Real>>        g_ptr = getVector(g);
    ROL::Ptr<const Thyra::VectorBase<Real>> un_ptr = getVector(unew);
    const Thyra::ConstDetachedVectorView<Real> un(un_ptr);
    ROL::Ptr<Thyra::VectorBase<Real>> uDiff_ptr = un_ptr->clone_v();
    Thyra::put_scalar(0.0, uDiff_ptr.ptr());
    Real target(0);
    { // scope for uDiff eval
      Thyra::DetachedVectorView<Real> uDiff(uDiff_ptr);
      for (uint n = 0; n < nx_; n++) {
        target = evaluate_target(static_cast<Real>(n)*dx_);
        uDiff[n] = un[n] - target;
      }
    } // end scope for uDiff eval
    apply_mass(g_ptr, *uDiff_ptr);
    g.scale((one-theta_)*dt_);
  }

  void gradient_z(ROL::Vector<Real> &g,
                  const ROL::Vector<Real> &uold,
                  const ROL::Vector<Real> &unew,
                  const ROL::Vector<Real> &z,
                  const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       g_ptr = getVector(g);
    ROL::Ptr<const Thyra::VectorBase<Real>> z_ptr = getVector(z);
    apply_mass(g_ptr, *z_ptr);
    g.scale(dt_*alpha_);
  }

  void hessVec_uo_uo(ROL::Vector<Real> &hv,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &uold,
                     const ROL::Vector<Real> &unew,
                     const ROL::Vector<Real> &z,
                     const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       hv_ptr = getVector(hv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    apply_mass(hv_ptr, *v_ptr);
    hv.scale(theta_*dt_);
  }

  void hessVec_uo_un(ROL::Vector<Real> &hv,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &uold,
                     const ROL::Vector<Real> &unew,
                     const ROL::Vector<Real> &z,
                     const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_uo_z(ROL::Vector<Real> &hv,
                    const ROL::Vector<Real> &v,
                    const ROL::Vector<Real> &uold,
                    const ROL::Vector<Real> &unew,
                    const ROL::Vector<Real> &z,
                    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_un_uo(ROL::Vector<Real> &hv,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &uold,
                     const ROL::Vector<Real> &unew,
                     const ROL::Vector<Real> &z,
                     const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_un_un(ROL::Vector<Real> &hv,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &uold,
                     const ROL::Vector<Real> &unew,
                     const ROL::Vector<Real> &z,
                     const ROL::TimeStamp<Real> &ts ) const {
    const Real one(1);
    ROL::Ptr<Thyra::VectorBase<Real>>       hv_ptr = getVector(hv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    apply_mass(hv_ptr, *v_ptr);
    hv.scale((one-theta_)*dt_);
  }

  void hessVec_un_z(ROL::Vector<Real> &hv,
                    const ROL::Vector<Real> &v,
                    const ROL::Vector<Real> &uold,
                    const ROL::Vector<Real> &unew,
                    const ROL::Vector<Real> &z,
                    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_z_uo(ROL::Vector<Real> &hv,
                    const ROL::Vector<Real> &v,
                    const ROL::Vector<Real> &uold,
                    const ROL::Vector<Real> &unew,
                    const ROL::Vector<Real> &z,
                    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_z_un(ROL::Vector<Real> &hv,
                    const ROL::Vector<Real> &v,
                    const ROL::Vector<Real> &uold,
                    const ROL::Vector<Real> &unew,
                    const ROL::Vector<Real> &z,
                    const ROL::TimeStamp<Real> &ts ) const {
    hv.zero();
  }

  void hessVec_z_z(ROL::Vector<Real> &hv,
                   const ROL::Vector<Real> &v,
                   const ROL::Vector<Real> &uold,
                   const ROL::Vector<Real> &unew,
                   const ROL::Vector<Real> &z,
                   const ROL::TimeStamp<Real> &ts ) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       hv_ptr = getVector(hv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    apply_mass(hv_ptr, *v_ptr);
    hv.scale(dt_*alpha_);
  }
};
