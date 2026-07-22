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
********************  PARABOLIC MODEL EVALUATOR  ***********************
************************************************************************/
template<class Scalar>
class ParabolicModel : public Thyra::StateFuncModelEvaluatorBase<Scalar> {
public:

  ParabolicModel();
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_f_space() const;
  Thyra::ModelEvaluatorBase ::InArgs<Scalar>         getNominalValues() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar>          createInArgs() const;

private:

  void setupInOutArgs_() const;

  // methods overridden from base class
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

private: // ParabolicModel data members

  int dim_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> f_space_;
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;

};

/********************
** IMPLEMENTATION. **
********************/

template<class Scalar>
ParabolicModel<Scalar>::ParabolicModel() {
  isInitialized_ = false;
  dim_ = 2;
  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> ParabolicModel<Scalar>::get_x_space() const {
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> ParabolicModel<Scalar>::get_f_space() const {
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ParabolicModel<Scalar>::getNominalValues() const {
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ParabolicModel<Scalar>::createInArgs() const {
  setupInOutArgs_();
  return inArgs_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> ParabolicModel<Scalar>::createOutArgsImpl() const {
  setupInOutArgs_();
  return outArgs_;
}


template<class Scalar>
void ParabolicModel<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const {

  using Teuchos::RCP;
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");

  const Thyra::ConstDetachedVectorView<Scalar> x_in_view(inArgs.get_x());
  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();

  // Evaluate the Explicit ODE f(x,t) [= 0]
  if (!is_null(f_out)) {
    Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
    f_out_view[0] =  x_in_view[1];
    f_out_view[1] = -x_in_view[0];
  }

}


template<class Scalar>
void ParabolicModel<Scalar>::setupInOutArgs_() const {
  if (isInitialized_) {
    return;
  }

  { // Set up prototypical InArgs.
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_t );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x );
    inArgs_ = inArgs;
  }

  { // Set up prototypical OutArgs,
    Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f );
    outArgs_ = outArgs;
  }

  // Set up nominal values.
  nominalValues_ = inArgs_;
  const Teuchos::RCP<Thyra::VectorBase<Scalar>> x_ic = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
    x_ic_view[0] = 1.0;
    x_ic_view[1] = 1.0;
  }
  nominalValues_.set_x(x_ic);

  isInitialized_ = true;
}



/***********************************************************************
****************************  CONSTRAINT  *******************************
************************************************************************/
template<class Real>
class Constraint_ParabolicControl : public ROL::DynamicConstraint<Real> {

  typedef typename std::vector<Real>::size_type uint;

private:
  Real eps1_;
  Real eps2_;
  uint nx_;
  Real dx_;
  Real dt_;
  int  type_;

  std::vector<Real> pts_;
  std::vector<Real> wts_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  void apply_mass(ROL::Ptr<Thyra::VectorBase<Real>> &Mu_ptr, ROL::Ptr<const Thyra::VectorBase<Real>> &u_ptr ) const {
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

  void compute_pde_jacobian(std::vector<Real> &d, std::vector<Real> &o, const Thyra::VectorBase<Real> &u_vb) const {
    const Real half(0.5), two(2), three(3), four(4), six(6);
    const Thyra::ConstDetachedVectorView<Real> u(u_vb);
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(nx_,four*dx_/six + dt_*eps1_*two/dx_);
    d[0]     = dx_/three + dt_*eps1_/dx_;
    d[nx_-1] = dx_/three + dt_*eps1_/dx_;
    o.clear();
    o.resize(nx_-1,dx_/six - dt_*eps1_/dx_);
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
          d[i]+= dt_*w*f*phi1*phi1;
          // Off diagonal contribution
          phi2 = (x-static_cast<Real>(i)*dx_)/dx_;
          o[i]+= dt_*w*f*phi1*phi2;
        }
      }
      if (i>0) {
        for (uint j=0; j<4; j++) {
          x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*i-1);
          w = half*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = (x-static_cast<Real>(i-1)*dx_)/dx_;
          d[i]+= dt_*w*f*phi1*phi1;
        }
      }
    }
  }

  void apply_pde_jacobian(ROL::Ptr<Thyra::VectorBase<Real>> &jv_ptr,
                          ROL::Ptr<const Thyra::VectorBase<Real>> &v_ptr,
                          ROL::Ptr<const Thyra::VectorBase<Real>> &u_ptr) const {
    const Real zero(0), half(0.5), two(2), six(6);
    const Thyra::ConstDetachedVectorView<Real> v(v_ptr);
    const Thyra::ConstDetachedVectorView<Real> u(u_ptr);
    Thyra::put_scalar(zero, jv_ptr.ptr());
    Real phi1(0), phi2(0), f(0), x(0), w(0);
    { // scope for jv eval
      Thyra::DetachedVectorView<Real> jv(jv_ptr);
      for (uint n = 0; n < nx_; n++) {
        if ( n < nx_-1 ) {
          jv[n] += dx_/six*(two*v[n]+v[n+1]) + dt_*eps1_/dx_*(v[n]-v[n+1]); // Mass & stiffness
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n+1);
            w = half*dx_*wts_[j];
            f = evaluate_nonlinearity(x,u,1);
            // Diagonal contribution
            phi1 = (static_cast<Real>(n+1)*dx_-x)/dx_;
            jv[n]+= dt_*w*f*phi1*phi1*v[n];
            // Off diagonal contribution
            phi2 = (x-static_cast<Real>(n)*dx_)/dx_;
            jv[n]+= dt_*w*f*phi1*phi2*v[n+1];
          }
        }
        if ( n > 0 ) {
          jv[n] += dx_/six*(v[n-1]+two*v[n]) + dt_*eps1_/dx_*(v[n]-v[n-1]); // Mass & stiffness
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n-1);
            w = half*dx_*wts_[j];
            f = evaluate_nonlinearity(x,u,1);
            // Diagonal contribution
            phi1 = (x-static_cast<Real>(n-1)*dx_)/dx_;
            jv[n]+= dt_*w*f*phi1*phi1*v[n];
            // Off diagonal contribution
            phi2 = (static_cast<Real>(n)*dx_-x)/dx_;
            jv[n]+= dt_*w*f*phi1*phi2*v[n-1];
          }
        }
      }
    } // end scope for jv eval
  }

  void apply_control_jacobian(ROL::Ptr<Thyra::VectorBase<Real>> &jv_ptr,
                              ROL::Ptr<const Thyra::VectorBase<Real>> &v_ptr,
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
            jv[n] = -dt_*dx_/six*v[n];
          }
          else if ( n == 1 ) {
            jv[n] = -dt_*dx_/six*(four*v[n-1]+v[n]);
          }
          else if ( n == nx_ ) {
            jv[n] = -dt_*dx_/six*(four*v[n-1]+v[n-2]);
          }
          else if ( n == nx_+1 ) {
            jv[n] = -dt_*dx_/six*v[n-2];
          }
          else {
            jv[n] = -dt_*dx_/six*(v[n-2]+four*v[n-1]+v[n]);
          }
        }
        else {
          jv[n] -= dt_*dx_/six*(v[n]+four*v[n+1]+v[n+2]);
        }
      }
    } // end scope for jv eval
  }

  void apply_pde_hessian(ROL::Ptr<Thyra::VectorBase<Real>> &r_ptr,
                         ROL::Ptr<const Thyra::VectorBase<Real>> &u_ptr,
                         ROL::Ptr<const Thyra::VectorBase<Real>> &p_ptr,
                         ROL::Ptr<const Thyra::VectorBase<Real>> &s_ptr) const {
    const Real zero(0), half(0.5);
    const Thyra::ConstDetachedVectorView<Real> u(u_ptr);
    const Thyra::ConstDetachedVectorView<Real> p(p_ptr);
    const Thyra::ConstDetachedVectorView<Real> s(s_ptr);
    Thyra::put_scalar(zero, r_ptr.ptr());
    // Contribution from nonlinearity
    Real phi(0), fx(0), px(0), sx(0), x(0), w(0);
    { // scope for r eval
      Thyra::DetachedVectorView<Real> r(r_ptr);
      for (uint n = 0; n < nx_; n++) {
        if (n < nx_-1) {
          for (uint j=0; j<4; j++) {
            x  = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n+1);
            w  = half*dx_*wts_[j];
            fx = evaluate_nonlinearity(x,u,2);
            px = evaluate_solution(x,p);
            sx = evaluate_solution(x,s);
            phi = (static_cast<Real>(n+1)*dx_-x)/dx_;
            r[n]+= dt_*w*fx*px*sx*phi;
          }
        }
        if (n > 0) {
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n-1);
            w = half*dx_*wts_[j];
            fx = evaluate_nonlinearity(x,u,2);
            px = evaluate_solution(x,p);
            sx = evaluate_solution(x,s);
            phi = (x-static_cast<Real>(n-1)*dx_)/dx_;
            r[n]+= dt_*w*fx*px*sx*phi;
          }
        }
      }
    } // end scope for r eval
  }

  Real evaluate_solution(const Real x, const Thyra::ConstDetachedVectorView<Real> &u) const {
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

  Real evaluate_nonlinearity(const Real x, const Thyra::ConstDetachedVectorView<Real> &u, const int deriv = 0) const {
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

  void compute_residual(ROL::Ptr<Thyra::VectorBase<Real>> &r_ptr,
                        const Thyra::VectorBase<Real> &uold_vb,
                        const Thyra::VectorBase<Real> &unew_vb,
                        const Thyra::VectorBase<Real> &z_vb) const {
    const Real zero(0), half(0.5), two(2), four(4), six(6);
    Thyra::put_scalar(zero, r_ptr.ptr());
    Real x(0), w(0);
    const Thyra::ConstDetachedVectorView<Real> up(uold_vb);
    const Thyra::ConstDetachedVectorView<Real> u(unew_vb);
    const Thyra::ConstDetachedVectorView<Real> z(z_vb);
    { // scope for r eval
      Thyra::DetachedVectorView<Real> r(r_ptr);
      for (uint n = 0; n < nx_; n++) {
        if ( n < nx_-1 ) {
          r[n] += dx_/six*(two*u[n]+u[n+1]) + dt_*eps1_/dx_*(u[n]-u[n+1]); // Mass & stiffness
          r[n] -= dx_/six*(two*up[n]+up[n+1]); // Previous time step
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n+1);
            w = half*dx_*wts_[j];
            r[n]+= dt_*w*evaluate_nonlinearity(x,u,0)*(static_cast<Real>(n+1)*dx_-x)/dx_;
          }
        }
        if ( n > 0 ) {
          r[n] += dx_/six*(u[n-1]+two*u[n]) + dt_*eps1_/dx_*(u[n]-u[n-1]); // Mass & stiffness
          r[n] -= dx_/six*(two*up[n]+up[n-1]); // Previous time step
          // Nonlinearity
          for (uint j=0; j<4; j++) {
            x = half*dx_*pts_[j] + half*dx_*static_cast<Real>(2*n-1);
            w = half*dx_*wts_[j];
            r[n]+= dt_*w*evaluate_nonlinearity(x,u,0)*(x-static_cast<Real>(n-1)*dx_)/dx_;
          }
        }
        // Control
        r[n] -= dt_*dx_/six*(z[n]+four*z[n+1]+z[n+2]);
      }
    }
  } // end scope for r eval

  Real compute_norm(const Thyra::VectorBase<Real> &v_vb) const {
    return Thyra::norm_2(v_vb);
  }

  /*void linear_solve(ROL::Ptr<Thyra::VectorBase<Real>> &u,
                    std::vector<Real> &d,
                    std::vector<Real> &o,
                    ROL::Ptr<const Thyra::VectorBase<Real>> &r) const {
    Thyra::copy(*r, u.ptr());
    Teuchos::Ptr<Teuchos::ArrayRCP<Real>> u_data;
    ROL::Ptr<Thyra::SpmdVectorBase<Real>> u_spmd = ROL::dynamicPtrCast(u);
    u_spmd->getNonconstLocalData(u_data);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int nx = static_cast<int>(nx_);
    int info;
    int ldb  = nx;
    int nhrs = 1;
    lp.PTTRF(nx,&d[0],&o[0],&info);
    //lp.PTTRS(nx,nhrs,&d[0],&o[0],&u[0],ldb,&info);
    lp.PTTRS(nx,nhrs,&d[0],&o[0],u_data->getRawPtr(),ldb,&info);
  }*/

  void linear_solve(ROL::Ptr<Thyra::VectorBase<Real>> &u_ptr,
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

  void run_newton(ROL::Ptr<Thyra::VectorBase<Real>> &u,
                  const Thyra::VectorBase<Real>    &up,
                  const Thyra::VectorBase<Real>     &z) const {
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

  Constraint_ParabolicControl(ROL::ParameterList &pl) {
    const Real one(1);
    uint nt   = pl.get("Temporal Discretization", 100);
    Real  T   = pl.get("End Time",                1.0);
    type_     = pl.get("Reaction Type",             1);
    eps1_     = pl.get("Diffusion Scale",         1.0);
    eps2_     = pl.get("Reaction Scale",          1.0);
    nx_       = pl.get("Spatial Discretization",  128);
    dx_       = one/(static_cast<Real>(nx_)-one);
    dt_       = T/(static_cast<Real>(nt)-one);

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
  }

  void value(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &uold,
             const ROL::Vector<Real> &unew,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       c_ptr    = getVector(c);
    ROL::Ptr<const Thyra::VectorBase<Real>> uold_ptr = getVector(uold);
    ROL::Ptr<const Thyra::VectorBase<Real>> unew_ptr = getVector(unew);
    ROL::Ptr<const Thyra::VectorBase<Real>> z_ptr    = getVector(z);
    compute_residual(c_ptr, *uold_ptr, *unew_ptr, *z_ptr);
  }

  void solve(ROL::Vector<Real>       &c,
             const ROL::Vector<Real> &uold,
             ROL::Vector<Real>       &unew,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts) {
    ROL::Ptr<const Thyra::VectorBase<Real>> uold_ptr = getVector(uold);
    ROL::Ptr<Thyra::VectorBase<Real>>       unew_ptr = getVector(unew);
    ROL::Ptr<const Thyra::VectorBase<Real>> z_ptr    = getVector(z);
    run_newton(unew_ptr, *uold_ptr, *z_ptr);
    value(c, uold, unew, z, ts);
  }

  void applyJacobian_uo(ROL::Vector<Real> &jv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &uold,
                        const ROL::Vector<Real> &unew,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       jv_ptr = getVector(jv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    apply_mass(jv_ptr, v_ptr);
    jv.scale(static_cast<Real>(-1));
  }

  void applyJacobian_un(ROL::Vector<Real> &jv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &uold,
                        const ROL::Vector<Real> &unew,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       jv_ptr = getVector(jv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    ROL::Ptr<const Thyra::VectorBase<Real>> un_ptr = getVector(unew);
    apply_pde_jacobian(jv_ptr, v_ptr, un_ptr);
  }

  void applyJacobian_z(ROL::Vector<Real> &jv,
                       const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &uold,
                       const ROL::Vector<Real> &unew,
                       const ROL::Vector<Real> &z,
                       const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       jv_ptr = getVector(jv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    apply_control_jacobian(jv_ptr, v_ptr);
  }

  void applyAdjointJacobian_uo(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_uo(jv, v, uold, unew, z, ts);
  }

  void applyAdjointJacobian_un(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_un(jv, v, uold, unew, z, ts);
  }

  void applyAdjointJacobian_z(ROL::Vector<Real> &jv,
                              const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &uold,
                              const ROL::Vector<Real> &unew,
                              const ROL::Vector<Real> &z,
                              const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       jv_ptr = getVector(jv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    apply_control_jacobian(jv_ptr, v_ptr, true);
  }

  void applyInverseJacobian_un(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       jv_ptr = getVector(jv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    ROL::Ptr<const Thyra::VectorBase<Real>> un_ptr = getVector(unew);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    compute_pde_jacobian(d, o, *un_ptr);
    linear_solve(jv_ptr, d, o, *v_ptr);
  }

  void applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,
                                      const ROL::Vector<Real> &v,
                                      const ROL::Vector<Real> &uold,
                                      const ROL::Vector<Real> &unew,
                                      const ROL::Vector<Real> &z,
                                      const ROL::TimeStamp<Real> &ts) const {
    applyInverseJacobian_un(jv, v, uold, unew, z, ts);
  }

  void applyAdjointHessian_un_un(ROL::Vector<Real> &ahwv,
                                 const ROL::Vector<Real> &w,
                                 const ROL::Vector<Real> &v,
                                 const ROL::Vector<Real> &uold,
                                 const ROL::Vector<Real> &unew,
                                 const ROL::Vector<Real> &z,
                                 const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Thyra::VectorBase<Real>>       ahwv_ptr = getVector(ahwv);
    ROL::Ptr<const Thyra::VectorBase<Real>>  w_ptr = getVector(w);
    ROL::Ptr<const Thyra::VectorBase<Real>>  v_ptr = getVector(v);
    ROL::Ptr<const Thyra::VectorBase<Real>> un_ptr = getVector(unew);
    apply_pde_hessian(ahwv_ptr, un_ptr, w_ptr, v_ptr);
  }

  void applyAdjointHessian_un_uo(ROL::Vector<Real> &ahwv,
                                 const ROL::Vector<Real> &w,
                                 const ROL::Vector<Real> &v,
                                 const ROL::Vector<Real> &uold,
                                 const ROL::Vector<Real> &unew,
                                 const ROL::Vector<Real> &z,
                                 const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_un_z(ROL::Vector<Real> &ahwv,
                                const ROL::Vector<Real> &w,
                                const ROL::Vector<Real> &v,
                                const ROL::Vector<Real> &uold,
                                const ROL::Vector<Real> &unew,
                                const ROL::Vector<Real> &z,
                                const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_uo_un(ROL::Vector<Real> &ahwv,
                                 const ROL::Vector<Real> &w,
                                 const ROL::Vector<Real> &v,
                                 const ROL::Vector<Real> &uold,
                                 const ROL::Vector<Real> &unew,
                                 const ROL::Vector<Real> &z,
                                 const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_uo_uo(ROL::Vector<Real> &ahwv,
                                 const ROL::Vector<Real> &w,
                                 const ROL::Vector<Real> &v,
                                 const ROL::Vector<Real> &uold,
                                 const ROL::Vector<Real> &unew,
                                 const ROL::Vector<Real> &z,
                                 const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_uo_z(ROL::Vector<Real> &ahwv,
                                const ROL::Vector<Real> &w,
                                const ROL::Vector<Real> &v,
                                const ROL::Vector<Real> &uold,
                                const ROL::Vector<Real> &unew,
                                const ROL::Vector<Real> &z,
                                const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_z_un(ROL::Vector<Real> &ahwv,
                                const ROL::Vector<Real> &w,
                                const ROL::Vector<Real> &v,
                                const ROL::Vector<Real> &uold,
                                const ROL::Vector<Real> &unew,
                                const ROL::Vector<Real> &z,
                                const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_z_uo(ROL::Vector<Real> &ahwv,
                                const ROL::Vector<Real> &w,
                                const ROL::Vector<Real> &v,
                                const ROL::Vector<Real> &uold,
                                const ROL::Vector<Real> &unew,
                                const ROL::Vector<Real> &z,
                                const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_z_z(ROL::Vector<Real> &ahwv,
                               const ROL::Vector<Real> &w,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }
};



/***********************************************************************
****************************  OBJECTIVE  *******************************
************************************************************************/

template<class Real>
class Objective_ParabolicControl : public ROL::DynamicObjective<Real> {

  typedef typename std::vector<Real>::size_type uint;

private:
  Real alpha_;
  Real theta_;
  Teuchos_Ordinal nx_;
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
    Real c = ((x.subDim()==nx_) ? four : two);
    for (Teuchos_Ordinal i=0; i<x.subDim(); i++) {
      if ( i == 0 ) {
        ip += dx_/six*(c*x[i] + x[i+1])*y[i];
      }
      else if ( i == x.subDim()-1 ) {
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
    Real c = ((u.subDim()==nx_) ? four : two);
    { // scope for Mu eval
      Thyra::DetachedVectorView<Real> Mu(Mu_ptr);
      for (Teuchos_Ordinal i=0; i<u.subDim(); i++) {
        if ( i == 0 ) {
          Mu[i] = dx_/six*(c*u[i] + u[i+1]);
        }
        else if ( i == u.subDim()-1 ) {
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
      for (Teuchos_Ordinal n = 0; n < nx_; n++) {
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
      for (Teuchos_Ordinal n = 0; n < nx_; n++) {
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
      for (Teuchos_Ordinal n = 0; n < nx_; n++) {
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
