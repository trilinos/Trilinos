// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_05.cpp
    \brief Shows how to solve a control problem governed by a semilinear 
           parabolic equation.  The nonlinearity is of Cahn-Hilliard-type.  
           The control problem is 
           \f[
              \min_{u,z} \;\frac{1}{2} \int_0^1 (u(T,x)-\bar{u}(x))^2\,\mathrm{d}x
                         \frac{\alpha}{2} \int_0^T z(t)^2\,\mathrm{d}t
           \f]
           subject to 
           \f[
               u_t(t,x) - u_{xx}(t,x) + u(t,x)^3 - u(t,x) = 0 \quad t\in (0,T], \; x\in (0,1)
           \f]
           with boundary conditions
           \f[
               u_x(t,0) = 0, \; u_x(t,1) = z(t)
           \f]
           and initial condition
           \f[
               u(0,x) = 0.
           \f]
*/

#include "ROL_Algorithm.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_Types.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>
#include <ctime>

template<class Real>
class Constraint_ParabolicControl : public ROL::Constraint_SimOpt<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;

  typedef typename vector::size_type uint;

private:
  std::vector<Real> u0_;
  Real eps1_;
  Real eps2_;
  uint nx_;
  uint nt_;
  Real T_;
  Real dx_;
  Real dt_;

  std::vector<Real> pts_;
  std::vector<Real> wts_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
    Mu.clear();
    Mu.resize(nx_,0.0);
    for (uint n = 0; n < nx_; n++) {
      if ( n < nx_-1 ) {
        Mu[n] += dx_/6.0*(2.0*u[n]+u[n+1]);
      }
      if ( n > 0 ) {
        Mu[n] += dx_/6.0*(u[n-1]+2.0*u[n]);
      }
    }
  }

  void compute_pde_jacobian(std::vector<Real> &d, std::vector<Real> &o, const std::vector<Real> &u) {
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(nx_,4.0*dx_/6.0 + dt_*eps1_*2.0/dx_);
    d[0]     = dx_/3.0 + dt_*eps1_/dx_;
    d[nx_-1] = dx_/3.0 + dt_*eps1_/dx_;
    o.clear();
    o.resize(nx_-1,dx_/6.0 - dt_*eps1_/dx_);
    // Contribution from nonlinearity
    Real phi1 = 0.0, phi2 = 0.0, f = 0.0, x = 0.0, w = 0.0;
    for (uint i=0; i<nx_; i++) {
      if (i<nx_-1) {
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*i+1);
          w = 0.5*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = ((Real)(i+1)*dx_-x)/dx_;
          d[i]+= dt_*w*f*phi1*phi1;
          // Off diagonal contribution
          phi2 = (x-(Real)(i)*dx_)/dx_;
          o[i]+= dt_*w*f*phi1*phi2;
        }
      }
      if (i>0) {
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*i-1);
          w = 0.5*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = (x-(Real)(i-1)*dx_)/dx_;
          d[i]+= dt_*w*f*phi1*phi1;
        }
      }
    }
  }

  void apply_pde_jacobian(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &u) {
    jv.clear();
    jv.resize(nx_,0.0);
    Real phi1 = 0.0, phi2 = 0.0, f = 0.0, x = 0.0, w = 0.0;
    for (uint n = 0; n < nx_; n++) {
      if ( n < nx_-1 ) {
        jv[n] += dx_/6.0*(2.0*v[n]+v[n+1]) + dt_*eps1_/dx_*(v[n]-v[n+1]); // Mass & stiffness
        // Nonlinearity
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*n+1);
          w = 0.5*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = ((Real)(n+1)*dx_-x)/dx_;
          jv[n]+= dt_*w*f*phi1*phi1*v[n];
          // Off diagonal contribution
          phi2 = (x-(Real)(n)*dx_)/dx_;
          jv[n]+= dt_*w*f*phi1*phi2*v[n+1];
        }
      }
      if ( n > 0 ) {
        jv[n] += dx_/6.0*(v[n-1]+2.0*v[n]) + dt_*eps1_/dx_*(v[n]-v[n-1]); // Mass & stiffness
        // Nonlinearity
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*n-1);
          w = 0.5*dx_*wts_[j];
          f = evaluate_nonlinearity(x,u,1);
          // Diagonal contribution
          phi1 = (x-(Real)(n-1)*dx_)/dx_;
          jv[n]+= dt_*w*f*phi1*phi1*v[n];
          // Off diagonal contribution
          phi2 = ((Real)(n)*dx_-x)/dx_;
          jv[n]+= dt_*w*f*phi1*phi2*v[n-1];
        }
      }
    }
  }

  void apply_pde_hessian(std::vector<Real> &r, const std::vector<Real> &u, const std::vector<Real> &p, 
                         const std::vector<Real> &s) {
    r.clear();
    r.resize(nx_,0.0);
    // Contribution from nonlinearity
    Real phi = 0.0, fx = 0.0, px = 0.0, sx = 0.0, x = 0.0, w = 0.0;
    for (uint n = 0; n < nx_; n++) {
      if (n < nx_-1) {
        for (uint j=0; j<4; j++) {
          x  = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*n+1);
          w  = 0.5*dx_*wts_[j];
          fx = evaluate_nonlinearity(x,u,2);
          px = evaluate_solution(x,p);
          sx = evaluate_solution(x,s);
          phi = ((Real)(n+1)*dx_-x)/dx_;
          r[n]+= dt_*w*fx*px*sx*phi;
        }
      }
      if (n > 0) {
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*n-1);
          w = 0.5*dx_*wts_[j];
          fx = evaluate_nonlinearity(x,u,2);
          px = evaluate_solution(x,p);
          sx = evaluate_solution(x,s);
          phi = (x-(Real)(n-1)*dx_)/dx_;
          r[n]+= dt_*w*fx*px*sx*phi;
        }
      }
    }
  }

  Real evaluate_solution(const Real x, const std::vector<Real> &u) {
    // Determine u(x)
    Real pt  = 0.0;
    Real val = 0.0;
    for (uint n = 0; n < nx_; n++) {
      if (x <= (Real)(n+1)*dx_ && x >= (Real)(n)*dx_) {
        pt  = (Real)(n+1)*dx_;
        val = u[n]*(pt-x)/dx_;
        pt  = (Real)(n)*dx_;
        val+= u[n+1]*(x-pt)/dx_;
        break;
      }
      else if (x <= (Real)(n)*dx_ && x >= (Real)(n-1)*dx_) {
        pt  = (Real)(n)*dx_;
        val = u[n-1]*(pt-x)/dx_;
        pt  = (Real)(n-1)*dx_;
        val+= u[n]*(x-pt)/dx_;
        break;
      }
    }
    return val;
  }

  Real evaluate_nonlinearity(const Real x, const std::vector<Real> &u, const int deriv = 0) {
    // Compute u(x)^3 - u(x) or its derivatives 3 u(x)^2 - 1 and 6 u(x)
    Real val = evaluate_solution(x,u);
    if (deriv == 0) {
      return eps2_*(std::pow(val,3.0)-val);
    }
    else if (deriv == 1) {
      return eps2_*(3.0*std::pow(val,2.0)-1.0);
    }
    else {
      return eps2_*(6.0*val);
    }
  }


  void compute_residual(std::vector<Real> &r, const std::vector<Real> &up, 
                        const std::vector<Real> &u, const Real z) {
    r.clear();
    r.resize(nx_,0.0);
    Real x = 0.0, w = 0.0;
    for (uint n = 0; n < nx_; n++) {
      if ( n < nx_-1 ) {
        r[n] += dx_/6.0*(2.0*u[n]+u[n+1]) + dt_*eps1_/dx_*(u[n]-u[n+1]); // Mass & stiffness
        r[n] -= dx_/6.0*(2.0*up[n]+up[n+1]); // Previous time step
        // Nonlinearity
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*n+1);
          w = 0.5*dx_*wts_[j];
          r[n]+= dt_*w*evaluate_nonlinearity(x,u,0)*((Real)(n+1)*dx_-x)/dx_;
        }
      }
      if ( n > 0 ) {
        r[n] += dx_/6.0*(u[n-1]+2.0*u[n]) + dt_*eps1_/dx_*(u[n]-u[n-1]); // Mass & stiffness
        r[n] -= dx_/6.0*(2.0*up[n]+up[n-1]); // Previous time step
        // Nonlinearity
        for (uint j=0; j<4; j++) {
          x = 0.5*dx_*pts_[j] + 0.5*dx_*(Real)(2*n-1);
          w = 0.5*dx_*wts_[j];
          r[n]+= dt_*w*evaluate_nonlinearity(x,u,0)*(x-(Real)(n-1)*dx_)/dx_;
        }
      }
    }
    r[nx_-1] -= dt_*z; // Control
  }

  Real compute_norm(const std::vector<Real> &r) {
    Real norm = 0.0;
    for (uint i = 0; i < r.size(); i++) {
      norm += r[i]*r[i];
    }
    return std::sqrt(norm);
  }

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    for (uint i = 0; i < u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &d, std::vector<Real> &o, 
              const std::vector<Real> &r) {
    u.assign(r.begin(),r.end());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int nx = static_cast<int>(nx_);
    int info;
    int ldb  = nx;
    int nhrs = 1;
    lp.PTTRF(nx,&d[0],&o[0],&info);
    lp.PTTRS(nx,nhrs,&d[0],&o[0],&u[0],ldb,&info);
  }

  void run_newton(std::vector<Real> &u, const std::vector<Real> &up, const Real z) {
    // Set initial guess
    u.assign(up.begin(),up.end());
    // Compute residual and residual norm
    std::vector<Real> r(u.size(),0.0);
    compute_residual(r,up,u,z);
    Real rnorm = compute_norm(r);
    // Define tolerances
    Real tol   = 1.e2*ROL::ROL_EPSILON<Real>();
    Real maxit = 100;
    // Initialize Jacobian storage
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    // Iterate Newton's method
    Real alpha = 1.0, tmp = 0.0;
    std::vector<Real> s(nx_,0.0);
    std::vector<Real> utmp(nx_,0.0);
    for (uint i = 0; i < maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      compute_pde_jacobian(d,o,u);
      // Solve Newton system
      linear_solve(s,d,o,r);
      // Perform line search
      tmp = rnorm;
      alpha = 1.0;
      utmp.assign(u.begin(),u.end());
      update(utmp,s,-alpha);
      compute_residual(r,up,utmp,z);
      rnorm = compute_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON<Real>()) ) {
        alpha /= 2.0;
        utmp.assign(u.begin(),u.end());
        update(utmp,s,-alpha);
        compute_residual(r,up,utmp,z);
        rnorm = compute_norm(r); 
      }
      // Update iterate
      u.assign(utmp.begin(),utmp.end());
      if ( rnorm < tol ) {
        break;
      }
    }
  }

  ROL::Ptr<const std::vector<Real> > getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real> > getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Constraint_ParabolicControl(Real eps = 1.0, uint nx = 128, uint nt = 100, Real T = 1) 
    : eps1_(eps*eps), eps2_(1.0), nx_(nx), nt_(nt), T_(T) {
    u0_.resize(nx_,0.0);
    dx_ = 1.0/((Real)nx-1.0);
    dt_ = T/((Real)nt-1.0);

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

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
      

    ROL::Ptr<vector> cp = getVector(c);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
    std::vector<Real> C(nx_,0.0);
    std::vector<Real> uold(u0_);
    std::vector<Real> unew(u0_);
    for (uint t = 0; t < nt_; t++) {
      // Copy state and control at t time step
      for (uint n = 0; n < nx_; n++) {
        unew[n] = (*up)[t*nx_+n];
      }
      // Evaluate residual at t time step
      compute_residual(C,uold,unew,(*zp)[t]);
      // Copy residual at t time step
      for (uint n = 0; n < nx_; n++) {
        (*cp)[t*nx_+n] = C[n];
      }
      uold.assign(unew.begin(),unew.end());
    }
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
      
    ROL::Ptr<vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Initialize State Storage
    std::vector<Real> uold(u0_);
    std::vector<Real> unew(u0_);
    // Time Step Using Implicit Euler
    for ( uint t = 0; t < nt_; t++ ) {
      run_newton(unew,uold,(*zp)[t]);
      for( uint n = 0; n < nx_; n++) {
        (*up)[t*nx_+n] = unew[n];
      }
      uold.assign(unew.begin(),unew.end());
    }
    /* TEST SOLVE -- ||c(u,z)|| SHOULD BE ZERO */
    value(c,u,z,tol);
    //std::cout << c->norm() << "\n";
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {

       

    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    std::vector<Real> J(u0_.size(),0.0);
    std::vector<Real> M(u0_.size(),0.0);
    std::vector<Real> vold(u0_);
    std::vector<Real> unew(u0_);
    std::vector<Real> vnew(u0_);
    for (uint t = 0; t < nt_; t++) {
      for (uint n = 0; n < nx_; n++) {
        unew[n] = (*up)[t*nx_+n];
        vnew[n] = (*vp)[t*nx_+n];
      }
      apply_pde_jacobian(J,vnew,unew);
      if ( t > 0 ) {
        apply_mass(M,vold);
      }
      for (uint n = 0; n < nx_; n++) {
        (*jvp)[t*nx_+n] = J[n] - M[n];
      }
      vold.assign(vnew.begin(),vnew.end());
    }
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
      
    
    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Initialize State Storage
    std::vector<Real> M(u0_);
    std::vector<Real> sold(u0_);
    std::vector<Real> unew(u0_);
    std::vector<Real> vnew(u0_);
    std::vector<Real> snew(u0_);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    // Time Step Using Implicit Euler
    for (uint t = 0; t < nt_; t++) {
      for (uint n = 0; n < nx_; n++) {
        unew[n] = (*up)[t*nx_+n];
        vnew[n] = (*vp)[t*nx_+n];
      }
      // Get PDE Jacobian
      compute_pde_jacobian(d,o,unew);
      // Get Right Hand Side
      if ( t > 0 ) {
        apply_mass(M,sold); 
        update(vnew,M);
      }
      // Solve solve adjoint system at current time step
      linear_solve(snew,d,o,vnew);
      for(uint n = 0; n < nx_; n++) {
        (*jvp)[t*nx_+n] = snew[n];
      }
      sold.assign(snew.begin(),snew.end());
    }
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {

      
    
    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    std::vector<Real> J(u0_.size(),0.0);
    std::vector<Real> M(u0_.size(),0.0);
    std::vector<Real> vold(u0_);
    std::vector<Real> unew(u0_);
    std::vector<Real> vnew(u0_);
    for (uint t = nt_; t > 0; t--) {
      for (uint n = 0; n < nx_; n++) {
        unew[n] = (*up)[(t-1)*nx_+n];
        vnew[n] = (*vp)[(t-1)*nx_+n];
      }
      apply_pde_jacobian(J,vnew,unew);
      if ( t < nt_ ) {
        apply_mass(M,vold);
      }
      for (uint n = 0; n < nx_; n++) {
        (*jvp)[(t-1)*nx_+n] = J[n] - M[n];
      }
      vold.assign(vnew.begin(),vnew.end());
    }
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, 
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
      

    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Initialize State Storage
    std::vector<Real> M(u0_);
    std::vector<Real> sold(u0_);
    std::vector<Real> unew(u0_);
    std::vector<Real> vnew(u0_);
    std::vector<Real> snew(u0_);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    // Time Step Using Implicit Euler
    for (uint t = nt_; t > 0; t--) {
      for( uint n = 0; n < nx_; n++) {
        unew[n] = (*up)[(t-1)*nx_+n];
        vnew[n] = (*vp)[(t-1)*nx_+n];
      }
      // Get PDE Jacobian
      compute_pde_jacobian(d,o,unew);
      // Get Right Hand Side
      if ( t < nt_ ) {
        apply_mass(M,sold); 
        update(vnew,M);
      }
      // Solve solve adjoint system at current time step
      linear_solve(snew,d,o,vnew);
      for ( uint n = 0; n < nx_; n++) {
        (*jvp)[(t-1)*nx_+n] = snew[n];
      }
      sold.assign(snew.begin(),snew.end());
    }
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {

      

    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    for (uint t = 0; t < nt_; t++) {
      (*jvp)[t*nx_+(nx_-1)] = -dt_*(*vp)[t];
    }
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {

        

    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    for (uint t = 0; t < nt_; t++) {
      (*jvp)[t] = -dt_*(*vp)[t*nx_+(nx_-1)];
    }
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &hwv, const ROL::Vector<Real> &w, 
                              const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {

        
    ROL::Ptr<vector> hwvp = getVector(hwv);
    ROL::Ptr<const vector> wp = getVector(w);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Initialize State Storage
    std::vector<Real> unew(u0_);
    std::vector<Real> wnew(u0_);
    std::vector<Real> vnew(u0_);
    std::vector<Real> snew(u0_);

    // Time Step Using Implicit Euler
    for (uint t = nt_; t > 0; t--) {
      for (uint n = 0; n < nx_; n++) {
        unew[n] = (*up)[(t-1)*nx_+n];
        vnew[n] = (*vp)[(t-1)*nx_+n];
        wnew[n] = (*wp)[(t-1)*nx_+n];
      }
      // Get PDE Hessian
      apply_pde_hessian(snew,unew,wnew,vnew);
      for(uint n = 0; n < nx_; n++) {
        (*hwvp)[(t-1)*nx_+n] = snew[n];
      }
    }
  }

  void applyAdjointHessian_12(ROL::Vector<Real> &hwv, const ROL::Vector<Real> &w, 
                              const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    hwv.zero();
  }

  void applyAdjointHessian_21(ROL::Vector<Real> &hwv, const ROL::Vector<Real> &w, 
                              const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    hwv.zero();
  }

  void applyAdjointHessian_22(ROL::Vector<Real> &hwv, const ROL::Vector<Real> &w, 
                              const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    hwv.zero();
  }
};

template<class Real>
class Objective_ParabolicControl : public ROL::Objective_SimOpt<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
 
  typedef typename vector::size_type uint;

private:
  Real alpha_;
  uint nx_;
  uint nt_;
  Real dx_;
  Real dt_;
  Real T_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
    Mu.resize(nx_,0.0);
    for (uint n = 0; n < nx_; n++) {
      if ( n < nx_-1 ) {
        Mu[n] += dx_/6.0*(2.0*u[n] + u[n+1]);
      }
      if ( n > 0 ) {
        Mu[n] += dx_/6.0*(u[n-1] + 2.0*u[n]);
      }
    }
  }

  Real evaluate_target(Real x) {
    Real val = 0.0;
    int example = 2;
    switch (example) {
      case 1:  val = ((x<0.5) ? 0.5 : 0.0); break;
      case 2:  val = 0.5; break;
      case 3:  val = 0.5*std::abs(std::sin(8.0*M_PI*x)); break;
      case 4:  val = 0.5*std::exp(-0.5*(x-0.5)*(x-0.5)); break;
    }
    return val;
  }

  Real compute_dot(const std::vector<Real> &r, const std::vector<Real> &s) {
    Real ip = 0.0; 
    for (uint i = 0; i < r.size(); i++) {
      ip += r[i]*s[i];
    }
    return ip;
  }

  Real compute_norm(const std::vector<Real> &r) {
    return std::sqrt(compute_dot(r,r));
  }

  Real compute_weighted_dot(const std::vector<Real> &r, const std::vector<Real> &s) {
    std::vector<Real> Mr(nx_,0.0);
    apply_mass(Mr,r);
    return compute_dot(Mr,s);
  }

  Real compute_weighted_norm(const std::vector<Real> &r) {
    return std::sqrt(compute_weighted_dot(r,r));
  }

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    for (uint i = 0; i < u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  ROL::Ptr<const std::vector<Real> > getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real> > getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_ParabolicControl(Real alpha = 1.e-4, uint nx = 128, uint nt = 100, Real T = 1) 
    : alpha_(alpha), nx_(nx), nt_(nt), T_(T) {
    dx_ = 1.0/((Real)nx-1.0);
    dt_ = T/((Real)nt-1.0);
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

      
 
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Compute Norm of State
    std::vector<Real> uT(nx_,0.0);
    for (uint n = 0; n < nx_; n++) {
      uT[n] = (*up)[(nt_-1)*nx_ + n] - evaluate_target((Real)n*dx_);
    } 
    Real val = 0.5*compute_weighted_dot(uT,uT); 
    // Add Norm of Control
    for (uint t = 0; t < nt_; t++) {
      val += 0.5*alpha_*dt_*(*zp)[t]*(*zp)[t];
    }
    return val;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
    
      

    ROL::Ptr<vector> gp = getVector(g);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
 
    std::vector<Real> uT(nx_,0.0);
    for (unsigned n = 0; n < nx_; n++) {
      uT[n] = (*up)[(nt_-1)*nx_ + n] - evaluate_target((Real)n*dx_);
    } 
    std::vector<Real> M(nx_,0.0);
    apply_mass(M,uT);
    for (unsigned n = 0; n < nx_; n++) {
      (*gp)[(nt_-1)*nx_ + n] = M[n];
    }
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();

      

    ROL::Ptr<vector> gp = getVector(g);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Compute gradient
    for (uint n = 0; n < nt_; n++) {
      (*gp)[n] = dt_*alpha_*(*zp)[n];
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                   const ROL::Vector<Real> &z, Real &tol ) {

      

    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<vector> hvp = getVector(hv);

    // Compute HessVec
    std::vector<Real> vT(nx_,0.0);
    for (uint n = 0; n < nx_; n++) {
      vT[n] = (*vp)[(nt_-1)*nx_ + n];
    } 
    std::vector<Real> M(nx_,0.0);
    apply_mass(M,vT);
    for (uint n = 0; n < nx_; n++) {
      (*hvp)[(nt_-1)*nx_ + n] = M[n];
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                   const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                   const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                   const ROL::Vector<Real> &z, Real &tol ) {

      

    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<vector> hvp = getVector(hv);

    // Compute HessVec
    for (uint n = 0; n < nt_; n++) {
      (*hvp)[n] = dt_*alpha_*(*vp)[n];
    }
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>     vector;
  typedef ROL::Vector<RealT>     V;
  typedef ROL::StdVector<RealT>  SV;
  
  typedef typename vector::size_type uint;
  
      

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    // Initialize objective function.
    uint nx     = 64;   // Set spatial discretization.
    uint nt     = 10;   // Set temporal discretization.
    RealT T     = 1.0;   // Set end time.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT eps   = 5.e-1; // Set conductivity 
    Objective_ParabolicControl<RealT> obj(alpha,nx,nt,T);
    Constraint_ParabolicControl<RealT> con(eps,nx,nt,T);
    // Initialize iteration vectors.
    ROL::Ptr<vector> xz_ptr = ROL::makePtr<vector>(nt, 1.0);
    ROL::Ptr<vector> xu_ptr = ROL::makePtr<vector>(nx*nt, 1.0);
    ROL::Ptr<vector> gz_ptr = ROL::makePtr<vector>(nt, 1.0);
    ROL::Ptr<vector> gu_ptr = ROL::makePtr<vector>(nx*nt, 1.0);
    ROL::Ptr<vector> yz_ptr = ROL::makePtr<vector>(nt, 1.0);
    ROL::Ptr<vector> yu_ptr = ROL::makePtr<vector>(nx*nt, 1.0);

    for (uint i=0; i<nt; i++) {
      (*xz_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      for (uint n=0; n<nx; n++) {
        (*xu_ptr)[i*nx + n] = (RealT)rand()/(RealT)RAND_MAX;
        (*yu_ptr)[i*nx + n] = (RealT)rand()/(RealT)RAND_MAX;
      }
    }

    SV xz(xz_ptr);
    SV xu(xu_ptr);
    SV gz(gz_ptr);
    SV gu(gu_ptr);
    SV yz(yz_ptr);
    SV yu(yu_ptr);

    ROL::Ptr<V> xzp = ROL::makePtrFromRef(xz);
    ROL::Ptr<V> xup = ROL::makePtrFromRef(xu);
    ROL::Ptr<V> gzp = ROL::makePtrFromRef(gz);
    ROL::Ptr<V> gup = ROL::makePtrFromRef(gu);

    ROL::Ptr<V> yzp = ROL::makePtrFromRef(yz);
    ROL::Ptr<V> yup = ROL::makePtrFromRef(yu);

    ROL::Vector_SimOpt<RealT> x(xup,xzp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);

    ROL::Ptr<vector> c_ptr  = ROL::makePtr<vector>(nt*nx, 0.0);
    ROL::Ptr<vector> l_ptr  = ROL::makePtr<vector>(nt*nx, 0.0);

    SV c(c_ptr);
    SV l(l_ptr);

    ROL::Ptr<V> cp = ROL::makePtrFromRef(c);

    // Initialize reduced objective function
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj  = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = ROL::makePtrFromRef(con);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,xup,xzp,cp);

    // Check deriatives.
    obj.checkGradient(x,y,true);
    obj.checkHessVec(x,y,true);
    con.checkApplyJacobian(x,y,c,true);

    //con.checkApplyAdjointJacobian(x,yu,c,x,true);
    con.checkApplyAdjointHessian(x,yu,y,x,true);
    robj.checkGradient(xz,yz,true);
    robj.checkHessVec(xz,yz,true);

    // Initialize constraints -- these are set to -infinity and infinity.
    ROL::Ptr<vector> lo_ptr = ROL::makePtr<vector>(nt,ROL::ROL_NINF<RealT>());
    ROL::Ptr<vector> up_ptr = ROL::makePtr<vector>(nt,ROL::ROL_INF<RealT>());
 
    ROL::Ptr<V> lo = ROL::makePtr<SV>(lo_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(up_ptr);

    ROL::Bounds<RealT> icon(lo,up);

    // Primal dual active set.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    // Krylov parameters.
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-8);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);
    // PDAS parameters.
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-6);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 10);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);
    parlist->sublist("General").sublist("Secant").set("Use as Hessian",true);
    // Status test parameters.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    // Define algorithm.
    ROL::Ptr<ROL::Step<RealT>>       step   = ROL::makePtr<ROL::PrimalDualActiveSetStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Ptr<ROL::Algorithm<RealT>>  algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
    // Run algorithm.
    xz.zero();
    std::clock_t timer_pdas = std::clock();
    algo->run(xz, robj, icon, true, *outStream);
    *outStream << "Primal Dual Active Set required " << (std::clock()-timer_pdas)/(RealT)CLOCKS_PER_SEC 
               << " seconds.\n";

    // Projected Newton.
    // re-load parameters
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // Set algorithm.
    step   = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
    // Run Algorithm
    xz.zero();
    std::clock_t timer_tr = std::clock();
    algo->run(xz, robj, icon, true, *outStream);
    *outStream << "Projected Newton required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC 
               << " seconds.\n";

    // Composite step.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-10);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    // Set algorithm.
    step   = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
    x.zero();
    std::clock_t timer_cs = std::clock();
    algo->run(x, g, l, c, obj, con, true, *outStream);
    *outStream << "Composite Step required " << (std::clock()-timer_cs)/(RealT)CLOCKS_PER_SEC 
               << " seconds.\n";
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

