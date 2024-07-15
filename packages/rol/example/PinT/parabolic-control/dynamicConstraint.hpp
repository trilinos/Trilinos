// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  dynamicConstraint.cpp
    \brief Shows how to solve a 1D semilinear parabolic equation,
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
*/

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_DynamicConstraint.hpp"
#include "Teuchos_LAPACK.hpp"

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
  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) const {
    const Real zero(0), two(2), six(6);
    Mu.clear();
    Mu.resize(nx_,zero);
    for (uint n = 0; n < nx_; n++) {
      if ( n < nx_-1 ) {
        Mu[n] += dx_/six*(two*u[n]+u[n+1]);
      }
      if ( n > 0 ) {
        Mu[n] += dx_/six*(u[n-1]+two*u[n]);
      }
    }
  }

  void compute_pde_jacobian(std::vector<Real> &d, std::vector<Real> &o, const std::vector<Real> &u) const {
    const Real half(0.5), two(2), three(3), four(4), six(6);
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

  void apply_pde_jacobian(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &u) const {
    const Real zero(0), half(0.5), two(2), six(6);
    jv.clear();
    jv.resize(nx_,zero);
    Real phi1(0), phi2(0), f(0), x(0), w(0);
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
  }

  void apply_control_jacobian(std::vector<Real> &jv, const std::vector<Real> &v, bool adjoint = false) const {
    const Real zero(0), four(4), six(6);
    jv.clear();
    uint dim = ((adjoint == true) ? nx_+2 : nx_);
    jv.resize(dim,zero);
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
  }

  void apply_pde_hessian(std::vector<Real> &r, const std::vector<Real> &u, const std::vector<Real> &p, 
                         const std::vector<Real> &s) const {
    const Real zero(0), half(0.5);
    r.clear();
    r.resize(nx_,zero);
    // Contribution from nonlinearity
    Real phi(0), fx(0), px(0), sx(0), x(0), w(0);
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
  }

  Real evaluate_solution(const Real x, const std::vector<Real> &u) const {
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

  Real evaluate_nonlinearity(const Real x, const std::vector<Real> &u, const int deriv = 0) const {
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


  void compute_residual(std::vector<Real> &r, const std::vector<Real> &up, 
                        const std::vector<Real> &u, const std::vector<Real> &z) const {
    const Real zero(0), half(0.5), two(2), four(4), six(6);
    r.clear();
    r.resize(nx_,zero);
    Real x(0), w(0);
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

  Real compute_norm(const std::vector<Real> &r) const {
    Real norm(0);
    for (uint i = 0; i < r.size(); i++) {
      norm += r[i]*r[i];
    }
    return std::sqrt(norm);
  }

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) const {
    for (uint i = 0; i < u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &d, std::vector<Real> &o, 
              const std::vector<Real> &r) const {
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

  void run_newton(std::vector<Real> &u, const std::vector<Real> &up, const std::vector<Real> &z) const {
    const Real half(0.5), one(1), lstol(1e-4);
    // Set initial guess
    u.assign(up.begin(),up.end());
    // Compute residual and residual norm
    std::vector<Real> r(u.size(),0.0);
    compute_residual(r,up,u,z);
    Real rnorm = compute_norm(r);
    // Define tolerances
    Real tol   = static_cast<Real>(1e2)*ROL::ROL_EPSILON<Real>();
    Real maxit = 100;
    // Initialize Jacobian storage
    std::vector<Real> d(nx_);
    std::vector<Real> o(nx_-1);
    // Iterate Newton's method
    Real alpha(1), tmp(0);
    std::vector<Real> s(nx_);
    std::vector<Real> utmp(nx_);
    for (uint i = 0; i < maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      compute_pde_jacobian(d,o,u);
      // Solve Newton system
      linear_solve(s,d,o,r);
      // Perform line search
      tmp = rnorm;
      alpha = one;
      utmp.assign(u.begin(),u.end());
      update(utmp,s,-alpha);
      compute_residual(r,up,utmp,z);
      rnorm = compute_norm(r); 
      while ( rnorm > (one-lstol*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON<Real>()) ) {
        alpha *= half;
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

  ROL::Ptr<const std::vector<Real>> getVector( const ROL::Vector<Real>& x ) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real>> getVector( ROL::Vector<Real>& x ) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
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

  void value(ROL::Vector<Real>    &c,    const ROL::Vector<Real> &uold,
       const ROL::Vector<Real>    &unew, const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>        cp = getVector(c);
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    compute_residual(*cp, *uop, *unp, *zp);
  }

  void solve(ROL::Vector<Real>    &c,    const ROL::Vector<Real> &uold,
             ROL::Vector<Real>    &unew, const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) {
    ROL::Ptr<const std::vector<Real>> uop = getVector(uold);
    ROL::Ptr<std::vector<Real>>       unp = getVector(unew);
    ROL::Ptr<const std::vector<Real>>  zp = getVector(z);
    run_newton(*unp, *uop, *zp);
    value(c, uold, unew, z, ts);
  }

  void applyJacobian_uo(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                  const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    apply_mass(*jvp, *vp);
    jv.scale(static_cast<Real>(-1));
  }

  void applyJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                  const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    apply_pde_jacobian(*jvp, *vp, *unp);
  }

  void applyJacobian_z(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                 const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                 const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    apply_control_jacobian(*jvp, *vp);
  }

  void applyAdjointJacobian_uo(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                         const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_uo(jv, v, uold, unew, z, ts);
  }

  void applyAdjointJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                         const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyJacobian_un(jv, v, uold, unew, z, ts);
  }

  void applyAdjointJacobian_z(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                        const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                        const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    apply_control_jacobian(*jvp, *vp, true);
  }

  void applyInverseJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                         const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>       jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    compute_pde_jacobian(d, o, *unp);
    linear_solve(*jvp, d, o, *vp);
  }

  void applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,   const ROL::Vector<Real>    &v,
                                const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                                const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    applyInverseJacobian_un(jv, v, uold, unew, z, ts);
  }

  void applyAdjointHessian_un_un(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>     ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real>>  wp = getVector(w);
    ROL::Ptr<const std::vector<Real>>  vp = getVector(v);
    ROL::Ptr<const std::vector<Real>> unp = getVector(unew);
    apply_pde_hessian(*ahwvp, *unp, *wp, *vp);
  }

  void applyAdjointHessian_un_uo(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_un_z(ROL::Vector<Real> &ahwv,
                          const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                          const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_uo_un(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_uo_uo(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_uo_z(ROL::Vector<Real> &ahwv,
                          const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                          const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_z_un(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_z_uo(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                           const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }

  void applyAdjointHessian_z_z(ROL::Vector<Real> &ahwv,
                          const ROL::Vector<Real> &w,    const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real> &uold, const ROL::Vector<Real>    &unew,
                          const ROL::Vector<Real> &z,    const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }
};
