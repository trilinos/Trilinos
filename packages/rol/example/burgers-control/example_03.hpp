// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_04.hpp
    \brief Provides definitions of equality constraint and objective for
           example_04.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_Types.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>
#include <ctime>

template<class Real>
class Constraint_BurgersControl : public ROL::Constraint_SimOpt<Real> {
private:
  unsigned nx_;
  unsigned nt_;

  Real dx_;
  Real T_;
  Real dt_;

  Real nu_;
  Real u0_;
  Real u1_;
  Real f_;
  std::vector<Real> u_init_;

private:
  Real compute_norm(const std::vector<Real> &r) {
    return std::sqrt(dot(r,r));
  }

  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) {
    Real ip = 0.0;
    Real c = ((x.size()==nx_) ? 4.0 : 2.0);
    for (unsigned i = 0; i < x.size(); i++) {
      if ( i == 0 ) {
        ip += dx_/6.0*(c*x[i] + x[i+1])*y[i];
      }
      else if ( i == x.size()-1 ) {
        ip += dx_/6.0*(x[i-1] + c*x[i])*y[i];
      }
      else {
        ip += dx_/6.0*(x[i-1] + 4.0*x[i] + x[i+1])*y[i];
      }
    }
    return ip;
  }

  using ROL::Constraint_SimOpt<Real>::update;

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    for (unsigned i = 0; i < u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void scale(std::vector<Real> &u, const Real alpha=0.0) {
    for (unsigned i = 0; i < u.size(); i++) {
      u[i] *= alpha;
    }
  }

  void compute_residual(std::vector<Real> &r, const std::vector<Real> &uold, const std::vector<Real> &zold, 
                                              const std::vector<Real> &unew, const std::vector<Real> &znew) {
    r.clear();
    r.resize(nx_,0.0);
    for (unsigned n = 0; n < nx_; n++) {
      // Contribution from mass & stiffness term at time step t and t-1
      r[n] += (4.0*dx_/6.0 + 0.5*dt_*2.0*nu_/dx_)*unew[n];
      r[n] += (-4.0*dx_/6.0 + 0.5*dt_*2.0*nu_/dx_)*uold[n];
      if ( n > 0 ) {
        r[n] += (dx_/6.0 - 0.5*dt_*nu_/dx_)*unew[n-1];
        r[n] -= (dx_/6.0 + 0.5*dt_*nu_/dx_)*uold[n-1];
      }
      if ( n < nx_-1 ) {
        r[n] += (dx_/6.0 - 0.5*dt_*nu_/dx_)*unew[n+1];
        r[n] -= (dx_/6.0 + 0.5*dt_*nu_/dx_)*uold[n+1];
      }
      // Contribution from nonlinear term
      if ( n > 0 ) {
        r[n] -= 0.5*dt_*unew[n-1]*(unew[n-1]+unew[n])/6.0;
        r[n] -= 0.5*dt_*uold[n-1]*(uold[n-1]+uold[n])/6.0;
      }
      if ( n < nx_-1 ){
        r[n] += 0.5*dt_*unew[n+1]*(unew[n]+unew[n+1])/6.0;
        r[n] += 0.5*dt_*uold[n+1]*(uold[n]+uold[n+1])/6.0;
      }
      // Contribution from control
      r[n] -= 0.5*dt_*dx_/6.0*(znew[n]+4.0*znew[n+1]+znew[n+2]);
      r[n] -= 0.5*dt_*dx_/6.0*(zold[n]+4.0*zold[n+1]+zold[n+2]);
      // Contribution from right-hand side
      r[n] -= dt_*dx_*f_;
    }
    // Contribution from Dirichlet boundary terms
    r[0]     -= dt_*(0.5*u0_*(unew[    0] + uold[    0])/6.0 + u0_*u0_/6.0 + nu_*u0_/dx_);
    r[nx_-1] += dt_*(0.5*u1_*(unew[nx_-1] + uold[nx_-1])/6.0 + u1_*u1_/6.0 - nu_*u1_/dx_);
  }

  void compute_pde_jacobian(std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
                      const std::vector<Real> &u) {
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(nx_,4.0*dx_/6.0 + 0.5*dt_*nu_*2.0/dx_);
    dl.clear();
    dl.resize(nx_-1,dx_/6.0 - 0.5*dt_*nu_/dx_);
    du.clear();
    du.resize(nx_-1,dx_/6.0 - 0.5*dt_*nu_/dx_);
    // Contribution from nonlinearity
    for (unsigned n = 0; n < nx_; n++) {
      if ( n < nx_-1 ) {
        dl[n] += 0.5*dt_*(-2.0*u[n]-u[n+1])/6.0;
        d[n]  += 0.5*dt_*u[n+1]/6.0;
      }
      if ( n > 0 ) {
        d[n]    -= 0.5*dt_*u[n-1]/6.0;
        du[n-1] += 0.5*dt_*(u[n-1]+2.0*u[n])/6.0;
      }
    }
    // Contribution from Dirichlet boundary conditions
    d[0]     -= 0.5*dt_*u0_/6.0;
    d[nx_-1] += 0.5*dt_*u1_/6.0;
  }

  void apply_pde_jacobian_new(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &u,
                              bool adjoint = false) {
    jv.clear();
    jv.resize(nx_,0.0);
    // Fill Jacobian times a vector
    for (unsigned n = 0; n < nx_; n++) {
      jv[n] = (4.0*dx_/6.0 + 0.5*dt_*nu_/dx_*2.0)*v[n]; // Mass & Stiffness
      if ( n > 0 ) {
        jv[n] += dx_/6.0*v[n-1]-0.5*dt_*nu_/dx_*v[n-1]; // Mass & Stiffness
        if ( adjoint ) {
          jv[n] -= 0.5*dt_*(u[n-1]/6.0*v[n]-(u[n-1]+2.0*u[n])/6.0*v[n-1]); // Nonlinearity
        } 
        else {
          jv[n] -= 0.5*dt_*(u[n-1]/6.0*v[n]+(u[n]+2.0*u[n-1])/6.0*v[n-1]); // Nonlinearity
        }
      }
      if ( n < nx_-1 ) {
        jv[n] += dx_/6.0*v[n+1]-0.5*dt_*nu_/dx_*v[n+1]; // Mass & Stiffness
        if ( adjoint ) {
          jv[n] += 0.5*dt_*(u[n+1]/6.0*v[n]-(u[n+1]+2.0*u[n])/6.0*v[n+1]); // Nonlinearity
        } 
        else {
          jv[n] += 0.5*dt_*(u[n+1]/6.0*v[n]+(u[n]+2.0*u[n+1])/6.0*v[n+1]); // Nonlinearity
        }
      }
    }
    jv[0]     -= 0.5*dt_*u0_/6.0*v[0]; // Nonlinearity
    jv[nx_-1] += 0.5*dt_*u1_/6.0*v[nx_-1]; // Nonlinearity
  }

  void apply_pde_jacobian_old(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &u,
                              bool adjoint = false) {
    jv.clear();
    jv.resize(nx_,0.0);
    // Fill Jacobian times a vector
    for (unsigned n = 0; n < nx_; n++) {
      jv[n] = (-4.0*dx_/6.0 + 0.5*dt_*nu_/dx_*2.0)*v[n]; // Mass & Stiffness
      if ( n > 0 ) {
        jv[n] += -dx_/6.0*v[n-1]-0.5*dt_*nu_/dx_*v[n-1]; // Mass & Stiffness
        if ( adjoint ) {
          jv[n] -= 0.5*dt_*(u[n-1]/6.0*v[n]-(u[n-1]+2.0*u[n])/6.0*v[n-1]); // Nonlinearity
        } 
        else {
          jv[n] -= 0.5*dt_*(u[n-1]/6.0*v[n]+(u[n]+2.0*u[n-1])/6.0*v[n-1]); // Nonlinearity
        }
      }
      if ( n < nx_-1 ) {
        jv[n] += -dx_/6.0*v[n+1]-0.5*dt_*nu_/dx_*v[n+1]; // Mass & Stiffness
        if ( adjoint ) {
          jv[n] += 0.5*dt_*(u[n+1]/6.0*v[n]-(u[n+1]+2.0*u[n])/6.0*v[n+1]); // Nonlinearity
        } 
        else {
          jv[n] += 0.5*dt_*(u[n+1]/6.0*v[n]+(u[n]+2.0*u[n+1])/6.0*v[n+1]); // Nonlinearity
        }
      }
    }
    jv[0]     -= 0.5*dt_*u0_/6.0*v[0]; // Nonlinearity
    jv[nx_-1] += 0.5*dt_*u1_/6.0*v[nx_-1]; // Nonlinearity
  }

  void apply_pde_jacobian(std::vector<Real> &jv, const std::vector<Real> &vold, const std::vector<Real> &uold,
                          const std::vector<Real> &vnew, const std::vector<Real> unew, bool adjoint = false) {
    jv.clear();
    jv.resize(nx_,0.0);
    // Fill Jacobian times a vector
    for (unsigned n = 0; n < nx_; n++) {
      jv[n] += (4.0*dx_/6.0+0.5*dt_*nu_/dx_*2.0)*vnew[n]; // Mass & Stiffness
      jv[n] -= (4.0*dx_/6.0-0.5*dt_*nu_/dx_*2.0)*vold[n]; // Mass & Stiffness
      if ( n > 0 ) {
        jv[n] += dx_/6.0*vnew[n-1]-0.5*dt_*nu_/dx_*vnew[n-1]; // Mass & Stiffness
        jv[n] -= dx_/6.0*vold[n-1]+0.5*dt_*nu_/dx_*vold[n-1]; // Mass & Stiffness
        if ( adjoint ) {
          jv[n] -= 0.5*dt_*(unew[n-1]/6.0*vnew[n]-(unew[n-1]+2.0*unew[n])/6.0*vnew[n-1]); // Nonlinearity
          jv[n] -= 0.5*dt_*(uold[n-1]/6.0*vold[n]-(uold[n-1]+2.0*uold[n])/6.0*vold[n-1]); // Nonlinearity
        } 
        else {
          jv[n] -= 0.5*dt_*(unew[n-1]/6.0*vnew[n]+(unew[n]+2.0*unew[n-1])/6.0*vnew[n-1]); // Nonlinearity
          jv[n] -= 0.5*dt_*(uold[n-1]/6.0*vold[n]+(uold[n]+2.0*uold[n-1])/6.0*vold[n-1]); // Nonlinearity
        }
      }
      if ( n < nx_-1 ) {
        jv[n] += dx_/6.0*vnew[n+1]-0.5*dt_*nu_/dx_*vnew[n+1]; // Mass & Stiffness
        jv[n] -= dx_/6.0*vold[n+1]+0.5*dt_*nu_/dx_*vold[n+1]; // Mass & Stiffness
        if ( adjoint ) {
          jv[n] += 0.5*dt_*(unew[n+1]/6.0*vnew[n]-(unew[n+1]+2.0*unew[n])/6.0*vnew[n+1]); // Nonlinearity
          jv[n] += 0.5*dt_*(uold[n+1]/6.0*vold[n]-(uold[n+1]+2.0*uold[n])/6.0*vold[n+1]); // Nonlinearity
        } 
        else {
          jv[n] += 0.5*dt_*(unew[n+1]/6.0*vnew[n]+(unew[n]+2.0*unew[n+1])/6.0*vnew[n+1]); // Nonlinearity
          jv[n] += 0.5*dt_*(uold[n+1]/6.0*vold[n]+(uold[n]+2.0*uold[n+1])/6.0*vold[n+1]); // Nonlinearity
        }
      }
    }
    jv[0]     -= 0.5*dt_*u0_/6.0*vnew[0]; // Nonlinearity
    jv[0]     -= 0.5*dt_*u0_/6.0*vold[0]; // Nonlinearity
    jv[nx_-1] += 0.5*dt_*u1_/6.0*vnew[nx_-1]; // Nonlinearity
    jv[nx_-1] += 0.5*dt_*u1_/6.0*vold[nx_-1]; // Nonlinearity
  }

  void apply_pde_hessian(std::vector<Real> &hv, const std::vector<Real> &wold, const std::vector<Real> &vold,
                                                const std::vector<Real> &wnew, const std::vector<Real> &vnew) {
    hv.clear();
    hv.resize(nx_,0.0);
    for (unsigned n = 0; n < nx_; n++) {
      if ( n > 0 ) {
        hv[n] += 0.5*dt_*((wnew[n-1]*(vnew[n-1]+2.0*vnew[n]) - wnew[n]*vnew[n-1])/6.0);
        hv[n] += 0.5*dt_*((wold[n-1]*(vold[n-1]+2.0*vold[n]) - wold[n]*vold[n-1])/6.0);
      }
      if ( n < nx_-1 ){
        hv[n] += 0.5*dt_*((wnew[n]*vnew[n+1] - wnew[n+1]*(2.0*vnew[n]+vnew[n+1]))/6.0);
        hv[n] += 0.5*dt_*((wold[n]*vold[n+1] - wold[n+1]*(2.0*vold[n]+vold[n+1]))/6.0);
      }
    }
  }

  void apply_control_jacobian(std::vector<Real> &jv, const std::vector<Real> &v, bool adjoint = false) {
    jv.clear();
    unsigned dim = ((adjoint == true) ? nx_+2 : nx_);
    jv.resize(dim,0.0);
    for (unsigned n = 0; n < dim; n++) {
      if ( adjoint ) {
        if ( n == 0 ) {
          jv[n] = -0.5*dt_*dx_/6.0*v[n];
        }
        else if ( n == 1 ) {
          jv[n] = -0.5*dt_*dx_/6.0*(4.0*v[n-1]+v[n]);
        }
        else if ( n == nx_ ) {
          jv[n] = -0.5*dt_*dx_/6.0*(4.0*v[n-1]+v[n-2]);
        }
        else if ( n == nx_+1 ) {
          jv[n] = -0.5*dt_*dx_/6.0*v[n-2];
        }
        else {
          jv[n] = -0.5*dt_*dx_/6.0*(v[n-2]+4.0*v[n-1]+v[n]);
        }
      }
      else {
        jv[n] -= 0.5*dt_*dx_/6.0*(v[n]+4.0*v[n+1]+v[n+2]);
      }
    }
  }

  void run_newton(std::vector<Real> &u, const std::vector<Real> &znew,
            const std::vector<Real> &uold, const std::vector<Real> &zold) { 
    u.clear();
    u.resize(nx_,0.0);
    // Compute residual and residual norm
    std::vector<Real> r(nx_,0.0);
    compute_residual(r,uold,zold,u,znew);
    Real rnorm = compute_norm(r);
    // Define tolerances
    Real rtol  = 1.e2*ROL::ROL_EPSILON<Real>();
    unsigned maxit = 500;
    // Initialize Jacobian storage
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    std::vector<Real> du(nx_-1,0.0);
    // Iterate Newton's method
    Real alpha = 1.0, tmp = 0.0;
    std::vector<Real> s(nx_,0.0);
    std::vector<Real> utmp(nx_,0.0);
    for (unsigned i = 0; i < maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      compute_pde_jacobian(dl,d,du,u);
      // Solve Newton system
      linear_solve(s,dl,d,du,r);
      // Perform line search
      tmp = rnorm;
      alpha = 1.0;
      utmp.assign(u.begin(),u.end());
      update(utmp,s,-alpha);
      compute_residual(r,uold,zold,utmp,znew);
      rnorm = compute_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON<Real>()) ) {
        alpha *= 0.5;
        utmp.assign(u.begin(),u.end());
        update(utmp,s,-alpha);
        compute_residual(r,uold,zold,utmp,znew);
        rnorm = compute_norm(r); 
      }
      // Update iterate
      u.assign(utmp.begin(),utmp.end());
      if ( rnorm < rtol ) {
        break;
      }
    }
  }

  void linear_solve(std::vector<Real> &u, 
              const std::vector<Real> &dl, const std::vector<Real> &d, const std::vector<Real> &du, 
              const std::vector<Real> &r, const bool transpose = false) {
    bool useLAPACK = true;
    if ( useLAPACK ) { // DIRECT SOLVE: USE LAPACK
      u.assign(r.begin(),r.end());
      // Store matrix diagonal & off-diagonals.
      std::vector<Real> Dl(dl);
      std::vector<Real> Du(du);
      std::vector<Real> D(d);
      // Perform LDL factorization
      Teuchos::LAPACK<int,Real> lp;
      std::vector<Real> Du2(nx_-2,0.0);
      std::vector<int> ipiv(nx_,0);
      int info;
      int ldb  = nx_;
      int nhrs = 1;
      lp.GTTRF(nx_,&Dl[0],&D[0],&Du[0],&Du2[0],&ipiv[0],&info);
      char trans = ((transpose == true) ? 'T' : 'N');
      lp.GTTRS(trans,nx_,nhrs,&Dl[0],&D[0],&Du[0],&Du2[0],&ipiv[0],&u[0],ldb,&info);
    }
    else { // ITERATIVE SOLVE: USE GAUSS-SEIDEL
      u.clear();
      u.resize(nx_,0.0);
      unsigned maxit = 100;
      Real rtol  = std::min(1.e-12,1.e-4*std::sqrt(dot(r,r)));
      //Real rtol  = 1e-6;
      Real resid = 0.0;
      Real rnorm = 10.0*rtol;
      for (unsigned i = 0; i < maxit; i++) {
        for (unsigned n = 0; n < nx_; n++) {
          u[n] = r[n]/d[n];
          if ( n < nx_-1 ) {
            u[n] -= ((transpose == false) ? du[n] : dl[n])*u[n+1]/d[n];
          }
          if ( n > 0 ) {
            u[n] -= ((transpose == false) ? dl[n-1] : du[n-1])*u[n-1]/d[n];
          }
        }
        // Compute Residual
        rnorm = 0.0;
        for (unsigned n = 0; n < nx_; n++) {
          resid = r[n] - d[n]*u[n];
          if ( n < nx_-1 ) {
            resid -= ((transpose == false) ? du[n] : dl[n])*u[n+1];
          }
          if ( n > 0 ) {
            resid -= ((transpose == false) ? dl[n-1] : du[n-1])*u[n-1];
          }
          rnorm += resid*resid;
        }
        rnorm = std::sqrt(rnorm);
        if ( rnorm < rtol ) {
          //std::cout << "\ni = " << i+1 << "  rnorm = " << rnorm << "\n";
          break;
        }
      }
    }
  }

public:

  Constraint_BurgersControl(int nx = 128, int nt = 100, Real T = 1, 
                                    Real nu = 1.e-2, Real u0 = 0.0, Real u1 = 0.0, Real f = 0.0) 
    : nx_((unsigned)nx), nt_((unsigned)nt), T_(T), nu_(nu), u0_(u0), u1_(u1), f_(f) {
    dx_ = 1.0/((Real)nx+1.0);
    dt_ = T/((Real)nt);
    u_init_.clear();
    u_init_.resize(nx_,0.0);
    Real x = 0.0;
    for (unsigned n = 0; n < nx_; n++) {
      x = (Real)(n+1)*dx_;
      u_init_[n] = ((x <= 0.5) ? 1.0 : 0.0);
    }
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp =
      dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Initialize storage
    std::vector<Real> C(nx_,0.0);
    std::vector<Real> uold(u_init_);
    std::vector<Real> unew(u_init_);
    std::vector<Real> zold(nx_+2,0.0);
    std::vector<Real> znew(nx_+2,0.0);
    // Copy initial control
    for (unsigned n = 0; n < nx_+2; n++) {
      zold[n] = (*zp)[n];
    }
    // Evaluate residual
    for (unsigned t = 0; t < nt_; t++) {
      // Copy current state at time step t
      for (unsigned n = 0; n < nx_; n++) {
        unew[n] = (*up)[t*nx_+n];
      }
      // Copy current control at time step t
      for (unsigned n = 0; n < nx_+2; n++) {
        znew[n] = (*zp)[(t+1)*(nx_+2)+n];
      }
      // Compute residual at time step t
      compute_residual(C,uold,zold,unew,znew);
      // Store residual at time step t
      for (unsigned n = 0; n < nx_; n++) {
        (*cp)[t*nx_+n] = C[n];
      }
      // Copy previous state and control at time step t+1
      uold.assign(unew.begin(),unew.end());
      zold.assign(znew.begin(),znew.end());
    }
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > up =
      dynamic_cast<ROL::StdVector<Real>&>(u).getVector();
    up->assign(up->size(),z.norm()/up->size());
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Initialize storage
    std::vector<Real> uold(u_init_);
    std::vector<Real> unew(u_init_);
    std::vector<Real> zold(nx_+2,0.0);
    std::vector<Real> znew(nx_+2,0.0);
    // Copy initial control
    for (unsigned n = 0; n < nx_+2; n++) {
      zold[n] = (*zp)[n];
    }
    // Time step using the trapezoidal rule
    for (unsigned t = 0; t < nt_; t++) {
      // Copy current control at time step t
      for (unsigned n = 0; n < nx_+2; n++) {
        znew[n] = (*zp)[(t+1)*(nx_+2)+n];
      }
      // Solve nonlinear system at time step t
      run_newton(unew,znew,uold,zold);
      // store state at time step t
      for (unsigned n = 0; n < nx_; n++) {
        (*up)[t*nx_+n] = unew[n];
      }
      // Copy previous state and control at time step t+1
      uold.assign(unew.begin(),unew.end());
      zold.assign(znew.begin(),znew.end());
    }
    value(c,u,z,tol);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    std::vector<Real> J(nx_,0.0);
    std::vector<Real> vold(nx_,0.0);
    std::vector<Real> vnew(nx_,0.0);
    std::vector<Real> uold(nx_,0.0);
    std::vector<Real> unew(nx_,0.0);
    for (unsigned t = 0; t < nt_; t++) {
      for (unsigned n = 0; n < nx_; n++) {
        unew[n] = (*up)[t*nx_+n];
        vnew[n] = (*vp)[t*nx_+n];
      }
      apply_pde_jacobian(J,vold,uold,vnew,unew);
      for (unsigned n = 0; n < nx_; n++) {
        (*jvp)[t*nx_+n] = J[n];
      }  
      vold.assign(vnew.begin(),vnew.end());
      uold.assign(unew.begin(),unew.end());
    }
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    std::vector<Real> vnew(nx_+2,0.0);
    std::vector<Real> vold(nx_+2,0.0);
    std::vector<Real> jnew(nx_,0.0);
    std::vector<Real> jold(nx_,0.0);
    for (unsigned n = 0; n < nx_+2; n++) {
      vold[n] = (*vp)[n];
    }
    apply_control_jacobian(jold,vold);
    for (unsigned t = 0; t < nt_; t++) {
      for (unsigned n = 0; n < nx_+2; n++) {
        vnew[n] = (*vp)[(t+1)*(nx_+2)+n];
      }
      apply_control_jacobian(jnew,vnew);
      for (unsigned n = 0; n < nx_; n++) {
        // Contribution from control
        (*jvp)[t*nx_+n] = jnew[n] + jold[n];
      }
      jold.assign(jnew.begin(),jnew.end());
    }
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ijvp =
      dynamic_cast<ROL::StdVector<Real>&>(ijv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    std::vector<Real> J(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> s(nx_,0.0);
    std::vector<Real> uold(nx_,0.0);
    std::vector<Real> unew(nx_,0.0);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    std::vector<Real> du(nx_-1,0.0);
    for (unsigned t = 0; t < nt_; t++) {
      // Build rhs.
      for (unsigned n = 0; n < nx_; n++) {
        r[n] = (*vp)[t*nx_+n];
        unew[n] = (*up)[t*nx_+n];
      }
      apply_pde_jacobian_old(J,s,uold);
      update(r,J,-1.0);
      // Build Jacobian.
      compute_pde_jacobian(dl,d,du,unew);
      // Solve linear system.
      linear_solve(s,dl,d,du,r);
      // Copy solution.
      for (unsigned n = 0; n < nx_; n++) {
        (*ijvp)[t*nx_+n] = s[n];
      }
      uold.assign(unew.begin(),unew.end());
    }
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    std::vector<Real> J(nx_,0.0);
    std::vector<Real> vold(nx_,0.0);
    std::vector<Real> vnew(nx_,0.0);
    std::vector<Real> unew(nx_,0.0);
    for (unsigned t = nt_; t > 0; t--) {
      for (unsigned n = 0; n < nx_; n++) {
        unew[n] = (*up)[(t-1)*nx_+n];
        vnew[n] = (*vp)[(t-1)*nx_+n];
      }
      apply_pde_jacobian(J,vold,unew,vnew,unew,true);
      for (unsigned n = 0; n < nx_; n++) {
        (*jvp)[(t-1)*nx_+n] = J[n];
      }  
      vold.assign(vnew.begin(),vnew.end());
    }
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    std::vector<Real> vnew(nx_,0.0);
    std::vector<Real> vold(nx_,0.0);
    std::vector<Real> jnew(nx_+2,0.0);
    std::vector<Real> jold(nx_+2,0.0);
    for (unsigned t = nt_+1; t > 0; t--) {
      for (unsigned n = 0; n < nx_; n++) {
        if ( t > 1 ) {
          vnew[n] = (*vp)[(t-2)*nx_+n];
        }
        else {
          vnew[n] = 0.0;
        }
      }
      apply_control_jacobian(jnew,vnew,true);
      for (unsigned n = 0; n < nx_+2; n++) {
        // Contribution from control
        (*jvp)[(t-1)*(nx_+2)+n] = jnew[n] + jold[n];
      }
      jold.assign(jnew.begin(),jnew.end());
    }
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ijvp =
      dynamic_cast<ROL::StdVector<Real>&>(ijv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    std::vector<Real> J(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> s(nx_,0.0);
    std::vector<Real> unew(nx_,0.0);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    std::vector<Real> du(nx_-1,0.0);
    for (unsigned t = nt_; t > 0; t--) {
      // Build rhs.
      for (unsigned n = 0; n < nx_; n++) {
        r[n] = (*vp)[(t-1)*nx_+n];
        unew[n] = (*up)[(t-1)*nx_+n];
      }
      apply_pde_jacobian_old(J,s,unew,true);
      update(r,J,-1.0);
      // Build Jacobian.
      compute_pde_jacobian(dl,d,du,unew);
      // Solve linear system.
      linear_solve(s,dl,d,du,r,true);
      // Copy solution.
      for (unsigned n = 0; n < nx_; n++) {
        (*ijvp)[(t-1)*nx_+n] = s[n];
      }
    }
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &hwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > hwvp =
      dynamic_cast<ROL::StdVector<Real>&>(hwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp =
      dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    std::vector<Real> snew(nx_,0.0);
    std::vector<Real> wnew(nx_,0.0);
    std::vector<Real> wold(nx_,0.0); 
    std::vector<Real> vnew(nx_,0.0);
    for (unsigned t = nt_; t > 0; t--) {
      for (unsigned n = 0; n < nx_; n++) {
        vnew[n] = (*vp)[(t-1)*nx_+n];
        wnew[n] = (*wp)[(t-1)*nx_+n];
      }
      apply_pde_hessian(snew,wold,vnew,wnew,vnew);
      for (unsigned n = 0; n < nx_; n++) {
        (*hwvp)[(t-1)*nx_+n] = snew[n];
      }
      wold.assign(wnew.begin(),wnew.end());
    }
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class Objective_BurgersControl : public ROL::Objective_SimOpt<Real> {
private:
  Real alpha_;     // Penalty Parameter

  unsigned nx_;    // Number of interior nodes
  Real     dx_;    // Mesh spacing (i.e. 1/(nx+1))
  unsigned nt_;    // Number of time steps
  Real     dt_;    // Time step size
  Real     T_;     // Final time

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  Real evaluate_target(Real x) {
    Real val = 0.0;
    int example = 1;
    switch (example) {
      case 1:  val = ((x<=0.5) ? 1.0 : 0.0);         break;
      case 2:  val = 1.0;                            break;
      case 3:  val = std::abs(std::sin(8.0*M_PI*x)); break;
      case 4:  val = std::exp(-0.5*(x-0.5)*(x-0.5)); break;
    }
    return val;
  }

  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) {
    Real ip = 0.0;
    Real c = ((x.size()==nx_) ? 4.0 : 2.0);
    for (unsigned i=0; i<x.size(); i++) {
      if ( i == 0 ) {
        ip += dx_/6.0*(c*x[i] + x[i+1])*y[i];
      }
      else if ( i == x.size()-1 ) {
        ip += dx_/6.0*(x[i-1] + c*x[i])*y[i];
      }
      else {
        ip += dx_/6.0*(x[i-1] + 4.0*x[i] + x[i+1])*y[i];
      }
    }
    return ip;
  }

  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
    Mu.resize(u.size(),0.0);
    Real c = ((u.size()==nx_) ? 4.0 : 2.0);
    for (unsigned i=0; i<u.size(); i++) {
      if ( i == 0 ) {
        Mu[i] = dx_/6.0*(c*u[i] + u[i+1]);
      }
      else if ( i == u.size()-1 ) {
        Mu[i] = dx_/6.0*(u[i-1] + c*u[i]);
      }
      else {
        Mu[i] = dx_/6.0*(u[i-1] + 4.0*u[i] + u[i+1]);
      }
    }
  }
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_BurgersControl(Real alpha = 1.e-4, int nx = 128, int nt = 100, Real T = 1.0) 
    : alpha_(alpha), nx_((unsigned)nx), nt_((unsigned)nt), T_(T) {
    dx_ = 1.0/((Real)nx+1.0);
    dt_ = T/((Real)nt);
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE RESIDUAL
    std::vector<Real> U(nx_,0.0);
    std::vector<Real> G(nx_,0.0);
    std::vector<Real> Z(nx_+2,0.0);
    for (unsigned n = 0; n < nx_+2; n++) {
      Z[n] = (*zp)[n];
    }
    Real ss  = 0.5*dt_;
    Real val = 0.5*ss*alpha_*dot(Z,Z);
    for (unsigned t = 0; t < nt_; t++) {
      ss = ((t < nt_-1) ? dt_ : 0.5*dt_);
      for (unsigned n = 0; n < nx_; n++) {
        U[n] = (*up)[t*nx_+n]-evaluate_target((Real)(n+1)*dx_);
        G[n] = evaluate_target((Real)(n+1)*dx_);
      }
      val += 0.5*ss*dot(U,U);
      val -= 0.5*ss*dot(G,G); // subtract constant term
      for (unsigned n = 0; n < nx_+2; n++) {
        Z[n] = (*zp)[(t+1)*(nx_+2)+n];
      }
      val += 0.5*ss*alpha_*dot(Z,Z);
    }
    return val;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp =
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    // COMPUTE GRADIENT WRT U
    std::vector<Real> U(nx_,0.0);
    std::vector<Real> M(nx_,0.0);
    Real ss = dt_;
    for (unsigned t = 0; t < nt_; t++) {
      ss = ((t < nt_-1) ? dt_ : 0.5*dt_);
      for (unsigned n = 0; n < nx_; n++) {
        U[n] = (*up)[t*nx_+n]-evaluate_target((Real)(n+1)*dx_);
      }
      apply_mass(M,U);
      for (unsigned n = 0; n < nx_; n++) {
        (*gp)[t*nx_+n] = ss*M[n];
      }
    }
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp =
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE GRADIENT WRT Z
    std::vector<Real> Z(nx_+2,0.0);
    std::vector<Real> M(nx_+2,0.0);
    Real ss = dt_;
    for (unsigned t = 0; t < nt_+1; t++) {
      ss = ((t < nt_ && t > 0) ? dt_ : 0.5*dt_);
      for (unsigned n = 0; n < nx_+2; n++) {
        Z[n] = (*zp)[t*(nx_+2)+n];
      }
      apply_mass(M,Z);
      for (unsigned n = 0; n < nx_+2; n++) {
        (*gp)[t*(nx_+2)+n] = ss*alpha_*M[n];
      }
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp =
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE GRADIENT WRT U
    std::vector<Real> V(nx_,0.0);
    std::vector<Real> M(nx_,0.0);
    Real ss  = 0.5*dt_;
    for (unsigned t = 0; t < nt_; t++) {
      ss = ((t < nt_-1) ? dt_ : 0.5*dt_);
      for (unsigned n = 0; n < nx_; n++) {
        V[n] = (*vp)[t*nx_+n];
      }
      apply_mass(M,V);
      for (unsigned n = 0; n < nx_; n++) {
        (*hvp)[t*nx_+n] = ss*M[n];
      }
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp = ROL::constPtrCast<std::vector<Real> >(
      (dynamic_cast<const ROL::StdVector<Real>&>(hv)).getVector());
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE GRADIENT WRT Z
    std::vector<Real> V(nx_+2,0.0);
    std::vector<Real> M(nx_+2,0.0);
    Real ss = 0.0;
    for (unsigned t = 0; t < nt_+1; t++) {
      ss = ((t < nt_ && t > 0) ? dt_ : 0.5*dt_);
      for (unsigned n = 0; n < nx_+2; n++) {
        V[n] = (*vp)[t*(nx_+2)+n];
      }
      apply_mass(M,V);
      for (unsigned n = 0; n < nx_+2; n++) {
        (*hvp)[t*(nx_+2)+n] = ss*alpha_*M[n];
      }
    }
  }
};
