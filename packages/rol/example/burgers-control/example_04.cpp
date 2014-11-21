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

/*! \file  example_04.cpp
    \brief Shows how to solve an optimal control problem constrained by 
           Burgers' equation with the SimOpt interface.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStepSQP.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>
#include <ctime>

template<class Real>
class EqualityConstraint_BurgersControl : public ROL::EqualityConstraint_SimOpt<Real> {
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
    Real rtol  = 1.e2*ROL::ROL_EPSILON;
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
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON) ) {
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
    bool useLAPACK = false;
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

  EqualityConstraint_BurgersControl(int nx = 128, int nt = 100, Real T = 1, 
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
    Teuchos::RCP<std::vector<Real> > cp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(c)).getVector());
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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

  void solve(ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > up =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(u)).getVector());
    up->assign(up->size(),z.norm()/up->size());
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
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
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > ijvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(ijv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
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
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(ajv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
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
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > ijvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(ijv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
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
    Teuchos::RCP<std::vector<Real> > hwvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hwv)).getVector());
    Teuchos::RCP<const std::vector<Real> > wp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(w))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
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
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE RESIDUAL
    std::vector<Real> U(nx_,0.0);
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
      }
      val += 0.5*ss*dot(U,U); 
      for (unsigned n = 0; n < nx_+2; n++) {
        Z[n] = (*zp)[(t+1)*(nx_+2)+n];
      }
      val += 0.5*ss*alpha_*dot(Z,Z);
    }
    return val;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<std::vector<Real> > gp = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(g)).getVector());
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
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
    Teuchos::RCP<std::vector<Real> > gp = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(g)).getVector());
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > hvp = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(hv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
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
    Teuchos::RCP<std::vector<Real> > hvp = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(hv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
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

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {
    // Initialize full objective function.
    int nx      = 20;    // Set spatial discretization.
    int nt      = 20;    // Set temporal discretization.
    RealT T     = 1.0;   // Set end time.
    RealT alpha = 0.05;  // Set penalty parameter.
    RealT nu    = 1.e-2; // Set viscosity parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx,nt,T);
    // Initialize equality constraints
    EqualityConstraint_BurgersControl<RealT> con(nx, nt, T, nu);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT> ((nx+2)*(nt+1), 1.0) );
    Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp( new std::vector<RealT> ((nx+2)*(nt+1), 1.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp( new std::vector<RealT> ((nx+2)*(nt+1), 1.0) );
    for (int i=0; i<(nx+2)*(nt+1); i++) {
      (*z_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> z(z_rcp);
    ROL::StdVector<RealT> gz(gz_rcp);
    ROL::StdVector<RealT> yz(yz_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > zp  = Teuchos::rcp(&z,false);
    Teuchos::RCP<ROL::Vector<RealT> > gzp = Teuchos::rcp(&gz,false);
    Teuchos::RCP<ROL::Vector<RealT> > yzp = Teuchos::rcp(&yz,false);

    Teuchos::RCP<std::vector<RealT> > u_rcp  = Teuchos::rcp( new std::vector<RealT> (nx*nt, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gu_rcp = Teuchos::rcp( new std::vector<RealT> (nx*nt, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yu_rcp = Teuchos::rcp( new std::vector<RealT> (nx*nt, 1.0) );
    for (int i=0; i<nx*nt; i++) {
      (*u_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yu_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> u(u_rcp);
    ROL::StdVector<RealT> gu(gu_rcp);
    ROL::StdVector<RealT> yu(yu_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > up  = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > gup = Teuchos::rcp(&gu,false);
    Teuchos::RCP<ROL::Vector<RealT> > yup = Teuchos::rcp(&yu,false);

    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (nx*nt, 1.0) );
    Teuchos::RCP<std::vector<RealT> > l_rcp = Teuchos::rcp( new std::vector<RealT> (nx*nt, 1.0) );
    ROL::StdVector<RealT> c(c_rcp);
    ROL::StdVector<RealT> l(l_rcp);

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);
    // Check derivatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);
    con.checkApplyJacobian(x,y,c,true,*outStream);
    con.checkApplyAdjointJacobian(x,yu,c,x,true,*outStream);
    con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);
    // Check consistency of Jacobians and adjoint Jacobians.
    con.checkJacobian_1(c,yu,u,z,true,*outStream);
    con.checkJacobian_2(c,yz,u,z,true,*outStream);
    // Check consistency of solves.
    con.checkSolve(u,z,c,true,*outStream);
    con.checkInverseJacobian_1(c,yu,u,z,true,*outStream);
    con.checkInverseAdjointJacobian_1(yu,c,u,z,true,*outStream);

    // Initialize reduced objective function.
    Teuchos::RCP<std::vector<RealT> > p_rcp  = Teuchos::rcp( new std::vector<RealT> (nx*nt, 1.0) );
    ROL::StdVector<RealT> p(p_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > pp  = Teuchos::rcp(&p,false);
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > pobj = Teuchos::rcp(&obj,false);
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > pcon = Teuchos::rcp(&con,false);
    ROL::Reduced_Objective_SimOpt<RealT> robj(pobj,pcon,up,pp);
    // Check derivatives.
    robj.checkGradient(z,z,yz,true,*outStream);
    robj.checkHessVec(z,z,yz,true,*outStream);
    // Optimization 
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist_tr = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist_tr) );

    // Trust Region Newton.
    RealT gtol  = 1e-14;  // norm of gradient tolerance
    RealT stol  = 1e-16;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    ROL::StatusTest<RealT> status_tr(gtol, stol, maxit);    
    ROL::TrustRegionStep<RealT> step_tr(*parlist_tr);
    ROL::DefaultAlgorithm<RealT> algo_tr(step_tr,status_tr,false);
    z.zero();
    std::clock_t timer_tr = std::clock();
    algo_tr.run(z,robj,true,*outStream);
    *outStream << "Trust-Region Newton required " << (std::clock()-timer_tr)/(RealT)CLOCKS_PER_SEC
               << " seconds.\n";

    // SQP.
    RealT ctol = 1.e-14;
    ROL::StatusTestSQP<RealT> status_sqp(gtol,ctol,stol,maxit);
    ROL::CompositeStepSQP<RealT> step_sqp(*parlist_tr);
    ROL::DefaultAlgorithm<RealT> algo_sqp(step_sqp,status_sqp,false);
    x.zero();
    std::clock_t timer_sqp = std::clock();
    algo_sqp.run(x,g,l,c,obj,con,true,*outStream);
    *outStream << "Composite-Step SQP required " << (std::clock()-timer_sqp)/(RealT)CLOCKS_PER_SEC
               << " seconds.\n";
 
    std::ofstream control;
    control.open("control.txt");
    for (int t = 0; t < nt+1; t++) {
      for (int n = 0; n < nx+2; n++) {
        control << (RealT)t/(RealT)nt       << "  " 
                << (RealT)n/((RealT)(nx+1)) << "  " 
                << (*z_rcp)[t*(nx+2)+n]     << "\n";
      }
    } 
    control.close();

    std::ofstream state;
    state.open("state.txt");
    for (int t = 0; t < nt; t++) {
      for (int n = 0; n < nx; n++) {
        state << (RealT)(t+1)/(RealT)nt       << "  " 
              << (RealT)(n+1)/((RealT)(nx+1)) << "  " 
              << (*u_rcp)[t*nx+n]             << "\n";
      }
    } 
    state.close();
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

