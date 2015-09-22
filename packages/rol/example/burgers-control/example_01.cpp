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

/*! \file  example_01.cpp
    \brief Shows how to solve an optimal control problem constrained by 
           steady Burgers' equation with bound constraints.
*/

#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

template<class Real>
class Objective_BurgersControl : public ROL::Objective<Real> {
private:
  Real alpha_; // Penalty Parameter

  int  nx_;    // Number of interior nodes
  Real dx_;    // Mesh spacing (i.e. 1/(nx+1))

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  Real evaluate_target(Real x) {
    Real val = 0.0;
    int example = 2;
    switch (example) {
      case 1:  val = ((x<0.5) ? 1.0 : 0.0);          break;
      case 2:  val = 1.0;                            break;
      case 3:  val = std::abs(std::sin(8.0*M_PI*x)); break;
      case 4:  val = std::exp(-0.5*(x-0.5)*(x-0.5)); break;
    }
    return val;
  }

  Real compute_norm(const std::vector<Real> &r) {
    return std::sqrt(this->dot(r,r));
  }

  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) {
    Real ip = 0.0;
    Real c = (((int)x.size()==this->nx_) ? 4.0 : 2.0);
    for (unsigned i=0; i<x.size(); i++) {
      if ( i == 0 ) {
        ip += this->dx_/6.0*(c*x[i] + x[i+1])*y[i];
      }
      else if ( i == x.size()-1 ) {
        ip += this->dx_/6.0*(x[i-1] + c*x[i])*y[i];
      }
      else {
        ip += this->dx_/6.0*(x[i-1] + 4.0*x[i] + x[i+1])*y[i];
      }
    }
    return ip;
  }

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    for (unsigned i=0; i<u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void scale(std::vector<Real> &u, const Real alpha=0.0) {
    for (unsigned i=0; i<u.size(); i++) {
      u[i] *= alpha;
    }
  }

  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
    Mu.resize(u.size(),0.0);
    Real c = (((int)u.size()==this->nx_) ? 4.0 : 2.0);
    for (unsigned i=0; i<u.size(); i++) {
      if ( i == 0 ) {
        Mu[i] = this->dx_/6.0*(c*u[i] + u[i+1]);
      }
      else if ( i == u.size()-1 ) {
        Mu[i] = this->dx_/6.0*(u[i-1] + c*u[i]);
      }
      else {
        Mu[i] = this->dx_/6.0*(u[i-1] + 4.0*u[i] + u[i+1]);
      }
    }
  }

  void compute_residual(std::vector<Real> &r, const std::vector<Real> &u, 
                  const std::vector<Real> &z, const std::vector<Real> &param) {
    r.clear();
    r.resize(this->nx_,0.0);
    Real nu = std::pow(10.0,param[0]-2.0);
    Real f  = param[1]/100.0;
    Real u0 = 1.0+param[2]/1000.0;
    Real u1 = param[3]/1000.0;
    for (int i=0; i<this->nx_; i++) {
      // Contribution from stiffness term
      if ( i==0 ) {
        r[i] = nu/this->dx_*(2.0*u[i]-u[i+1]);
      }
      else if (i==this->nx_-1) {
        r[i] = nu/this->dx_*(2.0*u[i]-u[i-1]);
      }
      else {
        r[i] = nu/this->dx_*(2.0*u[i]-u[i-1]-u[i+1]);
      }
      // Contribution from nonlinear term
      if (i<this->nx_-1){
        r[i] += u[i+1]*(u[i]+u[i+1])/6.0;
      }
      if (i>0) {
        r[i] -= u[i-1]*(u[i-1]+u[i])/6.0;
      }
      // Contribution from control
      r[i] -= this->dx_/6.0*(z[i]+4.0*z[i+1]+z[i+2]);
      // Contribution from right-hand side
      r[i] -= this->dx_*f;
    }
    // Contribution from Dirichlet boundary terms
    r[0]           -= u0*u[          0]/6.0 + u0*u0/6.0 + nu*u0/this->dx_;
    r[this->nx_-1] += u1*u[this->nx_-1]/6.0 + u1*u1/6.0 - nu*u1/this->dx_;
  }

  void compute_pde_jacobian(std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
                      const std::vector<Real> &u, const std::vector<Real> &param) {
    Real nu = std::pow(10.0,param[0]-2.0);
    Real u0 = 1.0+param[2]/1000.0;
    Real u1 = param[3]/1000.0;
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(this->nx_,nu*2.0/this->dx_);
    dl.clear();
    dl.resize(this->nx_-1,-nu/this->dx_);
    du.clear();
    du.resize(this->nx_-1,-nu/this->dx_);
    // Contribution from nonlinearity
    for (int i=0; i<this->nx_; i++) {
      if (i<this->nx_-1) {
        dl[i] += (-2.0*u[i]-u[i+1])/6.0;
        d[i]  += u[i+1]/6.0;
      }
      if (i>0) {
        d[i]    += -u[i-1]/6.0;
        du[i-1] += (u[i-1]+2.0*u[i])/6.0;
      }
    }
    // Contribution from Dirichlet boundary conditions
    d[0]           -= u0/6.0;
    d[this->nx_-1] += u1/6.0;
  }

  void add_pde_hessian(std::vector<Real> &r, const std::vector<Real> &u, const std::vector<Real> &p, 
                 const std::vector<Real> &s, Real alpha = 1.0) {
    for (int i=0; i<this->nx_; i++) {
      // Contribution from nonlinear term
      if (i<this->nx_-1){
        r[i] += alpha*(p[i]*s[i+1] - p[i+1]*(2.0*s[i]+s[i+1]))/6.0;
      }
      if (i>0) {
        r[i] += alpha*(p[i-1]*(s[i-1]+2.0*s[i]) - p[i]*s[i-1])/6.0;
      }
    }
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
              const std::vector<Real> &r, const bool transpose = false) {
    u.assign(r.begin(),r.end());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    std::vector<Real> du2(this->nx_-2,0.0);
    std::vector<int> ipiv(this->nx_,0);
    int info;
    int ldb  = this->nx_;
    int nhrs = 1;
    lp.GTTRF(this->nx_,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&info);
    char trans = 'N';
    if ( transpose ) { 
      trans = 'T';
    }
    lp.GTTRS(trans,this->nx_,nhrs,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&u[0],ldb,&info);
  }

  void run_newton(std::vector<Real> &u, const std::vector<Real> &z, const std::vector<Real> &param) {
    // Compute residual and residual norm
    std::vector<Real> r(u.size(),0.0);
    this->compute_residual(r,u,z,param);
    Real rnorm = this->compute_norm(r);
    // Define tolerances
    Real tol   = 1.e2*ROL::ROL_EPSILON;
    Real maxit = 500;
    // Initialize Jacobian storage
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> dl(this->nx_-1,0.0);
    std::vector<Real> du(this->nx_-1,0.0);
    // Iterate Newton's method
    Real alpha = 1.0, tmp = 0.0;
    std::vector<Real> s(this->nx_,0.0);
    std::vector<Real> utmp(this->nx_,0.0);
    for (int i=0; i<maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      this->compute_pde_jacobian(dl,d,du,u,param);
      // Solve Newton system
      this->linear_solve(s,dl,d,du,r);
      // Perform line search
      tmp = rnorm;
      alpha = 1.0;
      utmp.assign(u.begin(),u.end());
      this->update(utmp,s,-alpha);
      this->compute_residual(r,utmp,z,param);
      rnorm = this->compute_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON) ) {
        alpha /= 2.0;
        utmp.assign(u.begin(),u.end());
        this->update(utmp,s,-alpha);
        this->compute_residual(r,utmp,z,param);
        rnorm = this->compute_norm(r); 
      }
      // Update iterate
      u.assign(utmp.begin(),utmp.end());
      if ( rnorm < tol ) {
        break;
      }
    }
  }
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_BurgersControl(Real alpha = 1.e-4, int nx = 128) : alpha_(alpha), nx_(nx) {
    dx_ = 1.0/((Real)nx+1.0);
  }

  void solve_state(std::vector<Real> &u, const std::vector<Real> &z, const std::vector<Real> &param) {
    u.clear();
    u.resize(this->nx_,1.0);
    this->run_newton(u,z,param);
  }

  void solve_adjoint(std::vector<Real> &p, const std::vector<Real> &u, const std::vector<Real> &param) {
    // Initialize State Storage
    p.clear();
    p.resize(this->nx_);
    // Get PDE Jacobian
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> du(this->nx_-1,0.0);
    std::vector<Real> dl(this->nx_-1,0.0);
    this->compute_pde_jacobian(dl,d,du,u,param);
    // Get Right Hand Side
    std::vector<Real> r(this->nx_,0.0);
    std::vector<Real> diff(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      diff[i] = -(u[i]-this->evaluate_target((Real)(i+1)*this->dx_));
    }
    this->apply_mass(r,diff);
    // Solve solve adjoint system at current time step
    this->linear_solve(p,dl,d,du,r,true);
  }

  void solve_state_sensitivity(std::vector<Real> &v, const std::vector<Real> &u, 
                         const std::vector<Real> &z, const std::vector<Real> &param) {
    // Initialize State Storage
    v.clear();
    v.resize(this->nx_);
    // Get PDE Jacobian
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> dl(this->nx_-1,0.0);
    std::vector<Real> du(this->nx_-1,0.0);
    this->compute_pde_jacobian(dl,d,du,u,param);
    // Get Right Hand Side
    std::vector<Real> r(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      r[i] = this->dx_/6.0*(z[i]+4.0*z[i+1]+z[i+2]);
    }
    // Solve solve state sensitivity system at current time step
    this->linear_solve(v,dl,d,du,r);
  }

  void solve_adjoint_sensitivity(std::vector<Real> &q, const std::vector<Real> &u, 
                           const std::vector<Real> &p, const std::vector<Real> &v, 
                           const std::vector<Real> &z, const std::vector<Real> &param) {
    // Initialize State Storage
    q.clear();
    q.resize(this->nx_);
    // Get PDE Jacobian
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> dl(this->nx_-1,0.0);
    std::vector<Real> du(this->nx_-1,0.0);
    this->compute_pde_jacobian(dl,d,du,u,param);
    // Get Right Hand Side
    std::vector<Real> r(this->nx_,0.0);
    this->apply_mass(r,v);
    this->scale(r,-1.0);
    this->add_pde_hessian(r,u,p,v,-1.0);
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    this->linear_solve(q,dl,d,du,r,true);
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // SOLVE STATE EQUATION
    std::vector<Real> param(4,0.0);
    std::vector<Real> u;
    this->solve_state(u,*zp,param);
    // COMPUTE RESIDUAL
    Real val  = this->alpha_*0.5*this->dot(*zp,*zp);
    Real res  = 0.0, res1 = 0.0, res2 = 0.0, res3 = 0.0;
    for (int i=0; i<this->nx_; i++) {
      if ( i == 0 ) {
        res1 = u[i]-evaluate_target((Real)(i+1)*this->dx_);
        res2 = u[i+1]-evaluate_target((Real)(i+2)*this->dx_);
        res  = this->dx_/6.0*(4.0*res1 + res2)*res1;
      }
      else if ( i == this->nx_-1 ) {
        res1 = u[i-1]-evaluate_target((Real)i*this->dx_);
        res2 = u[i]-evaluate_target((Real)(i+1)*this->dx_);
        res  = this->dx_/6.0*(res1 + 4.0*res2)*res2;
      }
      else {
        res1 = u[i-1]-evaluate_target((Real)i*this->dx_);
        res2 = u[i]-evaluate_target((Real)(i+1)*this->dx_);
        res3 = u[i+1]-evaluate_target((Real)(i+2)*this->dx_);
        res  = this->dx_/6.0*(res1 + 4.0*res2 + res3)*res2;
      }
      val += 0.5*res;
    }
    return val;
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    Teuchos::RCP<std::vector<Real> > gp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());
    // SOLVE STATE EQUATION
    std::vector<Real> param(4,0.0);
    std::vector<Real> u;
    this->solve_state(u,*zp,param);
    // SOLVE ADJOINT EQUATION
    std::vector<Real> p;
    this->solve_adjoint(p,u,param);
    // COMPUTE GRADIENT
    for (int i=0; i<this->nx_+2; i++) {
      if (i==0) {
        (*gp)[i]  = this->alpha_*this->dx_/6.0*(2.0*(*zp)[i]+(*zp)[i+1]);
        (*gp)[i] -= this->dx_/6.0*(p[i]);
      }
      else if (i==this->nx_+1) {
        (*gp)[i]  = this->alpha_*this->dx_/6.0*(2.0*(*zp)[i]+(*zp)[i-1]);
        (*gp)[i] -= this->dx_/6.0*(p[i-2]);
      }
      else {
        (*gp)[i]  = this->alpha_*this->dx_/6.0*((*zp)[i-1]+4.0*(*zp)[i]+(*zp)[i+1]);
        if (i==1) {
          (*gp)[i] -= this->dx_/6.0*(4.0*p[i-1]+p[i]);
        }
        else if (i==this->nx_) {
          (*gp)[i] -= this->dx_/6.0*(4.0*p[i-1]+p[i-2]);
        }
        else {
          (*gp)[i] -= this->dx_/6.0*(p[i-2]+4.0*p[i-1]+p[i]);
        }
      }
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<std::vector<Real> > hvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
    // SOLVE STATE EQUATION
    std::vector<Real> param(4,0.0);
    std::vector<Real> u;
    this->solve_state(u,*zp,param);
    // SOLVE ADJOINT EQUATION
    std::vector<Real> p;
    this->solve_adjoint(p,u,param);
    // SOLVE STATE SENSITIVITY EQUATION
    std::vector<Real> s;
    this->solve_state_sensitivity(s,u,*vp,param);
    // SOLVE ADJOINT SENSITIVITY EQUATION
    std::vector<Real> q;
    this->solve_adjoint_sensitivity(q,u,p,s,*vp,param);
    // COMPUTE HESSVEC
    for (int i=0; i<this->nx_+2; i++) {
      if (i==0) {
        (*hvp)[i]  = this->alpha_*this->dx_/6.0*(2.0*(*vp)[i]+(*vp)[i+1]);
        (*hvp)[i] -= this->dx_/6.0*(q[i]);
      }
      else if (i==this->nx_+1) {
        (*hvp)[i]  = this->alpha_*this->dx_/6.0*(2.0*(*vp)[i]+(*vp)[i-1]);
        (*hvp)[i] -= this->dx_/6.0*(q[i-2]);
      }
      else {
        (*hvp)[i]  = this->alpha_*this->dx_/6.0*((*vp)[i-1]+4.0*(*vp)[i]+(*vp)[i+1]);
        if (i==1) {
          (*hvp)[i] -= this->dx_/6.0*(4.0*q[i-1]+q[i]);
        }
        else if (i==this->nx_) {
          (*hvp)[i] -= this->dx_/6.0*(4.0*q[i-1]+q[i-2]);
        }
        else {
          (*hvp)[i] -= this->dx_/6.0*(q[i-2]+4.0*q[i-1]+q[i]);
        }
      }
    }
  }
};

template<class Real>
class BoundConstraint_BurgersControl : public ROL::BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
public:
  BoundConstraint_BurgersControl(int dim) : dim_(dim), min_diff_(0.5) {
    x_lo_.resize(dim_,0.0);
    x_up_.resize(dim_,1.0);
  }
  bool isFeasible( const ROL::Vector<Real> &x ) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    bool val = true;
    int  cnt = 1;
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( (*ex)[i] >= this->x_lo_[i] && (*ex)[i] <= this->x_up_[i] ) { cnt *= 1; }
      else                                                            { cnt *= 0; }
    }
    if ( cnt == 0 ) { val = false; }
    return val;
  }
  void project( ROL::Vector<Real> &x ) {
    Teuchos::RCP<std::vector<Real> > ex =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(x)).getVector());
    for ( int i = 0; i < this->dim_; i++ ) {
      (*ex)[i] = std::max(this->x_lo_[i],std::min(this->x_up_[i],(*ex)[i]));
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ){
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void setVectorToUpperBound( ROL::Vector<Real> &u ) {
    Teuchos::RCP<std::vector<Real> > us = Teuchos::rcp( new std::vector<Real>(this->dim_,0.0) );
    us->assign(this->x_up_.begin(),this->x_up_.end()); 
    Teuchos::RCP<ROL::Vector<Real> > up = Teuchos::rcp( new ROL::StdVector<Real>(us) );
    u.set(*up);
  }
  void setVectorToLowerBound( ROL::Vector<Real> &l ) {
    Teuchos::RCP<std::vector<Real> > ls = Teuchos::rcp( new std::vector<Real>(this->dim_,0.0) );
    ls->assign(this->x_lo_.begin(),this->x_lo_.end()); 
    Teuchos::RCP<ROL::Vector<Real> > lp = Teuchos::rcp( new ROL::StdVector<Real>(ls) );
    l.set(*lp);
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
    // Initialize objective function.
    int nx      = 1028;  // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 0.0) );
    for (int i=0; i<nx+2; i++) {
      (*x_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*y_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);
    // Check deriatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);
    // Initialize Constraints
    BoundConstraint_BurgersControl<RealT> icon(nx+2);

    // Primal dual active set.
    Teuchos::ParameterList parlist;
    // Krylov parameters.
    parlist.set("Absolute Krylov Tolerance",              1.e-8);
    parlist.set("Relative Krylov Tolerance",              1.e-4);
    parlist.set("Maximum Number of Krylov Iterations",    50);
    // PDAS parameters.
    parlist.set("PDAS Relative Step Tolerance",           1.e-10);
    parlist.set("PDAS Relative Gradient Tolerance",       1.e-8);
    parlist.set("PDAS Maximum Number of Iterations",      10);
    parlist.set("PDAS Dual Scaling",                      (alpha>0.0) ? alpha : 1.e-4 );      
    // Define step.
    ROL::PrimalDualActiveSetStep<RealT> step_pdas(parlist);
    // Define status test.
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-16;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    ROL::StatusTest<RealT> status(gtol, stol, maxit);    
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo_pdas(step_pdas,status,false);
    // Run algorithm.
    x.zero();
    algo_pdas.run(x, obj, icon, true, *outStream);
    // Output control to file.
    std::ofstream file_pdas;
    file_pdas.open("control_PDAS.txt");
    for ( unsigned i = 0; i < (unsigned)nx+2; i++ ) {
      file_pdas << (*x_rcp)[i] << "\n";
    }
    file_pdas.close();

    // Projected Newtion.
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist_tr = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist_tr) );
    // Define step.
    ROL::TrustRegionStep<RealT> step_tr(*parlist_tr);
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo_tr(step_tr,status,false);
    // Run Algorithm
    y.zero();
    algo_tr.run(y,obj,icon,true,*outStream);
    // Output control to file.
    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( unsigned i = 0; i < (unsigned)nx+2; i++ ) {
      file_tr << (*y_rcp)[i] << "\n";
    }
    file_tr.close();
    // Output state to file.
    std::vector<RealT> u(nx,0.0);
    std::vector<RealT> param(4,0.0);
    obj.solve_state(u,*x_rcp,param);
    std::ofstream file;
    file.open("state.txt");
    for (unsigned i=0; i<(unsigned)nx; i++) {
      file << i/((RealT)(nx+1)) << "  " << u[i] << "\n";
    }
    file.close();
    // Compute error 
    Teuchos::RCP<ROL::Vector<RealT> > diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm();
    *outStream << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1e2*std::sqrt(ROL::ROL_EPSILON)) ? 1 : 0);
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

