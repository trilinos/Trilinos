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

/*! \file  example_02.cpp
    \brief Shows how to solve a linear-quadratic parabolic control problem 
           with bound constraints.
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
class Objective_ParabolicControl : public ROL::Objective<Real> {
private:
  std::vector<Real> u0_;
  Real alpha_;
  int nx_;
  int nt_;
  Real T_;
  Real dx_;
  Real dt_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
    Mu.resize(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      if ( i == 0 ) {
        Mu[i] = this->dx_/6.0*(2.0*u[i] + u[i+1]);
      }
      else if ( i == this->nx_-1 ) {
        Mu[i] = this->dx_/6.0*(u[i-1] + 2.0*u[i]);
      }
      else {
        Mu[i] = this->dx_/6.0*(u[i-1] + 4.0*u[i] + u[i+1]);
      }
    }
  }

  void compute_pde_jacobian(std::vector<Real> &d, std::vector<Real> &o, const std::vector<Real> &u) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    d.clear();
    d.resize(this->nx_,4.0*this->dx_/6.0 + this->dt_*2.0/this->dx_);
    d[0]           = this->dx_/3.0 + this->dt_/this->dx_;
    d[this->nx_-1] = this->dx_/3.0 + this->dt_/this->dx_ + this->dt_*4.0*std::pow(u[this->nx_-1],3.0);
    o.clear();
    o.resize(this->nx_-1,this->dx_/6.0 - this->dt_/this->dx_);
  }

  void compute_residual(std::vector<Real> &r, 
                  const std::vector<Real> &up, const std::vector<Real> &u, const Real z) {
    r.clear();
    r.resize(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      if ( i==0 ) {
        r[i] = this->dx_/6.0*(2.0*u[i]+u[i+1]) + this->dt_/this->dx_*(u[i]-u[i+1]);
        r[i]-= this->dx_/6.0*(2.0*up[i]+up[i+1]); // Contribution from previous state
      }
      else if ( i==this->nx_-1 ) {
        r[i] = this->dx_/6.0*(u[i-1]+2.0*u[i]) + this->dt_/this->dx_*(u[i]-u[i-1]);
        r[i]+= this->dt_*std::pow(u[i],4.0); // Stefan-Boltzmann boundary condition
        r[i]-= this->dx_/6.0*(2.0*up[i]+up[i-1]); // Contribution from previous state
        r[i]-= this->dt_*z; // Contribution from control
      }
      else {
        r[i] = this->dx_/6.0*(u[i-1]+4.0*u[i]+u[i+1]) + this->dt_/this->dx_*(2.0*u[i]-u[i-1]-u[i+1]);
        r[i]-= this->dx_/6.0*(up[i-1]+4.0*up[i]+up[i+1]); // Contribution from previous state
      }
    }
  }

  Real compute_norm(const std::vector<Real> &r) {
    Real norm = 0.0;
    for (unsigned i=0; i<r.size(); i++) {
      norm += r[i]*r[i];
    }
    return std::sqrt(norm);
  }

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    for (unsigned i=0; i<u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &d, std::vector<Real> &o, 
              const std::vector<Real> &r) {
    u.assign(r.begin(),r.end());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nx_;
    int nhrs = 1;
    lp.PTTRF(this->nx_,&d[0],&o[0],&info);
    lp.PTTRS(this->nx_,nhrs,&d[0],&o[0],&u[0],ldb,&info);
  }

  void run_newton(std::vector<Real> &u, const std::vector<Real> &up, const Real z) {
    // Set initial guess
    u.assign(up.begin(),up.end());
    // Compute residual and residual norm
    std::vector<Real> r(u.size(),0.0);
    this->compute_residual(r,up,u,z);
    Real rnorm = this->compute_norm(r);
    // Define tolerances
    Real tol   = 1.e2*ROL::ROL_EPSILON;
    Real maxit = 100;
    // Initialize Jacobian storage
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> o(this->nx_-1,0.0);
    // Iterate Newton's method
    Real alpha = 1.0, tmp = 0.0;
    std::vector<Real> s(this->nx_,0.0);
    std::vector<Real> utmp(this->nx_,0.0);
    for (int i=0; i<maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      this->compute_pde_jacobian(d,o,u);
      // Solve Newton system
      this->linear_solve(s,d,o,r);
      // Perform line search
      tmp = rnorm;
      alpha = 1.0;
      utmp.assign(u.begin(),u.end());
      this->update(utmp,s,-alpha);
      this->compute_residual(r,up,utmp,z);
      rnorm = this->compute_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON) ) {
        alpha /= 2.0;
        utmp.assign(u.begin(),u.end());
        this->update(utmp,s,-alpha);
        this->compute_residual(r,up,utmp,z);
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

  Objective_ParabolicControl(Real alpha = 1.e-4, int nx = 128, int nt = 100, Real T = 1) 
    : alpha_(alpha), nx_(nx), nt_(nt), T_(T) {
    u0_.resize(this->nx_,0.0);
    dx_ = 1.0/((Real)nx-1.0);
    dt_ = T/((Real)nt-1.0);
  }

  void solve_state(std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Initialize State Storage
    U.clear();
    U.resize(this->nt_+1);
    (U[0]).assign(this->u0_.begin(),this->u0_.end());
    std::vector<Real> up(this->u0_);
    std::vector<Real> u(this->u0_);
    // Time Step Using Implicit Euler
    for ( int t = 0; t < this->nt_; t++ ) {
      this->run_newton(u,up,z[t]);
      (U[t+1]).assign(u.begin(),u.end());
      up.assign(u.begin(),u.end());
    }
  }

  void solve_adjoint(std::vector<std::vector<Real> > &P, const std::vector<std::vector<Real> > &U) {
    // Initialize State Storage
    P.clear();
    P.resize(this->nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> p(this->nx_,0.0);
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> r(this->nx_,0.0);
    std::vector<Real> o(this->nx_-1,0.0);
    for ( int t = this->nt_; t > 0; t-- ) {
      // Get PDE Jacobian
      this->compute_pde_jacobian(d,o,U[t]);
      // Get Right Hand Side
      r.assign(this->nx_,0.0);
      if ( t==this->nt_ ) {
        std::vector<Real> diff(this->nx_,0.0);
        for (int i=0; i<this->nx_; i++) {
          diff[i] = -((U[t])[i]-this->evaluate_target((Real)i*this->dx_));
        }
        this->apply_mass(r,diff);
      }
      else {
        this->apply_mass(r,P[t]);
      } 
      // Solve solve adjoint system at current time step
      this->linear_solve(p,d,o,r);
      // Update State Storage
      (P[t-1]).assign(p.begin(),p.end());
    }
  }

  void solve_state_sensitivity(std::vector<std::vector<Real> > &V, 
                         const std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Initialize State Storage
    V.clear();
    V.resize(this->nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> v(this->nx_,0.0);
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> r(this->nx_,0.0);
    std::vector<Real> o(this->nx_-1,0.0);
    for ( int t = 0; t < this->nt_; t++ ) {
      // Get PDE Jacobian
      this->compute_pde_jacobian(d,o,U[t+1]);
      // Get Right Hand Side
      if( t == 0 ) {
        r.assign(this->nx_,0.0);
        r[this->nx_-1] = this->dt_*z[t];
      }
      else {
        this->apply_mass(r,V[t-1]);
        r[this->nx_-1] += this->dt_*z[t];
      }
      // Solve solve adjoint system at current time step
      this->linear_solve(v,d,o,r);
      // Update State Storage
      (V[t]).assign(v.begin(),v.end());
    }
  }

  void solve_adjoint_sensitivity(std::vector<std::vector<Real> > &Q, 
                           const std::vector<std::vector<Real> > &U, const std::vector<std::vector<Real> > &P,
                           const std::vector<std::vector<Real> > &V, const std::vector<Real> &z) {
    // Initialize State Storage
    Q.clear();
    Q.resize(this->nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> q(this->nx_,0.0);
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> r(this->nx_,0.0);
    std::vector<Real> o(this->nx_-1,0.0);
    for ( int t = this->nt_; t > 0; t-- ) {
      // Get PDE Jacobian
      this->compute_pde_jacobian(d,o,U[t]);
      // Get Right Hand Side
      if ( t == this->nt_ ) {
        std::vector<Real> tmp(this->nx_,0.0);
        r.assign(this->nx_,0.0);
        this->apply_mass(tmp,V[t-1]);
        this->update(r,tmp,-1.0);
      }
      else {
        this->apply_mass(r,Q[t]);
      }
      r[this->nx_-1]-=this->dt_*12.0*std::pow(U[t][this->nx_-1],2.0)*P[t-1][this->nx_-1]*V[t-1][this->nx_-1];
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      this->linear_solve(q,d,o,r);
      // Update State Storage
      (Q[t-1]).assign(q.begin(),q.end());
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

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // SOLVE STATE EQUATION
    std::vector<std::vector<Real> > U;
    this->solve_state(U,*zp);
    // COMPUTE RESIDUAL
    Real val  = 0.0;
    Real res  = 0.0;
    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (int t=0; t<this->nt_; t++) {
      val += (*zp)[t]*(*zp)[t];
    }
    val *= 0.5*this->alpha_*this->dt_;

    for (int i=0; i<this->nx_; i++) {
      if ( i == 0 ) {
        res1 = (U[this->nt_])[i]-evaluate_target((Real)i*this->dx_);
        res2 = (U[this->nt_])[i+1]-evaluate_target((Real)(i+1)*this->dx_);
        res  = this->dx_/6.0*(2.0*res1 + res2)*res1;
      }
      else if ( i == this->nx_-1 ) {
        res1 = (U[this->nt_])[i-1]-evaluate_target((Real)(i-1)*this->dx_);
        res2 = (U[this->nt_])[i]-evaluate_target((Real)i*this->dx_);
        res  = this->dx_/6.0*(res1 + 2.0*res2)*res2;
      }
      else {
        res1 = (U[this->nt_])[i-1]-evaluate_target((Real)(i-1)*this->dx_);
        res2 = (U[this->nt_])[i]-evaluate_target((Real)i*this->dx_);
        res3 = (U[this->nt_])[i+1]-evaluate_target((Real)(i+1)*this->dx_);
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
    std::vector<std::vector<Real> > U;
    this->solve_state(U,*zp);
    // SOLVE ADJOINT EQUATION
    std::vector<std::vector<Real> > P;
    this->solve_adjoint(P,U);
    // COMPUTE GRADIENT
    for (int t=0; t<this->nt_; t++) {
      (*gp)[t] = this->dt_*(this->alpha_*(*zp)[t] - (P[t])[this->nx_-1]);
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
    std::vector<std::vector<Real> > U;
    this->solve_state(U,*zp);
    // SOLVE ADJOINT EQUATION
    std::vector<std::vector<Real> > P;
    this->solve_adjoint(P,U);
    // SOLVE STATE SENSITIVITY EQUATION
    std::vector<std::vector<Real> > V;
    this->solve_state_sensitivity(V,U,*vp);
    // SOLVE ADJOINT SENSITIVITY EQUATION
    std::vector<std::vector<Real> > Q;
    this->solve_adjoint_sensitivity(Q,U,P,V,*vp);
    // COMPUTE HESSVEC
    for (int t=0; t<this->nt_; t++) {
      (*hvp)[t] = this->dt_*(this->alpha_*(*vp)[t] - (Q[t])[this->nx_-1]);
    }
  }
};

template<class Real>
class BoundConstraint_ParabolicControl : public ROL::BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
public:
  BoundConstraint_ParabolicControl(int dim) : dim_(dim), min_diff_(0.5) {
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
    int nx      = 100;   // Set spatial discretization.
    int nt      = 300;   // Set temporal discretization.
    RealT T     = 1.0;   // Set end time.
    RealT alpha = 1.e-3; // Set penalty parameter.
    Objective_ParabolicControl<RealT> obj(alpha,nx,nt,T);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (nt, 1.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (nt, 0.0) );
    for (int i=0; i<nt; i++) {
      (*x_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*y_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);
    // Check deriatives.
    obj.checkGradient(x,y,true);
    obj.checkHessVec(x,y,true);
    // Initialize Constraints
    BoundConstraint_ParabolicControl<RealT> icon(nt);

    // Primal dual active set.
    Teuchos::ParameterList parlist;
    // Krylov parameters.
    parlist.set("Absolute Krylov Tolerance",              1.e-8);
    parlist.set("Relative Krylov Tolerance",              1.e-4);
    parlist.set("Maximum Number of Krylov Iterations",    50);
    // PDAS parameters.
    parlist.set("PDAS Relative Step Tolerance",           1.e-8);
    parlist.set("PDAS Relative Gradient Tolerance",       1.e-6);
    parlist.set("PDAS Maximum Number of Iterations",      10);
    parlist.set("PDAS Dual Scaling",                      (alpha>0.0) ? alpha : 1.e-4 );      
    // Define step.
    ROL::PrimalDualActiveSetStep<RealT> step_pdas(parlist);
    // Define status test.
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-14;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    ROL::StatusTest<RealT> status(gtol, stol, maxit);    
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo_pdas(step_pdas,status,false);
    // Run algorithm.
    x.zero();
    algo_pdas.run(x, obj, icon, true);
    // Output control to file.
    std::ofstream file;
    file.open("control_PDAS.txt");
    for ( unsigned i = 0; i < (unsigned)nt; i++ ) {
      file << (*x_rcp)[i] << "\n";
    }
    file.close();

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
    algo_tr.run(y,obj,icon,true);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( unsigned i = 0; i < (unsigned)nt; i++ ) {
      file_tr << (*y_rcp)[i] << "\n";
    }
    file_tr.close();
   
    Teuchos::RCP<ROL::Vector<RealT> > diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm()/std::sqrt((RealT)nt-1.0);
    std::cout << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > std::sqrt(ROL::ROL_EPSILON)) ? 1 : 0);
#if 1
    // Output state to file.
    std::vector<std::vector<RealT> > U(nt);
    obj.solve_state(U,*y_rcp);
    std::ofstream file1;
    file1.open("state_tx.txt");
    for (unsigned t=0; t<(unsigned)nt; t++) {
      file1 << t*(T/((RealT)nt-1.0)) << "  ";
      for (unsigned i=0; i<(unsigned)nx; i++) {
        file1 << (U[t])[i] << "  ";
      }
      file1 << "\n";
    }
    file1.close();
    std::ofstream file2;
    file2.open("state_xt.txt");
    for (unsigned i=0; i<(unsigned)nx; i++) {
      file2 << i*(1.0/((RealT)nx-1.0)) << "  ";
      for (unsigned t=0; t<(unsigned)nt; t++) {
        file2 << (U[t])[i] << "  ";
      }
      file2 << "\n";
    }
    file2.close();
#endif
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

