// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.hpp
    \brief Provides definitions of equality constraint and objective for
           example_01.
*/

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

  using ROL::Objective<Real>::update;

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
    Real tol   = 1.e2*ROL::ROL_EPSILON<Real>();
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
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON<Real>()) ) {
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
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
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
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    ROL::Ptr<std::vector<Real> > gp =
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
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
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<std::vector<Real> > hvp =
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
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
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
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
    ROL::Ptr<std::vector<Real> > ex =
      dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    for ( int i = 0; i < this->dim_; i++ ) {
      (*ex)[i] = std::max(this->x_lo_[i],std::min(this->x_up_[i],(*ex)[i]));
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps = Real(0)) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps = Real(0)) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<const std::vector<Real> > eg =
      dynamic_cast<const ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(xeps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > geps) ){
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) {
    ROL::Ptr<const std::vector<Real> > ex =
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    ROL::Ptr<const std::vector<Real> > eg =
      dynamic_cast<const ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<std::vector<Real> > ev =
      dynamic_cast<ROL::StdVector<Real>&>(v).getVector();
    Real epsn = std::min(xeps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < -geps) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void setVectorToUpperBound( ROL::Vector<Real> &u ) {
    ROL::Ptr<std::vector<Real> > us = ROL::makePtr<std::vector<Real>>(this->dim_,0.0);
    us->assign(this->x_up_.begin(),this->x_up_.end()); 
    ROL::Ptr<ROL::Vector<Real> > up = ROL::makePtr<ROL::StdVector<Real>>(us);
    u.set(*up);
  }
  void setVectorToLowerBound( ROL::Vector<Real> &l ) {
    ROL::Ptr<std::vector<Real> > ls = ROL::makePtr<std::vector<Real>>(this->dim_,0.0);
    ls->assign(this->x_lo_.begin(),this->x_lo_.end()); 
    ROL::Ptr<ROL::Vector<Real> > lp = ROL::makePtr<ROL::StdVector<Real>>(ls);
    l.set(*lp);
  }
};
