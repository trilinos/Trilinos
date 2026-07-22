// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.hpp
    \brief Provides definitions of equality constraint and objective for
           example_02.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"

template<class Real>
class Constraint_BurgersControl : public ROL::Constraint_SimOpt<Real> {
private:
  int nx_;
  Real dx_;
  Real nu_;
  Real u0_;
  Real u1_;
  Real f_;

private:
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

  using ROL::Constraint_SimOpt<Real>::update;

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

  void compute_residual(std::vector<Real> &r, const std::vector<Real> &u, 
                  const std::vector<Real> &z) {
    r.clear();
    r.resize(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      // Contribution from stiffness term
      if ( i==0 ) {
        r[i] = this->nu_/this->dx_*(2.0*u[i]-u[i+1]);
      }
      else if (i==this->nx_-1) {
        r[i] = this->nu_/this->dx_*(2.0*u[i]-u[i-1]);
      }
      else {
        r[i] = this->nu_/this->dx_*(2.0*u[i]-u[i-1]-u[i+1]);
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
      r[i] -= this->dx_*this->f_;
    }
    // Contribution from Dirichlet boundary terms
    r[0]           -= this->u0_*u[          0]/6.0 + this->u0_*this->u0_/6.0 + this->nu_*this->u0_/this->dx_;
    r[this->nx_-1] += this->u1_*u[this->nx_-1]/6.0 + this->u1_*this->u1_/6.0 - this->nu_*this->u1_/this->dx_;
  }

  void compute_pde_jacobian(std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
                      const std::vector<Real> &u) {
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(this->nx_,this->nu_*2.0/this->dx_);
    dl.clear();
    dl.resize(this->nx_-1,-this->nu_/this->dx_);
    du.clear();
    du.resize(this->nx_-1,-this->nu_/this->dx_);
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
    d[0]           -= this->u0_/6.0;
    d[this->nx_-1] += this->u1_/6.0;
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

public:

  Constraint_BurgersControl(int nx = 128, Real nu = 1.e-2, Real u0 = 1.0, Real u1 = 0.0, Real f = 0.0) 
    : nx_(nx), nu_(nu), u0_(u0), u1_(u1), f_(f) {
    dx_ = 1.0/((Real)nx+1.0);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp =
      dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    this->compute_residual(*cp,*up,*zp);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Fill jvp
    for (int i = 0; i < this->nx_; i++) {
      (*jvp)[i] = this->nu_/this->dx_*2.0*(*vp)[i];
      if ( i > 0 ) {
        (*jvp)[i] += -this->nu_/this->dx_*(*vp)[i-1]
                     -(*up)[i-1]/6.0*(*vp)[i] 
                     -((*up)[i]+2.0*(*up)[i-1])/6.0*(*vp)[i-1];
      }
      if ( i < this->nx_-1 ) {
        (*jvp)[i] += -this->nu_/this->dx_*(*vp)[i+1]
                     +(*up)[i+1]/6.0*(*vp)[i] 
                     +((*up)[i]+2.0*(*up)[i+1])/6.0*(*vp)[i+1];
      }
    }
    (*jvp)[0]           -= this->u0_/6.0*(*vp)[0];
    (*jvp)[this->nx_-1] += this->u1_/6.0*(*vp)[this->nx_-1];
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0; i<this->nx_; i++) {
      // Contribution from control
      (*jvp)[i] = -this->dx_/6.0*((*vp)[i]+4.0*(*vp)[i+1]+(*vp)[i+2]);
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
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Get PDE Jacobian
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> dl(this->nx_-1,0.0);
    std::vector<Real> du(this->nx_-1,0.0);
    this->compute_pde_jacobian(dl,d,du,*up);
    // Solve solve state sensitivity system at current time step
    this->linear_solve(*ijvp,dl,d,du,*vp);
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Fill jvp
    for (int i = 0; i < this->nx_; i++) {
      (*jvp)[i] = this->nu_/this->dx_*2.0*(*vp)[i];
      if ( i > 0 ) {
        (*jvp)[i] += -this->nu_/this->dx_*(*vp)[i-1] 
                     -(*up)[i-1]/6.0*(*vp)[i] 
                     +((*up)[i-1]+2.0*(*up)[i])/6.0*(*vp)[i-1];
      }
      if ( i < this->nx_-1 ) {
        (*jvp)[i] += -this->nu_/this->dx_*(*vp)[i+1] 
                     +(*up)[i+1]/6.0*(*vp)[i]
                     -((*up)[i+1]+2.0*(*up)[i])/6.0*(*vp)[i+1];
      }
    }
    (*jvp)[0]           -= this->u0_/6.0*(*vp)[0];
    (*jvp)[this->nx_-1] += this->u1_/6.0*(*vp)[this->nx_-1];
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0; i<this->nx_+2; i++) {
      if ( i == 0 ) {
        (*jvp)[i] = -this->dx_/6.0*(*vp)[i];
      }
      else if ( i == 1 ) {
        (*jvp)[i] = -this->dx_/6.0*(4.0*(*vp)[i-1]+(*vp)[i]);
      }
      else if ( i == this->nx_ ) {
        (*jvp)[i] = -this->dx_/6.0*(4.0*(*vp)[i-1]+(*vp)[i-2]);
      }
      else if ( i == this->nx_+1 ) {
        (*jvp)[i] = -this->dx_/6.0*(*vp)[i-2];
      }
      else {
        (*jvp)[i] = -this->dx_/6.0*((*vp)[i-2]+4.0*(*vp)[i-1]+(*vp)[i]);
      }
    }
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > iajvp =
      dynamic_cast<ROL::StdVector<Real>&>(iajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    // Get PDE Jacobian
    std::vector<Real> d(this->nx_,0.0);
    std::vector<Real> du(this->nx_-1,0.0);
    std::vector<Real> dl(this->nx_-1,0.0);
    this->compute_pde_jacobian(dl,d,du,*up);
    // Solve solve adjoint system at current time step
    this->linear_solve(*iajvp,dl,d,du,*vp,true);
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp =
      dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp =
      dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    for (int i=0; i<this->nx_; i++) {
      // Contribution from nonlinear term
      (*ahwvp)[i] = 0.0;
      if (i<this->nx_-1){
        (*ahwvp)[i] += ((*wp)[i]*(*vp)[i+1] - (*wp)[i+1]*(2.0*(*vp)[i]+(*vp)[i+1]))/6.0;
      }
      if (i>0) {
        (*ahwvp)[i] += ((*wp)[i-1]*((*vp)[i-1]+2.0*(*vp)[i]) - (*wp)[i]*(*vp)[i-1])/6.0;
      }
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

//  void solveAugmentedSystem(ROL::Vector<Real> &v1, ROL::Vector<Real> &v2, const ROL::Vector<Real> &b1,
//                            const ROL::Vector<Real> &b2, const ROL::Vector<Real> &x, Real &tol) {}
};

template<class Real>
class Objective_BurgersControl : public ROL::Objective_SimOpt<Real> {
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
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_BurgersControl(Real alpha = 1.e-4, int nx = 128) : alpha_(alpha), nx_(nx) {
    dx_ = 1.0/((Real)nx+1.0);
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE RESIDUAL
    Real res1 = 0.0, res2 = 0.0, res3 = 0.0;
    Real valu = 0.0, valz = this->dot(*zp,*zp);
    for (int i=0; i<this->nx_; i++) {
      if ( i == 0 ) {
        res1  = (*up)[i]-evaluate_target((Real)(i+1)*this->dx_);
        res2  = (*up)[i+1]-evaluate_target((Real)(i+2)*this->dx_);
        valu += this->dx_/6.0*(4.0*res1 + res2)*res1;
      }
      else if ( i == this->nx_-1 ) {
        res1  = (*up)[i-1]-evaluate_target((Real)i*this->dx_);
        res2  = (*up)[i]-evaluate_target((Real)(i+1)*this->dx_);
        valu += this->dx_/6.0*(res1 + 4.0*res2)*res2;
      }
      else {
        res1  = (*up)[i-1]-evaluate_target((Real)i*this->dx_);
        res2  = (*up)[i]-evaluate_target((Real)(i+1)*this->dx_);
        res3  = (*up)[i+1]-evaluate_target((Real)(i+2)*this->dx_);
        valu += this->dx_/6.0*(res1 + 4.0*res2 + res3)*res2;
      }
    }
    return 0.5*(valu + this->alpha_*valz);
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<std::vector<Real> > gup = ROL::constPtrCast<std::vector<Real> >(
      (dynamic_cast<const ROL::StdVector<Real>&>(g)).getVector());
    // Unwrap x
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE GRADIENT WRT U
    std::vector<Real> diff(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      diff[i] = ((*up)[i]-this->evaluate_target((Real)(i+1)*this->dx_));
    }
    this->apply_mass(*gup,diff);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    ROL::Ptr<std::vector<Real> > gzp = ROL::constPtrCast<std::vector<Real> >(
      (dynamic_cast<const ROL::StdVector<Real>&>(g)).getVector());
    // Unwrap x
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE GRADIENT WRT Z
    for (int i=0; i<this->nx_+2; i++) {
      if (i==0) {
        (*gzp)[i] = this->alpha_*this->dx_/6.0*(2.0*(*zp)[i]+(*zp)[i+1]);
      }
      else if (i==this->nx_+1) {
        (*gzp)[i] = this->alpha_*this->dx_/6.0*(2.0*(*zp)[i]+(*zp)[i-1]);
      }
      else {
        (*gzp)[i] = this->alpha_*this->dx_/6.0*((*zp)[i-1]+4.0*(*zp)[i]+(*zp)[i+1]);
      }
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvup = 
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    // Unwrap v
    ROL::Ptr<const std::vector<Real> > vup =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE GRADIENT WRT U
    this->apply_mass(*hvup,*vup);
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
    ROL::Ptr<std::vector<Real> > hvzp = 
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    // Unwrap v
    ROL::Ptr<const std::vector<Real> > vzp =
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE GRADIENT WRT Z
    for (int i=0; i<this->nx_+2; i++) {
      if (i==0) {
        (*hvzp)[i] = this->alpha_*this->dx_/6.0*(2.0*(*vzp)[i]+(*vzp)[i+1]);
      }
      else if (i==this->nx_+1) {
        (*hvzp)[i] = this->alpha_*this->dx_/6.0*(2.0*(*vzp)[i]+(*vzp)[i-1]);
      }
      else {
        (*hvzp)[i] = this->alpha_*this->dx_/6.0*((*vzp)[i-1]+4.0*(*vzp)[i]+(*vzp)[i+1]);
      }
    }
  }
};
