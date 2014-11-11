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

/*! \file  example_03.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStepSQP.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

template<class Real>
class EqualityConstraint_BurgersControl : public ROL::EqualityConstraint_SimOpt<Real> {
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

  EqualityConstraint_BurgersControl(int nx = 128, Real nu = 1.e-2, Real u0 = 1.0, Real u1 = 0.0, Real f = 0.0) 
    : nx_(nx), nu_(nu), u0_(u0), u1_(u1), f_(f) {
    dx_ = 1.0/((Real)nx+1.0);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > cp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(c)).getVector());
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    this->compute_residual(*cp,*up,*zp);
  }

  void solve(ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > up =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(u)).getVector());
    up->assign(up->size(),z.norm()/up->size());
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // Compute residual and residual norm
    std::vector<Real> r(up->size(),0.0);
    this->compute_residual(r,*up,*zp);
    Real rnorm = this->compute_norm(r);
    // Define tolerances
    Real rtol  = 1.e2*ROL::ROL_EPSILON;
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
      this->compute_pde_jacobian(dl,d,du,*up);
      // Solve Newton system
      this->linear_solve(s,dl,d,du,r);
      // Perform line search
      tmp = rnorm;
      alpha = 1.0;
      utmp.assign(up->begin(),up->end());
      this->update(utmp,s,-alpha);
      this->compute_residual(r,utmp,*zp);
      rnorm = this->compute_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON) ) {
        alpha /= 2.0;
        utmp.assign(up->begin(),up->end());
        this->update(utmp,s,-alpha);
        this->compute_residual(r,utmp,*zp);
        rnorm = this->compute_norm(r); 
      }
      // Update iterate
      up->assign(utmp.begin(),utmp.end());
      if ( rnorm < rtol ) {
        break;
      }
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
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    for (int i=0; i<this->nx_; i++) {
      // Contribution from control
      (*jvp)[i] = -this->dx_/6.0*((*vp)[i]+4.0*(*vp)[i+1]+(*vp)[i+2]);
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
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(ajv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > iajvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(iajv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
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
    Teuchos::RCP<std::vector<Real> > ahwvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(ahwv)).getVector());
    Teuchos::RCP<const std::vector<Real> > wp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(w))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > gup = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT U
    std::vector<Real> diff(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      diff[i] = ((*up)[i]-this->evaluate_target((Real)(i+1)*this->dx_));
    }
    this->apply_mass(*gup,diff);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap g
    Teuchos::RCP<std::vector<Real> > gzp = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
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
    Teuchos::RCP<std::vector<Real> > hvup = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(hv)).getVector());
    // Unwrap v
    Teuchos::RCP<const std::vector<Real> > vup =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
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
    Teuchos::RCP<std::vector<Real> > hvzp = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(hv)).getVector());
    // Unwrap v
    Teuchos::RCP<const std::vector<Real> > vzp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
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
    int nx      = 128;  // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nu = 1e-2; // Viscosity parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx);
    // Initialize equality constraints
    EqualityConstraint_BurgersControl<RealT> con(nx,nu);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    for (int i=0; i<nx+2; i++) {
      (*z_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> z(z_rcp);
    ROL::StdVector<RealT> gz(gz_rcp);
    ROL::StdVector<RealT> yz(yz_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > zp  = Teuchos::rcp(&z,false);
    Teuchos::RCP<ROL::Vector<RealT> > gzp = Teuchos::rcp(&z,false);
    Teuchos::RCP<ROL::Vector<RealT> > yzp = Teuchos::rcp(&yz,false);

    Teuchos::RCP<std::vector<RealT> > u_rcp  = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gu_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yu_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    for (int i=0; i<nx; i++) {
      (*u_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yu_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    ROL::StdVector<RealT> u(u_rcp);
    ROL::StdVector<RealT> gu(gu_rcp);
    ROL::StdVector<RealT> yu(yu_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > up  = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > gup = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > yup = Teuchos::rcp(&yu,false);

    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > l_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    ROL::StdVector<RealT> c(c_rcp);
    ROL::StdVector<RealT> l(l_rcp);

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);

    // Check derivatives.
    obj.checkGradient(x,y,true);
    obj.checkHessVec(x,y,true);
    con.checkApplyJacobian(x,y,c,true);
    con.checkApplyAdjointJacobian(x,yu,c,true);
    con.checkApplyAdjointHessian(x,yu,y,true);

    // Optimization 
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist_tr = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist_tr) );
    // Define status test.
    RealT gtol  = 1e-14;  // norm of gradient of Lagrangian tolerance
    RealT ctol  = 1e-14;  // norm of constraint tolerance
    RealT stol  = 1e-16;  // norm of step tolerance
    int   maxit = 100;   // maximum number of iterations
    ROL::StatusTestSQP<RealT> status(gtol, ctol, stol, maxit);    
    // Define step.
    ROL::CompositeStepSQP<RealT> step(*parlist_tr);
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);
    // Run Algorithm
    RealT zerotol = 0.0;
    z.scale(50.0);
    con.solve(u,z,zerotol);
    //u_rcp->assign(u_rcp->size(),z.norm()/u_rcp->size());
    //u_rcp->assign(u_rcp->size(),1.0);
    //u.zero();
    c.zero();
    l.zero();
    std::vector<std::string> output = algo.run(x, g, l, c, obj, con, true);
    //for ( unsigned i = 0; i < output.size(); i++ ) {
    //  std::cout << output[i];
    //}
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

