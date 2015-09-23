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
#include "ROL_StatusTestSQP.hpp"
#include "ROL_AugmentedLagrangianStep.hpp"
#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Vector_SimOpt.hpp"
#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

template<class Real>
class BurgersFEM {
private:
  int nx_;
  Real dx_;
  Real nu_;
  Real u0_;
  Real u1_;
  Real f_;

private:
  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) const {
    for (unsigned i=0; i<u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void axpy(std::vector<Real> &out, const Real a, const std::vector<Real> &x, const std::vector<Real> &y) const {
    for (unsigned i=0; i < x.size(); i++) {
      out[i] = a*x[i] + y[i];
    }
  }

  void scale(std::vector<Real> &u, const Real alpha=0.0) const {
    for (unsigned i=0; i<u.size(); i++) {
      u[i] *= alpha;
    }
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
              const std::vector<Real> &r, const bool transpose = false) const {
    if ( r.size() == 1 ) {
      u.resize(1,r[0]/d[0]);
    }
    else if ( r.size() == 2 ) {
      u.resize(2,0.0);
      Real det = d[0]*d[1] - dl[0]*du[0];
      u[0] = (d[1]*r[0] - du[0]*r[1])/det;
      u[1] = (d[0]*r[1] - dl[0]*r[0])/det;
    }
    else {
      u.assign(r.begin(),r.end());
      // Perform LDL factorization
      Teuchos::LAPACK<int,Real> lp;
      std::vector<Real> du2(r.size()-2,0.0);
      std::vector<int> ipiv(r.size(),0);
      int info;
      int dim  = r.size();
      int ldb  = r.size();
      int nhrs = 1;
      lp.GTTRF(dim,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&info);
      char trans = 'N';
      if ( transpose ) { 
        trans = 'T';
      }
      lp.GTTRS(trans,dim,nhrs,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&u[0],ldb,&info);
    }
  }

public:
  BurgersFEM(int nx = 128, Real nu = 1.e-2, Real u0 = 1.0, Real u1 = 0.0, Real f = 0.0) 
    : nx_(nx), dx_(1.0/((Real)nx+1.0)), nu_(nu), u0_(u0), u1_(u1), f_(f) {}

  int num_dof(void) const {
    return nx_;
  }

  Real mesh_spacing(void) const {
    return dx_;
  }

  /***************************************************************************/
  /*********************** L2 INNER PRODUCT FUNCTIONS ************************/
  /***************************************************************************/
  // Compute L2 inner product
  Real compute_L2_dot(const std::vector<Real> &x, const std::vector<Real> &y) const {
    Real ip = 0.0;
    Real c = (((int)x.size()==nx_) ? 4.0 : 2.0);
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

  // compute L2 norm
  Real compute_L2_norm(const std::vector<Real> &r) const {
    return std::sqrt(compute_L2_dot(r,r));
  }

  // Apply L2 Reisz operator
  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) const {
    Mu.resize(u.size(),0.0);
    Real c = (((int)u.size()==nx_) ? 4.0 : 2.0);
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

  // Apply L2 inverse Reisz operator
  void apply_inverse_mass(std::vector<Real> &Mu, const std::vector<Real> &u) const {
    unsigned nx = u.size();
    // Build mass matrix
    std::vector<Real> dl(nx-1,dx_/6.0);
    std::vector<Real> d(nx,2.0*dx_/3.0);
    std::vector<Real> du(nx-1,dx_/6.0);
    if ( (int)nx != nx_ ) {
      d[   0] = dx_/3.0;
      d[nx-1] = dx_/3.0;
    }
    linear_solve(Mu,dl,d,du,u);
  }

  void test_inverse_mass(std::ostream &outStream = std::cout) {
    // State Mass Matrix
    std::vector<Real> u(nx_,0.0), Mu(nx_,0.0), iMMu(nx_,0.0), diff(nx_,0.0);
    for (int i = 0; i < nx_; i++) {
      u[i] = 2.0*(Real)rand()/(Real)RAND_MAX - 1.0;
    }
    apply_mass(Mu,u);
    apply_inverse_mass(iMMu,Mu);
    axpy(diff,-1.0,iMMu,u);
    Real error = compute_L2_norm(diff);
    Real normu = compute_L2_norm(u);
    outStream << "Test Inverse State Mass Matrix\n";
    outStream << "  ||u - inv(M)Mu|| = " << error << "\n";
    outStream << "  ||u||            = " << normu << "\n";
    outStream << "  Relative Error   = " << error/normu << "\n";
    outStream << "\n";
    // Control Mass Matrix
    u.resize(nx_+2,0.0); Mu.resize(nx_+2,0.0); iMMu.resize(nx_+2,0.0); diff.resize(nx_+2,0.0);
    for (int i = 0; i < nx_+2; i++) {
      u[i] = 2.0*(Real)rand()/(Real)RAND_MAX - 1.0;
    }
    apply_mass(Mu,u);
    apply_inverse_mass(iMMu,Mu);
    axpy(diff,-1.0,iMMu,u);
    error = compute_L2_norm(diff);
    normu = compute_L2_norm(u);
    outStream << "Test Inverse Control Mass Matrix\n";
    outStream << "  ||z - inv(M)Mz|| = " << error << "\n";
    outStream << "  ||z||            = " << normu << "\n";
    outStream << "  Relative Error   = " << error/normu << "\n";
    outStream << "\n";
  }

  // Apply L2 inverse Reisz operator
  void apply_inverse_mass(std::vector<Real> &Mu, const std::vector<Real> &u,
                    const std::vector<unsigned> &index) const {
    unsigned nx = index.size();
    if ( nx != 0 ) {
      std::vector<Real> ui(nx,0.0);
      for (unsigned i = 0; i < nx; i++) {
        ui[i] = u[index[i]];
      }
      // Build mass matrix
      std::vector<Real> dl(nx-1,dx_/6.0);
      std::vector<Real> d(nx,2.0*dx_/3.0);
      std::vector<Real> du(nx-1,dx_/6.0);
      if ( (int)nx != nx_ ) {
        if ( index[0] == 0 ) {
          d[   0] = dx_/3.0;
        }
        if ( index[nx-1] == u.size()-1 ) {
          d[nx-1] = dx_/3.0;
        }
      }
      linear_solve(Mu,dl,d,du,ui);
    }
  }

  /***************************************************************************/
  /*********************** PDE RESIDUAL AND SOLVE ****************************/
  /***************************************************************************/
  // Compute current PDE residual
  void compute_residual(std::vector<Real> &r, const std::vector<Real> &u, 
                  const std::vector<Real> &z) const {
    r.clear();
    r.resize(nx_,0.0);
    for (int i=0; i<nx_; i++) {
      // Contribution from stiffness term
      if ( i==0 ) {
        r[i] = nu_/dx_*(2.0*u[i]-u[i+1]);
      }
      else if (i==nx_-1) {
        r[i] = nu_/dx_*(2.0*u[i]-u[i-1]);
      }
      else {
        r[i] = nu_/dx_*(2.0*u[i]-u[i-1]-u[i+1]);
      }
      // Contribution from nonlinear term
      if (i<nx_-1){
        r[i] += u[i+1]*(u[i]+u[i+1])/6.0;
      }
      if (i>0) {
        r[i] -= u[i-1]*(u[i-1]+u[i])/6.0;
      }
      // Contribution from control
      r[i] -= dx_/6.0*(z[i]+4.0*z[i+1]+z[i+2]);
      // Contribution from right-hand side
      r[i] -= dx_*f_;
    }
    // Contribution from Dirichlet boundary terms
    r[0]     -= u0_*u[    0]/6.0 + u0_*u0_/6.0 + nu_*u0_/dx_;
    r[nx_-1] += u1_*u[nx_-1]/6.0 + u1_*u1_/6.0 - nu_*u1_/dx_;
  }

  // Solve PDE
  void solve(std::vector<Real> &u, const std::vector<Real> &z) const {
    // Compute residual and residual norm
    std::vector<Real> r(u.size(),0.0);
    compute_residual(r,u,z);
    Real rnorm = compute_L2_norm(r);
    // Define tolerances
    Real rtol  = 1.e2*ROL::ROL_EPSILON;
    Real maxit = 500;
    // Initialize Jacobian storage
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    std::vector<Real> du(nx_-1,0.0);
    // Iterate Newton's method
    Real alpha = 1.0, tmp = 0.0;
    std::vector<Real> s(nx_,0.0);
    std::vector<Real> utmp(nx_,0.0);
    for (int i=0; i<maxit; i++) {
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
      compute_residual(r,utmp,z);
      rnorm = compute_L2_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON) ) {
        alpha /= 2.0;
        utmp.assign(u.begin(),u.end());
        update(utmp,s,-alpha);
        compute_residual(r,utmp,z);
        rnorm = compute_L2_norm(r); 
      }
      // Update iterate
      u.assign(utmp.begin(),utmp.end());
      if ( rnorm < rtol ) {
        break;
      }
    }
  }

  /***************************************************************************/
  /*********************** PDE JACOBIAN FUNCTIONS ****************************/
  /***************************************************************************/
  // Build PDE Jacobian trigiagonal matrix
  void compute_pde_jacobian(std::vector<Real> &dl, // Lower diagonal
                            std::vector<Real> &d,  // Diagonal
                            std::vector<Real> &du, // Upper diagonal
                      const std::vector<Real> &u) const { // State variable
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();
    d.resize(nx_,nu_*2.0/dx_);
    dl.clear();
    dl.resize(nx_-1,-nu_/dx_);
    du.clear();
    du.resize(nx_-1,-nu_/dx_);
    // Contribution from nonlinearity
    for (int i=0; i<nx_; i++) {
      if (i<nx_-1) {
        dl[i] += (-2.0*u[i]-u[i+1])/6.0;
        d[i]  += u[i+1]/6.0;
      }
      if (i>0) {
        d[i]    += -u[i-1]/6.0;
        du[i-1] += (u[i-1]+2.0*u[i])/6.0;
      }
    }
    // Contribution from Dirichlet boundary conditions
    d[    0] -= u0_/6.0;
    d[nx_-1] += u1_/6.0;
  }

  // Apply PDE Jacobian to a vector
  void apply_pde_jacobian(std::vector<Real> &jv,
                    const std::vector<Real> &v,
                    const std::vector<Real> &u,
                    const std::vector<Real> &z) const {
    // Fill jv
    for (int i = 0; i < nx_; i++) {
      jv[i] = nu_/dx_*2.0*v[i];
      if ( i > 0 ) {
        jv[i] += -nu_/dx_*v[i-1]-u[i-1]/6.0*v[i]-(u[i]+2.0*u[i-1])/6.0*v[i-1];
      }
      if ( i < nx_-1 ) {
        jv[i] += -nu_/dx_*v[i+1]+u[i+1]/6.0*v[i]+(u[i]+2.0*u[i+1])/6.0*v[i+1];
      }
    }
    jv[    0] -= u0_/6.0*v[0];
    jv[nx_-1] += u1_/6.0*v[nx_-1];
  }

  // Apply inverse PDE Jacobian to a vector
  void apply_inverse_pde_jacobian(std::vector<Real> &ijv,
                            const std::vector<Real> &v,
                            const std::vector<Real> &u,
                            const std::vector<Real> &z) const {
    // Get PDE Jacobian
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    std::vector<Real> du(nx_-1,0.0);
    compute_pde_jacobian(dl,d,du,u);
    // Solve solve state sensitivity system at current time step
    linear_solve(ijv,dl,d,du,v);
  }

  // Apply adjoint PDE Jacobian to a vector
  void apply_adjoint_pde_jacobian(std::vector<Real> &ajv,
                            const std::vector<Real> &v,
                            const std::vector<Real> &u,
                            const std::vector<Real> &z) const {
    // Fill jvp
    std::vector<Real> Mv(v.size(),0.0);
    std::vector<Real> tmp(v.size(),0.0);
    apply_mass(Mv,v);
    for (int i = 0; i < nx_; i++) {
      tmp[i] = nu_/dx_*2.0*Mv[i];
      if ( i > 0 ) {
        tmp[i] += -nu_/dx_*Mv[i-1]-u[i-1]/6.0*Mv[i]
                  +(u[i-1]+2.0*u[i])/6.0*Mv[i-1];
      }
      if ( i < nx_-1 ) {
        tmp[i] += -nu_/dx_*Mv[i+1]+u[i+1]/6.0*Mv[i]
                  -(u[i+1]+2.0*u[i])/6.0*Mv[i+1];
      }
    }
    tmp[    0] -= u0_/6.0*Mv[0];
    tmp[nx_-1] += u1_/6.0*Mv[nx_-1];
    apply_inverse_mass(ajv,tmp);
  }

  // Apply inverse adjoint PDE Jacobian to a vector
  void apply_inverse_adjoint_pde_jacobian(std::vector<Real> &iajv,
                                    const std::vector<Real> &v,
                                    const std::vector<Real> &u,
                                    const std::vector<Real> &z) const {
    std::vector<Real> Mv, tmp;
    apply_mass(Mv,v);
    // Get PDE Jacobian
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> du(nx_-1,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    compute_pde_jacobian(dl,d,du,u);
    // Solve solve adjoint system at current time step
    linear_solve(tmp,dl,d,du,Mv,true);
    apply_inverse_mass(iajv,tmp);
  }

  /***************************************************************************/
  /*********************** CONTROL JACOBIAN FUNCTIONS ************************/
  /***************************************************************************/
  // Apply control Jacobian to a vector
  void apply_control_jacobian(std::vector<Real> &jv,
                        const std::vector<Real> &v,
                        const std::vector<Real> &u,
                        const std::vector<Real> &z) const {
    for (int i=0; i<nx_; i++) {
      // Contribution from control
      jv[i] = -dx_/6.0*(v[i]+4.0*v[i+1]+v[i+2]);
    }
  }

  // Apply adjoint control Jacobian to a vector
  void apply_adjoint_control_jacobian(std::vector<Real> &jv,
                                const std::vector<Real> &v,
                                const std::vector<Real> &u,
                                const std::vector<Real> &z) const {
    std::vector<Real> Mv(nx_,0.0);
    std::vector<Real> tmp(nx_+2,0.0);
    apply_mass(Mv,v);
    for (int i=0; i<nx_+2; i++) {
      if ( i == 0 ) {
        tmp[i] = -dx_/6.0*Mv[i];
      }
      else if ( i == 1 ) {
        tmp[i] = -dx_/6.0*(4.0*Mv[i-1]+Mv[i]);
      }
      else if ( i == nx_ ) {
        tmp[i] = -dx_/6.0*(4.0*Mv[i-1]+Mv[i-2]);
      }
      else if ( i == nx_+1 ) {
        tmp[i] = -dx_/6.0*Mv[i-2];
      }
      else {
        tmp[i] = -dx_/6.0*(Mv[i-2]+4.0*Mv[i-1]+Mv[i]);
      }
    }
    apply_inverse_mass(jv,tmp);
  }

  /***************************************************************************/
  /*********************** AJDOINT HESSIANS **********************************/
  /***************************************************************************/
  void apply_adjoint_pde_hessian(std::vector<Real> &ahwv,
                           const std::vector<Real> &w,
                           const std::vector<Real> &v,
                           const std::vector<Real> &u,
                           const std::vector<Real> &z) const {
    std::vector<Real> Mw(nx_,0.0), tmp(nx_,0.0);
    apply_mass(Mw,w);
    for (int i=0; i<nx_; i++) {
      // Contribution from nonlinear term
      tmp[i] = 0.0;
      if (i<nx_-1){
        tmp[i] += (Mw[i]*v[i+1] - Mw[i+1]*(2.0*v[i]+v[i+1]))/6.0;
      }
      if (i>0) {
        tmp[i] += (Mw[i-1]*(v[i-1]+2.0*v[i]) - Mw[i]*v[i-1])/6.0;
      }
    }
    apply_inverse_mass(ahwv,tmp);
  }
  void apply_adjoint_pde_control_hessian(std::vector<Real> &ahwv,
                                   const std::vector<Real> &w,
                                   const std::vector<Real> &v,
                                   const std::vector<Real> &u,
                                   const std::vector<Real> &z) {
    ahwv.assign(u.size(),0.0);
  }
  void apply_adjoint_control_pde_hessian(std::vector<Real> &ahwv,
                                   const std::vector<Real> &w,
                                   const std::vector<Real> &v,
                                   const std::vector<Real> &u,
                                   const std::vector<Real> &z) {
    ahwv.assign(z.size(),0.0);
  }
  void apply_adjoint_control_hessian(std::vector<Real> &ahwv,
                               const std::vector<Real> &w,
                               const std::vector<Real> &v,
                               const std::vector<Real> &u,
                               const std::vector<Real> &z) {
    ahwv.assign(z.size(),0.0);
  }
};

template<class Real>
class L2Vector : public ROL::Vector<Real> {
private:
  Teuchos::RCP<std::vector<Real> > vec_;
  Teuchos::RCP<BurgersFEM<Real> > fem_;

public:
  L2Vector(const Teuchos::RCP<std::vector<Real> > & vec,
           const Teuchos::RCP<BurgersFEM<Real> > &fem) : vec_(vec), fem_(fem) {}

  void set( const ROL::Vector<Real> &x ) {
    const L2Vector &ex = Teuchos::dyn_cast<const L2Vector>(x);
    const std::vector<Real>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),vec_->begin());
  }

  void plus( const ROL::Vector<Real> &x ) {
    const L2Vector &ex = Teuchos::dyn_cast<const L2Vector>(x);
    const std::vector<Real>& xval = *ex.getVector();
    unsigned dimension  = vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*vec_)[i] += xval[i];
    }
  }

  void scale( const Real alpha ) {
    unsigned dimension = vec_->size();
    for (unsigned i=0; i<dimension; i++) {
      (*vec_)[i] *= alpha;
    }
  }

  Real dot( const ROL::Vector<Real> &x ) const {
    const L2Vector & ex = Teuchos::dyn_cast<const L2Vector>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return fem_->compute_L2_dot(xval,*vec_);
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  Teuchos::RCP<ROL::Vector<Real> > clone() const {
    return Teuchos::rcp( new L2Vector( Teuchos::rcp(new std::vector<Real>(vec_->size(),0.0)),fem_));
  }

  Teuchos::RCP<const std::vector<Real> > getVector() const {
    return vec_;
  }

  Teuchos::RCP<std::vector<Real> > getVector() {
    return vec_;
  }

  Teuchos::RCP<ROL::Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<L2Vector> e
      = Teuchos::rcp( new L2Vector( Teuchos::rcp(new std::vector<Real>(vec_->size(),0.0)),fem_));
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return vec_->size();
  }
};

template<class Real>
class EqualityConstraint_BurgersControl : public ROL::EqualityConstraint_SimOpt<Real> {
private:
  Teuchos::RCP<BurgersFEM<Real> > fem_;

public:
  EqualityConstraint_BurgersControl(Teuchos::RCP<BurgersFEM<Real> > &fem) : fem_(fem) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > cp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(c)).getVector());
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->compute_residual(*cp,*up,*zp);
  }

  void solve(ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > up =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(u)).getVector());
    up->assign(up->size(),z.norm()/up->size());
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->solve(*up,*zp);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_pde_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_control_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > ijvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(ijv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_inverse_pde_jacobian(*ijvp,*vp,*up,*zp);
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(ajv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_adjoint_pde_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > jvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_adjoint_control_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > iajvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(iajv)).getVector());
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_inverse_adjoint_pde_jacobian(*iajvp,*vp,*up,*zp);
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > ahwvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(ahwv)).getVector());
    Teuchos::RCP<const std::vector<Real> > wp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(w))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_adjoint_pde_hessian(*ahwvp,*wp,*vp,*up,*zp);
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > ahwvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(ahwv)).getVector());
    Teuchos::RCP<const std::vector<Real> > wp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(w))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_adjoint_pde_control_hessian(*ahwvp,*wp,*vp,*up,*zp);
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > ahwvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(ahwv)).getVector());
    Teuchos::RCP<const std::vector<Real> > wp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(w))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_adjoint_control_pde_hessian(*ahwvp,*wp,*vp,*up,*zp);
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > ahwvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(ahwv)).getVector());
    Teuchos::RCP<const std::vector<Real> > wp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(w))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    fem_->apply_adjoint_control_hessian(*ahwvp,*wp,*vp,*up,*zp);
  }
};

template<class Real>
class Objective_BurgersControl : public ROL::Objective_SimOpt<Real> {
private:
  Real alpha_; // Penalty Parameter
  Teuchos::RCP<BurgersFEM<Real> > fem_;
  Teuchos::RCP<ROL::Vector<Real> > ud_;
  Teuchos::RCP<ROL::Vector<Real> > diff_;

public:
  Objective_BurgersControl(const Teuchos::RCP<BurgersFEM<Real> > &fem, 
                           const Teuchos::RCP<ROL::Vector<Real> > &ud,
                           Real alpha = 1.e-4) : alpha_(alpha), fem_(fem), ud_(ud) {
    diff_ = ud_->clone();
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    diff_->set(u);
    diff_->axpy(-1.0,*ud_);
    return 0.5*(diff_->dot(*diff_) + alpha_*z.dot(z)); 
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    g.set(u);
    g.axpy(-1.0,*ud_);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    g.set(z);
    g.scale(alpha_);
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.set(v);
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
    hv.set(v);
    hv.scale(alpha_);
  }
};

template<class Real>
class L2BoundConstraint : public ROL::BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
  Real scale_;
  Teuchos::RCP<BurgersFEM<Real> > fem_;

  void axpy(std::vector<Real> &out, const Real a,
      const std::vector<Real> &x, const std::vector<Real> &y) const{
    out.resize(dim_,0.0);
    for (unsigned i = 0; i < dim_; i++) {
      out[i] = a*x[i] + y[i];
    }
  }

  void projection(std::vector<Real> &x) {
    for ( int i = 0; i < dim_; i++ ) {
      x[i] = std::max(x_lo_[i],std::min(x_up_[i],x[i]));
    }
  }

public:
  L2BoundConstraint(std::vector<Real> &l, std::vector<Real> &u,
              const Teuchos::RCP<BurgersFEM<Real> > &fem, Real scale = 1.0)
    : x_lo_(l), x_up_(u), scale_(scale), fem_(fem) {
    dim_ = x_lo_.size();
    for ( int i = 0; i < dim_; i++ ) {
      if ( i == 0 ) {
        min_diff_ = x_up_[i] - x_lo_[i];
      }
      else {
        min_diff_ = ( (min_diff_ < (x_up_[i] - x_lo_[i])) ? min_diff_ : (x_up_[i] - x_lo_[i]) );
      }
    }
    min_diff_ *= 0.5;
  }

  bool isFeasible( const ROL::Vector<Real> &x ) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    bool val = true;
    int  cnt = 1;
    for ( int i = 0; i < dim_; i++ ) {
      if ( (*ex)[i] >= x_lo_[i] && (*ex)[i] <= x_up_[i] ) { cnt *= 1; }
      else                                                { cnt *= 0; }
    }
    if ( cnt == 0 ) { val = false; }
    return val;
  }

  void project( ROL::Vector<Real> &x ) {
    Teuchos::RCP<std::vector<Real> > ex =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(x)).getVector());
    projection(*ex);
  }

//  void project( ROL::Vector<Real> &x ) {
//    Teuchos::RCP<std::vector<Real> > ex =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(x)).getVector());
//    // Store current x value
//    std::vector<Real> xold(*ex);
//    // Compute initial guess
//    projection(*ex);
//    // Run semismooth Newton
//    Real tol = 1.e-10;
//    std::vector<Real> diff(dim_,0.0), Mdiff(dim_,0.0), PMdiff(dim_,0.0);
//    std::vector<Real> res(dim_,0.0), sA(dim_,0.0), MsA(dim_,0.0), rhs(dim_,0.0);
//    std::vector<Real> sI;
//    std::vector<unsigned> inactiveSet;
//    Real rnorm = tol+1.0;
//    unsigned iter = 0;
//    while (rnorm > tol) {
//      // Compute Residual
//      axpy(diff,-1.0,xold,*ex);
//      fem_->apply_mass(Mdiff,diff);
//      axpy(PMdiff,-1.0,Mdiff,*ex);
//      projection(PMdiff);
//      axpy(res,-1.0,PMdiff,*ex);
//      rnorm = fem_->compute_L2_norm(res);
//std::cout << "iter = " << iter << "  rnorm = " << rnorm << "\n";
//      // Solve semismooth Newton system
//      inactiveSet.clear(); sI.clear();
//      sA.assign(dim_,0.0);
//      for (unsigned i = 0; i < dim_; i++) {
//        if (PMdiff[i] == x_up_[i]) {
//          sA[i] = x_up_[i]-(*ex)[i];
//
//          (*ex)[i] = x_up_[i];
//        }
//        else if (PMdiff[i] == x_lo_[i]) {
//          sA[i] = x_lo_[i]-(*ex)[i];
//
//          (*ex)[i] = x_lo_[i];
//        }
//        else {
//          inactiveSet.push_back(i);
//        }
//      }
//      fem_->apply_mass(MsA,sA);
//      axpy(rhs,-1.0,MsA,sA);
//      axpy(rhs,-1.0,res,rhs);
//      fem_->apply_inverse_mass(sI,rhs,inactiveSet);
//      // Update iterate
//      for (unsigned i = 0; i < inactiveSet.size(); i++) {
//        (*ex)[inactiveSet[i]] += sI[i];
//      }
//    }
//  }

  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(v)).getVector());
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; i++ ) {
      if ( ((*ex)[i] <= x_lo_[i]+epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(v)).getVector());
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; i++ ) {
      if ( ((*ex)[i] >= x_up_[i]-epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(v)).getVector());
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; i++ ) {
      if ( ((*ex)[i] <= x_lo_[i]+epsn) ||
           ((*ex)[i] >= x_up_[i]-epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(v)).getVector());
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; i++ ) {
      if ( ((*ex)[i] <= x_lo_[i]+epsn && (*eg)[i] > 0.0) ) {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(v)).getVector());
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; i++ ) {
      if ( ((*ex)[i] >= x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<L2Vector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<L2Vector<Real> >(v)).getVector());
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; i++ ) {
      if ( ((*ex)[i] <= x_lo_[i]+epsn && (*eg)[i] > 0.0) ||
           ((*ex)[i] >= x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
        (*ev)[i] = 0.0;
      }
    }
  }

  void setVectorToUpperBound( ROL::Vector<Real> &u ) {
    Teuchos::RCP<std::vector<Real> > us = Teuchos::rcp( new std::vector<Real>(dim_,0.0) );
    us->assign(x_up_.begin(),x_up_.end());
    Teuchos::RCP<ROL::Vector<Real> > up = Teuchos::rcp( new L2Vector<Real>(us,fem_) );
    u.set(*up);
  }

  void setVectorToLowerBound( ROL::Vector<Real> &l ) {
    Teuchos::RCP<std::vector<Real> > ls = Teuchos::rcp( new std::vector<Real>(dim_,0.0) );
    ls->assign(x_lo_.begin(),x_lo_.end());
    Teuchos::RCP<ROL::Vector<Real> > lp = Teuchos::rcp( new L2Vector<Real>(ls,fem_) );
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
    // Initialize finite element class.
    int nx      = 128;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nu    = 1e-2;  // Viscosity parameter.
    Teuchos::RCP<BurgersFEM<RealT> > fem
      = Teuchos::rcp(new BurgersFEM<RealT>(nx,nu));
    fem->test_inverse_mass(*outStream);
    // Initialize objective function.
    Teuchos::RCP<std::vector<RealT> > ud_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<ROL::Vector<RealT> > ud = Teuchos::rcp(new L2Vector<RealT>(ud_rcp,fem));
    Objective_BurgersControl<RealT> obj(fem,ud,alpha);
    // Initialize equality constraints
    EqualityConstraint_BurgersControl<RealT> con(fem);
    // Initialize equality constraints
    std::vector<RealT> lo1(nx, 0.0), hi1(nx, 1.0);
    //std::vector<RealT> lo1(nx, -1.e6), hi1(nx, 1.e6);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd1 = Teuchos::rcp(new L2BoundConstraint<RealT>(lo1,hi1,fem));
    //bnd1.deactivate();
    //std::vector<RealT> lo2(nx+2, -0.1*ROL::ROL_OVERFLOW), hi2(nx+2, 0.1*ROL::ROL_OVERFLOW);
    std::vector<RealT> lo2(nx+2,0.0), hi2(nx+2,2.0);
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd2 = Teuchos::rcp(new L2BoundConstraint<RealT>(lo2,hi2,fem));
    //bnd2.deactivate();
    ROL::BoundConstraint_SimOpt<RealT> bnd(bnd1,bnd2);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    for (int i=0; i<nx+2; i++) {
      (*z_rcp)[i]  = (RealT)rand()/(RealT)RAND_MAX;
      (*yz_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    L2Vector<RealT> z(z_rcp,fem);
    L2Vector<RealT> gz(gz_rcp,fem);
    L2Vector<RealT> yz(yz_rcp,fem);
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
    L2Vector<RealT> u(u_rcp,fem);
    L2Vector<RealT> gu(gu_rcp,fem);
    L2Vector<RealT> yu(yu_rcp,fem);
    Teuchos::RCP<ROL::Vector<RealT> > up  = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > gup = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > yup = Teuchos::rcp(&yu,false);

    Teuchos::RCP<std::vector<RealT> > c_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > l_rcp = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    L2Vector<RealT> c(c_rcp,fem);
    L2Vector<RealT> l(l_rcp,fem);
    for (int i=0; i<nx; i++) {
      (*l_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);

    // Check derivatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);

    con.checkApplyJacobian(x,y,c,true,*outStream);
    con.checkApplyAdjointHessian(x,yu,y,x,true,*outStream);

    con.checkSolve(u,z,c,true,*outStream);
    con.checkAdjointConsistencyJacobian_1(c,yu,u,z,true,*outStream);
    con.checkAdjointConsistencyJacobian_2(c,yz,u,z,true,*outStream);
    con.checkInverseJacobian_1(yu,l,u,z,true,*outStream);
    con.checkInverseAdjointJacobian_1(yu,l,u,z,true,*outStream);
    *outStream << "\n";

    ROL::AugmentedLagrangian<RealT> augLag(obj,con,x,c);
    augLag.updateMultipliers(l,1.0);
    augLag.checkGradient(x, y, true, *outStream);
    augLag.checkHessVec(x, y, true, *outStream);

    // Optimization 
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist) );
    // Define status test.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-6);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-10);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    ROL::StatusTestSQP<RealT> status(*parlist);
    // Define step.
    ROL::AugmentedLagrangianStep<RealT> step(*parlist);
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);
    // Run Algorithm
    RealT zerotol = 0.0;
    //z.scale(50.0);
    con.solve(u,z,zerotol);
    obj.gradient_1(gu,u,z,zerotol);
    gu.scale(-1.0);
    con.applyInverseAdjointJacobian_1(l,gu,u,z,zerotol);
    gu.zero();
    c.zero();
    algo.run(x, g, l, c, obj, con, bnd, true, *outStream);

//    for ( int i = 0; i < nx+2; i++ ) {
//      std::cout << std::scientific << std::setprecision(10);
//      std::cout << std::setw(20) << std::left << (*z_rcp)[i];
//      if ( i == 0 ) {
//        std::cout << std::setw(20) << std::left << 1.0;
//      }
//      if ( i != 0 && i != nx+1 ) {
//        std::cout << std::setw(20) << std::left << (*u_rcp)[i-1];
//      }
//      if ( i == nx+1 ) {
//        std::cout << std::setw(20) << std::left << 0.0;
//      }
//      std::cout << "\n";
//    }
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
