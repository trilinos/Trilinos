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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

// ROL_Types contains predefined constants and objects
#include "ROL_Types.hpp"
// ROL algorithmic information
#include "ROL_Algorithm.hpp"
// ROL vectors
#include "ROL_StdVector.hpp"
// ROL objective functions and constraints
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_StochasticProblem.hpp"
// ROL sample generators
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"

template<class Real>
class EqualityConstraint_BurgersControl : public ROL::EqualityConstraint_SimOpt<Real> {
private:
  int nx_;
  Real dx_;

  Real compute_norm(const std::vector<Real> &r) {
    return std::sqrt(dot(r,r));
  }

  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) {
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

  using ROL::EqualityConstraint_SimOpt<Real>::update;

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
    r.clear(); r.resize(nx_,0.0);
    const std::vector<Real> param =
      ROL::EqualityConstraint_SimOpt<Real>::getParameter();
    Real nu = std::pow(10.0,param[0]-2.0);
    Real f  = param[1]/100.0;
    Real u0 = 1.0+param[2]/1000.0;
    Real u1 = param[3]/1000.0;
    for (int i=0; i<nx_; i++) {
      // Contribution from stiffness term
      if ( i==0 ) {
        r[i] = nu/dx_*(2.0*u[i]-u[i+1]);
      }
      else if (i==nx_-1) {
        r[i] = nu/dx_*(2.0*u[i]-u[i-1]);
      }
      else {
        r[i] = nu/dx_*(2.0*u[i]-u[i-1]-u[i+1]);
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
      r[i] -= dx_*f;
    }
    // Contribution from Dirichlet boundary terms
    r[    0] -= u0*u[    0]/6.0 + u0*u0/6.0 + nu*u0/dx_;
    r[nx_-1] += u1*u[nx_-1]/6.0 + u1*u1/6.0 - nu*u1/dx_;
  }

  void compute_pde_jacobian(std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
                      const std::vector<Real> &u) {
    const std::vector<Real> param =
      ROL::EqualityConstraint_SimOpt<Real>::getParameter();
    Real nu = std::pow(10.0,param[0]-2.0);
    Real u0 = 1.0+param[2]/1000.0;
    Real u1 = param[3]/1000.0;
    // Get Diagonal and Off-Diagonal Entries of linear PDE Jacobian
    d.clear();  d.resize(nx_,nu*2.0/dx_);
    dl.clear(); dl.resize(nx_-1,-nu/dx_);
    du.clear(); du.resize(nx_-1,-nu/dx_);
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
    d[    0] -= u0/6.0;
    d[nx_-1] += u1/6.0;
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
              const std::vector<Real> &r, const bool transpose = false) {
    u.assign(r.begin(),r.end());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    std::vector<Real> du2(nx_-2,0.0);
    std::vector<int> ipiv(nx_,0);
    int info;
    int ldb  = nx_;
    int nhrs = 1;
    lp.GTTRF(nx_,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&info);
    char trans = 'N';
    if ( transpose ) { 
      trans = 'T';
    }
    lp.GTTRS(trans,nx_,nhrs,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&u[0],ldb,&info);
  }

public:

  EqualityConstraint_BurgersControl(int nx = 128) : nx_(nx), dx_(1.0/((Real)nx+1.0)) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > cp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(c)).getVector());
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    compute_residual(*cp,*up,*zp);
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > up =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(u)).getVector());
    up->assign(up->size(),static_cast<Real>(1));
    ROL::EqualityConstraint_SimOpt<Real>::solve(c,u,z,tol);
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
    const std::vector<Real> param =
      ROL::EqualityConstraint_SimOpt<Real>::getParameter();
    Real nu = std::pow(10.0,param[0]-2.0);
    Real u0 = 1.0+param[2]/1000.0;
    Real u1 = param[3]/1000.0;
    // Fill jvp
    for (int i = 0; i < nx_; i++) {
      (*jvp)[i] = nu/dx_*2.0*(*vp)[i];
      if ( i > 0 ) {
        (*jvp)[i] += -nu/dx_*(*vp)[i-1]
                     -(*up)[i-1]/6.0*(*vp)[i] 
                     -((*up)[i]+2.0*(*up)[i-1])/6.0*(*vp)[i-1];
      }
      if ( i < nx_-1 ) {
        (*jvp)[i] += -nu/dx_*(*vp)[i+1]
                     +(*up)[i+1]/6.0*(*vp)[i] 
                     +((*up)[i]+2.0*(*up)[i+1])/6.0*(*vp)[i+1];
      }
    }
    (*jvp)[    0] -= u0/6.0*(*vp)[0];
    (*jvp)[nx_-1] += u1/6.0*(*vp)[nx_-1];
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
    for (int i=0; i<nx_; i++) {
      // Contribution from control
      (*jvp)[i] = -dx_/6.0*((*vp)[i]+4.0*(*vp)[i+1]+(*vp)[i+2]);
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
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    std::vector<Real> du(nx_-1,0.0);
    compute_pde_jacobian(dl,d,du,*up);
    // Solve solve state sensitivity system at current time step
    linear_solve(*ijvp,dl,d,du,*vp);
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
    const std::vector<Real> param =
      ROL::EqualityConstraint_SimOpt<Real>::getParameter();
    Real nu = std::pow(10.0,param[0]-2.0);
    Real u0 = 1.0+param[2]/1000.0;
    Real u1 = param[3]/1000.0;
    // Fill jvp
    for (int i = 0; i < nx_; i++) {
      (*jvp)[i] = nu/dx_*2.0*(*vp)[i];
      if ( i > 0 ) {
        (*jvp)[i] += -nu/dx_*(*vp)[i-1] 
                     -(*up)[i-1]/6.0*(*vp)[i] 
                     +((*up)[i-1]+2.0*(*up)[i])/6.0*(*vp)[i-1];
      }
      if ( i < nx_-1 ) {
        (*jvp)[i] += -nu/dx_*(*vp)[i+1] 
                     +(*up)[i+1]/6.0*(*vp)[i]
                     -((*up)[i+1]+2.0*(*up)[i])/6.0*(*vp)[i+1];
      }
    }
    (*jvp)[    0] -= u0/6.0*(*vp)[0];
    (*jvp)[nx_-1] += u1/6.0*(*vp)[nx_-1];
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
    for (int i=0; i<nx_+2; i++) {
      if ( i == 0 ) {
        (*jvp)[i] = -dx_/6.0*(*vp)[i];
      }
      else if ( i == 1 ) {
        (*jvp)[i] = -dx_/6.0*(4.0*(*vp)[i-1]+(*vp)[i]);
      }
      else if ( i == nx_ ) {
        (*jvp)[i] = -dx_/6.0*(4.0*(*vp)[i-1]+(*vp)[i-2]);
      }
      else if ( i == nx_+1 ) {
        (*jvp)[i] = -dx_/6.0*(*vp)[i-2];
      }
      else {
        (*jvp)[i] = -dx_/6.0*((*vp)[i-2]+4.0*(*vp)[i-1]+(*vp)[i]);
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
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> du(nx_-1,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    compute_pde_jacobian(dl,d,du,*up);
    // Solve solve adjoint system at current time step
    linear_solve(*iajvp,dl,d,du,*vp,true);
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
    for (int i=0; i<nx_; i++) {
      // Contribution from nonlinear term
      (*ahwvp)[i] = 0.0;
      if (i<nx_-1){
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

  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
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
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:

  Objective_BurgersControl(Real alpha = 1.e-4, int nx = 128) : alpha_(alpha), nx_(nx) {
    dx_ = 1.0/((Real)nx+1.0);
  }

  using ROL::Objective_SimOpt<Real>::value;

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > up =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE RESIDUAL
    Real res1 = 0.0, res2 = 0.0, res3 = 0.0;
    Real valu = 0.0, valz = dot(*zp,*zp);
    for (int i=0; i<nx_; i++) {
      if ( i == 0 ) {
        res1  = (*up)[i]-evaluate_target((Real)(i+1)*dx_);
        res2  = (*up)[i+1]-evaluate_target((Real)(i+2)*dx_);
        valu += dx_/6.0*(4.0*res1 + res2)*res1;
      }
      else if ( i == nx_-1 ) {
        res1  = (*up)[i-1]-evaluate_target((Real)i*dx_);
        res2  = (*up)[i]-evaluate_target((Real)(i+1)*dx_);
        valu += dx_/6.0*(res1 + 4.0*res2)*res2;
      }
      else {
        res1  = (*up)[i-1]-evaluate_target((Real)i*dx_);
        res2  = (*up)[i]-evaluate_target((Real)(i+1)*dx_);
        res3  = (*up)[i+1]-evaluate_target((Real)(i+2)*dx_);
        valu += dx_/6.0*(res1 + 4.0*res2 + res3)*res2;
      }
    }
    return 0.5*(valu + alpha_*valz);
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
    std::vector<Real> diff(nx_,0.0);
    for (int i=0; i<nx_; i++) {
      diff[i] = ((*up)[i]-evaluate_target((Real)(i+1)*dx_));
    }
    apply_mass(*gup,diff);
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
    for (int i=0; i<nx_+2; i++) {
      if (i==0) {
        (*gzp)[i] = alpha_*dx_/6.0*(2.0*(*zp)[i]+(*zp)[i+1]);
      }
      else if (i==nx_+1) {
        (*gzp)[i] = alpha_*dx_/6.0*(2.0*(*zp)[i]+(*zp)[i-1]);
      }
      else {
        (*gzp)[i] = alpha_*dx_/6.0*((*zp)[i-1]+4.0*(*zp)[i]+(*zp)[i+1]);
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
    apply_mass(*hvup,*vup);
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
    for (int i=0; i<nx_+2; i++) {
      if (i==0) {
        (*hvzp)[i] = alpha_*dx_/6.0*(2.0*(*vzp)[i]+(*vzp)[i+1]);
      }
      else if (i==nx_+1) {
        (*hvzp)[i] = alpha_*dx_/6.0*(2.0*(*vzp)[i]+(*vzp)[i-1]);
      }
      else {
        (*hvzp)[i] = alpha_*dx_/6.0*((*vzp)[i-1]+4.0*(*vzp)[i]+(*vzp)[i+1]);
      }
    }
  }
};
