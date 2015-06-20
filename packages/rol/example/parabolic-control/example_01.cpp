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
class Objective_PoissonControl : public ROL::Objective<Real> {
private:
  std::vector<Real> u0_;
  Real alpha_;
  int nx_;
  int nt_;
  Real T_;
  Real dx_;
  Real dt_;

public:

  Objective_PoissonControl(std::vector<Real> &u0, Real alpha = 1.e-4, int nx = 128, int nt = 100, Real T = 1) 
    : u0_(u0), alpha_(alpha), nx_(nx), nt_(nt), T_(T) {
    dx_ = 1.0/((Real)nx-1.0);
    dt_ = T/((Real)nt-1);
  }

  void apply_mass(std::vector<Real> &Mz, const std::vector<Real> &z ) {
    Mz.resize(this->nx_,0.0);
    for (int i=0; i<this->nx_; i++) {
      if ( i == 0 ) {
        Mz[i] = this->dx_/6.0*(2.0*z[i] + z[i+1]);
      }
      else if ( i == this->nx_-1 ) {
        Mz[i] = this->dx_/6.0*(z[i-1] + 2.0*z[i]);
      }
      else {
        Mz[i] = this->dx_/6.0*(z[i-1] + 4.0*z[i] + z[i+1]);
      }
    }
  }

  void solve_state(std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nx_,4.0*this->dx_/6.0 + this->dt_*2.0/this->dx_);
    d[0]           = this->dx_/3.0 + this->dt_/this->dx_;
    d[this->nx_-1] = this->dx_/3.0 + this->dt_/this->dx_;
    std::vector<Real> o(this->nx_-1,this->dx_/6.0 - this->dt_/this->dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nx_;
    int nhrs = 1;
    lp.PTTRF(this->nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    U.clear();
    U.resize(this->nt_+1);
    (U[0]).assign(this->u0_.begin(),this->u0_.end());
    // Time Step Using Implicit Euler
    std::vector<Real> b(this->nx_,0.0);
    for ( int t = 0; t < this->nt_; t++ ) {
      // Get Right Hand Side
      this->apply_mass(b,U[t]);
      b[this->nx_-1] += this->dt_*z[t];
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(this->nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (U[t+1]).assign(b.begin(),b.end());
    }
  }

  void solve_adjoint(std::vector<std::vector<Real> > &P, const std::vector<std::vector<Real> > &U) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nx_,4.0*this->dx_/6.0 + this->dt_*2.0/this->dx_);
    d[0]           = this->dx_/3.0 + this->dt_/this->dx_;
    d[this->nx_-1] = this->dx_/3.0 + this->dt_/this->dx_;
    std::vector<Real> o(this->nx_-1,this->dx_/6.0 - this->dt_/this->dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nx_;
    int nhrs = 1;
    lp.PTTRF(this->nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    P.clear();
    P.resize(this->nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> b(this->nx_,0.0);
    for ( int t = this->nt_; t > 0; t-- ) {
      // Get Right Hand Side
      if ( t == this->nt_ ) {
        std::vector<Real> res(this->nx_,0.0);
        for (int i=0; i<this->nx_; i++) {
          res[i] = -((U[t])[i]-this->evaluate_target((Real)i*this->dx_));
        }
        this->apply_mass(b,res);
      }
      else {
        this->apply_mass(b,P[t]);
      }
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(this->nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (P[t-1]).assign(b.begin(),b.end());
    }
  }

  void solve_state_sensitivity(std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nx_,4.0*this->dx_/6.0 + this->dt_*2.0/this->dx_);
    d[0]           = this->dx_/3.0 + this->dt_/this->dx_;
    d[this->nx_-1] = this->dx_/3.0 + this->dt_/this->dx_;
    std::vector<Real> o(this->nx_-1,this->dx_/6.0 - this->dt_/this->dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nx_;
    int nhrs = 1;
    lp.PTTRF(this->nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    U.clear();
    U.resize(this->nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> b(this->nx_,0.0);
    for ( int t = 0; t < this->nt_; t++ ) {
      // Get Right Hand Side
      if( t == 0 ) {
        b.resize(this->nx_,0.0);
        b[this->nx_-1] = -this->dt_*z[t];
      }
      else {
        this->apply_mass(b,U[t-1]);
        b[this->nx_-1] -= this->dt_*z[t];
      }
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(this->nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (U[t]).assign(b.begin(),b.end());
    }
  }

  void solve_adjoint_sensitivity(std::vector<std::vector<Real> > &P, const std::vector<std::vector<Real> > &U) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nx_,4.0*this->dx_/6.0 + this->dt_*2.0/this->dx_);
    d[0]           = this->dx_/3.0 + this->dt_/this->dx_;
    d[this->nx_-1] = this->dx_/3.0 + this->dt_/this->dx_;
    std::vector<Real> o(this->nx_-1,this->dx_/6.0 - this->dt_/this->dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nx_;
    int nhrs = 1;
    lp.PTTRF(this->nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    P.clear();
    P.resize(this->nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> b(this->nx_,0.0);
    for ( int t = this->nt_; t > 0; t-- ) {
      // Get Right Hand Side
      if ( t == this->nt_ ) {
        this->apply_mass(b,U[t-1]);
      }
      else {
        this->apply_mass(b,P[t]);
      }
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(this->nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (P[t-1]).assign(b.begin(),b.end());
    }
  }

  Real evaluate_target(Real x) {
    Real val = 0.0;
    int example = 1;
    switch (example) {
      case 1:  val = ((x<0.5) ? 1.0 : 0.0); break;
      case 2:  val = 1.0; break;
      default: val = 1.0/3.0*std::pow(x,4.0) - 2.0/3.0*std::pow(x,3.0) + 1.0/3.0*x + 8.0*this->alpha_; break;
    }
    return val;
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // SOLVE STATE EQUATION
    std::vector<std::vector<Real> > U;
    this->solve_state(U,*zp);

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

    for (int t=0; t<this->nt_; t++) {
      (*gp)[t] = this->dt_*(this->alpha_*(*zp)[t] - (P[t])[this->nx_-1]);
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<std::vector<Real> > hvp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

    // SOLVE STATE SENSITIVITY EQUATION
    std::vector<std::vector<Real> > U;
    this->solve_state_sensitivity(U,*vp);

    // SOLVE ADJOINT SENSITIVITY EQUATION
    std::vector<std::vector<Real> > P;
    this->solve_adjoint_sensitivity(P,U);

    for (int t=0; t<this->nt_; t++) {
      (*hvp)[t] = this->dt_*(this->alpha_*(*vp)[t] - (P[t])[this->nx_-1]);
    }
  }
};

template<class Real>
class BoundConstraint_PoissonControl : public ROL::BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
public:
  BoundConstraint_PoissonControl(std::vector<Real> &lo, std::vector<Real> &up) {
    dim_ = lo.size();
    x_lo_.clear();
    x_lo_.assign(lo.begin(),lo.end());
    x_up_.clear();
    x_up_.assign(up.begin(),up.end());
    for ( unsigned i = 0; i < (unsigned)dim_; i++ ) {
      if ( i == 0 ) {
        min_diff_ = x_up_[i]-x_lo_[i];
      }
      else {
        min_diff_ = std::min(min_diff_,x_up_[i]-x_lo_[i]);
      }
    }
    min_diff_ *= 0.5;
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
//  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
//    Teuchos::RCP<const std::vector<Real> > ex =
//      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
//    Teuchos::RCP<std::vector<Real> > ev =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
//    Real epsn = std::min(eps,this->min_diff_);
//    for ( int i = 0; i < this->dim_; i++ ) {
//      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ||
//           ((*ex)[i] >= this->x_up_[i]-epsn) ) {
//        (*ev)[i] = 0.0;
//      }
//    }
//  }
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
//  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
//    Teuchos::RCP<const std::vector<Real> > ex =
//      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
//    Teuchos::RCP<const std::vector<Real> > eg =
//      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
//    Teuchos::RCP<std::vector<Real> > ev =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
//    Real epsn = std::min(eps,this->min_diff_);
//    for ( int i = 0; i < this->dim_; i++ ) {
//      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ||
//           ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
//        (*ev)[i] = 0.0;
//      }
//    }
//  }
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
    RealT alpha = 1.e-2; // Set penalty parameter.
    std::vector<RealT> u0(nx,0.0); // Set initial conditions
    Objective_PoissonControl<RealT> obj(u0,alpha,nx,nt,T);
    // Initialize iteration vectors.
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (nt, 0.0) );
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
    std::vector<RealT> lo(nt,-1.0);
    std::vector<RealT> up(nt,1.0);
    BoundConstraint_PoissonControl<RealT> icon(lo,up);

    // Primal dual active set.
    Teuchos::ParameterList parlist;
    // Krylov parameters.
    parlist.set("Absolute Krylov Tolerance",              1.e-8);
    parlist.set("Relative Krylov Tolerance",              1.e-4);
    parlist.set("Maximum Number of Krylov Iterations",    50);
    // PDAS parameters.
    parlist.set("PDAS Relative Step Tolerance",           1.e-8);
    parlist.set("PDAS Relative Gradient Tolerance",       1.e-6);
    parlist.set("PDAS Maximum Number of Iterations",      1);
    parlist.set("PDAS Dual Scaling",                      alpha);      
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

