// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve a linear-quadratic parabolic control problem 
           with bound constraints.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Bounds.hpp"

template<class Real>
class Objective_ParabolicControl : public ROL::Objective<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;

  typedef typename vector::size_type  uint;  

private:
  std::vector<Real> u0_;
  Real alpha_;
  uint nx_;
  uint nt_;
  Real T_;
  Real dx_;
  Real dt_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/

  ROL::Ptr<const vector> getVector( const V& x ) {
     
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

  void apply_mass(std::vector<Real> &Mu, const std::vector<Real> &u ) {
    Mu.resize(nx_,0.0);
    for (uint i=0; i<nx_; i++) {
      if ( i == 0 ) {
        Mu[i] = dx_/6.0*(2.0*u[i] + u[i+1]);
      }
      else if ( i == nx_-1 ) {
        Mu[i] = dx_/6.0*(u[i-1] + 2.0*u[i]);
      }
      else {
        Mu[i] = dx_/6.0*(u[i-1] + 4.0*u[i] + u[i+1]);
      }
    }
  }

  void compute_pde_jacobian(std::vector<Real> &d, std::vector<Real> &o, const std::vector<Real> &u) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    d.clear();
    d.resize(nx_,4.0*dx_/6.0 + dt_*2.0/dx_);
    d[0]           = dx_/3.0 + dt_/dx_;
    d[nx_-1] = dx_/3.0 + dt_/dx_ + dt_*4.0*std::pow(u[nx_-1],3.0);
    o.clear();
    o.resize(nx_-1,dx_/6.0 - dt_/dx_);
  }

  void compute_residual(std::vector<Real> &r, 
                  const std::vector<Real> &up, const std::vector<Real> &u, const Real z) {
    r.clear();
    r.resize(nx_,0.0);
    for (uint i=0; i<nx_; i++) {
      if ( i==0 ) {
        r[i] = dx_/6.0*(2.0*u[i]+u[i+1]) + dt_/dx_*(u[i]-u[i+1]);
        r[i]-= dx_/6.0*(2.0*up[i]+up[i+1]); // Contribution from previous state
      }
      else if ( i==nx_-1 ) {
        r[i] = dx_/6.0*(u[i-1]+2.0*u[i]) + dt_/dx_*(u[i]-u[i-1]);
        r[i]+= dt_*std::pow(u[i],4.0); // Stefan-Boltzmann boundary condition
        r[i]-= dx_/6.0*(2.0*up[i]+up[i-1]); // Contribution from previous state
        r[i]-= dt_*z; // Contribution from control
      }
      else {
        r[i] = dx_/6.0*(u[i-1]+4.0*u[i]+u[i+1]) + dt_/dx_*(2.0*u[i]-u[i-1]-u[i+1]);
        r[i]-= dx_/6.0*(up[i-1]+4.0*up[i]+up[i+1]); // Contribution from previous state
      }
    }
  }

  Real compute_norm(const std::vector<Real> &r) {
    Real norm = 0.0;
    for (uint i=0; i<r.size(); i++) {
      norm += r[i]*r[i];
    }
    return std::sqrt(norm);
  }

  void update(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    for (uint i=0; i<u.size(); i++) {
      u[i] += alpha*s[i];
    }
  }

  void linear_solve(std::vector<Real> &u, std::vector<Real> &d, std::vector<Real> &o, 
              const std::vector<Real> &r) {
    u.assign(r.begin(),r.end());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int nx   = static_cast<int>(nx_);
    int ldb  = nx;
    int nhrs = 1;
    lp.PTTRF(nx,&d[0],&o[0],&info);
    lp.PTTRS(nx,nhrs,&d[0],&o[0],&u[0],ldb,&info);
  }

  void run_newton(std::vector<Real> &u, const std::vector<Real> &up, const Real z) {
    // Set initial guess
    u.assign(up.begin(),up.end());
    // Compute residual and residual norm
    std::vector<Real> r(u.size(),0.0);
    compute_residual(r,up,u,z);
    Real rnorm = compute_norm(r);
    // Define tolerances
    Real tol   = 1.e2*ROL::ROL_EPSILON<Real>();
    Real maxit = 100;
    // Initialize Jacobian storage
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    // Iterate Newton's method
    Real alpha = 1.0, tmp = 0.0;
    std::vector<Real> s(nx_,0.0);
    std::vector<Real> utmp(nx_,0.0);
    for (uint i=0; i<maxit; i++) {
      //std::cout << i << "  " << rnorm << "\n";
      // Get Jacobian
      compute_pde_jacobian(d,o,u);
      // Solve Newton system
      linear_solve(s,d,o,r);
      // Perform line search
      tmp = rnorm;
      alpha = 1.0;
      utmp.assign(u.begin(),u.end());
      update(utmp,s,-alpha);
      compute_residual(r,up,utmp,z);
      rnorm = compute_norm(r); 
      while ( rnorm > (1.0-1.e-4*alpha)*tmp && alpha > std::sqrt(ROL::ROL_EPSILON<Real>()) ) {
        alpha /= 2.0;
        utmp.assign(u.begin(),u.end());
        update(utmp,s,-alpha);
        compute_residual(r,up,utmp,z);
        rnorm = compute_norm(r); 
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

  Objective_ParabolicControl(Real alpha = 1.e-4, uint nx = 128, uint nt = 100, Real T = 1) 
    : alpha_(alpha), nx_(nx), nt_(nt), T_(T) {
    u0_.resize(nx_,0.0);
    dx_ = 1.0/((Real)nx-1.0);
    dt_ = T/((Real)nt-1.0);
  }

  void solve_state(std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Initialize State Storage
    U.clear();
    U.resize(nt_+1);
    (U[0]).assign(u0_.begin(),u0_.end());
    std::vector<Real> up(u0_);
    std::vector<Real> u(u0_);
    // Time Step Using Implicit Euler
    for ( uint t = 0; t < nt_; t++ ) {
      this->run_newton(u,up,z[t]);
      (U[t+1]).assign(u.begin(),u.end());
      up.assign(u.begin(),u.end());
    }
  }

  void solve_adjoint(std::vector<std::vector<Real> > &P, const std::vector<std::vector<Real> > &U) {
    // Initialize State Storage
    P.clear();
    P.resize(nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> p(nx_,0.0);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    for ( uint t = nt_; t > 0; t-- ) {
      // Get PDE Jacobian
      compute_pde_jacobian(d,o,U[t]);
      // Get Right Hand Side
      r.assign(nx_,0.0);
      if ( t==nt_ ) {
        std::vector<Real> diff(nx_,0.0);
        for (uint i=0; i<nx_; i++) {
          diff[i] = -((U[t])[i]-evaluate_target((Real)i*dx_));
        }
        apply_mass(r,diff);
      }
      else {
        apply_mass(r,P[t]);
      } 
      // Solve solve adjoint system at current time step
      linear_solve(p,d,o,r);
      // Update State Storage
      (P[t-1]).assign(p.begin(),p.end());
    }
  }

  void solve_state_sensitivity(std::vector<std::vector<Real> > &V, 
                         const std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Initialize State Storage
    V.clear();
    V.resize(nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> v(nx_,0.0);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    for ( uint t = 0; t < nt_; t++ ) {
      // Get PDE Jacobian
      compute_pde_jacobian(d,o,U[t+1]);
      // Get Right Hand Side
      if( t == 0 ) {
        r.assign(nx_,0.0);
        r[nx_-1] = dt_*z[t];
      }
      else {
        apply_mass(r,V[t-1]);
        r[nx_-1] += dt_*z[t];
      }
      // Solve solve adjoint system at current time step
      linear_solve(v,d,o,r);
      // Update State Storage
      (V[t]).assign(v.begin(),v.end());
    }
  }

  void solve_adjoint_sensitivity(std::vector<std::vector<Real> > &Q, 
                           const std::vector<std::vector<Real> > &U, const std::vector<std::vector<Real> > &P,
                           const std::vector<std::vector<Real> > &V, const std::vector<Real> &z) {
    // Initialize State Storage
    Q.clear();
    Q.resize(nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> q(nx_,0.0);
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> r(nx_,0.0);
    std::vector<Real> o(nx_-1,0.0);
    for ( uint t = nt_; t > 0; t-- ) {
      // Get PDE Jacobian
      compute_pde_jacobian(d,o,U[t]);
      // Get Right Hand Side
      if ( t == nt_ ) {
        std::vector<Real> tmp(nx_,0.0);
        r.assign(nx_,0.0);
        apply_mass(tmp,V[t-1]);
        update(r,tmp,-1.0);
      }
      else {
        apply_mass(r,Q[t]);
      }
      r[nx_-1]-=dt_*12.0*std::pow(U[t][nx_-1],2.0)*P[t-1][nx_-1]*V[t-1][nx_-1];
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      linear_solve(q,d,o,r);
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

    
    ROL::Ptr<const vector> zp = getVector(z);

    // SOLVE STATE EQUATION
    std::vector<std::vector<Real> > U;
    solve_state(U,*zp);
    // COMPUTE RESIDUAL
    Real val  = 0.0;
    Real res  = 0.0;
    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (uint t=0; t<nt_; t++) {
      val += (*zp)[t]*(*zp)[t];
    }
    val *= 0.5*alpha_*dt_;

    for (uint i=0; i<nx_; i++) {
      if ( i == 0 ) {
        res1 = (U[nt_])[i]-evaluate_target((Real)i*dx_);
        res2 = (U[nt_])[i+1]-evaluate_target((Real)(i+1)*dx_);
        res  = dx_/6.0*(2.0*res1 + res2)*res1;
      }
      else if ( i == nx_-1 ) {
        res1 = (U[nt_])[i-1]-evaluate_target((Real)(i-1)*dx_);
        res2 = (U[nt_])[i]-evaluate_target((Real)i*dx_);
        res  = dx_/6.0*(res1 + 2.0*res2)*res2;
      }
      else {
        res1 = (U[nt_])[i-1]-evaluate_target((Real)(i-1)*dx_);
        res2 = (U[nt_])[i]-evaluate_target((Real)i*dx_);
        res3 = (U[nt_])[i+1]-evaluate_target((Real)(i+1)*dx_);
        res  = dx_/6.0*(res1 + 4.0*res2 + res3)*res2;
      }
      val += 0.5*res;
    }
    return val;
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<vector> gp = getVector(g);

    // SOLVE STATE EQUATION
    std::vector<std::vector<Real> > U;
    solve_state(U,*zp);
    // SOLVE ADJOINT EQUATION
    std::vector<std::vector<Real> > P;
    solve_adjoint(P,U);
    // COMPUTE GRADIENT
    for (uint t=0; t<nt_; t++) {
      (*gp)[t] = dt_*(alpha_*(*zp)[t] - (P[t])[nx_-1]);
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<vector> hvp = getVector(hv);

    // SOLVE STATE EQUATION
    std::vector<std::vector<Real> > U;
    solve_state(U,*zp);
    // SOLVE ADJOINT EQUATION
    std::vector<std::vector<Real> > P;
    solve_adjoint(P,U);
    // SOLVE STATE SENSITIVITY EQUATION
    std::vector<std::vector<Real> > V;
    solve_state_sensitivity(V,U,*vp);
    // SOLVE ADJOINT SENSITIVITY EQUATION
    std::vector<std::vector<Real> > Q;
    solve_adjoint_sensitivity(Q,U,P,V,*vp);
    // COMPUTE HESSVEC
    for (uint t=0; t<nt_; t++) {
      (*hvp)[t] = dt_*(alpha_*(*vp)[t] - (Q[t])[nx_-1]);
    }
  }
};


typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>         vector;
  typedef typename vector::size_type uint;

    
  

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    // Initialize objective function.
    uint nx      = 100;   // Set spatial discretization.
    uint nt      = 300;   // Set temporal discretization.
    RealT T     = 1.0;   // Set end time.
    RealT alpha = 1.e-3; // Set penalty parameter.
    Objective_ParabolicControl<RealT> obj(alpha,nx,nt,T);
    // Initialize iteration vectors.
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(nt, 1.0);
    ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(nt, 0.0);

    for (uint i=0; i<nt; i++) {
      (*x_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*y_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }

    ROL::StdVector<RealT> x(x_ptr);
    ROL::StdVector<RealT> y(y_ptr);
    // Check deriatives.
    obj.checkGradient(x,y,true,*outStream);
    obj.checkHessVec(x,y,true,*outStream);

    // Initialize Constraints
    ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(nt,0.0);
    ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(nt,1.0);

    ROL::Ptr<ROL::Vector<RealT> > lo = ROL::makePtr<ROL::StdVector<RealT>>(l_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtr<ROL::StdVector<RealT>>(u_ptr);

    ROL::Bounds<RealT> icon(lo,up);
//    BoundConstraint_ParabolicControl<RealT> icon(nt);

    // Primal dual active set.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    // Krylov parameters.
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-8);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);
    // PDAS parameters.
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-6);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 1);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);
    // Status test parameters.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    // Define algorithm.
    ROL::Ptr<ROL::Step<RealT>>       step   = ROL::makePtr<ROL::PrimalDualActiveSetStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Ptr<ROL::Algorithm<RealT>>  algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
    // Run algorithm.
    x.zero();
    algo->run(x, obj, icon, true, *outStream);
    // Output control to file.
    std::ofstream file;
    file.open("control_PDAS.txt");
    for ( uint i = 0; i < nt; i++ ) {
      file << (*x_ptr)[i] << "\n";
    }
    file.close();

    // Projected Newton.
    // re-load parameters
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // Set algorithm.
    step   = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
    // Run Algorithm
    y.zero();
    algo->run(y, obj, icon, true, *outStream);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( uint i = 0; i < nt; i++ ) {
      file_tr << (*y_ptr)[i] << "\n";
    }
    file_tr.close();
   
    ROL::Ptr<ROL::Vector<RealT> > diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm()/std::sqrt((RealT)nt-1.0);
    *outStream << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);

    // Output state to file.
    std::vector<std::vector<RealT> > U(nt);
    obj.solve_state(U,*y_ptr);
    std::ofstream file1;
    file1.open("state_tx.txt");
    for (uint t=0; t<nt; t++) {
      file1 << t*(T/((RealT)nt-1.0)) << "  ";
      for (uint i=0; i<nx; i++) {
        file1 << (U[t])[i] << "  ";
      }
      file1 << "\n";
    }
    file1.close();
    std::ofstream file2;
    file2.open("state_xt.txt");
    for (uint i=0; i<nx; i++) {
      file2 << i*(1.0/((RealT)nx-1.0)) << "  ";
      for (uint t=0; t<nt; t++) {
        file2 << (U[t])[i] << "  ";
      }
      file2 << "\n";
    }
    file2.close();
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

