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
class Objective_PoissonControl : public ROL::Objective<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;

  typedef typename vector::size_type uint;

private:
  std::vector<Real> u0_;
  Real alpha_;
  uint nx_;
  uint nt_;
  Real T_;
  Real dx_;
  Real dt_;

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector(); 
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();  
  }

public:

  Objective_PoissonControl(std::vector<Real> &u0, Real alpha = 1.e-4, uint nx = 128, uint nt = 100, Real T = 1) 
    : u0_(u0), alpha_(alpha), nx_(nx), nt_(nt), T_(T) {
    dx_ = 1.0/((Real)nx-1.0);
    dt_ = T/((Real)nt-1);
  }

  void apply_mass(std::vector<Real> &Mz, const std::vector<Real> &z ) {
    Mz.resize(nx_,0.0);
    for (uint i=0; i<nx_; i++) {
      if ( i == 0 ) {
        Mz[i] = dx_/6.0*(2.0*z[i] + z[i+1]);
      }
      else if ( i == nx_-1 ) {
        Mz[i] = dx_/6.0*(z[i-1] + 2.0*z[i]);
      }
      else {
        Mz[i] = dx_/6.0*(z[i-1] + 4.0*z[i] + z[i+1]);
      }
    }
  }

  void solve_state(std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(nx_,4.0*dx_/6.0 + dt_*2.0/dx_);
    d[0]           = dx_/3.0 + dt_/dx_;
    d[nx_-1] = dx_/3.0 + dt_/dx_;
    std::vector<Real> o(nx_-1,dx_/6.0 - dt_/dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nx_;
    int nhrs = 1;
    lp.PTTRF(nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    U.clear();
    U.resize(nt_+1);
    (U[0]).assign(u0_.begin(),u0_.end());
    // Time Step Using Implicit Euler
    std::vector<Real> b(nx_,0.0);
    for ( uint t = 0; t < nt_; t++ ) {
      // Get Right Hand Side
      apply_mass(b,U[t]);
      b[nx_-1] += dt_*z[t];
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (U[t+1]).assign(b.begin(),b.end());
    }
  }

  void solve_adjoint(std::vector<std::vector<Real> > &P, const std::vector<std::vector<Real> > &U) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(nx_,4.0*dx_/6.0 + dt_*2.0/dx_);
    d[0]           = dx_/3.0 + dt_/dx_;
    d[nx_-1] = dx_/3.0 + dt_/dx_;
    std::vector<Real> o(nx_-1,dx_/6.0 - dt_/dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nx_;
    int nhrs = 1;
    lp.PTTRF(nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    P.clear();
    P.resize(nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> b(nx_,0.0);
    for ( uint t = nt_; t > 0; t-- ) {
      // Get Right Hand Side
      if ( t == nt_ ) {
        std::vector<Real> res(nx_,0.0);
        for (uint i=0; i<nx_; i++) {
          res[i] = -((U[t])[i]-evaluate_target((Real)i*dx_));
        }
        apply_mass(b,res);
      }
      else {
        apply_mass(b,P[t]);
      }
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (P[t-1]).assign(b.begin(),b.end());
    }
  }

  void solve_state_sensitivity(std::vector<std::vector<Real> > &U, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(nx_,4.0*dx_/6.0 + dt_*2.0/dx_);
    d[0]           = dx_/3.0 + dt_/dx_;
    d[nx_-1] = dx_/3.0 + dt_/dx_;
    std::vector<Real> o(nx_-1,dx_/6.0 - dt_/dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nx_;
    int nhrs = 1;
    lp.PTTRF(nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    U.clear();
    U.resize(nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> b(nx_,0.0);
    for ( uint t = 0; t < nt_; t++ ) {
      // Get Right Hand Side
      if( t == 0 ) {
        b.resize(nx_,0.0);
        b[nx_-1] = -dt_*z[t];
      }
      else {
        apply_mass(b,U[t-1]);
        b[nx_-1] -= dt_*z[t];
      }
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
      // Update State Storage
      (U[t]).assign(b.begin(),b.end());
    }
  }

  void solve_adjoint_sensitivity(std::vector<std::vector<Real> > &P, const std::vector<std::vector<Real> > &U) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(nx_,4.0*dx_/6.0 + dt_*2.0/dx_);
    d[0]           = dx_/3.0 + dt_/dx_;
    d[nx_-1] = dx_/3.0 + dt_/dx_;
    std::vector<Real> o(nx_-1,dx_/6.0 - dt_/dx_);
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nx_;
    int nhrs = 1;
    lp.PTTRF(nx_,&d[0],&o[0],&info);
    // Initialize State Storage
    P.clear();
    P.resize(nt_);
    // Time Step Using Implicit Euler
    std::vector<Real> b(nx_,0.0);
    for ( uint t = nt_; t > 0; t-- ) {
      // Get Right Hand Side
      if ( t == nt_ ) {
        apply_mass(b,U[t-1]);
      }
      else {
        apply_mass(b,P[t]);
      }
      // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
      lp.PTTRS(nx_,nhrs,&d[0],&o[0],&b[0],ldb,&info);
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
      default: val = 1.0/3.0*std::pow(x,4.0) - 2.0/3.0*std::pow(x,3.0) + 1.0/3.0*x + 8.0*alpha_; break;
    }
    return val;
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<const vector> zp = getVector(z);

    // SOLVE STATE EQUATION
    std::vector<std::vector<Real> > U;
    this->solve_state(U,*zp);

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

    for (uint t=0; t<nt_; t++) {
      (*gp)[t] = dt_*(alpha_*(*zp)[t] - (P[t])[nx_-1]);
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<vector> hvp = getVector(hv);

    // SOLVE STATE SENSITIVITY EQUATION
    std::vector<std::vector<Real> > U;
    solve_state_sensitivity(U,*vp);

    // SOLVE ADJOINT SENSITIVITY EQUATION
    std::vector<std::vector<Real> > P;
    solve_adjoint_sensitivity(P,U);

    for (uint t=0; t<nt_; t++) {
      (*hvp)[t] = dt_*(alpha_*(*vp)[t] - (P[t])[nx_-1]);
    }
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>    vector;
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

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
    uint nx     = 100;   // Set spatial discretization.
    uint nt     = 300;   // Set temporal discretization.
    RealT T     = 1.0;   // Set end time.
    RealT alpha = 1.e-2; // Set penalty parameter.
    vector u0(nx,0.0); // Set initial conditions
    Objective_PoissonControl<RealT> obj(u0,alpha,nx,nt,T);

    // Initialize iteration vectors.
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(nt, 0.0);
    ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(nt, 0.0);

    for (uint i=0; i<nt; i++) {
      (*x_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*y_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    SV x(x_ptr);
    SV y(y_ptr);

    // Check deriatives.
    obj.checkGradient(x,y,true,*outStream);
    obj.checkHessVec(x,y,true,*outStream);

    // Initialize Constraints
    ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(nt,-1.0);
    ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(nt, 1.0);
    ROL::Ptr<V> lo = ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(u_ptr);

    ROL::Bounds<RealT> icon(lo,up);

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
   
    ROL::Ptr<V> diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm()/std::sqrt((RealT)nt-1.0);
    *outStream << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);
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

