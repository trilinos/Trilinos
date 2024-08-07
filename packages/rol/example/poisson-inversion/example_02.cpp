// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the inverse Poisson problem using trust-region
           methods with dense Hessian diagnostics.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"
#include "ROL_Stream.hpp"


#include "Teuchos_GlobalMPISession.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

template<class Real>
class Objective_PoissonInversion : public ROL::Objective<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
  
  typedef typename vector::size_type uint;  


private:
  uint nu_;
  uint nz_;

  Real hu_;
  Real hz_;

  Real u0_;
  Real u1_;

  Real alpha_;

  bool useCorrection_;
  Teuchos::SerialDenseMatrix<int, Real> H_;

  ROL::Ptr<const vector> getVector( const V& x ) {
     
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  /* CONSTRUCTOR */
  Objective_PoissonInversion(int nz = 32, Real alpha = 1.e-4) 
    : nz_(nz), u0_(1.0), u1_(0.0), alpha_(alpha), useCorrection_(false) {
    nu_ = nz_-1;
    hu_ = 1.0/((Real)nu_+1.0);
    hz_ = hu_; 
  }

  void apply_mass(std::vector<Real> &Mz, const std::vector<Real> &z ) {
    Mz.resize(nu_,0.0);
    for (uint i=0; i<nu_; i++) {
      if ( i == 0 ) {
        Mz[i] = hu_/6.0*(2.0*z[i] + z[i+1]);
      }
      else if ( i == nu_-1 ) {
        Mz[i] = hu_/6.0*(z[i-1] + 2.0*z[i]);
      }
      else {
        Mz[i] = hu_/6.0*(z[i-1] + 4.0*z[i] + z[i+1]);
      }
    }
  }

  Real evaluate_target(Real x) {
    return (x <= 0.5) ? 1.0 : 0.0;
  }

  void apply_linearized_control_operator( std::vector<Real> &Bd, const std::vector<Real> &z, 
                                    const std::vector<Real> &d,  const std::vector<Real> &u,
                                          bool addBC = true ) {
    Bd.clear();
    Bd.resize(nu_,0.0);
    for (uint i = 0; i < nu_; i++) {
      if ( i == 0 ) {
        Bd[i] = 1.0/hu_*( u[i]*d[i] + (u[i]-u[i+1])*d[i+1] );
      }
      else if ( i == nu_-1 ) {
        Bd[i] = 1.0/hu_*( (u[i]-u[i-1])*d[i] + u[i]*d[i+1] );
      }
      else {
        Bd[i] = 1.0/hu_*( (u[i]-u[i-1])*d[i] + (u[i]-u[i+1])*d[i+1] );
      }
    }
    if ( addBC ) {
      Bd[    0] -= u0_*d[    0]/hu_;
      Bd[nu_-1] -= u1_*d[nz_-1]/hu_;
    }
  }

  void apply_transposed_linearized_control_operator( std::vector<Real> &Bd, const std::vector<Real> &z,
                                               const std::vector<Real> &d,  const std::vector<Real> &u,
                                                     bool addBC = true ) {
    Bd.clear();
    Bd.resize(nz_,0.0);
    for (uint i = 0; i < nz_; i++) {
      if ( i == 0 ) {
        Bd[i] = 1.0/hu_*u[i]*d[i];
      }
      else if ( i == nz_-1 ) {
        Bd[i] = 1.0/hu_*u[i-1]*d[i-1];
      }
      else {
        Bd[i] = 1.0/hu_*( (u[i]-u[i-1])*(d[i]-d[i-1]) );
      }
    }
    if ( addBC ) {
      Bd[    0] -= u0_*d[    0]/hu_;
      Bd[nz_-1] -= u1_*d[nu_-1]/hu_;
    }
  }

  /* STATE AND ADJOINT EQUATION DEFINTIONS */
  void solve_state_equation(std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(nu_,1.0);
    std::vector<Real> o(nu_-1,1.0);
    for ( uint i = 0; i < nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/hu_;
      if ( i < nu_-1 ) {
        o[i] *= -z[i+1]/hu_;
      }
    }
    // Set right hand side
    u.clear();
    u.resize(nu_,0.0);
    u[    0] = z[    0]/hu_ * u0_;
    u[nu_-1] = z[nz_-1]/hu_ * u1_;
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nu_;
    int nhrs = 1;
    lp.PTTRF(nu_,&d[0],&o[0],&info);
    lp.PTTRS(nu_,nhrs,&d[0],&o[0],&u[0],ldb,&info);
  }

  void solve_adjoint_equation(std::vector<Real> &p, const std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    vector d(nu_,1.0);
    vector o(nu_-1,1.0);
    for ( uint i = 0; i < nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/hu_;
      if ( i < nu_-1 ) {
        o[i] *= -z[i+1]/hu_;
      }
    }
    // Set right hand side
    vector r(nu_,0.0);
    for (uint i = 0; i < nu_; i++) {
      r[i] = -(u[i]-evaluate_target((Real)(i+1)*hu_));
    }
    p.clear();
    p.resize(nu_,0.0);
    apply_mass(p,r);    
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nu_;
    int nhrs = 1;
    lp.PTTRF(nu_,&d[0],&o[0],&info);
    lp.PTTRS(nu_,nhrs,&d[0],&o[0],&p[0],ldb,&info);
  }

  void solve_state_sensitivity_equation(std::vector<Real> &w, const std::vector<Real> &v, 
                                  const std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    vector d(nu_,1.0);
    vector o(nu_-1,1.0);
    for ( uint i = 0; i < nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/hu_;
      if ( i < nu_-1 ) {
        o[i] *= -z[i+1]/hu_;
      }
    }
    // Set right hand side
    w.clear();
    w.resize(nu_,0.0);
    apply_linearized_control_operator(w,z,v,u);
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nu_;
    int nhrs = 1;
    lp.PTTRF(nu_,&d[0],&o[0],&info);
    lp.PTTRS(nu_,nhrs,&d[0],&o[0],&w[0],ldb,&info);
  }

  void solve_adjoint_sensitivity_equation(std::vector<Real> &q, const std::vector<Real> &w, 
                                    const std::vector<Real> &v, const std::vector<Real> &p, 
                                    const std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    vector d(nu_,1.0);
    vector o(nu_-1,1.0);
    for ( uint i = 0; i < nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/hu_;
      if ( i < nu_-1 ) {
        o[i] *= -z[i+1]/hu_;
      }
    }
    // Set right hand side
    q.clear();
    q.resize(nu_,0.0);
    apply_mass(q,w);
    std::vector<Real> res(nu_,0.0);
    apply_linearized_control_operator(res,z,v,p,false);
    for (uint i = 0; i < nu_; i++) {
      q[i] -= res[i];
    }
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = nu_;
    int nhrs = 1;
    lp.PTTRF(nu_,&d[0],&o[0],&info);
    lp.PTTRS(nu_,nhrs,&d[0],&o[0],&q[0],ldb,&info);
  }

  void update(const ROL::Vector<Real> &z, bool flag, int iter) {

    

    if ( flag && useCorrection_ ) {
      Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
      H_.shape(nz_,nz_); 
      ROL::Ptr<V> e = z.clone();
      ROL::Ptr<V> h = z.clone();
      for ( uint i = 0; i < nz_; i++ ) {
        e = z.basis(i);
        hessVec_true(*h,*e,z,tol);
        for ( uint j = 0; j < nz_; j++ ) {
          e = z.basis(j);
          (H_)(j,i) = e->dot(*h);
        }
      }
      std::vector<vector> eigenvals = ROL::computeEigenvalues<Real>(H_);
      std::sort((eigenvals[0]).begin(), (eigenvals[0]).end());
      Real inertia = (eigenvals[0])[0];
      Real correction = 0.0;
      if ( inertia <= 0.0 ) {
        correction = (1.0+std::sqrt(ROL::ROL_EPSILON<Real>()))*std::abs(inertia);
        if ( inertia == 0.0 ) {
          uint cnt = 0;
          while ( eigenvals[0][cnt] == 0.0 ) {
            cnt++;
          }
          correction = std::sqrt(ROL::ROL_EPSILON<Real>())*eigenvals[0][cnt];
          if ( cnt == nz_-1 ) {
            correction = 1.0;
          }
        }
        for ( uint i = 0; i < nz_; i++ ) {
          (H_)(i,i) += correction;
        }
      }  
    }
  }

  /* OBJECTIVE FUNCTION DEFINITIONS */
  Real value( const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<const vector> zp = getVector(z); 

    // SOLVE STATE EQUATION
    vector u(nu_,0.0);
    solve_state_equation(u,*zp);

    // EVALUATE OBJECTIVE
    Real val  = 0.0;
    for (uint i=0; i<nz_;i++) {
      val += hz_*alpha_*0.5*(*zp)[i]*(*zp)[i];
    }
    Real res  = 0.0;
    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (uint i=0; i<nu_; i++) {
      if ( i == 0 ) {
        res1 = u[i]-evaluate_target((Real)(i+1)*hu_);
        res2 = u[i+1]-evaluate_target((Real)(i+2)*hu_);
        res  = hu_/6.0*(2.0*res1 + res2)*res1;
      }
      else if ( i == nu_-1 ) {
        res1 = u[i-1]-evaluate_target((Real)i*hu_);
        res2 = u[i]-evaluate_target((Real)(i+1)*hu_);
        res  = hu_/6.0*(res1 + 2.0*res2)*res2;
      }
      else {
        res1 = u[i-1]-evaluate_target((Real)i*hu_);
        res2 = u[i]-evaluate_target((Real)(i+1)*hu_);
        res3 = u[i+1]-evaluate_target((Real)(i+2)*hu_);
        res  = hu_/6.0*(res1 + 4.0*res2 + res3)*res2;
      }
      val += 0.5*res;
    }
    return val;
  } 

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
 
    

    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<vector> gp = getVector(g);

    // SOLVE STATE EQUATION
    vector u(nu_,0.0);
    solve_state_equation(u,*zp);

    // SOLVE ADJOINT EQUATION
    vector p(nu_,0.0);
    solve_adjoint_equation(p,u,*zp);

    // Apply Transpose of Linearized Control Operator
    apply_transposed_linearized_control_operator(*gp,*zp,p,u);
    // Build Gradient
    for ( uint i = 0; i < nz_; i++ ) {
      (*gp)[i] += hz_*alpha_*(*zp)[i];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    if ( useCorrection_ ) {
      hessVec_inertia(hv,v,z,tol);
    }
    else {
      hessVec_true(hv,v,z,tol);
    }
  }

  void activateInertia(void) {
    useCorrection_ = true;
  }

  void deactivateInertia(void) {
    useCorrection_ = false;
  }

  void hessVec_true( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {

    

    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<vector> hvp = getVector(hv);

    // SOLVE STATE EQUATION
    vector u(nu_,0.0);
    solve_state_equation(u,*zp);

    // SOLVE ADJOINT EQUATION
    vector p(nu_,0.0);
    solve_adjoint_equation(p,u,*zp);

    // SOLVE STATE SENSITIVITY EQUATION
    vector w(nu_,0.0);
    solve_state_sensitivity_equation(w,*vp,u,*zp);
    // SOLVE ADJOINT SENSITIVITY EQUATION
    vector q(nu_,0.0);
    solve_adjoint_sensitivity_equation(q,w,*vp,p,u,*zp);

    // Apply Transpose of Linearized Control Operator
    apply_transposed_linearized_control_operator(*hvp,*zp,q,u);

    // Apply Transpose of Linearized Control Operator
    std::vector<Real> tmp(nz_,0.0);
    apply_transposed_linearized_control_operator(tmp,*zp,w,p,false);
    for (uint i=0; i < nz_; i++) {
      (*hvp)[i] -= tmp[i];
    }
    // Regularization hessVec
    for (uint i=0; i < nz_; i++) {
      (*hvp)[i] += hz_*alpha_*(*vp)[i];
    }
  }

  void hessVec_inertia( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {

     
    using ROL::constPtrCast;

    ROL::Ptr<vector> hvp = getVector(hv);

    
    ROL::Ptr<vector> vp  = ROL::constPtrCast<vector>(getVector(v));

    Teuchos::SerialDenseVector<int, Real> hv_teuchos(Teuchos::View, &((*hvp)[0]), static_cast<int>(nz_));
    Teuchos::SerialDenseVector<int, Real>  v_teuchos(Teuchos::View, &(( *vp)[0]), static_cast<int>(nz_));
    hv_teuchos.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, H_, v_teuchos, 0.0);
  }

};



typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>     vector;
  typedef ROL::Vector<RealT>     V;
  typedef ROL::StdVector<RealT>  SV;
  
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

    uint dim = 128; // Set problem dimension.
    RealT alpha = 1.e-6;
    Objective_PoissonInversion<RealT> obj(dim, alpha);

    // Iteration vector.
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim, 0.0);
    ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(dim, 0.0);

    // Set initial guess.
    for (uint i=0; i<dim; i++) {
      (*x_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX + 1.e2;
      (*y_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX + 1.e2;
    }

    SV x(x_ptr);
    SV y(y_ptr);

    obj.checkGradient(x,y,true);
    obj.checkHessVec(x,y,true);

    ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(dim,1.0);
    ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(dim,10.0);

    ROL::Ptr<V> lo = ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(u_ptr);

    ROL::Bounds<RealT> icon(lo,up);

    ROL::ParameterList parlist;

    // Krylov parameters.
    parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-8);
    parlist.sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-4);
    parlist.sublist("General").sublist("Krylov").set("Iteration Limit",static_cast<int>(dim));
    // PDAS parameters.
    parlist.sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-8);
    parlist.sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-6);
    parlist.sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 10);
    parlist.sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);
    parlist.sublist("General").sublist("Secant").set("Use as Hessian",true);
    // Status test parameters.
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);

    ROL::Ptr<ROL::Step<RealT>>       step;
    ROL::Ptr<ROL::StatusTest<RealT>> status;

    // Define algorithm.
    step = ROL::makePtr<ROL::PrimalDualActiveSetStep<RealT>>(parlist);
    status = ROL::makePtr<ROL::StatusTest<RealT>>(parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    x.zero();
    obj.deactivateInertia();
    algo.run(x,obj,icon,true,*outStream);

    // Output control to file.
    std::ofstream file;
    file.open("control_PDAS.txt");
    for ( uint i = 0; i < dim; i++ ) {
      file << (*x_ptr)[i] << "\n";
    }
    file.close();

    // Projected Newtion.
    // Define step.
    parlist.sublist("General").sublist("Secant").set("Use as Hessian",false);
    parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver", "Truncated CG");
    parlist.sublist("Step").sublist("Trust Region").set("Initial Radius", 1e3);
    parlist.sublist("Step").sublist("Trust Region").set("Maximum Radius", 1e8);
    step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(parlist);
    status = ROL::makePtr<ROL::StatusTest<RealT>>(parlist);
    ROL::Algorithm<RealT> algo_tr(step,status,false);
    // Run Algorithm
    y.zero();
    obj.deactivateInertia();
    algo_tr.run(y,obj,icon,true,*outStream);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( uint i = 0; i < dim; i++ ) {
      file_tr << (*y_ptr)[i] << "\n";
    }
    file_tr.close();

    std::vector<RealT> u;
    obj.solve_state_equation(u,*y_ptr);
    std::ofstream file_u;
    file_u.open("state.txt");
    for ( uint i = 0; i < (dim-1); i++ ) {
      file_u << u[i] << "\n";
    }
    file_u.close();
   
    ROL::Ptr<V> diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm()/std::sqrt((RealT)dim-1.0);
    std::cout << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1.e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);

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

