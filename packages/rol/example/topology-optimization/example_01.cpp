// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve a topology optimization problem using Moreau-Yoshida 
           regularization for the volume constraint.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_LAPACK.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_ParameterList.hpp"

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

template<class Real>
class FEM {

  typedef typename std::vector<Real>::size_type uint;

private:
  uint nx_;
  uint ny_;
  int p_;
  int prob_;
  Real E_;
  Real nu_;
  Teuchos::SerialDenseMatrix<int,Real> KE_;

public:

  FEM(uint nx, uint ny, int p = 3, int prob = 1, Real E = 1.0, Real nu = 0.3) 
  : nx_(nx), ny_(ny), p_(p), prob_(prob), E_(E), nu_(nu) {
    std::vector<Real> k(8,0.0);
    k[0] =  1.0/2.0-nu_/6.0; 
    k[1] =  1.0/8.0+nu_/8.0;
    k[2] = -1.0/4.0-nu_/12.0;
    k[3] = -1.0/8.0+3.0*nu_/8.0;
    k[4] = -1.0/4.0+nu_/12.0;
    k[5] = -1.0/8.0-nu_/8.0;
    k[6] = nu_/6.0;
    k[7] =  1.0/8.0-3.0*nu_/8.0;
    KE_.shape(8,8);
    // Row 1
    KE_(0,0) = E_/(1.0-nu_*nu_)*k[0];
    KE_(0,1) = E_/(1.0-nu_*nu_)*k[1];
    KE_(0,2) = E_/(1.0-nu_*nu_)*k[2];
    KE_(0,3) = E_/(1.0-nu_*nu_)*k[3];
    KE_(0,4) = E_/(1.0-nu_*nu_)*k[4];
    KE_(0,5) = E_/(1.0-nu_*nu_)*k[5];
    KE_(0,6) = E_/(1.0-nu_*nu_)*k[6];
    KE_(0,7) = E_/(1.0-nu_*nu_)*k[7];
    // Row 2
    KE_(1,0) = E_/(1.0-nu_*nu_)*k[1];
    KE_(1,1) = E_/(1.0-nu_*nu_)*k[0];
    KE_(1,2) = E_/(1.0-nu_*nu_)*k[7];
    KE_(1,3) = E_/(1.0-nu_*nu_)*k[6];
    KE_(1,4) = E_/(1.0-nu_*nu_)*k[5];
    KE_(1,5) = E_/(1.0-nu_*nu_)*k[4];
    KE_(1,6) = E_/(1.0-nu_*nu_)*k[3];
    KE_(1,7) = E_/(1.0-nu_*nu_)*k[2];
    // Row 3
    KE_(2,0) = E_/(1.0-nu_*nu_)*k[2];
    KE_(2,1) = E_/(1.0-nu_*nu_)*k[7];
    KE_(2,2) = E_/(1.0-nu_*nu_)*k[0];
    KE_(2,3) = E_/(1.0-nu_*nu_)*k[5];
    KE_(2,4) = E_/(1.0-nu_*nu_)*k[6];
    KE_(2,5) = E_/(1.0-nu_*nu_)*k[3];
    KE_(2,6) = E_/(1.0-nu_*nu_)*k[4];
    KE_(2,7) = E_/(1.0-nu_*nu_)*k[1];
    // Row 4
    KE_(3,0) = E_/(1.0-nu_*nu_)*k[3];
    KE_(3,1) = E_/(1.0-nu_*nu_)*k[6];
    KE_(3,2) = E_/(1.0-nu_*nu_)*k[5];
    KE_(3,3) = E_/(1.0-nu_*nu_)*k[0];
    KE_(3,4) = E_/(1.0-nu_*nu_)*k[7];
    KE_(3,5) = E_/(1.0-nu_*nu_)*k[2];
    KE_(3,6) = E_/(1.0-nu_*nu_)*k[1];
    KE_(3,7) = E_/(1.0-nu_*nu_)*k[4];
    // Row 5
    KE_(4,0) = E_/(1.0-nu_*nu_)*k[4];
    KE_(4,1) = E_/(1.0-nu_*nu_)*k[5];
    KE_(4,2) = E_/(1.0-nu_*nu_)*k[6];
    KE_(4,3) = E_/(1.0-nu_*nu_)*k[7];
    KE_(4,4) = E_/(1.0-nu_*nu_)*k[0];
    KE_(4,5) = E_/(1.0-nu_*nu_)*k[1];
    KE_(4,6) = E_/(1.0-nu_*nu_)*k[2];
    KE_(4,7) = E_/(1.0-nu_*nu_)*k[3];
    // Row 6
    KE_(5,0) = E_/(1.0-nu_*nu_)*k[5];
    KE_(5,1) = E_/(1.0-nu_*nu_)*k[4];
    KE_(5,2) = E_/(1.0-nu_*nu_)*k[3];
    KE_(5,3) = E_/(1.0-nu_*nu_)*k[2];
    KE_(5,4) = E_/(1.0-nu_*nu_)*k[1];
    KE_(5,5) = E_/(1.0-nu_*nu_)*k[0];
    KE_(5,6) = E_/(1.0-nu_*nu_)*k[7];
    KE_(5,7) = E_/(1.0-nu_*nu_)*k[6];
    // Row 7
    KE_(6,0) = E_/(1.0-nu_*nu_)*k[6];
    KE_(6,1) = E_/(1.0-nu_*nu_)*k[3];
    KE_(6,2) = E_/(1.0-nu_*nu_)*k[4];
    KE_(6,3) = E_/(1.0-nu_*nu_)*k[1];
    KE_(6,4) = E_/(1.0-nu_*nu_)*k[2];
    KE_(6,5) = E_/(1.0-nu_*nu_)*k[7];
    KE_(6,6) = E_/(1.0-nu_*nu_)*k[0];
    KE_(6,7) = E_/(1.0-nu_*nu_)*k[5];
    // Row 8
    KE_(7,0) = E_/(1.0-nu_*nu_)*k[7];
    KE_(7,1) = E_/(1.0-nu_*nu_)*k[2];
    KE_(7,2) = E_/(1.0-nu_*nu_)*k[1];
    KE_(7,3) = E_/(1.0-nu_*nu_)*k[4];
    KE_(7,4) = E_/(1.0-nu_*nu_)*k[3];
    KE_(7,5) = E_/(1.0-nu_*nu_)*k[6];
    KE_(7,6) = E_/(1.0-nu_*nu_)*k[5];
    KE_(7,7) = E_/(1.0-nu_*nu_)*k[0];
  }

  uint numX(void) { return nx_; }
  uint numY(void) { return ny_; }
  uint numZ(void) { return nx_*ny_; }
  uint numU(void) { return 2*(nx_+1)*(ny_+1); }
  int power(void) { return p_; }

  int get_index(int r, int n1, int n2) {
    int ind = 0;
    switch(r) {
      case 0: ind = 2*n1-2; break;
      case 1: ind = 2*n1-1; break;
      case 2: ind = 2*n2-2; break;
      case 3: ind = 2*n2-1; break;
      case 4: ind = 2*n2;   break;
      case 5: ind = 2*n2+1; break;
      case 6: ind = 2*n1;   break;
      case 7: ind = 2*n1+1; break;
    }
    return ind;
  }

  bool check_on_boundary(uint r) {
    switch(prob_) {
      case 0: 
        if ( (r < 2*(ny_+1) && r%2 == 0) || ( r == 2*(nx_+1)*(ny_+1)-1 ) ) {
          return true;
        }
        break;
      case 1: 
        if ( r < 2*(ny_+1) ) {
          return true;
        }
        break;
    }
    return false;
  }

  void set_boundary_conditions(Teuchos::SerialDenseVector<int, Real> &U) {
    for (int i=0; i<U.length(); i++) {
      if ( check_on_boundary(i) ) {
        U(i) = 0.0;
      }
    }
  }

  void set_boundary_conditions(std::vector<Real> &u) {
    for (uint i=0; i<u.size(); i++) {
      if ( check_on_boundary(i) ) {
        u[i] = 0.0;
      }
    }
  }

  void build_force(std::vector<Real> &F) {
    F.assign(numU(),0.0);
    switch(prob_) {
      case 0: F[1] = -1;                   break;
      case 1: F[2*(nx_+1)*(ny_+1)-1] = -1; break;
    }
  }

  void build_force(Teuchos::SerialDenseVector<int, Real> &F) {
    F.resize(numU());
    F.putScalar(0.0);
    switch(prob_) {
      case 0: F(1) = -1;                   break;
      case 1: F(2*(nx_+1)*(ny_+1)-1) = -1; break;
    }
  }

  void build_jacobian(Teuchos::SerialDenseMatrix<int, Real> &K, const std::vector<Real> &z,
                      bool transpose = false) {
    // Fill jacobian
    K.shape(2*(nx_+1)*(ny_+1),2*(nx_+1)*(ny_+1));
    int n1 = 0, n2 = 0, row = 0, col = 0;
    Real Zp = 0.0, val = 0.0;
    for (uint i=0; i<nx_; i++) {
      for (uint j=0; j<ny_; j++) {
        n1 = (ny_+1)* i   +(j+1);
        n2 = (ny_+1)*(i+1)+(j+1); 
        Zp = std::pow(z[i+j*nx_],(Real)p_);
        for (uint r=0; r<8; r++) {
          row = get_index(r,n1,n2); 
          if ( check_on_boundary(row) ) {
            K(row,row) = 1.0;
          }
          else {
            for (uint c=0; c<8; c++) {
              col = get_index(c,n1,n2);
              if ( !check_on_boundary(col) ) {
                val = Zp*(KE_)(r,c);
                if (transpose) {
                  K(col,row) += val;
                }
                else {
                  K(row,col) += val;
                }
              }
            }
          }
        }
      }
    }
  }

  void build_jacobian(Teuchos::SerialDenseMatrix<int, Real> &K, const std::vector<Real> &z, 
                      const std::vector<Real> &v, bool transpose = false) {
    // Fill jacobian
    K.shape(2*(nx_+1)*(ny_+1),2*(nx_+1)*(ny_+1));
    uint n1 = 0, n2 = 0, row = 0, col = 0;
    Real Zp = 0.0, V = 0.0, val = 0.0;
    for (uint i=0; i<nx_; i++) {
      for (uint j=0; j<ny_; j++) {
        n1 = (ny_+1)* i   +(j+1);
        n2 = (ny_+1)*(i+1)+(j+1); 
        Zp = (p_ == 1) ? 1.0 : (Real)p_*std::pow(z[i+j*nx_],(Real)p_-1.0);
        V  = v[i+j*nx_];
        for (uint r=0; r<8; r++) {
          row = get_index(r,n1,n2); 
          if ( check_on_boundary(row) ) {
            K(row,row) = 1.0;
          }
          else {
            for (uint c=0; c<8; c++) {
              col = get_index(c,n1,n2);
              if ( !check_on_boundary(col) ) {
                val = Zp*V*(KE_)(r,c);
                if (transpose) {
                  K(col,row) += val;
                }
                else {
                  K(row,col) += val;
                }
              }
            }
          }
        }
      }
    }
  }

  void apply_jacobian(std::vector<Real> &ju, const std::vector<Real> &u, const std::vector<Real> &z) {
    uint n1 = 0, n2 = 0, row = 0, col = 0;
    Real Zp = 0.0;
    ju.assign(u.size(),0.0);
    for (uint i=0; i<nx_; i++) {
      for (uint j=0; j<ny_; j++) {
        n1 = (ny_+1)* i   +(j+1);
        n2 = (ny_+1)*(i+1)+(j+1); 
        Zp = std::pow(z[i+j*nx_],(Real)p_);
        for (uint r=0; r<8; r++) {
          row = get_index(r,n1,n2); 
          if ( check_on_boundary(row) ) {
            ju[row] = u[row];
          }
          else {
            for (uint c=0; c<8; c++) {
              col = get_index(c,n1,n2);
              if ( !check_on_boundary(col) ) {
                ju[row] += Zp*(KE_)(r,c)*u[col];
              }
            }
          }
        }
      }
    }
  } 

  void apply_jacobian(std::vector<Real> &ju, const std::vector<Real> &u, const std::vector<Real> &z, 
                      const std::vector<Real> &v) {
    // Fill jacobian
    uint n1 = 0, n2 = 0, row = 0, col = 0;
    Real Zp = 0.0, V = 0.0;
    ju.assign(u.size(),0.0);
    for (uint i=0; i<nx_; i++) {
      for (uint j=0; j<ny_; j++) {
        n1 = (ny_+1)* i   +(j+1);
        n2 = (ny_+1)*(i+1)+(j+1); 
        Zp = (p_ == 1) ? 1.0 : (Real)p_*std::pow(z[i+j*nx_],(Real)p_-1.0);
        V  = v[i+j*nx_];
        for (uint r=0; r<8; r++) {
          row = get_index(r,n1,n2); 
          if ( check_on_boundary(row) ) {
            ju[row] = u[row];
          }
          else {
            for (uint c=0; c<8; c++) {
              col = get_index(c,n1,n2);
              if ( !check_on_boundary(col) ) {
                ju[row] += V*Zp*(KE_)(r,c)*u[col];
              }
            }
          }
        }
      }
    }
  } 

  void apply_adjoint_jacobian(std::vector<Real> &jv, const std::vector<Real> &u, const std::vector<Real> &z, 
                              const std::vector<Real> &v) {
    // Fill jacobian
    uint n1 = 0, n2 = 0, row = 0, col = 0;
    Real Zp = 0.0, VKU = 0.0;
    for (uint i=0; i<nx_; i++) {
      for (uint j=0; j<ny_; j++) {
        n1 = (ny_+1)* i   +(j+1);
        n2 = (ny_+1)*(i+1)+(j+1); 
        Zp = (p_ == 1) ? 1.0 : (Real)p_*std::pow(z[i+j*nx_],(Real)p_-1.0);
        VKU = 0.0;
        for (uint r=0; r<8; r++) {
          row = get_index(r,n1,n2); 
          if ( check_on_boundary(row) ) {
            VKU += v[row]*u[row];
          }
          else {
            for (uint c=0; c<8; c++) {
              col = get_index(c,n1,n2);
              if ( !check_on_boundary(col) ) {
                VKU += Zp*(KE_)(r,c)*u[col]*v[row];
              }
            }
          }
        }
        jv[i+j*nx_] = VKU;
      }
    }
  }

  void apply_adjoint_jacobian(std::vector<Real> &jv, const std::vector<Real> &u, const std::vector<Real> &z, 
                              const std::vector<Real> &v, const std::vector<Real> &w) {
    // Fill jacobian
    uint n1 = 0, n2 = 0, row = 0, col = 0;
    Real Zp = 0.0, V = 0.0, VKU = 0.0;
    for (uint i=0; i<nx_; i++) {
      for (uint j=0; j<ny_; j++) {
        n1 = (ny_+1)* i   +(j+1);
        n2 = (ny_+1)*(i+1)+(j+1); 
        Zp = (p_ == 1) ? 0.0 : 
               ((p_ == 2) ? 2.0 : 
                 (Real)p_*((Real)p_-1.0)*std::pow(z[i+j*nx_],(Real)p_-2.0));
        V  = v[i+j*nx_];
        VKU = 0.0;
        for (uint r=0; r<8; r++) {
          row = get_index(r,n1,n2); 
          if ( check_on_boundary(row) ) {
            VKU += w[row]*u[row];
          }
          else {
            for (uint c=0; c<8; c++) {
              col = get_index(c,n1,n2);
              if ( !check_on_boundary(col) ) {
                VKU += Zp*V*(KE_)(r,c)*u[col]*w[row];
              }
            }
          }
        }
        jv[i+j*nx_] = VKU;
      }
    }
  }
};

template<class Real>
class Constraint_TopOpt : public ROL::Constraint_SimOpt<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
  
  typedef typename vector::size_type uint;
 
private:
  ROL::Ptr<FEM<Real> > FEM_;

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }
  
  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  Constraint_TopOpt(ROL::Ptr<FEM<Real> > & FEM) : FEM_(FEM) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {

    

    ROL::Ptr<vector> cp = getVector(c);
    applyJacobian_1(c, u, u, z, tol);
    vector f;
    FEM_->build_force(f);
    for (uint i = 0; i < f.size(); i++) {
      (*cp)[i] -= f[i];
    }
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {

    

    ROL::Ptr<vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Assemble Jacobian
    Teuchos::SerialDenseMatrix<int, Real> K;
    FEM_->build_jacobian(K,*zp);
    // Assemble RHS
    Teuchos::SerialDenseVector<int, Real> F(K.numRows());
    FEM_->build_force(F);
    // Solve
    Teuchos::SerialDenseVector<int, Real> U(K.numCols());
    Teuchos::SerialDenseSolver<int, Real> solver;
    solver.setMatrix( Teuchos::rcpFromRef(K) );
    solver.setVectors( Teuchos::rcpFromRef(U), Teuchos::rcpFromRef(F) );
    solver.factorWithEquilibration(true);
    solver.factor();
    solver.solve();
    FEM_->set_boundary_conditions(U);
    // Retrieve solution
    up->resize(U.length(),0.0);
    for (uint i=0; i<static_cast<uint>(U.length()); i++) {
      (*up)[i] = U(i);
    }
    // Compute residual
    this->value(c,u,z,tol);
    //std::cout << " IN SOLVE: ||c(u,z)|| = " << c->norm() << "\n";
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {

    

    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    vector V;
    V.assign(vp->begin(),vp->end());
    FEM_->set_boundary_conditions(V);
    FEM_->apply_jacobian(*jvp,V,*zp);
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {

    
    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z); 

    // Apply Jacobian
    vector U;
    U.assign(up->begin(),up->end());
    FEM_->set_boundary_conditions(U);
    FEM_->apply_jacobian(*jvp,U,*zp,*vp);
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {

    
    ROL::Ptr<vector> ijvp = getVector(ijv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Assemble Jacobian
    Teuchos::SerialDenseMatrix<int, Real> K;
    FEM_->build_jacobian(K,*zp);
    // Solve
    Teuchos::SerialDenseVector<int, Real> U(K.numCols());
    Teuchos::SerialDenseVector<int, Real> F(vp->size());
    for (uint i=0; i<vp->size(); i++) {
      F(i) = (*vp)[i];
    }
    Teuchos::SerialDenseSolver<int, Real> solver;
    solver.setMatrix(Teuchos::rcpFromRef(K));
    solver.setVectors(Teuchos::rcpFromRef(U),Teuchos::rcpFromRef(F));
    solver.factorWithEquilibration(true);
    solver.factor();
    solver.solve();
    FEM_->set_boundary_conditions(U);
    // Retrieve solution
    ijvp->resize(U.length(),0.0);
    for (uint i=0; i<static_cast<uint>(U.length()); i++) {
      (*ijvp)[i] = U(i);
    }
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    
 
    ROL::Ptr<vector> ajvp = getVector(ajv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // apply jacobian
    vector V;
    V.assign(vp->begin(),vp->end());
    FEM_->set_boundary_conditions(V);
    FEM_->apply_jacobian(*ajvp,V,*zp);
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {

    
    ROL::Ptr<vector> ajvp = getVector(ajv);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    vector U;
    U.assign(up->begin(),up->end());
    FEM_->set_boundary_conditions(U);
    std::vector<Real> V;
    V.assign(vp->begin(),vp->end());
    FEM_->set_boundary_conditions(V);
    FEM_->apply_adjoint_jacobian(*ajvp,U,*zp,V);
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
     
      
     ROL::Ptr<vector> iajvp = getVector(iajv);
     ROL::Ptr<const vector> vp = getVector(v);
     ROL::Ptr<const vector> up = getVector(u);
     ROL::Ptr<const vector> zp = getVector(z);

    // Assemble Jacobian
    Teuchos::SerialDenseMatrix<int, Real> K;
    FEM_->build_jacobian(K,*zp);
    // Solve
    Teuchos::SerialDenseVector<int, Real> U(K.numCols());
    Teuchos::SerialDenseVector<int, Real> F(vp->size());
    for (uint i=0; i<vp->size(); i++) {
      F(i) = (*vp)[i];
    }
    Teuchos::SerialDenseSolver<int, Real> solver;
    solver.setMatrix(Teuchos::rcpFromRef(K));
    solver.setVectors(Teuchos::rcpFromRef(U), Teuchos::rcpFromRef(F));
    solver.factorWithEquilibration(true);
    solver.factor();
    solver.solve();
    FEM_->set_boundary_conditions(U);
    // Retrieve solution
    iajvp->resize(U.length(),0.0);
    for (int i=0; i<U.length(); i++) {
      (*iajvp)[i] = U(i);
    }
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    this->applyAdjointJacobian_2(ahwv,w,v,z,tol);
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    this->applyJacobian_2(ahwv,v,w,z,tol);
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {

    
    ROL::Ptr<vector> ahwvp = getVector(ahwv);
    ROL::Ptr<const vector> wp = getVector(w);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    vector U;
    U.assign(up->begin(),up->end());
    FEM_->set_boundary_conditions(U);
    vector W;
    W.assign(wp->begin(),wp->end());
    FEM_->set_boundary_conditions(W);
    FEM_->apply_adjoint_jacobian(*ahwvp,U,*zp,*vp,W);
  }

  //void solveAugmentedSystem(ROL::Vector<Real> &v1, ROL::Vector<Real> &v2, const ROL::Vector<Real> &b1,
  //                          const ROL::Vector<Real> &b2, const ROL::Vector<Real> &x, Real &tol) {}
};

template<class Real>
class Objective_TopOpt : public ROL::Objective_SimOpt<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
  
  typedef typename vector::size_type uint;  

private:
  ROL::Ptr<FEM<Real> > FEM_;
  Real frac_; 
  Real reg_; 
  Real pen_;
  Real rmin_;

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  Objective_TopOpt(ROL::Ptr<FEM<Real> > FEM, 
    Real frac = 0.5, Real reg = 1.0, Real pen = 1.0, Real rmin = -1.0 )
    : FEM_(FEM), frac_(frac), reg_(reg), pen_(pen), rmin_(rmin), useLC_(true) {}

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
 
    // Apply Jacobian
    vector KU(up->size(),0.0);
    if ( useLC_ ) {
      FEM_->build_force(KU);
    }
    else {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp);
    }
    // Compliance
    Real c = 0.0;
    for (uint i=0; i<up->size(); i++) {
      c += (*up)[i]*KU[i];
    }
    // Compute Moreau-Yoshida term
    Real vol = 0.0, r = 0.0;
    for (uint i=0; i<zp->size(); i++) {
      vol += (*zp)[i]; 
    }
    vol  = (vol <= frac_*FEM_->numZ()) ? 0.0 : (vol - frac_*FEM_->numZ());
    r = (reg_)*std::pow(vol,3.0);
    // Compute 0-1 penalty
    Real val = 0.0, p = 0.0;
    if ( rmin_ <= 0.0 ) {
      for (uint i=0; i<zp->size(); i++) {
        val += (*zp)[i]*(1.0-(*zp)[i]);
      }
    }
    else {
      Real sum = 0.0, area = 0.0;
      uint i1 = 0, i2 = 0, j1 = 0, j2 = 0;
      for (uint i=0; i<FEM_->numX(); i++) {
        for (uint j=0; j<FEM_->numY(); j++) {
          sum = 0.0;
          i1 = std::max(i-(uint)floor(rmin_),(uint)1)-1;
          i2 = std::min(i+(uint)floor(rmin_),FEM_->numX());
          j1 = std::max(j-(uint)floor(rmin_),(uint)1)-1;
          j2 = std::min(j+(uint)floor(rmin_),FEM_->numY());
          area = (Real)(i2-i1)*(j2-j1); 
          for (uint ii=i1; ii<i2; ii++) {
            for (uint jj=j1; jj<j2; jj++) {
              sum += (*zp)[ii+FEM_->numX()*jj]/area;
            }
          }
          val += sum*(1.0-sum);
        }
      }
    }
    p = (pen_)*val;
    return c + r + p;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    

    // Unwrap g
    ROL::Ptr<vector> gp = getVector(g);
    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    vector KU(up->size(),0.0);
    if ( useLC_ ) {
      FEM_->build_force(KU);
      // Apply jacobian to u
      for (uint i=0; i<up->size(); i++) {
        (*gp)[i] = KU[i];
      }
    }
    else {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp);
      // Apply jacobian to u
      for (uint i=0; i<up->size(); i++) {
        (*gp)[i] = 2.0*KU[i];
      }
    }
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    

    // Unwrap g
    ROL::Ptr<vector> gp = getVector(g);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    g.zero();
    if ( !useLC_ ) {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_adjoint_jacobian(*gp,U,*zp,U);
    }
    // Compute Moreau-Yoshida term
    Real vol = 0.0;
    for (uint i=0; i<zp->size(); i++) {
      vol += (*zp)[i]; 
    }
    vol = (vol <= frac_*FEM_->numZ()) ? 0.0 : (vol - frac_*FEM_->numZ());
    for (uint i=0; i<zp->size(); i++) {
      (*gp)[i] += (reg_)*3.0*std::pow(vol,2.0);
      // Compute 0-1 penalty
      if ( rmin_ <= 0.0 ) {
        (*gp)[i] += (pen_)*(1.0-2.0*(*zp)[i]);
      }
    }
    // Compute 0-1 penalty
    vector Tz(zp->size(),0.0);
    Real area = 0.0;
    uint i1 = 0, i2 = 0, j1 = 0, j2 = 0;
    if (rmin_ > 0.0) {
      for (uint i=0; i<FEM_->numX(); i++) {
        for (uint j=0; j<FEM_->numY(); j++) {
          i1 = std::max(i-(uint)floor(rmin_),(uint)1)-1;
          i2 = std::min(i+(uint)floor(rmin_),FEM_->numX());
          j1 = std::max(j-(uint)floor(rmin_),(uint)1)-1;
          j2 = std::min(j+(uint)floor(rmin_),FEM_->numY());
          area = (Real)(i2-i1)*(j2-j1);
          for (uint ii=i1; ii<i2; ii++) {
            for (uint jj=j1; jj<j2; jj++) {
              Tz[i+j*FEM_->numX()] += (*zp)[ii+FEM_->numX()*jj]/area;
            }
          }
        }
      }
      for (uint i=0; i<FEM_->numX(); i++) {
        for (uint j=0; j<FEM_->numY(); j++) {
          i1 = std::max(i-(uint)floor(rmin_),(uint)1)-1;
          i2 = std::min(i+(uint)floor(rmin_),FEM_->numX());
          j1 = std::max(j-(uint)floor(rmin_),(uint)1)-1;
          j2 = std::min(j+(uint)floor(rmin_),FEM_->numY());
          area = (Real)(i2-i1)*(j2-j1);
          for (uint ii=i1; ii<i2; ii++) {
            for (uint jj=j1; jj<j2; jj++) {
              (*gp)[ii+jj*FEM_->numX()] += (pen_)*(1.0-2.0*Tz[i+FEM_->numX()*j])/area;
            }
          }
        }
      }
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      vector KV(vp->size(),0.0);
      vector V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_jacobian(KV,V,*zp);
      for (uint i=0; i<vp->size(); i++) {
        (*hvp)[i] = 2.0*KV[i];
      }
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    

    // Unwrap hv
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);

    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      vector KU(up->size(),0.0);
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      FEM_->apply_jacobian(KU,U,*zp,*vp);
      for (uint i=0; i<up->size(); i++) {
        (*hvp)[i] = 2.0*KU[i];
      }
    }
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    

    // Unwrap g
    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
 
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      std::vector<Real> U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      std::vector<Real> V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,V);
      for (uint i=0; i<hvp->size(); i++) {
        (*hvp)[i] *= 2.0;
      }
    }
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    

    ROL::Ptr<vector> hvp = getVector(hv);

    // Unwrap v
    ROL::Ptr<const vector> vp = getVector(v);

    // Unwrap x
    ROL::Ptr<const vector> up = getVector(u);
    ROL::Ptr<const vector> zp = getVector(z);
    
    // Apply Jacobian
    hv.zero();
    if ( !useLC_ ) {
      vector U;
      U.assign(up->begin(),up->end());
      FEM_->set_boundary_conditions(U);
      vector V;
      V.assign(vp->begin(),vp->end());
      FEM_->set_boundary_conditions(V);
      FEM_->apply_adjoint_jacobian(*hvp,U,*zp,*vp,U);
    }
    // Compute Moreau-Yoshida term
    Real vol = 0.0, vvol = 0.0;
    for (uint i=0; i<zp->size(); i++) {
      vol  += (*zp)[i]; 
      vvol += (*vp)[i];
    }
    vol  = (vol <= frac_*FEM_->numZ()) ? 0.0 : (vol - frac_*FEM_->numZ());
    for (uint i=0; i<zp->size(); i++) {
      (*hvp)[i] += (reg_)*6.0*vol*vvol; //(*vzp)[i];
      // Compute 0-1 penalty
      if ( rmin_ <= 0.0 ) {
        (*hvp)[i] -= (pen_)*2.0*(*vp)[i];
      }
    }
    // Compute 0-1 penalty
    vector Tv(zp->size(),0.0);
    uint i1 = 0, i2 = 0, j1 = 0, j2 = 0;
    Real area = 0.0;
    if (rmin_ > 0.0) {
      for (uint i=0; i<FEM_->numX(); i++) {
        for (uint j=0; j<FEM_->numY(); j++) {
          i1 = std::max(i-(uint)floor(rmin_),(uint)1)-1;
          i2 = std::min(i+(uint)floor(rmin_),FEM_->numX());
          j1 = std::max(j-(uint)floor(rmin_),(uint)1)-1;
          j2 = std::min(j+(uint)floor(rmin_),FEM_->numY());
          area = (Real)(i2-i1)*(j2-j1);
          for (uint ii=i1; ii<i2; ii++) {
            for (uint jj=j1; jj<j2; jj++) {
              Tv[i+j*FEM_->numX()] += (*vp)[ii+FEM_->numX()*jj]/area;
            }
          }
        }
      }
      for (uint i=0; i<FEM_->numX(); i++) {
        for (uint j=0; j<FEM_->numY(); j++) {
          i1 = std::max(i-(uint)floor(rmin_),(uint)1)-1;
          i2 = std::min(i+(uint)floor(rmin_),FEM_->numX());
          j1 = std::max(j-(uint)floor(rmin_),(uint)1)-1;
          j2 = std::min(j+(uint)floor(rmin_),FEM_->numY());
          area = (Real)(i2-i1)*(j2-j1);
          for (uint ii=i1; ii<i2; ii++) {
            for (uint jj=j1; jj<j2; jj++) {
              (*hvp)[ii+jj*FEM_->numX()] -= (pen_)*2.0*Tv[i+FEM_->numX()*j]/area;
            }
          }
        }
      }
    }
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  typedef typename std::vector<RealT>::size_type uint;
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
    // FEM problem description.
    int prob = 1;  // prob = 0 is the MBB beam example, prob = 1 is the cantilever beam example.
    uint nx  = 12; // Number of x-elements (60 for prob = 0, 32 for prob = 1).
    uint ny  = 8; // Number of y-elements (20 for prob = 0, 20 for prob = 1).
    int P    = 1;  // SIMP penalization power.
    ROL::Ptr<FEM<RealT> > pFEM = ROL::makePtr<FEM<RealT>>(nx,ny,P,prob);
    // Objective function description.
    int   nreg = 11;       // # of Moreau-Yoshida parameter updates (e.g., 21).
    int   npen = 2;        // # of penalty parameter updates (e.g., 10).
    RealT frac = 0.4;      // Volume fraction.
    RealT reg  = 1.0;      // Moreau-Yoshida regularization parameter.
    RealT pen  = 1.e-4;    // 0-1 penalty parameter. 
    RealT rmin = 1.2;      // Radius for spatial average.
    // Optimization parameters.
    bool derivCheck = true;  // Derivative check flag.
    bool useTR      = false; // Use trust-region or line-search.
    RealT gtol      = 1e-5;  // Norm of gradient tolerance.
    RealT stol      = 1e-8;  // Norm of step tolerance.
    int   maxit     = 100;   // Maximum number of iterations (e.g., 500).
    // Read optimization input parameter list.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    // Initialize ROL::Ptrs.
    ROL::Ptr<ROL::Objective_SimOpt<RealT> >         pobj;   // Full objective.
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj;   // Reduced objective.
    ROL::Ptr<ROL::Algorithm<RealT> >                algo;   // Optimization algorithm.
    ROL::Ptr<ROL::Step<RealT> >                     step;   // Globalized step.
    ROL::Ptr<ROL::StatusTest<RealT> >               status; // Status test.
    // Initialize equality constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon = 
      ROL::makePtr<Constraint_TopOpt<RealT>>(pFEM);
    // Initialize bound constraints.
    std::vector<RealT> lo(pFEM->numZ(),1.e-3);
    std::vector<RealT> hi(pFEM->numZ(),1.0);
    ROL::StdBoundConstraint<RealT> bound(lo,hi);
    // Initialize control vector.
    ROL::Ptr<std::vector<RealT> > z_ptr  = ROL::makePtr<std::vector<RealT>> (pFEM->numZ(), frac);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtrFromRef(z);
    // Initialize state vector.
    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(pFEM->numU(), 0.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtrFromRef(u);
    // Initialize adjoint vector.
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(pFEM->numU(), 0.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp  = ROL::makePtrFromRef(p);
    // Derivative check.
    if (derivCheck) {
      // Initialize control vectors.
      ROL::Ptr<std::vector<RealT> > yz_ptr = ROL::makePtr<std::vector<RealT>>(pFEM->numZ(), 0.0);
      for (uint i=0; i<pFEM->numZ(); i++) {
        (*yz_ptr)[i] = frac * (RealT)rand()/(RealT)RAND_MAX;
      }
      ROL::StdVector<RealT> yz(yz_ptr);
      ROL::Ptr<ROL::Vector<RealT> > yzp = ROL::makePtrFromRef(yz);
      // Initialize state vectors.
      ROL::Ptr<std::vector<RealT> > yu_ptr = ROL::makePtr<std::vector<RealT>>(pFEM->numU(), 0.0);
      for (uint i=0; i<pFEM->numU(); i++) {
      (*u_ptr)[i]  = (RealT)rand()/(RealT)RAND_MAX;
        (*yu_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      }
      ROL::StdVector<RealT> yu(yu_ptr);
      ROL::Ptr<ROL::Vector<RealT> > yup = ROL::makePtrFromRef(yu);
      // Initialize Jacobian vector.
      ROL::Ptr<std::vector<RealT> > jv_ptr  = ROL::makePtr<std::vector<RealT>>(pFEM->numU(), 0.0);
      ROL::StdVector<RealT> jv(jv_ptr);
      ROL::Ptr<ROL::Vector<RealT> > jvp = ROL::makePtrFromRef(jv);
      // Initialize SimOpt Vectors 
      ROL::Vector_SimOpt<RealT> x(up,zp);
      ROL::Vector_SimOpt<RealT> y(yup,yzp);
      // Test equality constraint.
      pcon->checkApplyJacobian(x,y,jv,true,*outStream);
      //pcon->checkApplyAdjointJacobian(x,yu,jv,x,true);
      pcon->checkApplyAdjointHessian(x,yu,y,x,true,*outStream);
      // Test full objective function.
      pobj = ROL::makePtr<Objective_TopOpt<RealT>>(pFEM,frac,reg,pen,rmin);
      pobj->checkGradient(x,y,true,*outStream);
      pobj->checkHessVec(x,y,true,*outStream);
      // Test reduced objective function.
      robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobj,pcon,up,zp,pp);
      robj->checkGradient(z,yz,true,*outStream);
      robj->checkHessVec(z,yz,true,*outStream);
    }
    // Run optimization.
    for ( int j=0; j<npen; j++ ) {
      *outStream << "\nPenalty parameter: " << pen << "\n";
      for ( int i=0; i<nreg; i++ ) {
        *outStream << "\nMoreau-Yoshida regularization parameter: " << reg << "\n";
        // Initialize full objective function.
        pobj = ROL::makePtr<Objective_TopOpt<RealT>>(pFEM,frac,reg,pen,rmin);
        // Initialize reduced objective function.
        robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobj,pcon,up,zp,pp);
        if ( !useTR ) {
          // Run line-search secant step.
          parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Quasi-Newton Method");
          parlist->sublist("General").sublist("Secant").set("Type", "Limited-Memory SR1");
          if ( maxit > 0 ) {
            step   = ROL::makePtr<ROL::LineSearchStep<RealT>>(*parlist);
            status = ROL::makePtr<ROL::StatusTest<RealT>>(gtol,stol,maxit);
            algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
            algo->run(z,*robj,bound,true,*outStream);
          }
          // Run line-search Newton-Krylov step.
          parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Newton-Krylov");
          step   = ROL::makePtr<ROL::LineSearchStep<RealT>>(*parlist);
          status = ROL::makePtr<ROL::StatusTest<RealT>>(gtol,stol,maxit);
          algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
          algo->run(z,*robj,bound,true,*outStream);
        }
        else {
          // Run trust-region step.
          step   = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
          status = ROL::makePtr<ROL::StatusTest<RealT>>(gtol,stol,maxit);
          algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
          algo->run(z,*robj,bound,true,*outStream);
        }
        // Compute volume.
        RealT vol = 0.0;
        for (uint i=0; i<nx; i++) {
          for (uint j=0; j<ny; j++) {
            vol += (*z_ptr)[i+j*nx];
          }
        }
        *outStream << "The volume fraction is " << vol/pFEM->numZ() << "\n";
        // Increase Moreau-Yoshida regularization parameter.
        reg *= 2.0;
      }
      // Print to file.
      std::stringstream name;
      name << "density" << j << ".txt";
      std::ofstream file;
      file.open(name.str().c_str());
      RealT val = 0.0;
      for (uint i=0; i<nx; i++) {
        for (uint j=0; j<ny; j++) {
          val = (*z_ptr)[i+j*nx];
          file << i << "  " << j << "  " << val << "\n"; 
        }
      }
      file.close();
      // Increase parameters.
      reg   = 1.0;
      pen  *= 10.0;
    }
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

