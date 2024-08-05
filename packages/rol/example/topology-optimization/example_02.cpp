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

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_OptimizationSolver.hpp"
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
};

template<class Real>
class Constraint_Volume : public ROL::Constraint<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
  
  typedef typename vector::size_type uint;
 
private:
  const Real frac_;

  ROL::Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const SV&>(x).getVector();
  }
  
  ROL::Ptr<vector> getVector( V& x ) {
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  Constraint_Volume(const Real frac) : frac_(frac) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<vector> cp = getVector(c);
    ROL::Ptr<const vector> zp = getVector(z);
    (*cp)[0] = 0;
    for (uint i = 0; i < zp->size(); i++) {
      (*cp)[0] += (*zp)[i];
    }
    (*cp)[0] -= frac_*static_cast<Real>(zp->size());
  }

  void applyJacobian(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<vector> jvp = getVector(jv);
    ROL::Ptr<const vector> vp = getVector(v);
    (*jvp)[0] = 0;
    for (uint i = 0; i < vp->size(); i++) {
      (*jvp)[0] += (*vp)[i];
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<vector> ajvp = getVector(ajv);
    ROL::Ptr<const vector> vp = getVector(v);
    for (uint i = 0; i < ajvp->size(); i++) {
      (*ajvp)[i] = (*vp)[0];
    }
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w,
                              const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class Objective_TopOpt : public ROL::Objective_SimOpt<Real> {

  typedef std::vector<Real>    vector;
  typedef ROL::Vector<Real>    V;
  typedef ROL::StdVector<Real> SV;
  
  typedef typename vector::size_type uint;  

private:
  ROL::Ptr<FEM<Real> > FEM_;

  bool useLC_; // Use linear form of compliance.  Otherwise use quadratic form.

  ROL::Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  Objective_TopOpt(ROL::Ptr<FEM<Real> > FEM) 
    : FEM_(FEM), useLC_(true) {}

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
    return c;
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
    uint nx  = 30; // Number of x-elements (60 for prob = 0, 32 for prob = 1).
    uint ny  = 10; // Number of y-elements (20 for prob = 0, 20 for prob = 1).
    int P    = 1;  // SIMP penalization power.
    RealT frac = 0.4;      // Volume fraction.
    ROL::Ptr<FEM<RealT> > pFEM = ROL::makePtr<FEM<RealT>>(nx,ny,P,prob);
    // Read optimization input parameter list.
    std::string filename = "input_ex02.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    // Initialize ROL::Ptrs.
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>         pobj;   // Full objective.
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> robj;   // Reduced objective.
    ROL::Ptr<ROL::BoundConstraint<RealT>>          bound;  // Bound constraint.
    // Initialize volume constraint.
    ROL::Ptr<ROL::Constraint<RealT>> vcon
      = ROL::makePtr<Constraint_Volume<RealT>>(frac);
    // Initialize equality constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> pcon
      = ROL::makePtr<Constraint_TopOpt<RealT>>(pFEM);
    // Initialize bound constraints.
    ROL::Ptr<std::vector<RealT>> lo_ptr = ROL::makePtr<std::vector<RealT>>(pFEM->numZ(),1.e-3);
    ROL::Ptr<std::vector<RealT>> hi_ptr = ROL::makePtr<std::vector<RealT>>(pFEM->numZ(),1.0);
    ROL::Ptr<ROL::Vector<RealT> > lop = ROL::makePtr<ROL::StdVector<RealT>>(lo_ptr);
    ROL::Ptr<ROL::Vector<RealT> > hip = ROL::makePtr<ROL::StdVector<RealT>>(hi_ptr);
    bound = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);
    // Initialize control vector.
    ROL::Ptr<std::vector<RealT> > z_ptr = ROL::makePtr<std::vector<RealT>> (pFEM->numZ(), frac);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp = ROL::makePtrFromRef(z);
    // Initialize state vector.
    ROL::Ptr<std::vector<RealT> > u_ptr = ROL::makePtr<std::vector<RealT>>(pFEM->numU(), 0.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtrFromRef(u);
    // Initialize adjoint vector.
    ROL::Ptr<std::vector<RealT> > p_ptr = ROL::makePtr<std::vector<RealT>>(pFEM->numU(), 0.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp = ROL::makePtrFromRef(p);
    // Initialize multiplier vector.
    ROL::Ptr<std::vector<RealT> > l_ptr = ROL::makePtr<std::vector<RealT>>(1, 0.0);
    ROL::StdVector<RealT> l(l_ptr);
    ROL::Ptr<ROL::Vector<RealT> > lp = ROL::makePtrFromRef(l);
    // Initialize objective function.
    pobj = ROL::makePtr<Objective_TopOpt<RealT>>(pFEM);
    // Initialize reduced objective function.
    robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobj,pcon,up,zp,pp);
    // Run optimization.
    ROL::OptimizationProblem<RealT> problem(robj, zp, bound, vcon, lp); 
    bool derivCheck = true;  // Derivative check flag.
    if (derivCheck) {
      problem.check(*outStream);
    }
    ROL::OptimizationSolver<RealT> solver(problem, *parlist);
    solver.solve(*outStream);
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

