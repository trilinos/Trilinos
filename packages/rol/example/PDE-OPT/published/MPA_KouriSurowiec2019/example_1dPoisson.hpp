// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include <cmath>
#include <limits>

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"


template<class Real>
class FEM {
private:
  int  nx_;    // Number of intervals
  Real kl_;    // Left diffusivity
  Real kr_;    // Right diffusivity

  std::vector<Real> pts_;
  std::vector<Real> wts_;
public:

  FEM(int nx = 128, Real kl = 0.1, Real kr = 10.0) : nx_(nx), kl_(kl), kr_(kr) {
    // Set quadrature points
    pts_.clear();
    pts_.resize(5,0.0);
    pts_[0] = -0.9061798459386640; 
    pts_[1] = -0.5384693101056831;
    pts_[2] =  0.0;
    pts_[3] =  0.5384693101056831;
    pts_[4] =  0.9061798459386640;
    wts_.clear();
    wts_.resize(5,0.0);
    wts_[0] = 0.2369268850561891;
    wts_[1] = 0.4786286704993665;
    wts_[2] = 0.5688888888888889;
    wts_[3] = 0.4786286704993665;
    wts_[4] = 0.2369268850561891;
    // Scale and normalize points and weights
    Real sum = 0.0;
    for (int i = 0; i < 5; i++) {
      pts_[i] = 0.5*pts_[i] + 0.5;
      sum += wts_[i];
    }
    for (int i = 0; i < 5; i++) {
      wts_[i] /= sum;
    }
  } 

  int nu(void) { return nx_-1; }
  int nz(void) { return nx_+1; }

  void build_mesh(std::vector<Real> &x, const std::vector<Real> &param) {
    // Partition mesh: 
    //     nx1 is the # of intervals on the left of xp -- n1 is the # of interior nodes
    //     nx2 is the # of intervals on the right of xp -- n2 is the # of interior nodes
    Real  xp = 0.1*param[0];
    int frac = nx_/2;
    int  rem = nx_%2;  
    int  nx1 = frac;
    int  nx2 = frac;
    if ( rem == 1 ) {
      if (xp > 0.0) {
        nx2++;
      }
      else {
        nx1++;
      }
    }
    int n1 = nx1-1;
    int n2 = nx2-1;
    // Compute mesh spacing
    Real dx1 = (xp+1.0)/(Real)nx1;
    Real dx2 = (1.0-xp)/(Real)nx2;
    // Build uniform meshes on [-1,xp] u [xp,1]
    x.clear();
    x.resize(n1+n2+3,0.0);
    for (int k = 0; k < n1+n2+3; k++) {
      x[k] = ((k < nx1) ? (Real)k*dx1-1.0 : ((k==nx1) ? xp : (Real)(k-nx1)*dx2+xp));
      //std::cout << x[k] << "\n";
    }
  }

  void build_mesh(std::vector<Real> &x) {
    int frac = nz()/16;
    int rem  = nz()%16;
    int nz1  = frac*4;
    int nz2  = frac;
    int nz3  = frac*9;
    int nz4  = frac*2;
    for ( int i = 0; i < rem; i++ ) {
      if ( i%4 == 0 )      { nz3++; }
      else if ( i%4 == 1 ) { nz1++; }
      else if ( i%4 == 2 ) { nz4++; } 
      else if ( i%4 == 3 ) { nz2++; }
    }
    x.clear();
    x.resize(nz(),0.0);
    for (int k = 0; k < nz(); k++) {
      if ( k < nz1 ) {
        x[k] = 0.25*(Real)k/(Real)(nz1-1) - 1.0;
      }
      if ( k >= nz1 && k < nz1+nz2 ) {
        x[k] = 0.5*(Real)(k-nz1+1)/(Real)nz2 - 0.75;
      }
      if ( k >= nz1+nz2 && k < nz1+nz2+nz3 ) {
        x[k] = 0.5*(Real)(k-nz1-nz2+1)/(Real)nz3 - 0.25;
      }
      if ( k >= nz1+nz2+nz3-1 ) {
        x[k] = 0.75*(Real)(k-nz1-nz2-nz3+1)/(Real)nz4 + 0.25;
      }
      //std::cout << x[k] << "\n";
    }
    //x.clear(); 
    //x.resize(nz(),0.0);
    //for (int k = 0; k < nz(); k++) {
    //  x[k] = 2.0*(Real)k/(Real)(nz()-1)-1.0;
    //  //std::cout << xz[k] << "\n";
    //}
  }

  // Build force.
  void build_force(std::vector<Real> &F, const std::vector<Real> &param) {
    // Build mesh
    std::vector<Real> x(nu()+2,0.0);
    build_mesh(x,param);
    // Build force term
    F.clear();
    F.resize(x.size()-2,0.0);
    Real pt = 0.0;
    int size = static_cast<int>(x.size());
    for (int i = 0; i < size-2; i++) {
      for (int j = 0; j < 5; j++) {
        // Integrate over [xl,x0]
        pt = (x[i+1]-x[i])*pts_[j] + x[i];
        F[i] += wts_[j]*(pt-x[i])*std::exp(-std::pow(pt-0.5*param[1],2.0));
        // Integrate over [x0,xu]
        pt = (x[i+2]-x[i+1])*pts_[j] + x[i+1];
        F[i] += wts_[j]*(x[i+2]-pt)*std::exp(-std::pow(pt-0.5*param[1],2.0));
      }
    }
  }

  // Build the PDE residual Jacobian.  
  void build_jacobian_1(std::vector<Real> &dl, std::vector<Real> &d, std::vector<Real> &du, 
                  const std::vector<Real> &param) {
    // Build mesh
    std::vector<Real> x(nu()+2,0.0);
    build_mesh(x,param);
    // Fill diagonal
    int xsize = static_cast<int>(x.size());
    d.clear();
    d.resize(xsize-2,0.0);
    for ( int i = 0; i < xsize-2; i++ ) {
      if ( x[i+1] < 0.1*param[0] ) {
        d[i] = kl_/(x[i+1]-x[i]) + kl_/(x[i+2]-x[i+1]);
      }
      else if ( x[i+1] > 0.1*param[0] ) {
        d[i] = kr_/(x[i+1]-x[i]) + kr_/(x[i+2]-x[i+1]);
      }
      else {
        d[i] = kl_/(x[i+1]-x[i]) + kr_/(x[i+2]-x[i+1]);
      }
      //std::cout << d[i] << "\n";
    }
    // Fill off-diagonal
    dl.clear();
    dl.resize(xsize-3,0.0);
    for ( int i = 0; i < xsize-3; i++ ) {
      if ( x[i+2] <= 0.1*param[0] ) {
        dl[i] = -kl_/(x[i+2]-x[i+1]);
      }
      else if ( x[i+2] > 0.1*param[0] ) {
        dl[i] = -kr_/(x[i+2]-x[i+1]);
      }
      //std::cout << dl[i] << "\n";
    }
    du.clear();
    du.assign(dl.begin(),dl.end());
  }

  // Apply the PDE residual Jacobian.  
  void apply_jacobian_1(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &param) {
    // Build Jacobian
    std::vector<Real> dl(nu()-1,0.0);
    std::vector<Real> d(nu(),0.0);
    std::vector<Real> du(nu()-1,0.0);
    build_jacobian_1(dl,d,du,param);
    // Apply Jacobian
    jv.clear();
    jv.resize(d.size(),0.0);
    int size = static_cast<int>(d.size());
    for ( int i = 0; i < size; ++i ) {
      jv[i]  = d[i]*v[i];
      jv[i] += ((i>0) ? dl[i-1]*v[i-1] : 0.0);
      jv[i] += ((i<size-1) ? du[i]*v[i+1] : 0.0);
    }
  }

  void apply_jacobian_2(std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &param) {
    // Build u mesh
    std::vector<Real> xu(nu()+2,0.0);
    build_mesh(xu,param);
    // Build z mesh 
    std::vector<Real> xz(nz(),0.0);
    build_mesh(xz);
    // Build union of u and z meshes
    std::vector<Real> x(xu.size()+xz.size(),0.0);
    typename std::vector<Real>::iterator it;
    it = std::set_union(xu.begin(),xu.end(),xz.begin(),xz.end(),x.begin());
    x.resize(it-x.begin());
    // Interpolate v onto x basis
    std::vector<Real> vi;
    int xzsize = static_cast<int>(xz.size());
    for (it=x.begin(); it!=x.end(); it++) {
      for (int i = 0; i < xzsize-1; i++) {
        if ( *it >= xz[i] && *it <= xz[i+1] ) {
          vi.push_back(v[i+1]*(*it-xz[i])/(xz[i+1]-xz[i]) + v[i]*(xz[i+1]-*it)/(xz[i+1]-xz[i]));
          break;
        }
      }
    }
    int xsize = static_cast<int>(x.size());
    // Apply x basis mass matrix to interpolated v
    std::vector<Real> Mv(xsize,0.0);
    for (int i = 0; i < xsize; i++) {
      Mv[i] += ((i>0) ? (x[i]-x[i-1])/6.0*vi[i-1] + (x[i]-x[i-1])/3.0*vi[i] : 0.0);
      Mv[i] += ((i<xsize-1) ? (x[i+1]-x[i])/3.0*vi[i] + (x[i+1]-x[i])/6.0*vi[i+1] : 0.0);
    }
    // Reduced mass times v to u basis
    int xusize = static_cast<int>(xu.size());
    jv.clear();
    jv.resize(xusize-2,0.0);
    for (int i = 0; i < xusize-2; i++) {
      for (int j = 0; j < xsize-1; j++) {
        jv[i] += (((x[j]>=xu[i  ]) && (x[j]< xu[i+1])) ? Mv[j]*(x[j]-xu[i  ])/(xu[i+1]-xu[i  ]) : 0.0);
        jv[i] += (((x[j]>=xu[i+1]) && (x[j]<=xu[i+2])) ? Mv[j]*(xu[i+2]-x[j])/(xu[i+2]-xu[i+1]) : 0.0);
        if ( x[j] > xu[i+2] ) { break; }
      }
    }
  }

  void apply_adjoint_jacobian_2(std::vector<Real> &jv, const std::vector<Real> &v, 
                          const std::vector<Real> &param){
    // Build u mesh
    std::vector<Real> xu(nu()+2,0.0);
    build_mesh(xu,param);
    // Build z mesh 
    std::vector<Real> xz(nz(),0.0);
    build_mesh(xz);
    // Build union of u and z meshes
    std::vector<Real> x(xu.size()+xz.size(),0.0);
    typename std::vector<Real>::iterator it;
    it = std::set_union(xu.begin(),xu.end(),xz.begin(),xz.end(),x.begin());
    x.resize(it-x.begin());
    // Interpolate v onto x basis
    std::vector<Real> vi;
    Real val = 0.0;
    int xusize = static_cast<int>(xu.size());
    for (it=x.begin(); it!=x.end(); it++) {
      for (int i = 0; i < xusize-1; i++) {
        if ( *it >= xu[i] && *it <= xu[i+1] ) {
          val = 0.0;
          val += ((i!=0       ) ? v[i-1]*(xu[i+1]-*it)/(xu[i+1]-xu[i]) : 0.0);
          val += ((i!=xusize-2) ? v[i  ]*(*it-xu[i  ])/(xu[i+1]-xu[i]) : 0.0);
          vi.push_back(val);
          break;
        }
      }
    }
    // Apply x basis mass matrix to interpolated v
    int xsize = static_cast<int>(x.size());
    std::vector<Real> Mv(xsize,0.0);
    for (int i = 0; i < xsize; i++) {
      Mv[i] += ((i>0) ? (x[i]-x[i-1])/6.0*vi[i-1] + (x[i]-x[i-1])/3.0*vi[i] : 0.0);
      Mv[i] += ((i<xsize-1) ? (x[i+1]-x[i])/3.0*vi[i] + (x[i+1]-x[i])/6.0*vi[i+1] : 0.0);
    }
    // Reduced mass times v to u basis
    jv.clear();
    jv.resize(nz(),0.0);
    int xzsize = static_cast<int>(xz.size());
    for (int i = 0; i < xzsize; i++) {
      for (int j = 0; j < xsize; j++) {
        if ( i==0 ) {
          jv[i] += (((x[j]>=xz[i  ]) && (x[j]<=xz[i+1])) ? Mv[j]*(xz[i+1]-x[j])/(xz[i+1]-xz[i  ]) : 0.0);
          if ( x[j] > xz[i+1] ) { break; }
        }
        else if ( i == xzsize-1 ) {
          jv[i] += (((x[j]>=xz[i-1]) && (x[j]<=xz[i  ])) ? Mv[j]*(x[j]-xz[i-1])/(xz[i  ]-xz[i-1]) : 0.0);
        }
        else { 
          jv[i] += (((x[j]>=xz[i-1]) && (x[j]< xz[i  ])) ? Mv[j]*(x[j]-xz[i-1])/(xz[i  ]-xz[i-1]) : 0.0);
          jv[i] += (((x[j]>=xz[i  ]) && (x[j]<=xz[i+1])) ? Mv[j]*(xz[i+1]-x[j])/(xz[i+1]-xz[i  ]) : 0.0);
          if ( x[j] > xz[i+1] ) { break; }
        }
      }
    }
  }

  void apply_mass_state(std::vector<Real> &Mv, const std::vector<Real> &v, const std::vector<Real> &param) {
    // Build u mesh
    std::vector<Real> x;
    build_mesh(x,param);
    // Apply mass matrix
    int size = static_cast<int>(x.size());
    Mv.clear();
    Mv.resize(size-2,0.0);
    for (int i = 0; i < size-2; i++) {
      Mv[i] = (x[i+2]-x[i])/3.0*v[i];
      if ( i > 0 ) {
        Mv[i] += (x[i+1]-x[i])/6.0*v[i-1];
      }
      if ( i < size-3 ) {
        Mv[i] += (x[i+2]-x[i+1])/6.0*v[i+1];
      }
    }
  }

  void apply_mass_control(std::vector<Real> &Mv, const std::vector<Real> &v) {
    // Build z mesh 
    std::vector<Real> x(nz(),0.0);
    build_mesh(x);
    // Apply mass matrix
    int xsize = static_cast<Real>(x.size());
    Mv.clear();
    Mv.resize(xsize,0.0);
    for (int i = 0; i < xsize; i++) {
      if ( i > 0 ) {
        Mv[i] += (x[i]-x[i-1])/6.0*v[i-1] + (x[i]-x[i-1])/3.0*v[i];
      }
      if ( i < xsize-1 ) {
        Mv[i] += (x[i+1]-x[i])/3.0*v[i] + (x[i+1]-x[i])/6.0*v[i+1];
      }
    }
  }
   
  void apply_inverse_mass_control(std::vector<Real> &Mv, const std::vector<Real> &v) {
    // Build z mesh
    std::vector<Real> x(nz(),0.0);
    build_mesh(x);
    // Build mass matrix
    std::vector<Real> d(nz(),0.0);
    std::vector<Real> dl(nz()-1,0.0);
    std::vector<Real> du(nz()-1,0.0);
    for (int i = 0; i < nz(); i++) {
      if ( i > 0 ) {
        dl[i-1] = (x[i]-x[i-1])/6.0;
        d[i]   += (x[i]-x[i-1])/3.0;
      }
      if ( i < nz()-1 ) {
        d[i] += (x[i+1]-x[i])/3.0;
        du[i] = (x[i+1]-x[i])/6.0;
      }
    }
    // Solve linear system
    linear_solve(Mv,dl,d,du,v,false);
  }

  // Solve linear systems in which the matrix is triadiagonal.
  // This function uses LAPACK routines for:
  //   1.) Computing the LU factorization of a tridiagonal matrix.
  //   2.) Use LU factors to solve the linear system.
  void linear_solve(std::vector<Real> &u, std::vector<Real> &dl, std::vector<Real> &d, 
                    std::vector<Real> &du, const std::vector<Real> &r, const bool transpose = false) {
    u.clear();
    u.assign(r.begin(),r.end());
    // Perform LDL factorization
    Teuchos::LAPACK<int,Real> lp;
    std::vector<Real> du2(r.size()-2,0.0);
    std::vector<int> ipiv(r.size(),0);
    int n    = r.size();
    int info;
    int ldb  = n;
    int nhrs = 1;
    lp.GTTRF(n,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&info);
    char trans = 'N';
    if ( transpose ) { 
      trans = 'T';
    }
    lp.GTTRS(trans,n,nhrs,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&u[0],ldb,&info);
  }

  // Evaluates the difference between the solution to the PDE.
  // and the desired profile
  Real evaluate_target(int i, const std::vector<Real> &param) {
    std::vector<Real> x;
    build_mesh(x,param);
    Real val = 0.0;
    int example = 2;
    switch (example) {
      case 1:  val = ((x[i]<0.5) ? 1.0 : 0.0);             break;
      case 2:  val = 1.0;                                  break;
      case 3:  val = std::abs(std::sin(8.0*M_PI*x[i]));    break;
      case 4:  val = std::exp(-0.5*(x[i]-0.5)*(x[i]-0.5)); break;
    }
    return val;
  }
};

template<class Real>
class DiffusionConstraint : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<FEM<Real> > FEM_;
  int num_solves_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  // Returns u <- (u + alpha*s) where u, s are vectors and 
  // alpha is a scalar.
  void plus(std::vector<Real> &u, const std::vector<Real> &s, const Real alpha=1.0) {
    int size = static_cast<int>(u.size());
    for (int i=0; i<size; i++) {
      u[i] += alpha*s[i];
    }
  }

  // Returns u <- alpha*u where u is a vector and alpha is a scalar
  void scale(std::vector<Real> &u, const Real alpha=0.0) {
    int size = static_cast<int>(u.size());
    for (int i=0; i<size; i++) {
      u[i] *= alpha;
    }
  }
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:
  DiffusionConstraint(const ROL::Ptr<FEM<Real> > &FEM) : FEM_(FEM), num_solves_(0) {}

  int getNumSolves(void) const {
    return num_solves_;
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp = dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Apply PDE Jacobian
    FEM_->apply_jacobian_1(*cp, *up, this->getParameter());
    // Get Right Hand Side
    std::vector<Real> F(FEM_->nu(),0.0);
    FEM_->build_force(F,this->getParameter());
    std::vector<Real> Bz(FEM_->nu(),0.0);
    FEM_->apply_jacobian_2(Bz,*zp,this->getParameter());
    plus(F,Bz);
    // Add Right Hand Side to PDE Term
    plus(*cp,F,-1.0);
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp = dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<std::vector<Real> > up = dynamic_cast<ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // Get PDE Jacobian
    std::vector<Real> dl(FEM_->nu()-1,0.0);
    std::vector<Real> d(FEM_->nu(),0.0);
    std::vector<Real> du(FEM_->nu()-1,0.0);
    FEM_->build_jacobian_1(dl,d,du,this->getParameter());
    // Get Right Hand Side
    std::vector<Real> F(FEM_->nu(),0.0);
    FEM_->build_force(F,this->getParameter());
    std::vector<Real> Bz(FEM_->nu(),0.0);
    FEM_->apply_jacobian_2(Bz,*zp,this->getParameter());
    plus(F,Bz);
    // Solve solve state
    FEM_->linear_solve(*up,dl,d,du,F);
    num_solves_++;
    // Compute residual
    value(c,u,z,tol);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    FEM_->apply_jacobian_1(*jvp, *vp, this->getParameter());
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    FEM_->apply_jacobian_2(*jvp, *vp, this->getParameter());
    scale(*jvp,-1.0);
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // Get PDE Jacobian
    std::vector<Real> dl(FEM_->nu()-1,0.0);
    std::vector<Real> d(FEM_->nu(),0.0);
    std::vector<Real> du(FEM_->nu()-1,0.0);
    FEM_->build_jacobian_1(dl,d,du,this->getParameter());
    // Solve solve state
    FEM_->linear_solve(*jvp,dl,d,du,*vp);
    num_solves_++;
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    applyJacobian_1(jv,v,u,z,tol);
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    FEM_->apply_adjoint_jacobian_2(*jvp, *vp, this->getParameter());
    scale(*jvp,-1.0);
  }
 
  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol) {
    applyInverseJacobian_1(jv,v,u,z,tol);
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
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

// Class BurgersObjective contains the necessary information
// to compute the value and gradient of the objective function.
template<class Real>
class DiffusionObjective : public ROL::Objective_SimOpt<Real> {
private:
  const ROL::Ptr<FEM<Real> > FEM_;
  const Real alpha_;

/***************************************************************/
/********** BEGIN PRIVATE MEMBER FUNCTION DECLARATION **********/
/***************************************************************/
  // Evaluates the discretized L2 inner product of two vectors.
  Real dot(const std::vector<Real> &x, const std::vector<Real> &y) {
    Real ip = 0.0;
    int size = static_cast<int>(x.size());
    for (int i=0; i<size; i++) {
      ip += x[i]*y[i];
    }
    return ip;
  }

  // Returns u <- alpha*u where u is a vector and alpha is a scalar
  void scale(std::vector<Real> &u, const Real alpha=0.0) {
    int size = static_cast<int>(u.size());
    for (int i=0; i<size; i++) {
      u[i] *= alpha;
    }
  }
/*************************************************************/
/********** END PRIVATE MEMBER FUNCTION DECLARATION **********/
/*************************************************************/

public:
  DiffusionObjective(const ROL::Ptr<FEM<Real> > fem, const Real alpha = 1.e-4)
    : FEM_(fem), alpha_(alpha) {}

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    int nz = FEM_->nz(), nu = FEM_->nu();
    // COMPUTE CONTROL PENALTY
    std::vector<Real> Mz(nz);
    FEM_->apply_mass_control(Mz,*zp);
    Real zval = alpha_*0.5*dot(Mz,*zp);   
    // COMPUTE STATE TRACKING TERM
    std::vector<Real> diff(nu);
    for (int i=0; i<nu; i++) {
      diff[i] = ((*up)[i]-FEM_->evaluate_target(i+1,this->getParameter()));
    }
    std::vector<Real> Mu(nu);
    FEM_->apply_mass_state(Mu,diff,this->getParameter());
    Real uval = 0.5*dot(Mu,diff);  
    return uval+zval;
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > up = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    int nu = FEM_->nu();
    // COMPUTE GRADIENT
    std::vector<Real> diff(nu);
    for (int i=0; i<nu; i++) {
      diff[i] = ((*up)[i]-FEM_->evaluate_target(i+1,this->getParameter()));
    }
    FEM_->apply_mass_state(*gp,diff,this->getParameter());
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp = dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > zp = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    // COMPUTE GRADIENT
    FEM_->apply_mass_control(*gp,*zp);
    scale(*gp,alpha_);
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
             const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp = dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE HESSVEC
    FEM_->apply_mass_state(*hvp,*vp,this->getParameter());
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
    ROL::Ptr<std::vector<Real> > hvp = dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    // COMPUTE HESSVEC
    FEM_->apply_mass_control(*hvp,*vp);
    scale(*hvp,alpha_);
  }

};
