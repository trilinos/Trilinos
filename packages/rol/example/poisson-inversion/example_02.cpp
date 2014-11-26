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

/*! \file  example_01.cpp
    \brief Shows how to solve the inverse Poisson problem using trust-region
           methods with dense Hessian diagnostics.
*/

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <algorithm>

template<class Real>
class Objective_PoissonInversion : public ROL::Objective<Real> {
private:
  int nu_;
  int nz_;

  Real hu_;
  Real hz_;

  Real u0_;
  Real u1_;

  Real alpha_;

  bool useCorrection_;
  Teuchos::SerialDenseMatrix<int, Real> H_;

public:

  /* CONSTRUCTOR */
  Objective_PoissonInversion(int nz = 32, Real alpha = 1.e-4) 
    : nz_(nz), u0_(1.0), u1_(0.0), alpha_(alpha), useCorrection_(false) {
    nu_ = nz_-1;
    hu_ = 1.0/((Real)nu_+1.0);
    hz_ = hu_; 
  }

  void apply_mass(std::vector<Real> &Mz, const std::vector<Real> &z ) {
    Mz.resize(this->nu_,0.0);
    for (int i=0; i<this->nu_; i++) {
      if ( i == 0 ) {
        Mz[i] = this->hu_/6.0*(2.0*z[i] + z[i+1]);
      }
      else if ( i == this->nu_-1 ) {
        Mz[i] = this->hu_/6.0*(z[i-1] + 2.0*z[i]);
      }
      else {
        Mz[i] = this->hu_/6.0*(z[i-1] + 4.0*z[i] + z[i+1]);
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
    Bd.resize(this->nu_,0.0);
    for (int i = 0; i < this->nu_; i++) {
      if ( i == 0 ) {
        Bd[i] = 1.0/this->hu_*( u[i]*d[i] + (u[i]-u[i+1])*d[i+1] );
      }
      else if ( i == this->nu_-1 ) {
        Bd[i] = 1.0/this->hu_*( (u[i]-u[i-1])*d[i] + u[i]*d[i+1] );
      }
      else {
        Bd[i] = 1.0/this->hu_*( (u[i]-u[i-1])*d[i] + (u[i]-u[i+1])*d[i+1] );
      }
    }
    if ( addBC ) {
      Bd[          0] -= this->u0_*d[           0]/this->hu_;
      Bd[this->nu_-1] -= this->u1_*d[this-> nz_-1]/this->hu_;
    }
  }

  void apply_transposed_linearized_control_operator( std::vector<Real> &Bd, const std::vector<Real> &z,
                                               const std::vector<Real> &d,  const std::vector<Real> &u,
                                                     bool addBC = true ) {
    Bd.clear();
    Bd.resize(this->nz_,0.0);
    for (int i = 0; i < this->nz_; i++) {
      if ( i == 0 ) {
        Bd[i] = 1.0/this->hu_*u[i]*d[i];
      }
      else if ( i == this->nz_-1 ) {
        Bd[i] = 1.0/this->hu_*u[i-1]*d[i-1];
      }
      else {
        Bd[i] = 1.0/this->hu_*( (u[i]-u[i-1])*(d[i]-d[i-1]) );
      }
    }
    if ( addBC ) {
      Bd[          0] -= this->u0_*d[           0]/this->hu_;
      Bd[this->nz_-1] -= this->u1_*d[this-> nu_-1]/this->hu_;
    }
  }

  /* STATE AND ADJOINT EQUATION DEFINTIONS */
  void solve_state_equation(std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nu_,1.0);
    std::vector<Real> o(this->nu_-1,1.0);
    for ( int i = 0; i < this->nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/this->hu_;
      if ( i < this->nu_-1 ) {
        o[i] *= -z[i+1]/this->hu_;
      }
    }
    // Set right hand side
    u.clear();
    u.resize(this->nu_,0.0);
    u[          0] = z[          0]/this->hu_ * this->u0_;
    u[this->nu_-1] = z[this->nz_-1]/this->hu_ * this->u1_;
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nu_;
    int nhrs = 1;
    lp.PTTRF(this->nu_,&d[0],&o[0],&info);
    lp.PTTRS(this->nu_,nhrs,&d[0],&o[0],&u[0],ldb,&info);
  }

  void solve_adjoint_equation(std::vector<Real> &p, const std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nu_,1.0);
    std::vector<Real> o(this->nu_-1,1.0);
    for ( int i = 0; i < this->nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/this->hu_;
      if ( i < this->nu_-1 ) {
        o[i] *= -z[i+1]/this->hu_;
      }
    }
    // Set right hand side
    std::vector<Real> r(this->nu_,0.0);
    for (int i = 0; i < this->nu_; i++) {
      r[i] = -(u[i]-this->evaluate_target((Real)(i+1)*this->hu_));
    }
    p.clear();
    p.resize(this->nu_,0.0);
    this->apply_mass(p,r);    
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nu_;
    int nhrs = 1;
    lp.PTTRF(this->nu_,&d[0],&o[0],&info);
    lp.PTTRS(this->nu_,nhrs,&d[0],&o[0],&p[0],ldb,&info);
  }

  void solve_state_sensitivity_equation(std::vector<Real> &w, const std::vector<Real> &v, 
                                  const std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nu_,1.0);
    std::vector<Real> o(this->nu_-1,1.0);
    for ( int i = 0; i < this->nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/this->hu_;
      if ( i < this->nu_-1 ) {
        o[i] *= -z[i+1]/this->hu_;
      }
    }
    // Set right hand side
    w.clear();
    w.resize(this->nu_,0.0);
    this->apply_linearized_control_operator(w,z,v,u);
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nu_;
    int nhrs = 1;
    lp.PTTRF(this->nu_,&d[0],&o[0],&info);
    lp.PTTRS(this->nu_,nhrs,&d[0],&o[0],&w[0],ldb,&info);
  }

  void solve_adjoint_sensitivity_equation(std::vector<Real> &q, const std::vector<Real> &w, 
                                    const std::vector<Real> &v, const std::vector<Real> &p, 
                                    const std::vector<Real> &u, const std::vector<Real> &z) {
    // Get Diagonal and Off-Diagonal Entries of PDE Jacobian
    std::vector<Real> d(this->nu_,1.0);
    std::vector<Real> o(this->nu_-1,1.0);
    for ( int i = 0; i < this->nu_; i++ ) {
      d[i] = (z[i] + z[i+1])/this->hu_;
      if ( i < this->nu_-1 ) {
        o[i] *= -z[i+1]/this->hu_;
      }
    }
    // Set right hand side
    q.clear();
    q.resize(this->nu_,0.0);
    this->apply_mass(q,w);
    std::vector<Real> res(this->nu_,0.0);
    this->apply_linearized_control_operator(res,z,v,p,false);
    for (int i = 0; i < this->nu_; i++) {
      q[i] -= res[i];
    }
    // Solve Tridiagonal System Using LAPACK's SPD Tridiagonal Solver
    Teuchos::LAPACK<int,Real> lp;
    int info;
    int ldb  = this->nu_;
    int nhrs = 1;
    lp.PTTRF(this->nu_,&d[0],&o[0],&info);
    lp.PTTRS(this->nu_,nhrs,&d[0],&o[0],&q[0],ldb,&info);
  }

  void update(const ROL::Vector<Real> &z, bool flag, int iter) {
    if ( flag && this->useCorrection_ ) {
      Real tol = std::sqrt(ROL::ROL_EPSILON);
      this->H_.shape(this->nz_,this->nz_); 
      Teuchos::RCP<ROL::Vector<Real> > e = z.clone();
      Teuchos::RCP<ROL::Vector<Real> > h = z.clone();
      for ( int i = 0; i < this->nz_; i++ ) {
        e = z.basis(i);
        this->hessVec_true(*h,*e,z,tol);
        for ( int j = 0; j < this->nz_; j++ ) {
          e = z.basis(j);
          (this->H_)(j,i) = e->dot(*h);
        }
      }
      std::vector<std::vector<Real> > eigenvals = ROL::computeEigenvalues<Real>(this->H_);
      std::sort((eigenvals[0]).begin(), (eigenvals[0]).end());
      Real inertia = (eigenvals[0])[0];
      Real correction = 0.0;
      if ( inertia <= 0.0 ) {
        correction = (1.0+std::sqrt(ROL::ROL_EPSILON))*std::abs(inertia);
        if ( inertia == 0.0 ) {
          int cnt = 0;
          while ( eigenvals[0][cnt] == 0.0 ) {
            cnt++;
          }
          correction = std::sqrt(ROL::ROL_EPSILON)*eigenvals[0][cnt];
          if ( cnt == this->nz_-1 ) {
            correction = 1.0;
          }
        }
        for ( int i = 0; i < this->nz_; i++ ) {
          (this->H_)(i,i) += correction;
        }
      }  
    }
  }

  /* OBJECTIVE FUNCTION DEFINITIONS */
  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    // SOLVE STATE EQUATION
    std::vector<Real> u(this->nu_,0.0);
    this->solve_state_equation(u,*zp);
    // EVALUATE OBJECTIVE
    Real val  = 0.0;
    for (int i=0; i<this->nz_;i++) {
      val += this->hz_*this->alpha_*0.5*(*zp)[i]*(*zp)[i];
    }
    Real res  = 0.0;
    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (int i=0; i<this->nu_; i++) {
      if ( i == 0 ) {
        res1 = u[i]-evaluate_target((Real)(i+1)*this->hu_);
        res2 = u[i+1]-evaluate_target((Real)(i+2)*this->hu_);
        res  = this->hu_/6.0*(2.0*res1 + res2)*res1;
      }
      else if ( i == this->nu_-1 ) {
        res1 = u[i-1]-evaluate_target((Real)i*this->hu_);
        res2 = u[i]-evaluate_target((Real)(i+1)*this->hu_);
        res  = this->hu_/6.0*(res1 + 2.0*res2)*res2;
      }
      else {
        res1 = u[i-1]-evaluate_target((Real)i*this->hu_);
        res2 = u[i]-evaluate_target((Real)(i+1)*this->hu_);
        res3 = u[i+1]-evaluate_target((Real)(i+2)*this->hu_);
        res  = this->hu_/6.0*(res1 + 4.0*res2 + res3)*res2;
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
    std::vector<Real> u(this->nu_,0.0);
    this->solve_state_equation(u,*zp);
    // SOLVE ADJOINT EQUATION
    std::vector<Real> p(this->nu_,0.0);
    this->solve_adjoint_equation(p,u,*zp);
    // Apply Transpose of Linearized Control Operator
    this->apply_transposed_linearized_control_operator(*gp,*zp,p,u);
    // Build Gradient
    for ( int i = 0; i < this->nz_; i++ ) {
      (*gp)[i] += this->hz_*this->alpha_*(*zp)[i];
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    if ( this->useCorrection_ ) {
      this->hessVec_inertia(hv,v,z,tol);
    }
    else {
      this->hessVec_true(hv,v,z,tol);
    }
  }

  void activateInertia(void) {
    this->useCorrection_ = true;
  }

  void deactivateInertia(void) {
    this->useCorrection_ = false;
  }

  void hessVec_true( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(z))).getVector();
    Teuchos::RCP<std::vector<Real> > hvp = 
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
    // SOLVE STATE EQUATION
    std::vector<Real> u(this->nu_,0.0);
    this->solve_state_equation(u,*zp);
    // SOLVE ADJOINT EQUATION
    std::vector<Real> p(this->nu_,0.0);
    this->solve_adjoint_equation(p,u,*zp);
    // SOLVE STATE SENSITIVITY EQUATION
    std::vector<Real> w(this->nu_,0.0);
    this->solve_state_sensitivity_equation(w,*vp,u,*zp);
    // SOLVE ADJOINT SENSITIVITY EQUATION
    std::vector<Real> q(this->nu_,0.0);
    this->solve_adjoint_sensitivity_equation(q,w,*vp,p,u,*zp);
    // Apply Transpose of Linearized Control Operator
    this->apply_transposed_linearized_control_operator(*hvp,*zp,q,u);
    // Apply Transpose of Linearized Control Operator
    std::vector<Real> tmp(this->nz_,0.0);
    this->apply_transposed_linearized_control_operator(tmp,*zp,w,p,false);
    for (int i=0; i < this->nz_; i++) {
      (*hvp)[i] -= tmp[i];
    }
    // Regularization hessVec
    for (int i=0; i < this->nz_; i++) {
      (*hvp)[i] += this->hz_*this->alpha_*(*vp)[i];
    }
  }

  void hessVec_inertia( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<std::vector<Real> >  vp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector());
    Teuchos::RCP<std::vector<Real> > hvp = 
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

    Teuchos::SerialDenseVector<int, Real> hv_teuchos(Teuchos::View, &((*hvp)[0]), this->nz_);
    Teuchos::SerialDenseVector<int, Real>  v_teuchos(Teuchos::View, &(( *vp)[0]), this->nz_);
    hv_teuchos.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, this->H_, v_teuchos, 0.0);
  }

};

template<class Real>
class BoundConstraint_PoissonInversion : public ROL::BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
public:
  BoundConstraint_PoissonInversion(std::vector<Real> &lo, std::vector<Real> &up) {
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

    int dim = 128; // Set problem dimension.
    RealT alpha = 1.e-6;
    Objective_PoissonInversion<RealT> obj(dim, alpha);

    // Iteration vector.
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    // Set initial guess.
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX + 1.e2;
      (*y_rcp)[i] = (RealT)rand()/(RealT)RAND_MAX + 1.e2;
    }
    ROL::StdVector<RealT> x(x_rcp);
    ROL::StdVector<RealT> y(y_rcp);
    obj.checkGradient(x,y,true);
    obj.checkHessVec(x,y,true);

    std::vector<RealT> lo(dim,1.0);
    std::vector<RealT> up(dim,10.0);
    BoundConstraint_PoissonInversion<RealT> icon(lo,up);

    Teuchos::ParameterList parlist;
    // Basic algorithm.
    parlist.set("Trust-Region Subproblem Solver Type",    "Truncated CG");
    parlist.set("Initial Trust-Region Radius",            100.0);
    // Secant parameters.
    parlist.set("Secant Type",                            "Limited-Memory BFGS");
    parlist.set("Maximum Secant Storage",                 100);
    // Krylov parameters.
    parlist.set("Absolute Krylov Tolerance",              1.e-8);
    parlist.set("Relative Krylov Tolerance",              1.e-4);
    parlist.set("Maximum Number of Krylov Iterations",    dim);
    // PDAS parameters.
    parlist.set("PDAS Relative Step Tolerance",           1.e-8);
    parlist.set("PDAS Relative Gradient Tolerance",       1.e-6);
    parlist.set("PDAS Maximum Number of Iterations",      10);
    parlist.set("PDAS Dual Scaling",                      alpha);      
    // Define step.
    parlist.set("Use Secant Hessian-Times-A-Vector",      true);
    ROL::PrimalDualActiveSetStep<RealT> step(parlist);

    // Define status test.
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-14;  // norm of step tolerance
    int   maxit = 20000;  // maximum number of iterations
    ROL::StatusTest<RealT> status(gtol, stol, maxit);    

    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo(step,status,false);

    x.zero();
    obj.deactivateInertia();
    algo.run(x,obj,icon,true);

    // Output control to file.
    std::ofstream file;
    file.open("control_PDAS.txt");
    for ( unsigned i = 0; i < (unsigned)dim; i++ ) {
      file << (*x_rcp)[i] << "\n";
    }
    file.close();

    // Projected Newtion.
    // Define step.
    parlist.set("Use Secant Hessian-Times-A-Vector",      false);
    ROL::TrustRegionStep<RealT> step_tr(parlist);
    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo_tr(step_tr,status,false);
    // Run Algorithm
    y.zero();
    obj.deactivateInertia();
    algo_tr.run(y,obj,icon,true);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( unsigned i = 0; i < (unsigned)dim; i++ ) {
      file_tr << (*y_rcp)[i] << "\n";
    }
    file_tr.close();

    std::vector<RealT> u;
    obj.solve_state_equation(u,*y_rcp);
    std::ofstream file_u;
    file_u.open("state.txt");
    for ( unsigned i = 0; i < (unsigned)(dim-1); i++ ) {
      file_u << u[i] << "\n";
    }
    file_u.close();
   
    Teuchos::RCP<ROL::Vector<RealT> > diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm()/std::sqrt((RealT)dim-1.0);
    std::cout << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1.e2*std::sqrt(ROL::ROL_EPSILON)) ? 1 : 0);

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

