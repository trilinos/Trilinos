// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_06.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_TeuchosBatchManager.hpp"

#include "Teuchos_LAPACK.hpp"

template<class Real>
class L2VectorPrimal;

template<class Real>
class L2VectorDual;

template<class Real>
class H1VectorPrimal;

template<class Real>
class H1VectorDual;

template<class Real>
class BurgersFEM {
private:
  int nx_;
  Real dx_;
  Real nu_;
  Real nl_;
  Real u0_;
  Real u1_;
  Real f_;
  Real cH1_;
  Real cL2_;

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
  BurgersFEM(int nx = 128, Real nl = 1.0, Real cH1 = 1.0, Real cL2 = 1.0) 
    : nx_(nx), dx_(1.0/((Real)nx+1.0)), nl_(nl), cH1_(cH1), cL2_(cL2) {}

  void set_problem_data(const Real nu, const Real f, const Real u0, const Real u1) {
    nu_ = std::pow(10.0,nu-2.0);
    f_  = 0.01*f;
    u0_ = 1.0+0.001*u0;
    u1_ = 0.001*u1;
  }

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

  /***************************************************************************/
  /*********************** H1 INNER PRODUCT FUNCTIONS ************************/
  /***************************************************************************/
  // Compute H1 inner product
  Real compute_H1_dot(const std::vector<Real> &x, const std::vector<Real> &y) const {
    Real ip = 0.0;
    for (int i=0; i<nx_; i++) {
      if ( i == 0 ) {
        ip += cL2_*dx_/6.0*(4.0*x[i] + x[i+1])*y[i]; // Mass term
        ip += cH1_*(2.0*x[i] - x[i+1])/dx_*y[i];     // Stiffness term
      }
      else if ( i == nx_-1 ) {
        ip += cL2_*dx_/6.0*(x[i-1] + 4.0*x[i])*y[i]; // Mass term
        ip += cH1_*(2.0*x[i] - x[i-1])/dx_*y[i];     // Stiffness term
      }
      else {
        ip += cL2_*dx_/6.0*(x[i-1] + 4.0*x[i] + x[i+1])*y[i]; // Mass term
        ip += cH1_*(2.0*x[i] - x[i-1] - x[i+1])/dx_*y[i];     // Stiffness term
      }
    }
    return ip;
  }

  // compute H1 norm
  Real compute_H1_norm(const std::vector<Real> &r) const {
    return std::sqrt(compute_H1_dot(r,r));
  }

  // Apply H2 Reisz operator
  void apply_H1(std::vector<Real> &Mu, const std::vector<Real> &u ) const {
    Mu.resize(nx_,0.0);
    for (int i=0; i<nx_; i++) {
      if ( i == 0 ) {
        Mu[i] = cL2_*dx_/6.0*(4.0*u[i] + u[i+1])
              + cH1_*(2.0*u[i] - u[i+1])/dx_;
      }
      else if ( i == nx_-1 ) {
        Mu[i] = cL2_*dx_/6.0*(u[i-1] + 4.0*u[i])
              + cH1_*(2.0*u[i] - u[i-1])/dx_;
      }
      else {
        Mu[i] = cL2_*dx_/6.0*(u[i-1] + 4.0*u[i] + u[i+1])
              + cH1_*(2.0*u[i] - u[i-1] - u[i+1])/dx_;
      }
    }
  }

  // Apply H1 inverse Reisz operator
  void apply_inverse_H1(std::vector<Real> &Mu, const std::vector<Real> &u) const {
    // Build mass matrix
    std::vector<Real> dl(nx_-1,cL2_*dx_/6.0   - cH1_/dx_);
    std::vector<Real> d(nx_,2.0*(cL2_*dx_/3.0 + cH1_/dx_));
    std::vector<Real> du(nx_-1,cL2_*dx_/6.0   - cH1_/dx_);
    linear_solve(Mu,dl,d,du,u);
  }

  void test_inverse_H1(std::ostream &outStream = std::cout) {
    std::vector<Real> u(nx_,0.0), Mu(nx_,0.0), iMMu(nx_,0.0), diff(nx_,0.0);
    for (int i = 0; i < nx_; i++) {
      u[i] = 2.0*(Real)rand()/(Real)RAND_MAX - 1.0;
    }
    apply_H1(Mu,u);
    apply_inverse_H1(iMMu,Mu);
    axpy(diff,-1.0,iMMu,u);
    Real error = compute_H1_norm(diff);
    Real normu = compute_H1_norm(u);
    outStream << "Test Inverse State H1 Matrix\n";
    outStream << "  ||u - inv(M)Mu|| = " << error << "\n";
    outStream << "  ||u||            = " << normu << "\n";
    outStream << "  Relative Error   = " << error/normu << "\n";
    outStream << "\n";
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
        r[i] += nl_*u[i+1]*(u[i]+u[i+1])/6.0;
      }
      if (i>0) {
        r[i] -= nl_*u[i-1]*(u[i-1]+u[i])/6.0;
      }
      // Contribution from control
      r[i] -= dx_/6.0*(z[i]+4.0*z[i+1]+z[i+2]);
      // Contribution from right-hand side
      r[i] -= dx_*f_;
    }
    // Contribution from Dirichlet boundary terms
    r[0]     -= nl_*(u0_*u[    0]/6.0 + u0_*u0_/6.0) + nu_*u0_/dx_;
    r[nx_-1] += nl_*(u1_*u[nx_-1]/6.0 + u1_*u1_/6.0) - nu_*u1_/dx_;
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
        dl[i] += nl_*(-2.0*u[i]-u[i+1])/6.0;
        d[i]  += nl_*u[i+1]/6.0;
      }
      if (i>0) {
        d[i]    -= nl_*u[i-1]/6.0;
        du[i-1] += nl_*(u[i-1]+2.0*u[i])/6.0;
      }
    }
    // Contribution from Dirichlet boundary conditions
    d[    0] -= nl_*u0_/6.0;
    d[nx_-1] += nl_*u1_/6.0;
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
        jv[i] += -nu_/dx_*v[i-1]-nl_*(u[i-1]/6.0*v[i]+(u[i]+2.0*u[i-1])/6.0*v[i-1]);
      }
      if ( i < nx_-1 ) {
        jv[i] += -nu_/dx_*v[i+1]+nl_*(u[i+1]/6.0*v[i]+(u[i]+2.0*u[i+1])/6.0*v[i+1]);
      }
    }
    jv[    0] -= nl_*u0_/6.0*v[0];
    jv[nx_-1] += nl_*u1_/6.0*v[nx_-1];
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
    for (int i = 0; i < nx_; i++) {
      ajv[i] = nu_/dx_*2.0*v[i];
      if ( i > 0 ) {
        ajv[i] += -nu_/dx_*v[i-1]-nl_*(u[i-1]/6.0*v[i]
                  -(u[i-1]+2.0*u[i])/6.0*v[i-1]);
      }
      if ( i < nx_-1 ) {
        ajv[i] += -nu_/dx_*v[i+1]+nl_*(u[i+1]/6.0*v[i]
                  -(u[i+1]+2.0*u[i])/6.0*v[i+1]);
      }
    }
    ajv[    0] -= nl_*u0_/6.0*v[0];
    ajv[nx_-1] += nl_*u1_/6.0*v[nx_-1];
  }

  // Apply inverse adjoint PDE Jacobian to a vector
  void apply_inverse_adjoint_pde_jacobian(std::vector<Real> &iajv,
                                    const std::vector<Real> &v,
                                    const std::vector<Real> &u,
                                    const std::vector<Real> &z) const {
    // Get PDE Jacobian
    std::vector<Real> d(nx_,0.0);
    std::vector<Real> du(nx_-1,0.0);
    std::vector<Real> dl(nx_-1,0.0);
    compute_pde_jacobian(dl,d,du,u);
    // Solve solve adjoint system at current time step
    linear_solve(iajv,dl,d,du,v,true);
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
    for (int i=0; i<nx_+2; i++) {
      if ( i == 0 ) {
        jv[i] = -dx_/6.0*v[i];
      }
      else if ( i == 1 ) {
        jv[i] = -dx_/6.0*(4.0*v[i-1]+v[i]);
      }
      else if ( i == nx_ ) {
        jv[i] = -dx_/6.0*(4.0*v[i-1]+v[i-2]);
      }
      else if ( i == nx_+1 ) {
        jv[i] = -dx_/6.0*v[i-2];
      }
      else {
        jv[i] = -dx_/6.0*(v[i-2]+4.0*v[i-1]+v[i]);
      }
    }
  }

  /***************************************************************************/
  /*********************** AJDOINT HESSIANS **********************************/
  /***************************************************************************/
  void apply_adjoint_pde_hessian(std::vector<Real> &ahwv,
                           const std::vector<Real> &w,
                           const std::vector<Real> &v,
                           const std::vector<Real> &u,
                           const std::vector<Real> &z) const {
    for (int i=0; i<nx_; i++) {
      // Contribution from nonlinear term
      ahwv[i] = 0.0;
      if (i<nx_-1){
        ahwv[i] += (w[i]*v[i+1] - w[i+1]*(2.0*v[i]+v[i+1]))/6.0;
      }
      if (i>0) {
        ahwv[i] += (w[i-1]*(v[i-1]+2.0*v[i]) - w[i]*v[i-1])/6.0;
      }
    }
    //ahwv.assign(u.size(),0.0);
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
class L2VectorPrimal : public ROL::Vector<Real> {
private:
  ROL::Ptr<std::vector<Real> > vec_;
  ROL::Ptr<BurgersFEM<Real> > fem_;

  mutable ROL::Ptr<L2VectorDual<Real> > dual_vec_;

public:
  L2VectorPrimal(const ROL::Ptr<std::vector<Real> > & vec,
                 const ROL::Ptr<BurgersFEM<Real> > &fem)
    : vec_(vec), fem_(fem), dual_vec_(ROL::nullPtr) {}

  void set( const ROL::Vector<Real> &x ) {
    const L2VectorPrimal &ex = dynamic_cast<const L2VectorPrimal&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),vec_->begin());
  }

  void plus( const ROL::Vector<Real> &x ) {
    const L2VectorPrimal &ex = dynamic_cast<const L2VectorPrimal&>(x);
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
    const L2VectorPrimal & ex = dynamic_cast<const L2VectorPrimal&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return fem_->compute_L2_dot(xval,*vec_);
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  ROL::Ptr<ROL::Vector<Real> > clone() const {
    return ROL::makePtr<L2VectorPrimal>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
  }

  ROL::Ptr<const std::vector<Real> > getVector() const {
    return vec_;
  }

  ROL::Ptr<std::vector<Real> > getVector() {
    return vec_;
  }

  ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    ROL::Ptr<L2VectorPrimal> e
      = ROL::makePtr<L2VectorPrimal>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return vec_->size();
  }

  const ROL::Vector<Real>& dual() const {
    dual_vec_ = ROL::makePtr<L2VectorDual<Real>>(
      ROL::makePtr<std::vector<Real>>(*vec_),fem_);

    fem_->apply_mass(*(ROL::constPtrCast<std::vector<Real> >(dual_vec_->getVector())),*vec_);
    return *dual_vec_;
  }

  Real apply(const ROL::Vector<Real> &x) const {
    const L2VectorDual<Real> &ex = dynamic_cast<const L2VectorDual<Real>&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return std::inner_product(vec_->begin(), vec_->end(), xval.begin(), Real(0));
  }

};

template<class Real>
class L2VectorDual : public ROL::Vector<Real> {
private:
  ROL::Ptr<std::vector<Real> > vec_;
  ROL::Ptr<BurgersFEM<Real> > fem_;

  mutable ROL::Ptr<L2VectorPrimal<Real> > dual_vec_;

public:
  L2VectorDual(const ROL::Ptr<std::vector<Real> > & vec,
               const ROL::Ptr<BurgersFEM<Real> > &fem)
    : vec_(vec), fem_(fem), dual_vec_(ROL::nullPtr) {}

  void set( const ROL::Vector<Real> &x ) {
    const L2VectorDual &ex = dynamic_cast<const L2VectorDual&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),vec_->begin());
  }

  void plus( const ROL::Vector<Real> &x ) {
    const L2VectorDual &ex = dynamic_cast<const L2VectorDual&>(x);
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
    const L2VectorDual & ex = dynamic_cast<const L2VectorDual&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    unsigned dimension = vec_->size();
    std::vector<Real> Mx(dimension,0.0);
    fem_->apply_inverse_mass(Mx,xval);
    Real val = 0.0;
    for (unsigned i = 0; i < dimension; i++) {
      val += Mx[i]*(*vec_)[i];
    }
    return val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  ROL::Ptr<ROL::Vector<Real> > clone() const {
    return ROL::makePtr<L2VectorDual>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
  }

  ROL::Ptr<const std::vector<Real> > getVector() const {
    return vec_;
  }

  ROL::Ptr<std::vector<Real> > getVector() {
    return vec_;
  }

  ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    ROL::Ptr<L2VectorDual> e
      = ROL::makePtr<L2VectorDual>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return vec_->size();
  }

  const ROL::Vector<Real>& dual() const {
    dual_vec_ = ROL::makePtr<L2VectorPrimal<Real>>(
      ROL::makePtr<std::vector<Real>>(*vec_),fem_);

    fem_->apply_inverse_mass(*(ROL::constPtrCast<std::vector<Real> >(dual_vec_->getVector())),*vec_);
    return *dual_vec_;
  }

  Real apply(const ROL::Vector<Real> &x) const {
    const L2VectorPrimal<Real> &ex = dynamic_cast<const L2VectorPrimal<Real>&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return std::inner_product(vec_->begin(), vec_->end(), xval.begin(), Real(0));
  }

};

template<class Real>
class H1VectorPrimal : public ROL::Vector<Real> {
private:
  ROL::Ptr<std::vector<Real> > vec_;
  ROL::Ptr<BurgersFEM<Real> > fem_;

  mutable ROL::Ptr<H1VectorDual<Real> > dual_vec_;

public:
  H1VectorPrimal(const ROL::Ptr<std::vector<Real> > & vec,
                 const ROL::Ptr<BurgersFEM<Real> > &fem)
    : vec_(vec), fem_(fem), dual_vec_(ROL::nullPtr) {}

  void set( const ROL::Vector<Real> &x ) {
    const H1VectorPrimal &ex = dynamic_cast<const H1VectorPrimal&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),vec_->begin());
  }

  void plus( const ROL::Vector<Real> &x ) {
    const H1VectorPrimal &ex = dynamic_cast<const H1VectorPrimal&>(x);
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
    const H1VectorPrimal & ex = dynamic_cast<const H1VectorPrimal&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return fem_->compute_H1_dot(xval,*vec_);
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  ROL::Ptr<ROL::Vector<Real> > clone() const {
    return ROL::makePtr<H1VectorPrimal>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
  }

  ROL::Ptr<const std::vector<Real> > getVector() const {
    return vec_;
  }

  ROL::Ptr<std::vector<Real> > getVector() {
    return vec_;
  }

  ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    ROL::Ptr<H1VectorPrimal> e
      = ROL::makePtr<H1VectorPrimal>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return vec_->size();
  }

  const ROL::Vector<Real>& dual() const {
    dual_vec_ = ROL::makePtr<H1VectorDual<Real>>(
      ROL::makePtr<std::vector<Real>>(*vec_),fem_);

    fem_->apply_H1(*(ROL::constPtrCast<std::vector<Real> >(dual_vec_->getVector())),*vec_);
    return *dual_vec_;
  }

  Real apply(const ROL::Vector<Real> &x) const {
    const H1VectorDual<Real> &ex = dynamic_cast<const H1VectorDual<Real>&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return std::inner_product(vec_->begin(), vec_->end(), xval.begin(), Real(0));
  }

};

template<class Real>
class H1VectorDual : public ROL::Vector<Real> {
private:
  ROL::Ptr<std::vector<Real> > vec_;
  ROL::Ptr<BurgersFEM<Real> > fem_;

  mutable ROL::Ptr<H1VectorPrimal<Real> > dual_vec_;

public:
  H1VectorDual(const ROL::Ptr<std::vector<Real> > & vec,
               const ROL::Ptr<BurgersFEM<Real> > &fem)
    : vec_(vec), fem_(fem), dual_vec_(ROL::nullPtr) {}

  void set( const ROL::Vector<Real> &x ) {
    const H1VectorDual &ex = dynamic_cast<const H1VectorDual&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),vec_->begin());
  }

  void plus( const ROL::Vector<Real> &x ) {
    const H1VectorDual &ex = dynamic_cast<const H1VectorDual&>(x);
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
    const H1VectorDual & ex = dynamic_cast<const H1VectorDual&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    unsigned dimension = vec_->size();
    std::vector<Real> Mx(dimension,0.0);
    fem_->apply_inverse_H1(Mx,xval);
    Real val = 0.0;
    for (unsigned i = 0; i < dimension; i++) {
      val += Mx[i]*(*vec_)[i];
    }
    return val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  ROL::Ptr<ROL::Vector<Real> > clone() const {
    return ROL::makePtr<H1VectorDual>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
  }

  ROL::Ptr<const std::vector<Real> > getVector() const {
    return vec_;
  }

  ROL::Ptr<std::vector<Real> > getVector() {
    return vec_;
  }

  ROL::Ptr<ROL::Vector<Real> > basis( const int i ) const {
    ROL::Ptr<H1VectorDual> e
      = ROL::makePtr<H1VectorDual>( ROL::makePtr<std::vector<Real>>(vec_->size(),0.0),fem_);
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return vec_->size();
  }

  const ROL::Vector<Real>& dual() const {
    dual_vec_ = ROL::makePtr<H1VectorPrimal<Real>>(
      ROL::makePtr<std::vector<Real>>(*vec_),fem_);

    fem_->apply_inverse_H1(*(ROL::constPtrCast<std::vector<Real> >(dual_vec_->getVector())),*vec_);
    return *dual_vec_;
  }

  Real apply(const ROL::Vector<Real> &x) const {
    const H1VectorPrimal<Real> &ex = dynamic_cast<const H1VectorPrimal<Real>&>(x);
    const std::vector<Real>& xval = *ex.getVector();
    return std::inner_product(vec_->begin(), vec_->end(), xval.begin(), Real(0));
  }

};

template<class Real>
class Constraint_BurgersControl : public ROL::Constraint_SimOpt<Real> {
private:

  typedef H1VectorPrimal<Real> PrimalStateVector;
  typedef H1VectorDual<Real> DualStateVector;
  
  typedef L2VectorPrimal<Real> PrimalControlVector;
  typedef L2VectorDual<Real> DualControlVector;
  
  typedef H1VectorDual<Real> PrimalConstraintVector;
  typedef H1VectorPrimal<Real> DualConstraintVector;

  ROL::Ptr<BurgersFEM<Real> > fem_;
  bool useHessian_;

public:
  Constraint_BurgersControl(ROL::Ptr<BurgersFEM<Real> > &fem, bool useHessian = true)
   : fem_(fem), useHessian_(useHessian) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp =
      dynamic_cast<PrimalConstraintVector&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->compute_residual(*cp,*up,*zp);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<PrimalConstraintVector&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const PrimalStateVector&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->apply_pde_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<PrimalConstraintVector&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const PrimalControlVector&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->apply_control_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ijvp =
      dynamic_cast<PrimalStateVector&>(ijv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const PrimalConstraintVector&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->apply_inverse_pde_jacobian(*ijvp,*vp,*up,*zp);
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<DualStateVector&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const DualConstraintVector&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->apply_adjoint_pde_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp =
      dynamic_cast<DualControlVector&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const DualConstraintVector&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->apply_adjoint_control_jacobian(*jvp,*vp,*up,*zp);
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > iajvp =
      dynamic_cast<DualConstraintVector&>(iajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const DualStateVector&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    const std::vector<Real> param
      = ROL::Constraint_SimOpt<Real>::getParameter();
    fem_->set_problem_data(param[0],param[1],param[2],param[3]);

    fem_->apply_inverse_adjoint_pde_jacobian(*iajvp,*vp,*up,*zp);
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    if ( useHessian_ ) {
      ROL::Ptr<std::vector<Real> > ahwvp =
        dynamic_cast<DualStateVector&>(ahwv).getVector();
      ROL::Ptr<const std::vector<Real> > wp =
        dynamic_cast<const DualConstraintVector&>(w).getVector();
      ROL::Ptr<const std::vector<Real> > vp =
        dynamic_cast<const PrimalStateVector&>(v).getVector();
      ROL::Ptr<const std::vector<Real> > up =
        dynamic_cast<const PrimalStateVector&>(u).getVector();
      ROL::Ptr<const std::vector<Real> > zp =
        dynamic_cast<const PrimalControlVector&>(z).getVector();

      const std::vector<Real> param
        = ROL::Constraint_SimOpt<Real>::getParameter();
      fem_->set_problem_data(param[0],param[1],param[2],param[3]);

      fem_->apply_adjoint_pde_hessian(*ahwvp,*wp,*vp,*up,*zp);
    }
    else {
      ahwv.zero();
    }
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    if ( useHessian_ ) {
      ROL::Ptr<std::vector<Real> > ahwvp =
        dynamic_cast<DualControlVector&>(ahwv).getVector();
      ROL::Ptr<const std::vector<Real> > wp =
        dynamic_cast<const DualConstraintVector&>(w).getVector();
      ROL::Ptr<const std::vector<Real> > vp =
        dynamic_cast<const PrimalStateVector&>(v).getVector();
      ROL::Ptr<const std::vector<Real> > up =
        dynamic_cast<const PrimalStateVector&>(u).getVector();
      ROL::Ptr<const std::vector<Real> > zp =
        dynamic_cast<const PrimalControlVector&>(z).getVector();

      const std::vector<Real> param
        = ROL::Constraint_SimOpt<Real>::getParameter();
      fem_->set_problem_data(param[0],param[1],param[2],param[3]);

      fem_->apply_adjoint_control_pde_hessian(*ahwvp,*wp,*vp,*up,*zp);
    }
    else {
      ahwv.zero();
    }
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    if ( useHessian_ ) {
      ROL::Ptr<std::vector<Real> > ahwvp =
        dynamic_cast<DualStateVector&>(ahwv).getVector();
      ROL::Ptr<const std::vector<Real> > wp =
        dynamic_cast<const DualConstraintVector&>(w).getVector();
      ROL::Ptr<const std::vector<Real> > vp =
        dynamic_cast<const PrimalControlVector&>(v).getVector();
      ROL::Ptr<const std::vector<Real> > up =
        dynamic_cast<const PrimalStateVector&>(u).getVector();
      ROL::Ptr<const std::vector<Real> > zp =
        dynamic_cast<const PrimalControlVector&>(z).getVector();

      const std::vector<Real> param
        = ROL::Constraint_SimOpt<Real>::getParameter();
      fem_->set_problem_data(param[0],param[1],param[2],param[3]);

      fem_->apply_adjoint_pde_control_hessian(*ahwvp,*wp,*vp,*up,*zp);
    }
    else {
      ahwv.zero();
    }
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    if ( useHessian_ ) {
      ROL::Ptr<std::vector<Real> > ahwvp =
        dynamic_cast<DualControlVector&>(ahwv).getVector();
      ROL::Ptr<const std::vector<Real> > wp =
        dynamic_cast<const DualConstraintVector&>(w).getVector();
      ROL::Ptr<const std::vector<Real> > vp =
        dynamic_cast<const PrimalControlVector&>(v).getVector();
      ROL::Ptr<const std::vector<Real> > up =
        dynamic_cast<const PrimalStateVector&>(u).getVector();
      ROL::Ptr<const std::vector<Real> > zp =
        dynamic_cast<const PrimalControlVector&>(z).getVector();

      const std::vector<Real> param
        = ROL::Constraint_SimOpt<Real>::getParameter();
      fem_->set_problem_data(param[0],param[1],param[2],param[3]);

      fem_->apply_adjoint_control_hessian(*ahwvp,*wp,*vp,*up,*zp);
    }
    else {
      ahwv.zero();
    }
  }
};

template<class Real>
class Objective_BurgersControl : public ROL::Objective_SimOpt<Real> {
private:

  typedef H1VectorPrimal<Real> PrimalStateVector;
  typedef H1VectorDual<Real> DualStateVector;
  
  typedef L2VectorPrimal<Real> PrimalControlVector;
  typedef L2VectorDual<Real> DualControlVector;

  Real alpha_; // Penalty Parameter
  ROL::Ptr<BurgersFEM<Real> > fem_;
  ROL::Ptr<ROL::Vector<Real> > ud_;
  ROL::Ptr<ROL::Vector<Real> > diff_;

public:
  Objective_BurgersControl(const ROL::Ptr<BurgersFEM<Real> > &fem, 
                           const ROL::Ptr<ROL::Vector<Real> > &ud,
                           Real alpha = 1.e-4) : alpha_(alpha), fem_(fem), ud_(ud) {
    diff_ = ud_->clone();
  }

  using ROL::Objective_SimOpt<Real>::value;

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();
    ROL::Ptr<const std::vector<Real> > udp =
      dynamic_cast<const L2VectorPrimal<Real>&>(*ud_).getVector();

    std::vector<Real> diff(udp->size(),0.0);
    for (unsigned i = 0; i < udp->size(); i++) {
      diff[i] = (*up)[i] - (*udp)[i];
    }
    return 0.5*(fem_->compute_L2_dot(diff,diff) + alpha_*fem_->compute_L2_dot(*zp,*zp));
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp =
      dynamic_cast<DualStateVector&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > up =
      dynamic_cast<const PrimalStateVector&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > udp =
      dynamic_cast<const L2VectorPrimal<Real>&>(*ud_).getVector();

    std::vector<Real> diff(udp->size(),0.0);
    for (unsigned i = 0; i < udp->size(); i++) {
      diff[i] = (*up)[i] - (*udp)[i];
    }
    fem_->apply_mass(*gp,diff);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp =
      dynamic_cast<DualControlVector&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      dynamic_cast<const PrimalControlVector&>(z).getVector();

    fem_->apply_mass(*gp,*zp);
    for (unsigned i = 0; i < zp->size(); i++) {
      (*gp)[i] *= alpha_;
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp =
      dynamic_cast<DualStateVector&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const PrimalStateVector&>(v).getVector();

    fem_->apply_mass(*hvp,*vp);
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
    ROL::Ptr<std::vector<Real> > hvp =
      dynamic_cast<DualControlVector&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp =
      dynamic_cast<const PrimalControlVector&>(v).getVector();

    fem_->apply_mass(*hvp,*vp);
    for (unsigned i = 0; i < vp->size(); i++) {
      (*hvp)[i] *= alpha_;
    }
  }
};

template<class Real, class Ordinal>
class L2VectorBatchManager : public ROL::TeuchosBatchManager<Real,Ordinal> {
private:
  void cast_vector(ROL::Ptr<std::vector<Real> > &xvec,
                   ROL::Vector<Real> &x) const {
    try {
      xvec = dynamic_cast<L2VectorPrimal<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      xvec = dynamic_cast<L2VectorDual<Real>&>(x).getVector();
    }
  }

public:
  L2VectorBatchManager(const ROL::Ptr<const Teuchos::Comm<int> > &comm)
    : ROL::TeuchosBatchManager<Real,Ordinal>(comm) {}
  void sumAll(ROL::Vector<Real> &input, ROL::Vector<Real> &output) {
    ROL::Ptr<std::vector<Real> > input_ptr;
    cast_vector(input_ptr,input);
    int dim_i = input_ptr->size();
    ROL::Ptr<std::vector<Real> > output_ptr;
    cast_vector(output_ptr,output);
    int dim_o = output_ptr->size();
    if ( dim_i != dim_o ) {
      std::cout << "L2VectorBatchManager: DIMENSION MISMATCH ON RANK "
                << ROL::TeuchosBatchManager<Real,Ordinal>::batchID() << "\n";
    }
    else {
      ROL::TeuchosBatchManager<Real,Ordinal>::sumAll(&(*input_ptr)[0],&(*output_ptr)[0],dim_i);
    }
  }
};

template<class Real, class Ordinal>
class H1VectorBatchManager : public ROL::TeuchosBatchManager<Real,Ordinal> {
private:
  void cast_vector(ROL::Ptr<std::vector<Real> > &xvec,
                   ROL::Vector<Real> &x) const {
    try {
      xvec = dynamic_cast<H1VectorPrimal<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      xvec = dynamic_cast<H1VectorDual<Real>&>(x).getVector();
    }
  }

public:
  H1VectorBatchManager(const ROL::Ptr<const Teuchos::Comm<int> > &comm)
    : ROL::TeuchosBatchManager<Real,Ordinal>(comm) {}
  void sumAll(ROL::Vector<Real> &input, ROL::Vector<Real> &output) {
    ROL::Ptr<std::vector<Real> > input_ptr;
    cast_vector(input_ptr,input);
    int dim_i = input_ptr->size();
    ROL::Ptr<std::vector<Real> > output_ptr;
    cast_vector(output_ptr,output);
    int dim_o = output_ptr->size();
    if ( dim_i != dim_o ) {
      std::cout << "H1VectorBatchManager: DIMENSION MISMATCH ON RANK "
                << ROL::TeuchosBatchManager<Real,Ordinal>::batchID() << "\n";
    }
    else {
      ROL::TeuchosBatchManager<Real,Ordinal>::sumAll(&(*input_ptr)[0],&(*output_ptr)[0],dim_i);
    }
  }
};

template<class Real>
Real random(const ROL::Ptr<const Teuchos::Comm<int> > &comm) {
  Real val = 0.0;
  if ( Teuchos::rank<int>(*comm)==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*comm,0,1,&val);
  return val;
}
