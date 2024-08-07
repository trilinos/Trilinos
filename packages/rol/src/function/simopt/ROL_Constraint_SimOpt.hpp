// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_SIMOPT_H
#define ROL_CONSTRAINT_SIMOPT_H

#include "ROL_Constraint.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Types.hpp"
#include <iostream>

namespace ROL {

template <class Real>
class Constraint_SimOpt;

}

#include "ROL_NonlinearLeastSquaresObjective.hpp"
#include "ROL_SimConstraint.hpp"
#include "ROL_Objective_FSsolver.hpp"
#include "ROL_TypeU_TrustRegionAlgorithm.hpp"
#include "ROL_TypeE_AugmentedLagrangianAlgorithm.hpp"

/** @ingroup func_group
    \class ROL::Constraint_SimOpt
    \brief Defines the constraint operator interface for simulation-based optimization.

    This constraint interface inherits from ROL_Constraint, for the
    use case when \f$\mathcal{X}=\mathcal{U}\times\mathcal{Z}\f$ where \f$\mathcal{U}\f$ and 
    \f$\mathcal{Z}\f$ are Banach spaces.  \f$\mathcal{U}\f$ denotes the "simulation space"
    and \f$\mathcal{Z}\f$ denotes the "optimization space" (of designs, controls, parameters).
    The simulation-based constraints are of the form
    \f[
      c(u,z) = 0 \,.
    \f]
    The basic operator interface, to be implemented by the user, requires:
    \li #value -- constraint evaluation.
    \li #applyJacobian_1        -- action of the partial constraint Jacobian --derivatives are 
                                   with respect to the first component \f$\mathcal{U}\f$;
    \li #applyJacobian_2        -- action of the partial constraint Jacobian --derivatives are 
                                   with respect to the second component \f$\mathcal{Z}\f$;
    \li #applyAdjointJacobian_1 -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                   with respect to the first component \f$\mathcal{U}\f$;
    \li #applyAdjointJacobian_2 -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                   with respect to the second component \f$\mathcal{Z}\f$;

    The user may also overload:
    \li #applyAdjointHessian_11  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the first component only;
    \li #applyAdjointHessian_12  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the first and second components;
    \li #applyAdjointHessian_21  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the second and first components;
    \li #applyAdjointHessian_22  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the second component only;
    \li #solveAugmentedSystem -- solution of the augmented system --the default is an iterative
                                 scheme based on the action of the Jacobian and its adjoint.
    \li #applyPreconditioner  -- action of a constraint preconditioner --the default is null-op.

    ---
*/


namespace ROL {

template <class Real>
class Constraint_SimOpt : public Constraint<Real> {
private:
  // Additional vector storage for solve
  Ptr<Vector<Real>> unew_;
  Ptr<Vector<Real>> jv_;

  // Default parameters for solve (backtracking Newton)
  const Real DEFAULT_atol_;
  const Real DEFAULT_rtol_;
  const Real DEFAULT_stol_;
  const Real DEFAULT_factor_;
  const Real DEFAULT_decr_;
  const int  DEFAULT_maxit_;
  const bool DEFAULT_print_;
  const bool DEFAULT_zero_;
  const int  DEFAULT_solverType_;

  // User-set parameters for solve (backtracking Newton)

protected:
  Real atol_;
  Real rtol_;
  Real stol_;
  Real factor_;
  Real decr_;
  int  maxit_;
  bool print_;
  bool zero_;
  int  solverType_;

  // Flag to initialize vector storage in solve
  bool firstSolve_;

public:
  Constraint_SimOpt()
    : Constraint<Real>(),
      unew_(nullPtr), jv_(nullPtr),
      DEFAULT_atol_(1.e-4*std::sqrt(ROL_EPSILON<Real>())),
      DEFAULT_rtol_(1.e0),
      DEFAULT_stol_(std::sqrt(ROL_EPSILON<Real>())),
      DEFAULT_factor_(0.5),
      DEFAULT_decr_(1.e-4),
      DEFAULT_maxit_(500),
      DEFAULT_print_(false),
      DEFAULT_zero_(false),
      DEFAULT_solverType_(0),
      atol_(DEFAULT_atol_), rtol_(DEFAULT_rtol_), stol_(DEFAULT_stol_), factor_(DEFAULT_factor_),
      decr_(DEFAULT_decr_), maxit_(DEFAULT_maxit_), print_(DEFAULT_print_), zero_(DEFAULT_zero_),
      solverType_(DEFAULT_solverType_), firstSolve_(true) {}

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);  
  }
  virtual void update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1 ) {
    update_1(u,type,iter);
    update_2(z,type,iter);  
  }

  /** \brief Update constraint functions with respect to Sim variable.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1( const Vector<Real> &u, bool flag = true, int iter = -1 ) {}
  virtual void update_1( const Vector<Real> &u, UpdateType type, int iter = -1 ) {}

  /** \brief Update constraint functions with respect to Opt variable.
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) {}
  virtual void update_2( const Vector<Real> &z, UpdateType type, int iter = -1 ) {}

  /** \brief Update SimOpt constraint during solve (disconnected from optimization updates).
  
                @param[in]    x is the optimization variable
                @param[in] type is the update type
                @param[in] iter is the solver iteration count
  */
  virtual void solve_update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1) {}

  /** \brief Evaluate the constraint operator \f$c:\mathcal{U}\times\mathcal{Z} \rightarrow \mathcal{C}\f$
             at \f$(u,z)\f$.

             @param[out]      c   is the result of evaluating the constraint operator at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{c} = c(u,z)\f$,
             where \f$\mathsf{c} \in \mathcal{C}\f$, \f$\mathsf{u} \in \mathcal{U}\f$, and $\f$\mathsf{z} \in\mathcal{Z}\f$.

             ---
  */
  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) = 0;

  /** \brief Given \f$z\f$, solve \f$c(u,z)=0\f$ for \f$u\f$.

             @param[out]      c   is the result of evaluating the constraint operator at @b \f$(u,z)\f$; a constraint-space vector
             @param[in,out]   u   is the solution vector; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             The defualt implementation is Newton's method globalized with a backtracking line search.

             ---
  */
  virtual void solve(Vector<Real> &c,
                     Vector<Real> &u, 
                     const Vector<Real> &z,
                     Real &tol) {
    if ( zero_ ) u.zero();
    Ptr<std::ostream> stream = makeStreamPtr(std::cout, print_);
    solve_update(u,z,UpdateType::Initial,0);
    value(c,u,z,tol);
    Real cnorm = c.norm();
    const Real ctol = std::min(atol_, rtol_*cnorm);
    if (solverType_==0 || solverType_==3 || solverType_==4) {
      if ( firstSolve_ ) {
        unew_ = u.clone();
        jv_   = u.clone();
        firstSolve_ = false;
      }
      const Real one(1);
      Real alpha(1), tmp(0);
      int cnt = 0;
      *stream << std::endl;
      *stream << "     Default Constraint_SimOpt::solve" << std::endl;
      *stream << "       ";
      *stream << std::setw(6)  << std::left << "iter";
      *stream << std::setw(15) << std::left << "rnorm";
      *stream << std::setw(15) << std::left << "alpha";
      *stream << std::endl;
      for (cnt = 0; cnt < maxit_; ++cnt) {
        // Compute Newton step
        applyInverseJacobian_1(*jv_,c,u,z,tol);
        unew_->set(u);
        unew_->axpy(-alpha, *jv_);
        solve_update(*unew_,z,UpdateType::Trial);
        value(c,*unew_,z,tol);
        tmp = c.norm();
        // Perform backtracking line search
        while ( tmp > (one-decr_*alpha)*cnorm &&
                alpha > stol_ ) {
          alpha *= factor_;
          unew_->set(u);
          unew_->axpy(-alpha,*jv_);
          solve_update(*unew_,z,UpdateType::Trial);
          value(c,*unew_,z,tol);
          tmp = c.norm();
        }
        *stream << "       ";
        *stream << std::setw(6)  << std::left << cnt;
        *stream << std::scientific << std::setprecision(6);
        *stream << std::setw(15) << std::left << tmp;
        *stream << std::scientific << std::setprecision(6);
        *stream << std::setw(15) << std::left << alpha;
        *stream << std::endl;
        // Update iterate
        cnorm = tmp;
        u.set(*unew_);
        solve_update(u,z,UpdateType::Accept,cnt);
        if (cnorm <= ctol) break; // = covers the case of identically zero residual
        alpha = one;
      }
    }
    if (solverType_==1 || (solverType_==3 && cnorm > ctol)) {
      Ptr<Constraint<Real>> con = makePtr<SimConstraint<Real>>(makePtrFromRef(*this),makePtrFromRef(z),true);
      Ptr<Objective<Real>>  obj = makePtr<NonlinearLeastSquaresObjective<Real>>(con,u,c,true);
      ParameterList parlist;
      parlist.sublist("General").set("Output Level",1);
      parlist.sublist("General").sublist("Krylov").set("Iteration Limit",100);
      parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      parlist.sublist("Status Test").set("Gradient Tolerance",ctol);
      parlist.sublist("Status Test").set("Step Tolerance",stol_);
      parlist.sublist("Status Test").set("Iteration Limit",maxit_);
      Ptr<TypeU::Algorithm<Real>> algo = makePtr<TypeU::TrustRegionAlgorithm<Real>>(parlist);
      algo->run(u,*obj,*stream);
      value(c,u,z,tol);
    }
    if (solverType_==2 || (solverType_==4 && cnorm > ctol)) {
      Ptr<Constraint<Real>> con = makePtr<SimConstraint<Real>>(makePtrFromRef(*this),makePtrFromRef(z),true);
      Ptr<Objective<Real>>  obj = makePtr<Objective_FSsolver<Real>>();
      Ptr<Vector<Real>>       l = c.dual().clone();
      ParameterList parlist;
      parlist.sublist("General").set("Output Level",1);
      parlist.sublist("Status Test").set("Constraint Tolerance",ctol);
      parlist.sublist("Status Test").set("Step Tolerance",stol_);
      parlist.sublist("Status Test").set("Iteration Limit",maxit_);
      Ptr<TypeE::Algorithm<Real>> algo = makePtr<TypeE::AugmentedLagrangianAlgorithm<Real>>(parlist);
      algo->run(u,*obj,*con,*l,*stream);
      value(c,u,z,tol);
    }
    if (solverType_ > 4 || solverType_ < 0) {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL:Constraint_SimOpt:solve): Invalid solver type!");
    }
  }

  /** \brief Set solve parameters.

             @param[in]       parlist   ParameterList containing solve parameters

             For the default implementation, parlist has two sublist ("SimOpt"
             and "Solve") and the "Solve" sublist has six input parameters.

                - "Residual Tolerance": Absolute tolerance for the norm of the residual (Real)
                - "Iteration Limit": Maximum number of Newton iterations (int)
                - "Sufficient Decrease Tolerance": Tolerance signifying sufficient decrease in the residual norm, between 0 and 1 (Real)
                - "Step Tolerance": Absolute tolerance for the step size parameter (Real)
                - "Backtracking Factor": Rate for decreasing step size during backtracking, between 0 and 1 (Real)
                - "Output Iteration History": Set to true in order to print solve iteration history (bool)
                - "Zero Initial Guess": Use a vector of zeros as an initial guess for the solve (bool)
                - "Solver Type": Determine which solver to use (0: Newton with line search, 1: Levenberg-Marquardt, 2: SQP) (int)

             These parameters are accessed as parlist.sublist("SimOpt").sublist("Solve").get(...).

             ---
  */
  virtual void setSolveParameters(ParameterList &parlist) {
    ParameterList & list = parlist.sublist("SimOpt").sublist("Solve");
    atol_       = list.get("Absolute Residual Tolerance",   DEFAULT_atol_);
    rtol_       = list.get("Relative Residual Tolerance",   DEFAULT_rtol_);
    maxit_      = list.get("Iteration Limit",               DEFAULT_maxit_);
    decr_       = list.get("Sufficient Decrease Tolerance", DEFAULT_decr_);
    stol_       = list.get("Step Tolerance",                DEFAULT_stol_);
    factor_     = list.get("Backtracking Factor",           DEFAULT_factor_);
    print_      = list.get("Output Iteration History",      DEFAULT_print_);
    zero_       = list.get("Zero Initial Guess",            DEFAULT_zero_);
    solverType_ = list.get("Solver Type",                   DEFAULT_solverType_);
  }

  /** \brief Apply the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z) \in L(\mathcal{U}, \mathcal{C})\f$,
             to the vector \f$v\f$.

             @param[out]      jv  is the result of applying the constraint Jacobian to @b v at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       v   is a simulation-space vector
             @param[in]       u   is the constraint argument; an simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{jv} = c_u(u,z)v\f$, where
             \f$v \in \mathcal{U}\f$, \f$\mathsf{jv} \in \mathcal{C}\f$.

             ---
  */
  virtual void applyJacobian_1(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Update state vector to u + hv
    Ptr<Vector<Real>> unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    // Compute new constraint value
    update(*unew,z,UpdateType::Temp);
    value(jv,*unew,z,ctol);
    // Compute current constraint value
    Ptr<Vector<Real>> cold = jv.clone();
    update(u,z,UpdateType::Temp);
    value(*cold,u,z,ctol);
    // Compute Newton quotient
    jv.axpy(-1.0,*cold);
    jv.scale(1.0/h);
  }


  /** \brief Apply the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_z(u,z) \in L(\mathcal{Z}, \mathcal{C})\f$,
             to the vector \f$v\f$.

             @param[out]      jv  is the result of applying the constraint Jacobian to @b v at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       v   is an optimization-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{jv} = c_z(u,z)v\f$, where
             \f$v \in \mathcal{Z}\f$, \f$\mathsf{jv} \in \mathcal{C}\f$. 

             ---
  */
  virtual void applyJacobian_2(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) { 
    Real ctol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Update state vector to u + hv
    Ptr<Vector<Real>> znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    // Compute new constraint value
    update(u,*znew,UpdateType::Temp);
    value(jv,u,*znew,ctol);
    // Compute current constraint value
    Ptr<Vector<Real>> cold = jv.clone();
    update(u,z,UpdateType::Temp);
    value(*cold,u,z,ctol);
    // Compute Newton quotient
    jv.axpy(-1.0,*cold);
    jv.scale(1.0/h);
  }

  /** \brief Apply the inverse partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^{-1} \in L(\mathcal{C}, \mathcal{U})\f$,
             to the vector \f$v\f$.

             @param[out]      ijv is the result of applying the inverse constraint Jacobian to @b v at @b \f$(u,z)\f$; a simulation-space vector
             @param[in]       v   is a constraint-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ijv} = c_u(u,z)^{-1}v\f$, where
             \f$v \in \mathcal{C}\f$, \f$\mathsf{ijv} \in \mathcal{U}\f$.

             ---
  */
  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseJacobian_1 is used but not implemented!\n");
  }

  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^* \in L(\mathcal{C}^*, \mathcal{U}^*)\f$,
             to the vector \f$v\f$.  This is the primary interface.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b (u,z); a dual simulation-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_u(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    applyAdjointJacobian_1(ajv, v, u, z, v.dual(), tol);
  }


  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^* \in L(\mathcal{C}^*, \mathcal{U}^*)\f$,
             to the vector \f$v\f$.  This is the secondary interface, for use
             with dual spaces where the user does not define the dual() operation.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b (u,z); a dual simulation-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in]       dualv  is a vector used for temporary variables; a constraint-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_u(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      const Vector<Real> &dualv,
                                      Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON<Real>());
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    Ptr<Vector<Real>> cold = dualv.clone();
    Ptr<Vector<Real>> cnew = dualv.clone();
    update(u,z,UpdateType::Temp);
    value(*cold,u,z,ctol);
    Ptr<Vector<Real>> unew = u.clone();
    ajv.zero();
    for (int i = 0; i < u.dimension(); i++) {
      unew->set(u);
      unew->axpy(h,*(u.basis(i)));
      update(*unew,z,UpdateType::Temp);
      value(*cnew,*unew,z,ctol);
      cnew->axpy(-1.0,*cold);
      cnew->scale(1.0/h);
      ajv.axpy(cnew->dot(v),*((u.dual()).basis(i)));
    }
    update(u,z,UpdateType::Temp);
  }


  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_z(u,z)^* \in L(\mathcal{C}^*, \mathcal{Z}^*)\f$,
             to vector \f$v\f$.  This is the primary interface.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b \f$(u,z)\f$; a dual optimization-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_z(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    applyAdjointJacobian_2(ajv, v, u, z, v.dual(), tol);
  }


  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_z(u,z)^* \in L(\mathcal{C}^*, \mathcal{Z}^*)\f$,
             to vector \f$v\f$.  This is the secondary interface, for use
             with dual spaces where the user does not define the dual() operation.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b \f$(u,z)\f$; a dual optimization-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in]       dualv  is a vector used for temporary variables; a constraint-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_z(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      const Vector<Real> &dualv,
                                      Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON<Real>());
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    Ptr<Vector<Real>> cold = dualv.clone();
    Ptr<Vector<Real>> cnew = dualv.clone();
    update(u,z,UpdateType::Temp);
    value(*cold,u,z,ctol);
    Ptr<Vector<Real>> znew = z.clone();
    ajv.zero();
    for (int i = 0; i < z.dimension(); i++) {
      znew->set(z);
      znew->axpy(h,*(z.basis(i)));
      update(u,*znew,UpdateType::Temp);
      value(*cnew,u,*znew,ctol);
      cnew->axpy(-1.0,*cold);
      cnew->scale(1.0/h);
      ajv.axpy(cnew->dot(v),*((z.dual()).basis(i)));
    }
    update(u,z,UpdateType::Temp);
  }

  /** \brief Apply the inverse of the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^{-*} \in L(\mathcal{U}^*, \mathcal{C}^*)\f$,
             to the vector \f$v\f$.

             @param[out]      iajv is the result of applying the inverse adjoint of the constraint Jacobian to @b v at @b (u,z); a dual constraint-space vector
             @param[in]       v   is a dual simulation-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{iajv} = c_u(u,z)^{-*}v\f$, where
             \f$v \in \mathcal{U}^*\f$, \f$\mathsf{iajv} \in \mathcal{C}^*\f$.

             ---
  */
  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) {
    ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseAdjointJacobian_1 is used but not implemented!\n");
  };

  /** \brief Apply the simulation-space derivative of the adjoint of the constraint
             simulation-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{uu}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the simulation-space derivative of the adjoint of the constraint simulation-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual simulation-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{uu}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_11(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new state
    Ptr<Vector<Real>> unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    update(*unew,z,UpdateType::Temp);
    applyAdjointJacobian_1(ahwv,w,*unew,z,jtol);
    // Evaluate Jacobian at old state
    Ptr<Vector<Real>> jv = ahwv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_1(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }


  /** \brief Apply the optimization-space derivative of the adjoint of the constraint
             simulation-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{uz}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the optimization-space derivative of the adjoint of the constraint simulation-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{uz}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_12(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new state
    Ptr<Vector<Real>> unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    update(*unew,z,UpdateType::Temp);
    applyAdjointJacobian_2(ahwv,w,*unew,z,jtol);
    // Evaluate Jacobian at old state
    Ptr<Vector<Real>> jv = ahwv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_2(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }


  /** \brief Apply the simulation-space derivative of the adjoint of the constraint
             optimization-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{zu}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the simulation-space derivative of the adjoint of the constraint optimization-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual simulation-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a optimization-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{zu}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_21(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new control
    Ptr<Vector<Real>> znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    update(u,*znew,UpdateType::Temp);
    applyAdjointJacobian_1(ahwv,w,u,*znew,jtol);
    // Evaluate Jacobian at old control
    Ptr<Vector<Real>> jv = ahwv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_1(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }

  /** \brief Apply the optimization-space derivative of the adjoint of the constraint
             optimization-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{zz}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the optimization-space derivative of the adjoint of the constraint optimization-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a optimization-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{zz}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_22(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new control
    Ptr<Vector<Real>> znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    update(u,*znew,UpdateType::Temp);
    applyAdjointJacobian_2(ahwv,w,u,*znew,jtol);
    // Evaluate Jacobian at old control
    Ptr<Vector<Real>> jv = ahwv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_2(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
}

  /** \brief Approximately solves the <em> augmented system </em>
             \f[
                 \begin{pmatrix}
                   I     & c'(x)^* \\
                   c'(x) & 0
                 \end{pmatrix}
                 \begin{pmatrix}
                   v_{1} \\
                   v_{2}
                 \end{pmatrix}
                 =
                 \begin{pmatrix}
                   b_{1} \\
                   b_{2}
                 \end{pmatrix}
             \f]
             where \f$v_{1} \in \mathcal{X}\f$, \f$v_{2} \in \mathcal{C}^*\f$,
             \f$b_{1} \in \mathcal{X}^*\f$, \f$b_{2} \in \mathcal{C}\f$,
             \f$I : \mathcal{X} \rightarrow \mathcal{X}^*\f$ is an identity
             operator, and \f$0 : \mathcal{C}^* \rightarrow \mathcal{C}\f$
             is a zero operator.

             @param[out]      v1  is the optimization-space component of the result
             @param[out]      v2  is the dual constraint-space component of the result
             @param[in]       b1  is the dual optimization-space component of the right-hand side
             @param[in]       b2  is the constraint-space component of the right-hand side
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is the nominal relative residual tolerance

             On return, \f$ [\mathsf{v1} \,\, \mathsf{v2}] \f$ approximately
             solves the augmented system, where the size of the residual is
             governed by special stopping conditions. \n\n
             The default implementation is the preconditioned generalized
             minimal residual (GMRES) method, which enables the use of
             nonsymmetric preconditioners.

             ---
  */
  virtual std::vector<Real> solveAugmentedSystem(Vector<Real> &v1,
                                                 Vector<Real> &v2,
                                                 const Vector<Real> &b1,
                                                 const Vector<Real> &b2,
                                                 const Vector<Real> &x,
                                                 Real &tol) {
    return Constraint<Real>::solveAugmentedSystem(v1,v2,b1,b2,x,tol);
  }


  /** \brief Apply a constraint preconditioner at \f$x\f$, \f$P(x) \in L(\mathcal{C}, \mathcal{C})\f$,
             to vector \f$v\f$.  In general, this preconditioner satisfies the following relationship:
             \f[
               c'(x) c'(x)^* P(x) v \approx v \,.
             \f]
             It is used by the #solveAugmentedSystem method.

             @param[out]      pv  is the result of applying the constraint preconditioner to @b v at @b x; a constraint-space vector
             @param[in]       v   is a constraint-space vector
             @param[in]       x   is the preconditioner argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations

             On return, \f$\mathsf{pv} = P(x)v\f$, where
             \f$v \in \mathcal{C}\f$, \f$\mathsf{pv} \in \mathcal{C}\f$. \n\n
             The default implementation is a null-op.

             ---
  */
  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   const Vector<Real> &g,
                                   Real &tol) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(x);
    Ptr<Vector<Real>> ijv = (xs.get_1())->clone();

    try {
      applyInverseJacobian_1(*ijv, v, *(xs.get_1()), *(xs.get_2()), tol);
    }
    catch (const std::logic_error &e) {
      Constraint<Real>::applyPreconditioner(pv, v, x, g, tol);
      return;
    }

    const Vector_SimOpt<Real> &gs = dynamic_cast<const Vector_SimOpt<Real>&>(g);
    Ptr<Vector<Real>> ijv_dual = (gs.get_1())->clone();
    ijv_dual->set(ijv->dual());

    try {
      applyInverseAdjointJacobian_1(pv, *ijv_dual, *(xs.get_1()), *(xs.get_2()), tol);
    }
    catch (const std::logic_error &e) {
      Constraint<Real>::applyPreconditioner(pv, v, x, g, tol);
      return;
    }

  }

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    update(*(xs.get_1()),*(xs.get_2()),flag,iter);
  }
  virtual void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    update(*(xs.get_1()),*(xs.get_2()),type,iter);
  }

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &x,
                     Real &tol) {
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    value(c,*(xs.get_1()),*(xs.get_2()),tol);
  }


  virtual void applyJacobian(Vector<Real> &jv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol) { 
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    const Vector_SimOpt<Real> &vs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(v));
    applyJacobian_1(jv,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    Ptr<Vector<Real>> jv2 = jv.clone();
    applyJacobian_2(*jv2,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    jv.plus(*jv2);
  }


  using Constraint<Real>::applyAdjointJacobian;
  virtual void applyAdjointJacobian(Vector<Real> &ajv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    Real &tol) { 
    Vector_SimOpt<Real> &ajvs = dynamic_cast<Vector_SimOpt<Real>&>(
      dynamic_cast<Vector<Real>&>(ajv));
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    Ptr<Vector<Real>> ajv1 = (ajvs.get_1())->clone();
    applyAdjointJacobian_1(*ajv1,v,*(xs.get_1()),*(xs.get_2()),tol);
    ajvs.set_1(*ajv1);
    Ptr<Vector<Real>> ajv2 = (ajvs.get_2())->clone();
    applyAdjointJacobian_2(*ajv2,v,*(xs.get_1()),*(xs.get_2()),tol);
    ajvs.set_2(*ajv2);
  }


  virtual void applyAdjointHessian(Vector<Real> &ahwv,
                                   const Vector<Real> &w,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    Vector_SimOpt<Real> &ahwvs = dynamic_cast<Vector_SimOpt<Real>&>(
      dynamic_cast<Vector<Real>&>(ahwv));
    const Vector_SimOpt<Real> &xs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(x));
    const Vector_SimOpt<Real> &vs = dynamic_cast<const Vector_SimOpt<Real>&>(
      dynamic_cast<const Vector<Real>&>(v));
    // Block-row 1
    Ptr<Vector<Real>> C11 = (ahwvs.get_1())->clone();
    Ptr<Vector<Real>> C21 = (ahwvs.get_1())->clone();
    applyAdjointHessian_11(*C11,w,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    applyAdjointHessian_21(*C21,w,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    C11->plus(*C21);
    ahwvs.set_1(*C11); 
    // Block-row 2
    Ptr<Vector<Real>> C12 = (ahwvs.get_2())->clone();
    Ptr<Vector<Real>> C22 = (ahwvs.get_2())->clone();
    applyAdjointHessian_12(*C12,w,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    applyAdjointHessian_22(*C22,w,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    C22->plus(*C12);
    ahwvs.set_2(*C22); 
  }



  virtual Real checkSolve(const Vector<Real> &u, 
                          const Vector<Real> &z, 
                          const Vector<Real> &c,
                          const bool printToStream = true,
                          std::ostream & outStream = std::cout) {
    // Solve constraint for u. 
    Real tol = ROL_EPSILON<Real>();
    Ptr<Vector<Real>> r = c.clone();
    Ptr<Vector<Real>> s = u.clone();
    solve(*r,*s,z,tol);
    // Evaluate constraint residual at (u,z).
    Ptr<Vector<Real>> cs = c.clone();
    update(*s,z,UpdateType::Temp);
    value(*cs,*s,z,tol);
    // Output norm of residual.
    Real rnorm = r->norm();
    Real cnorm = cs->norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt solve at feasible (u,z):\n";
      hist << "  Solver Residual = " << rnorm << "\n";
      hist << "       ||c(u,z)|| = " << cnorm << "\n";
      outStream << hist.str();
    }
    return cnorm;
  }


  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the primary interface.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is a simulation-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_1(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    return checkAdjointConsistencyJacobian_1(w, v, u, z, w.dual(), v.dual(), printToStream, outStream);
  }


  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the secondary interface, for use with dual spaces where
             the user does not define the dual() operation.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is a simulation-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       dualw          is a constraint-space vector 
             @param[in]       dualv          is a dual simulation-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_1(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const Vector<Real> &dualw,
                                                 const Vector<Real> &dualv,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON<Real>();
    Ptr<Vector<Real>> Jv = dualw.clone();
    update(u,z,UpdateType::Temp);
    applyJacobian_1(*Jv,v,u,z,tol);
    //Real wJv = w.dot(Jv->dual());
    Real wJv = w.apply(*Jv);
    Ptr<Vector<Real>> Jw = dualv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_1(*Jw,w,u,z,tol);
    //Real vJw = v.dot(Jw->dual());
    Real vJw = v.apply(*Jw);
    Real diff = std::abs(wJv-vJw);
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of Jacobian_1 and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = " 
           << diff << "\n";
      hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
      hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW<Real>()) << "\n";
      outStream << hist.str();
    }
    return diff;
  }


  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the primary interface.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is an optimization-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_2(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    return checkAdjointConsistencyJacobian_2(w, v, u, z, w.dual(), v.dual(), printToStream, outStream);
  }

  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the secondary interface, for use with dual spaces where
             the user does not define the dual() operation.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is an optimization-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       dualw          is a constraint-space vector 
             @param[in]       dualv          is a dual optimization-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_2(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const Vector<Real> &dualw,
                                                 const Vector<Real> &dualv,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON<Real>();
    Ptr<Vector<Real>> Jv = dualw.clone();
    update(u,z,UpdateType::Temp);
    applyJacobian_2(*Jv,v,u,z,tol);
    //Real wJv = w.dot(Jv->dual());
    Real wJv = w.apply(*Jv);
    Ptr<Vector<Real>> Jw = dualv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_2(*Jw,w,u,z,tol);
    //Real vJw = v.dot(Jw->dual());
    Real vJw = v.apply(*Jw);
    Real diff = std::abs(wJv-vJw);
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of Jacobian_2 and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = "
           << diff << "\n";
      hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
      hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW<Real>()) << "\n";
      outStream << hist.str();
    }
    return diff;
  }

  virtual Real checkInverseJacobian_1(const Vector<Real> &jv, 
                                      const Vector<Real> &v, 
                                      const Vector<Real> &u, 
                                      const Vector<Real> &z, 
                                      const bool printToStream = true,
                                      std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON<Real>();
    Ptr<Vector<Real>> Jv = jv.clone();
    update(u,z,UpdateType::Temp);
    applyJacobian_1(*Jv,v,u,z,tol);
    Ptr<Vector<Real>> iJJv = u.clone();
    //update(u,z); // Does this update do anything?
    applyInverseJacobian_1(*iJJv,*Jv,u,z,tol);
    Ptr<Vector<Real>> diff = v.clone();
    diff->set(v);
    diff->axpy(-1.0,*iJJv);
    Real dnorm = diff->norm();
    Real vnorm = v.norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of inverse Jacobian_1: \n  ||v-inv(J)Jv|| = " 
           << dnorm << "\n";
      hist << "  ||v||          = " << vnorm << "\n";
      hist << "  Relative Error = " << dnorm / (vnorm+ROL_UNDERFLOW<Real>()) << "\n";
      outStream << hist.str();
    }
    return dnorm;
  }

  virtual Real checkInverseAdjointJacobian_1(const Vector<Real> &jv, 
                                             const Vector<Real> &v, 
                                             const Vector<Real> &u, 
                                             const Vector<Real> &z, 
                                             const bool printToStream = true,
                                             std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON<Real>();
    Ptr<Vector<Real>> Jv = jv.clone();
    update(u,z,UpdateType::Temp);
    applyAdjointJacobian_1(*Jv,v,u,z,tol);
    Ptr<Vector<Real>> iJJv = v.clone();
    //update(u,z);
    applyInverseAdjointJacobian_1(*iJJv,*Jv,u,z,tol);
    Ptr<Vector<Real>> diff = v.clone();
    diff->set(v);
    diff->axpy(-1.0,*iJJv);
    Real dnorm = diff->norm();
    Real vnorm = v.norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of inverse adjoint Jacobian_1: \n  ||v-inv(adj(J))adj(J)v|| = "
           << dnorm << "\n";
      hist << "  ||v||                    = " << vnorm << "\n";
      hist << "  Relative Error           = " << dnorm / (vnorm+ROL_UNDERFLOW<Real>()) << "\n";
      outStream << hist.str();
    }
    return dnorm;
  }



  std::vector<std::vector<Real>> checkApplyJacobian_1(const Vector<Real> &u,
                                                       const Vector<Real> &z,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &jv,
                                                       const bool printToStream = true,
                                                       std::ostream & outStream = std::cout,
                                                       const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                       const int order = 1) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyJacobian_1(u,z,v,jv,steps,printToStream,outStream,order);
  }
  
  
  
  
  std::vector<std::vector<Real>> checkApplyJacobian_1(const Vector<Real> &u,
                                                       const Vector<Real> &z,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &jv,
                                                       const std::vector<Real> &steps, 
                                                       const bool printToStream = true,
                                                       std::ostream & outStream = std::cout,
                                                       const int order = 1) {
 
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                                "Error: finite difference order must be 1,2,3, or 4" );
 
    Real one(1.0);
 
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
 
    Real tol = std::sqrt(ROL_EPSILON<Real>());
 
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> jvCheck(numSteps, tmp);
 
    // Save the format state of the original outStream.
    nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
 
    // Compute constraint value at x.
    Ptr<Vector<Real>> c = jv.clone();
    this->update(u,z,UpdateType::Temp);
    this->value(*c, u, z, tol);
 
    // Compute (Jacobian at x) times (vector v).
    Ptr<Vector<Real>> Jv = jv.clone();
    this->applyJacobian_1(*Jv, v, u, z, tol);
    Real normJv = Jv->norm();
 
    // Temporary vectors.
    Ptr<Vector<Real>> cdif = jv.clone();
    Ptr<Vector<Real>> cnew = jv.clone();
    Ptr<Vector<Real>> unew = u.clone();
 
    for (int i=0; i<numSteps; i++) {
 
      Real eta = steps[i];
 
      unew->set(u);
 
      cdif->set(*c);
      cdif->scale(weights[order-1][0]);
 
      for(int j=0; j<order; ++j) {
 
         unew->axpy(eta*shifts[order-1][j], v);

         if( weights[order-1][j+1] != 0 ) {
             this->update(*unew,z,UpdateType::Temp);
             this->value(*cnew,*unew,z,tol);
             cdif->axpy(weights[order-1][j+1],*cnew);
         }

      }
 
      cdif->scale(one/eta);
 
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      jvCheck[i][0] = eta;
      jvCheck[i][1] = normJv;
      jvCheck[i][2] = cdif->norm();
      cdif->axpy(-one, *Jv);
      jvCheck[i][3] = cdif->norm();
 
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(Jac*vec)"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << jvCheck[i][0]
             << std::setw(20) << jvCheck[i][1]
             << std::setw(20) << jvCheck[i][2]
             << std::setw(20) << jvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
 
    }
 
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
 
    return jvCheck;
  } // checkApplyJacobian


  std::vector<std::vector<Real>> checkApplyJacobian_2(const Vector<Real> &u,
                                                       const Vector<Real> &z,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &jv,
                                                       const bool printToStream = true,
                                                       std::ostream & outStream = std::cout,
                                                       const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                       const int order = 1) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyJacobian_2(u,z,v,jv,steps,printToStream,outStream,order);
  }
  
  
  
  
  std::vector<std::vector<Real>> checkApplyJacobian_2(const Vector<Real> &u,
                                                       const Vector<Real> &z,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &jv,
                                                       const std::vector<Real> &steps, 
                                                       const bool printToStream = true,
                                                       std::ostream & outStream = std::cout,
                                                       const int order = 1) {
 
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                                "Error: finite difference order must be 1,2,3, or 4" );
 
    Real one(1.0);
 
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
 
    Real tol = std::sqrt(ROL_EPSILON<Real>());
 
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> jvCheck(numSteps, tmp);
 
    // Save the format state of the original outStream.
    nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
 
    // Compute constraint value at x.
    Ptr<Vector<Real>> c = jv.clone();
    this->update(u,z,UpdateType::Temp);
    this->value(*c, u, z, tol);
 
    // Compute (Jacobian at x) times (vector v).
    Ptr<Vector<Real>> Jv = jv.clone();
    this->applyJacobian_2(*Jv, v, u, z, tol);
    Real normJv = Jv->norm();
 
    // Temporary vectors.
    Ptr<Vector<Real>> cdif = jv.clone();
    Ptr<Vector<Real>> cnew = jv.clone();
    Ptr<Vector<Real>> znew = z.clone();
 
    for (int i=0; i<numSteps; i++) {
 
      Real eta = steps[i];
 
      znew->set(z);
 
      cdif->set(*c);
      cdif->scale(weights[order-1][0]);
 
      for(int j=0; j<order; ++j) {
 
         znew->axpy(eta*shifts[order-1][j], v);

         if( weights[order-1][j+1] != 0 ) {
             this->update(u,*znew,UpdateType::Temp);
             this->value(*cnew,u,*znew,tol);
             cdif->axpy(weights[order-1][j+1],*cnew);
         }

      }
 
      cdif->scale(one/eta);
 
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      jvCheck[i][0] = eta;
      jvCheck[i][1] = normJv;
      jvCheck[i][2] = cdif->norm();
      cdif->axpy(-one, *Jv);
      jvCheck[i][3] = cdif->norm();
 
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(Jac*vec)"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << jvCheck[i][0]
             << std::setw(20) << jvCheck[i][1]
             << std::setw(20) << jvCheck[i][2]
             << std::setw(20) << jvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
 
    }
 
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
 
    return jvCheck;
  } // checkApplyJacobian



  std::vector<std::vector<Real>> checkApplyAdjointHessian_11(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                              const int order = 1 ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyAdjointHessian_11(u,z,p,v,hv,steps,printToStream,outStream,order);
  
  }
  
  std::vector<std::vector<Real>> checkApplyAdjointHessian_11(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const std::vector<Real> &steps,  
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int order = 1 ) {
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
    Real one(1.0);
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> ahpvCheck(numSteps, tmp);
  
    // Temporary vectors.
    Ptr<Vector<Real>> AJdif = hv.clone();
    Ptr<Vector<Real>> AJp = hv.clone();
    Ptr<Vector<Real>> AHpv = hv.clone();
    Ptr<Vector<Real>> AJnew = hv.clone();
    Ptr<Vector<Real>> unew = u.clone();
  
    // Save the format state of the original outStream.
    nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Apply adjoint Jacobian to p.
    this->update(u,z,UpdateType::Temp);
    this->applyAdjointJacobian_1(*AJp, p, u, z, tol);
  
    // Apply adjoint Hessian at (u,z), in direction v, to p.
    this->applyAdjointHessian_11(*AHpv, p, v, u, z, tol);
    Real normAHpv = AHpv->norm();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Apply adjoint Jacobian to p at (u+eta*v,z).
      unew->set(u);
  
      AJdif->set(*AJp);
      AJdif->scale(weights[order-1][0]);     
  
      for(int j=0; j<order; ++j) {
  
          unew->axpy(eta*shifts[order-1][j],v); 
  
          if( weights[order-1][j+1] != 0 ) {    
              this->update(*unew,z,UpdateType::Temp);
              this->applyAdjointJacobian_1(*AJnew, p, *unew, z, tol);
              AJdif->axpy(weights[order-1][j+1],*AJnew);
          }
      }
  
      AJdif->scale(one/eta);
  
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      ahpvCheck[i][0] = eta;
      ahpvCheck[i][1] = normAHpv;
      ahpvCheck[i][2] = AJdif->norm();
      AJdif->axpy(-one, *AHpv);
      ahpvCheck[i][3] = AJdif->norm();
  
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(adj(H)(u,v))"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-----------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << ahpvCheck[i][0]
             << std::setw(20) << ahpvCheck[i][1]
             << std::setw(20) << ahpvCheck[i][2]
             << std::setw(20) << ahpvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return ahpvCheck;
  } // checkApplyAdjointHessian_11

  /** 
     \brief \f$ u\in U \f$, \f$ z\in Z \f$, \f$ p\in C^\ast \f$, \f$ v \in U \f$, \f$ hv \in U^\ast \f$ 
  */
  std::vector<std::vector<Real>> checkApplyAdjointHessian_21(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                              const int order = 1 ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyAdjointHessian_21(u,z,p,v,hv,steps,printToStream,outStream,order);
  
  }
  
  /** 
     \brief \f$ u\in U \f$, \f$ z\in Z \f$, \f$ p\in C^\ast \f$, \f$ v \in U \f$, \f$ hv \in U^\ast \f$ 
  */
  std::vector<std::vector<Real>> checkApplyAdjointHessian_21(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const std::vector<Real> &steps,  
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int order = 1 ) {
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
    Real one(1.0);
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> ahpvCheck(numSteps, tmp);
  
    // Temporary vectors.
    Ptr<Vector<Real>> AJdif = hv.clone();
    Ptr<Vector<Real>> AJp = hv.clone();
    Ptr<Vector<Real>> AHpv = hv.clone();
    Ptr<Vector<Real>> AJnew = hv.clone();
    Ptr<Vector<Real>> znew = z.clone();
  
    // Save the format state of the original outStream.
    nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Apply adjoint Jacobian to p.
    this->update(u,z,UpdateType::Temp);
    this->applyAdjointJacobian_1(*AJp, p, u, z, tol);
  
    // Apply adjoint Hessian at (u,z), in direction v, to p.
    this->applyAdjointHessian_21(*AHpv, p, v, u, z, tol);
    Real normAHpv = AHpv->norm();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Apply adjoint Jacobian to p at (u,z+eta*v).
      znew->set(z);
  
      AJdif->set(*AJp);
      AJdif->scale(weights[order-1][0]);     
  
      for(int j=0; j<order; ++j) {
  
          znew->axpy(eta*shifts[order-1][j],v); 
  
          if( weights[order-1][j+1] != 0 ) {    
              this->update(u,*znew,UpdateType::Temp);
              this->applyAdjointJacobian_1(*AJnew, p, u, *znew, tol);
              AJdif->axpy(weights[order-1][j+1],*AJnew);
          }
      }
  
      AJdif->scale(one/eta);
  
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      ahpvCheck[i][0] = eta;
      ahpvCheck[i][1] = normAHpv;
      ahpvCheck[i][2] = AJdif->norm();
      AJdif->axpy(-one, *AHpv);
      ahpvCheck[i][3] = AJdif->norm();
  
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(adj(H)(u,v))"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-----------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << ahpvCheck[i][0]
             << std::setw(20) << ahpvCheck[i][1]
             << std::setw(20) << ahpvCheck[i][2]
             << std::setw(20) << ahpvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return ahpvCheck;
  } // checkApplyAdjointHessian_21

  /** 
     \brief \f$ u\in U \f$, \f$ z\in Z \f$, \f$ p\in C^\ast \f$, \f$ v \in U \f$, \f$ hv \in U^\ast \f$ 
  */
  std::vector<std::vector<Real>> checkApplyAdjointHessian_12(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                              const int order = 1 ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyAdjointHessian_12(u,z,p,v,hv,steps,printToStream,outStream,order);
  
  }
  

  std::vector<std::vector<Real>> checkApplyAdjointHessian_12(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const std::vector<Real> &steps,  
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int order = 1 ) {
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
    Real one(1.0);
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> ahpvCheck(numSteps, tmp);
  
    // Temporary vectors.
    Ptr<Vector<Real>> AJdif = hv.clone();
    Ptr<Vector<Real>> AJp = hv.clone();
    Ptr<Vector<Real>> AHpv = hv.clone();
    Ptr<Vector<Real>> AJnew = hv.clone();
    Ptr<Vector<Real>> unew = u.clone();
  
    // Save the format state of the original outStream.
    nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Apply adjoint Jacobian to p.
    this->update(u,z,UpdateType::Temp);
    this->applyAdjointJacobian_2(*AJp, p, u, z, tol);
  
    // Apply adjoint Hessian at (u,z), in direction v, to p.
    this->applyAdjointHessian_12(*AHpv, p, v, u, z, tol);
    Real normAHpv = AHpv->norm();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Apply adjoint Jacobian to p at (u+eta*v,z).
      unew->set(u);
  
      AJdif->set(*AJp);
      AJdif->scale(weights[order-1][0]);     
  
      for(int j=0; j<order; ++j) {
  
          unew->axpy(eta*shifts[order-1][j],v); 
  
          if( weights[order-1][j+1] != 0 ) {    
              this->update(*unew,z,UpdateType::Temp);
              this->applyAdjointJacobian_2(*AJnew, p, *unew, z, tol);
              AJdif->axpy(weights[order-1][j+1],*AJnew);
          }
      }
  
      AJdif->scale(one/eta);
  
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      ahpvCheck[i][0] = eta;
      ahpvCheck[i][1] = normAHpv;
      ahpvCheck[i][2] = AJdif->norm();
      AJdif->axpy(-one, *AHpv);
      ahpvCheck[i][3] = AJdif->norm();
  
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(adj(H)(u,v))"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-----------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << ahpvCheck[i][0]
             << std::setw(20) << ahpvCheck[i][1]
             << std::setw(20) << ahpvCheck[i][2]
             << std::setw(20) << ahpvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return ahpvCheck;
  } // checkApplyAdjointHessian_12

  std::vector<std::vector<Real>> checkApplyAdjointHessian_22(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                              const int order = 1 ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyAdjointHessian_22(u,z,p,v,hv,steps,printToStream,outStream,order);
  
  }
  
  std::vector<std::vector<Real>> checkApplyAdjointHessian_22(const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              const Vector<Real> &p,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &hv,
                                                              const std::vector<Real> &steps,  
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int order = 1 ) {
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
    Real one(1.0);
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> ahpvCheck(numSteps, tmp);
  
    // Temporary vectors.
    Ptr<Vector<Real>> AJdif = hv.clone();
    Ptr<Vector<Real>> AJp = hv.clone();
    Ptr<Vector<Real>> AHpv = hv.clone();
    Ptr<Vector<Real>> AJnew = hv.clone();
    Ptr<Vector<Real>> znew = z.clone();
  
    // Save the format state of the original outStream.
    nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Apply adjoint Jacobian to p.
    this->update(u,z,UpdateType::Temp);
    this->applyAdjointJacobian_2(*AJp, p, u, z, tol);
  
    // Apply adjoint Hessian at (u,z), in direction v, to p.
    this->applyAdjointHessian_22(*AHpv, p, v, u, z, tol);
    Real normAHpv = AHpv->norm();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Apply adjoint Jacobian to p at (u,z+eta*v).
      znew->set(z);
  
      AJdif->set(*AJp);
      AJdif->scale(weights[order-1][0]);     
  
      for(int j=0; j<order; ++j) {
  
          znew->axpy(eta*shifts[order-1][j],v); 
  
          if( weights[order-1][j+1] != 0 ) {    
              this->update(u,*znew,UpdateType::Temp);
              this->applyAdjointJacobian_2(*AJnew, p, u, *znew, tol);
              AJdif->axpy(weights[order-1][j+1],*AJnew);
          }
      }
  
      AJdif->scale(one/eta);
  
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      ahpvCheck[i][0] = eta;
      ahpvCheck[i][1] = normAHpv;
      ahpvCheck[i][2] = AJdif->norm();
      AJdif->axpy(-one, *AHpv);
      ahpvCheck[i][3] = AJdif->norm();
  
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(adj(H)(u,v))"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-----------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << ahpvCheck[i][0]
             << std::setw(20) << ahpvCheck[i][1]
             << std::setw(20) << ahpvCheck[i][2]
             << std::setw(20) << ahpvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return ahpvCheck;
  } // checkApplyAdjointHessian_22

}; // class Constraint_SimOpt

} // namespace ROL

#endif
