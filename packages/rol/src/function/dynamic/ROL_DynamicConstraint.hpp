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


#pragma once
#ifndef ROL_DYNAMICCONSTRAINT_HPP
#define ROL_DYNAMICCONSTRAINT_HPP

#include "ROL_DynamicFunction.hpp"
#include "ROL_TimeStamp.hpp"

#include "ROL_NonlinearLeastSquaresObjective_Dynamic.hpp"
#include "ROL_Constraint_DynamicState.hpp"
#include "ROL_Objective_FSsolver.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Types.hpp"

/** @ingroup dynamic_group
    \class ROL::DynamicConstraint
    \brief Defines the time-dependent constraint operator interface for 
           simulation-based optimization.

    This constraint interface inherits from ROL_Constraint_SimOpt. Though 
    the interface takes two simulation space vectors from spaces
    \f$\mathcal{U_o}\times\mathcal{U_n}\f$. The space \f$\mathcal{U_o}\f$ is 
    ``old'' information that accounts for the initial condition on the time 
     interval. The space \f$\mathcal{U_n}\f$ is the ``new'' variables that can 
    be determined by satisfying constraints in the form
    \f[
      c(u_o,u_n,z,t_o,t_n) = 0 \,.
    \f]

*/


namespace ROL {

template<typename Real>
class DynamicConstraint : public DynamicFunction<Real> {
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

  using V  = Vector<Real>;
  using PV = PartitionedVector<Real>;
  using TS = TimeStamp<Real>;
 
  virtual ~DynamicConstraint() {}

  DynamicConstraint(void)
    : unew_               (                             nullPtr ),
      jv_                 (                             nullPtr ),
      DEFAULT_atol_       ( 1e-4*std::sqrt(ROL_EPSILON<Real>()) ),
      DEFAULT_rtol_       (                                 1e0 ),
      DEFAULT_stol_       (      std::sqrt(ROL_EPSILON<Real>()) ),
      DEFAULT_factor_     (                                 0.5 ),
      DEFAULT_decr_       (                                1e-4 ),
      DEFAULT_maxit_      (                                 500 ),
      DEFAULT_print_      (                               false ),
      DEFAULT_zero_       (                               false ),
      DEFAULT_solverType_ (                                   0 ),
      atol_               (                       DEFAULT_atol_ ),
      rtol_               (                       DEFAULT_rtol_ ),
      stol_               (                       DEFAULT_stol_ ),
      factor_             (                     DEFAULT_factor_ ),
      decr_               (                       DEFAULT_decr_ ),
      maxit_              (                      DEFAULT_maxit_ ),
      print_              (                      DEFAULT_print_ ),
      zero_               (                       DEFAULT_zero_ ),
      solverType_         (                 DEFAULT_solverType_ ),
      firstSolve_         (                                true ) {}

  virtual void update( const V& uo, const V& un, const V& z, const TS& ts ) {
    update_uo( uo, ts );
    update_un( un, ts );
    update_z( z, ts );
  }

  virtual void update_uo( const V& uo, const TS& ts ) { }
  virtual void update_un( const V& un, const TS& ts ) { }
  virtual void update_z( const V& z, const TS& ts ) { }

  virtual void value( V& c, const V& uo, const V& un, 
                      const V& z, const TS& ts ) const = 0;

  virtual void solve( V& c, const V& uo, V& un, 
                      const V& z, const TS& ts ) {
    if ( zero_ ) {
      un.zero();
    }
    update(uo,un,z,ts);
    value(c,uo,un,z,ts);
    Real cnorm = c.norm();
    const Real ctol = std::min(atol_, rtol_*cnorm);
    if (solverType_==0 || solverType_==3 || solverType_==4) {
      if ( firstSolve_ ) {
        unew_ = un.clone();
        jv_   = un.clone();
        firstSolve_ = false;
      }
      Real alpha(1), tmp(0);
      int cnt = 0;
      if ( print_ ) {
        std::cout << "\n     Default DynamicConstraint::solve\n";
        std::cout << "       ";
        std::cout << std::setw(6)  << std::left << "iter";
        std::cout << std::setw(15) << std::left << "rnorm";
        std::cout << std::setw(15) << std::left << "alpha";
        std::cout << "\n";
      }
      for (cnt = 0; cnt < maxit_; ++cnt) {
        // Compute Newton step
        applyInverseJacobian_un(*jv_,c,uo,un,z,ts);
        unew_->set(un);
        unew_->axpy(-alpha, *jv_);
        //update_un(*unew_,ts);
        update(uo,*unew_,z,ts);
        value(c,uo,*unew_,z,ts);
        tmp = c.norm();
        // Perform backtracking line search
        while ( tmp > (1.0-decr_*alpha)*cnorm &&
                alpha > stol_ ) {
          alpha *= factor_;
          unew_->set(un);
          unew_->axpy(-alpha,*jv_);
          //update_un(*unew_,ts);
          update(uo,*unew_,z,ts);
          value(c,uo,*unew_,z,ts);
          tmp = c.norm();
        }
        if ( print_ ) {
          std::cout << "       ";
          std::cout << std::setw(6)  << std::left << cnt;
          std::cout << std::scientific << std::setprecision(6);
          std::cout << std::setw(15) << std::left << tmp;
          std::cout << std::scientific << std::setprecision(6);
          std::cout << std::setw(15) << std::left << alpha;
          std::cout << "\n";
        }
        // Update iterate
        cnorm = tmp;
        un.set(*unew_);
        if (cnorm <= ctol) { // = covers the case of identically zero residual
          break;
        }
        update(uo,un,z,ts);
        alpha = 1.0;
      }
    }
    if (solverType_==1 || (solverType_==3 && cnorm > ctol)) {
      ROL::Ptr<DynamicConstraint<Real>> con_ptr = ROL::makePtrFromRef(*this);
      ROL::Ptr<const Vector<Real>>       uo_ptr = ROL::makePtrFromRef(uo);
      ROL::Ptr<const Vector<Real>>        z_ptr = ROL::makePtrFromRef(z);
      ROL::Ptr<const TimeStamp<Real>>    ts_ptr = ROL::makePtrFromRef(ts);
      ROL::Ptr<Objective<Real>>         obj_ptr
        = ROL::makePtr<NonlinearLeastSquaresObjective_Dynamic<Real>>(con_ptr,c,uo_ptr,z_ptr,ts_ptr,true);
      ROL::ParameterList parlist;
      parlist.sublist("Status Test").set("Gradient Tolerance",ctol);
      parlist.sublist("Status Test").set("Step Tolerance",stol_);
      parlist.sublist("Status Test").set("Iteration Limit",maxit_);
      parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      parlist.sublist("General").sublist("Krylov").set("Iteration Limit",100);
      ROL::Ptr<Algorithm<Real>> algo = ROL::makePtr<Algorithm<Real>>("Trust Region",parlist,false);
      algo->run(un,*obj_ptr,print_);
      value(c,uo,un,z,ts);
    }
    if (solverType_==2 || (solverType_==4 && cnorm > ctol)) {
      ROL::Ptr<DynamicConstraint<Real>> con_ptr = ROL::makePtrFromRef(*this);
      ROL::Ptr<const Vector<Real>>       uo_ptr = ROL::makePtrFromRef(uo);
      ROL::Ptr<const Vector<Real>>        z_ptr = ROL::makePtrFromRef(z);
      ROL::Ptr<const TimeStamp<Real>>    ts_ptr = ROL::makePtrFromRef(ts);
      ROL::Ptr<Constraint<Real>>         cn_ptr
        = ROL::makePtr<Constraint_DynamicState<Real>>(con_ptr,uo_ptr,z_ptr,ts_ptr);
      ROL::Ptr<Objective<Real>>         obj_ptr
        = ROL::makePtr<Objective_FSsolver<Real>>();
      ROL::ParameterList parlist;
      parlist.sublist("Status Test").set("Constraint Tolerance",ctol);
      parlist.sublist("Status Test").set("Step Tolerance",stol_);
      parlist.sublist("Status Test").set("Iteration Limit",maxit_);
      ROL::Ptr<Algorithm<Real>> algo = ROL::makePtr<Algorithm<Real>>("Composite Step",parlist,false);
      ROL::Ptr<Vector<Real>> l = c.dual().clone();
      algo->run(un,*l,*obj_ptr,*cn_ptr,print_);
      value(c,uo,un,z,ts);
    }
    if (solverType_ > 4 || solverType_ < 0) {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL:DynamicConstraint:solve): Invalid solver type!");
    }
  }

  /** \brief Set solve parameters.

             @param[in]       parlist   ROL::ParameterList containing solve parameters

             For the default implementation, parlist has two sublist ("Dynamic Constraint"
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
  virtual void setSolveParameters(ROL::ParameterList &parlist) {
    ROL::ParameterList & list = parlist.sublist("Dynamic Constraint").sublist("Solve");
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

  //----------------------------------------------------------------------------
  // Partial Jacobians
  virtual void applyJacobian_uo( V& jv, const V& vo, const V& uo, 
                                    const V& un, const V& z, 
                                    const TS& ts ) const {}

  virtual void applyJacobian_un( V& jv, const V& vn, const V& uo, 
                                    const V& un, const V& z, 
                                    const TS& ts ) const {}

  virtual void applyJacobian_z( V& jv, const V& vz, const V& uo, 
                                const V& un, const V& z, 
                                const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Adjoint partial Jacobians
  virtual void applyAdjointJacobian_uo( V& ajv, const V& dualv, const V& uo, 
                                           const V& un, const V& z, 
                                           const TS& ts ) const {}

  virtual void applyAdjointJacobian_un( V& ajv, const V& dualv, const V& uo, 
                                           const V& un, const V& z, 
                                           const TS& ts ) const {}

  virtual void applyAdjointJacobian_z( V& ajv, const V& dualv, const V& uo, 
                                       const V& un, const V& z, 
                                       const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Inverses
  virtual void applyInverseJacobian_un( V& ijv, const V& vn, const V& uo, 
                                           const V& un, const V& z, 
                                           const TS& ts ) const {}
    
  virtual void applyInverseAdjointJacobian_un( V& iajv, const V& vn, const V& uo, 
                                                  const V& un, const V& z, 
                                                  const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Adjoint Hessian components
  virtual void applyAdjointHessian_un_un( V& ahwv, const V& wn, const V& vn,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_un_uo( V& ahwv, const V& w, const V& vn,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_un_z( V& ahwv, const V& w, const V& vn,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_uo_un( V& ahwv, const V& w, const V& vo,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_uo_uo( V& ahwv, const V& w, const V& v,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_uo_z( V& ahwv, const V& w, const V& vo,
                                         const V& uo, const V& un, 
                                         const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_z_un( V& ahwv, const V& w, const V& vz,
                                         const V& uo, const V& un, 
                                         const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_z_uo( V& ahwv, const V& w, const V& vz,
                                         const V& uo, const V& un, 
                                         const V& z, const TS& ts ) const {
    ahwv.zero();
  }
  
  virtual void applyAdjointHessian_z_z( V& ahwv, const V& w, const V& vz,
                                        const V& uo, const V& un, 
                                        const V& z, const TS& ts ) const {
    ahwv.zero();
  }

}; // DynamicConstraint


} // namespace ROL


#endif // ROL_DYNAMICCONSTRAINT_HPP

