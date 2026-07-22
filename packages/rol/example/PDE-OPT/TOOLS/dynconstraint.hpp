// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  constraint.hpp
    \brief Defines the SimOpt constraint for the 'poisson' example.
*/

#ifndef ROL_PDEOPT_DYNCONSTRAINT_H
#define ROL_PDEOPT_DYNCONSTRAINT_H

#include "ROL_DynamicConstraint.hpp"
#include "dynpde.hpp"
#include "assembler.hpp"
#include "solver.hpp"
#include "pdevector.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> DynConstraintSolverConstruct_Jacobian_un        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Solver Construction Time for Jacobian un");
    ROL::Ptr<Teuchos::Time> DynConstraintSolverSolve_Jacobian_un            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Solver Solution Time for Jacobian un");
    ROL::Ptr<Teuchos::Time> DynConstraintSolverSolve_AdjointJacobian_un     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Solver Solution Time for Adjoint Jacobian un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyJacobian_uo                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Jacobian uo");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyJacobian_un                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Jacobian un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyJacobian_zf                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Jacobian zf");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyJacobian_zp                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Jacobian zp");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyAdjointJacobian_uo            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Adjoint Jacobian uo");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyAdjointJacobian_un            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Adjoint Jacobian un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyAdjointJacobian_zf            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Adjoint Jacobian zf");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyAdjointJacobian_zp            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Adjoint Jacobian zp");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_uo_uo                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_uo_uo");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_uo_un                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_uo_un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_un_uo                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_un_uo");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_un_un                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_un_un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_uo_zf                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_uo_zf");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_uo_zp                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_uo_zp");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zf_uo                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zf_uo");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zp_uo                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zp_uo");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_un_zf                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_un_zf");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_un_zp                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_un_zp");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zf_un                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zf_un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zp_un                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zp_un");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zp_zp                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zp_zp");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zp_zf                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zp_zf");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zf_zp                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zf_zp");
    ROL::Ptr<Teuchos::Time> DynConstraintApplyHessian_zf_zf                 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Dyn Constraint Apply Hessian_zf_zf");
  }
}
#endif


template<class Real>
class DynConstraint : public ROL::DynamicConstraint<Real> {
private:
  const ROL::Ptr<DynamicPDE<Real>>   pde_;
  ROL::Ptr<Solver<Real>>          solver_;
  ROL::Ptr<Assembler<Real>>    assembler_;

  mutable ROL::Ptr<Tpetra::MultiVector<>> cvec_;

  // Storage for spatial PDE residual
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecR_;
  mutable bool isResAssembled_;

  // Storage for spatial PDE jacobian at old time
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matJuo_;
  mutable bool isJuoAssembled_, isJuoZero_, isJuoNotImplemented_;

  // Storage for spatial PDE jacobian at new time
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matJun_;
  mutable bool isJunAssembled_, isJunZero_, isJunNotImplemented_;

  // Storage for field control jacobian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matJzf_;
  mutable bool isJzfAssembled_, isJzfZero_, isJzfNotImplemented_;

  // Storage for parameter control jacobian
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecJzp_;
  mutable bool isJzpAssembled_, isJzpZero_, isJzpNotImplemented_;

  // Storage for spatial PDE old/old hessian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHuo_uo_;
  mutable bool isHuo_uoAssembled_, isHuo_uoZero_, isHuo_uoNotImplemented_;

  // Storage for spatial PDE old/field hessian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHuo_zf_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHzf_uo_;
  mutable bool isHuo_zfAssembled_, isHuo_zfZero_, isHuo_zfNotImplemented_;
  mutable bool isHzf_uoAssembled_, isHzf_uoZero_, isHzf_uoNotImplemented_;

  // Storage for spatial PDE old/parameter hessian
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecHuo_zp_;
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecHzp_uo_;
  mutable bool isHuo_zpAssembled_, isHuo_zpZero_, isHuo_zpNotImplemented_;
  mutable bool isHzp_uoAssembled_, isHzp_uoZero_, isHzp_uoNotImplemented_;

  // Storage for spatial PDE new/old hessian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHuo_un_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHun_uo_;
  mutable bool isHuo_unAssembled_, isHuo_unZero_, isHuo_unNotImplemented_;
  mutable bool isHun_uoAssembled_, isHun_uoZero_, isHun_uoNotImplemented_;

  // Storage for spatial PDE new/new hessian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHun_un_;
  mutable bool isHun_unAssembled_, isHun_unZero_, isHun_unNotImplemented_;

  // Storage for spatial PDE new/field hessian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHun_zf_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHzf_un_;
  mutable bool isHun_zfAssembled_, isHun_zfZero_, isHun_zfNotImplemented_;
  mutable bool isHzf_unAssembled_, isHzf_unZero_, isHzf_unNotImplemented_;

  // Storage for spatial PDE new/parameter hessian
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecHun_zp_;
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecHzp_un_;
  mutable bool isHun_zpAssembled_, isHun_zpZero_, isHun_zpNotImplemented_;
  mutable bool isHzp_unAssembled_, isHzp_unZero_, isHzp_unNotImplemented_;

  // Storage for field/field hessian
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matHzf_zf_;
  mutable bool isHzf_zfAssembled_, isHzf_zfZero_, isHzf_zfNotImplemented_;

  // Storage for parameter/parameter hessian
  mutable ROL::Ptr<std::vector<std::vector<Real>>> matHzp_zp_;
  mutable bool isHzp_zpAssembled_, isHzp_zpZero_, isHzp_zpNotImplemented_;

  // Storage for field/parameter hessian
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecHzf_zp_;
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecHzp_zf_;
  mutable bool isHzf_zpAssembled_, isHzf_zpZero_, isHzf_zpNotImplemented_;
  mutable bool isHzp_zfAssembled_, isHzp_zfZero_, isHzp_zfNotImplemented_;

  void assembleResidual(const ROL::Vector<Real> &uo,
                        const ROL::Vector<Real> &un,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const {
    if (!isResAssembled_) {
      ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
      ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
      ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
      ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);

      assembler_->assembleDynPDEResidual(vecR_,pde_,ts,uof,unf,zf,zp);
      isResAssembled_ = true;
    }
  }

  void setSolver(void) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintSolverConstruct_Jacobian_un);
    #endif
    solver_->setA(matJun_);
  }

  // Assemble the PDE Jacobian.
  void assembleJuo(const ROL::Vector<Real> &uo,
                   const ROL::Vector<Real> &un,
                   const ROL::Vector<Real> &z,
                   const ROL::TimeStamp<Real> &ts) const {
    if ( !isJuoZero_ && !isJuoNotImplemented_ ) {
      try {
        if (!isJuoAssembled_) {
	  ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
	  ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
	  ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
	  ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);

	  assembler_->assembleDynPDEJacobian_uo(matJuo_,pde_,ts,uof,unf,zf,zp);
	  isJuoAssembled_ = true;
	}
      }
      catch ( Exception::Zero & ez ) {
	isJuoZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
	isJuoNotImplemented_ = true;
      }
    }
  }

  void assembleJun(const ROL::Vector<Real> &uo,
        	   const ROL::Vector<Real> &un,
        	   const ROL::Vector<Real> &z,
        	   const ROL::TimeStamp<Real> &ts) const {
    if ( !isJunZero_ && !isJunNotImplemented_ ) {
      try {
        if (!isJunAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEJacobian_un(matJun_,pde_,ts,uof,unf,zf,zp);
          isJunAssembled_ = true;
          setSolver();
        }
      }
      catch ( Exception::Zero & ez ) {
        isJunZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isJunNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE Jacobian.
  void assembleJzf(const ROL::Vector<Real> &uo,
        	   const ROL::Vector<Real> &un,
        	   const ROL::Vector<Real> &z,
        	   const ROL::TimeStamp<Real> &ts) const {
    if ( !isJzfZero_ && !isJzfNotImplemented_ ) {
      try {
        if (!isJzfAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEJacobian_zf(matJzf_,pde_,ts,uof,unf,zf,zp);
          isJzfAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isJzfZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isJzfNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE Jacobian.
  void assembleJzp(const ROL::Vector<Real> &uo,
        	   const ROL::Vector<Real> &un,
        	   const ROL::Vector<Real> &z,
        	   const ROL::TimeStamp<Real> &ts) const {
    if ( !isJzpZero_ && !isJzpNotImplemented_ ) {
      try {
        if (!isJzpAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEJacobian_zp(vecJzp_,pde_,ts,uof,unf,zf,zp);
          isJzpAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isJzpZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isJzpNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHuo_uo(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHuo_uoZero_ && !isHuo_uoNotImplemented_ ) {
      try {
        if (!isHuo_uoAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_uo_uo(matHuo_uo_,pde_,ts,wf,uof,unf,zf,zp);
          isHuo_uoAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHuo_uoZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHuo_uoNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHuo_un(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHuo_unZero_ && !isHuo_unNotImplemented_ ) {
      try {
        if (!isHuo_uoAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_uo_un(matHuo_uo_,pde_,ts,wf,uof,unf,zf,zp);
          isHuo_unAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHuo_unZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHuo_unNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHuo_zf(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHuo_zfZero_ && !isHuo_zfNotImplemented_ ) {
      try {
        if (!isHuo_zfAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_uo_zf(matHuo_zf_,pde_,ts,wf,uof,unf,zf,zp);
          isHuo_zfAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHuo_zfZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHuo_zfNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzf_uo(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzf_uoZero_ && !isHzf_uoNotImplemented_ ) {
      try {
        if (!isHzf_uoAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zf_uo(matHzf_uo_,pde_,ts,wf,uof,unf,zf,zp);
          isHzf_uoAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzf_uoZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzf_uoNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHuo_zp(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHuo_zpZero_ && !isHuo_zpNotImplemented_ ) {
      try {
        if (!isHuo_zpAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_uo_zp(vecHuo_zp_,pde_,ts,wf,uof,unf,zf,zp);
          isHuo_zpAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHuo_zpZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHuo_zpNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzp_uo(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzp_uoZero_ && !isHzp_uoNotImplemented_ ) {
      try {
        if (!isHzp_uoAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zp_uo(vecHzp_uo_,pde_,ts,wf,uof,unf,zf,zp);
          isHzp_uoAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzp_uoZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzp_uoNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHun_un(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHun_unZero_ && !isHun_unNotImplemented_ ) {
      try {
        if (!isHun_unAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_un_un(matHun_un_,pde_,ts,wf,uof,unf,zf,zp);
          isHun_unAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHun_unZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHun_unNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHun_uo(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHun_uoZero_ && !isHun_uoNotImplemented_ ) {
      try {
        if (!isHun_uoAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_un_uo(matHun_un_,pde_,ts,wf,uof,unf,zf,zp);
          isHun_uoAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHun_uoZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHun_uoNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHun_zf(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHun_zfZero_ && !isHun_zfNotImplemented_ ) {
      try {
        if (!isHun_zfAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_un_zf(matHun_zf_,pde_,ts,wf,uof,unf,zf,zp);
          isHun_zfAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHun_zfZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHun_zfNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzf_un(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzf_unZero_ && !isHzf_unNotImplemented_ ) {
      try {
        if (!isHzf_unAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zf_un(matHzf_un_,pde_,ts,wf,uof,unf,zf,zp);
          isHzf_unAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzf_unZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzf_unNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHun_zp(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHun_zpZero_ && !isHun_zpNotImplemented_ ) {
      try {
        if (!isHun_zpAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_un_zp(vecHun_zp_,pde_,ts,wf,uof,unf,zf,zp);
          isHun_zpAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHun_zpZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHun_zpNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzp_un(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzp_unZero_ && !isHzp_unNotImplemented_ ) {
      try {
        if (!isHzp_unAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>> wf  = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>> zf  = getConstField(z);
          ROL::Ptr<const std::vector<Real>>     zp  = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zp_un(vecHzp_un_,pde_,ts,wf,uof,unf,zf,zp);
          isHzp_unAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzp_unZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzp_unNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzf_zf(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzf_zfZero_ && !isHzf_zfNotImplemented_ ) {
      try {
        if (!isHzf_zfAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>>  wf = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zf_zf(matHzf_zf_,pde_,ts,wf,uof,unf,zf,zp);
          isHzf_zfAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzf_zfZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzf_zfNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzf_zp(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzf_zpZero_ && !isHzf_zpNotImplemented_ ) {
      try {
        if (!isHzf_zpAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>>  wf = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zf_zp(vecHzf_zp_,pde_,ts,wf,uof,unf,zf,zp);
          isHzf_zpAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzf_zpZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzf_zpNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzp_zf(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzp_zfZero_ && !isHzp_zfNotImplemented_ ) {
      try {
        if (!isHzp_zfAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>>  wf = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zp_zf(vecHzp_zf_,pde_,ts,wf,uof,unf,zf,zp);
          isHzp_zfAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzp_zfZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzp_zfNotImplemented_ = true;
      }
    }
  }
  
  // Assemble the PDE adjoint Hessian.
  void assembleHzp_zp(const ROL::Vector<Real> &w,
        	      const ROL::Vector<Real> &uo,
        	      const ROL::Vector<Real> &un,
        	      const ROL::Vector<Real> &z,
        	      const ROL::TimeStamp<Real> &ts) const {
    if ( !isHzp_zpZero_ && !isHzp_zpNotImplemented_ ) {
      try {
        if (!isHzp_zpAssembled_) {
          ROL::Ptr<const Tpetra::MultiVector<>>  wf = getConstField(w);
          ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
          ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
          ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
          ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);
  
          assembler_->assembleDynPDEHessian_zp_zp(matHzp_zp_,pde_,ts,wf,uof,unf,zf,zp);
          isHzp_zpAssembled_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isHzp_zpZero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isHzp_zpNotImplemented_ = true;
      }
    }
  }
  
  void solveForward(ROL::Ptr<Tpetra::MultiVector<>> &x,
        	    const ROL::Ptr<const Tpetra::MultiVector<>> &b) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintSolverSolve_Jacobian_un);
    #endif
    solver_->solve(x,b,false);
  }
  
  void solveAdjoint(ROL::Ptr<Tpetra::MultiVector<>> &x,
        	    const ROL::Ptr<const Tpetra::MultiVector<>> &b) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintSolverSolve_AdjointJacobian_un);
    #endif
    solver_->solve(x,b,true);
  }
  
  void applyJacobian_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
        		const ROL::Ptr<const Tpetra::MultiVector<>> &v,
        		const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyJacobian_uo);
      #endif
      matJuo_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyAdjointJacobian_uo);
      #endif
      matJuo_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }
  
  void applyJacobian_un(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
        		const ROL::Ptr<const Tpetra::MultiVector<>> &v,
        		const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyJacobian_un);
      #endif
      matJun_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyAdjointJacobian_un);
      #endif
      matJun_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }
  
  void applyJacobian_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
        		const ROL::Ptr<const Tpetra::MultiVector<>> &v,
        		const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyJacobian_zf);
      #endif
      matJzf_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyAdjointJacobian_zf);
      #endif
      matJzf_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }
  
  void applyJacobian_zp(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
        		const ROL::Ptr<const std::vector<Real>> &v,
        		const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyJacobian_zp);
    #endif
    if ( zeroOut ) {
      Jv->putScalar(static_cast<Real>(0));
    }
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Jv->update((*v)[i],*(vecJzp_->subView(col)),static_cast<Real>(1));
    }
  }
  
  void applyAdjointJacobian_zp(const ROL::Ptr<std::vector<Real>> &Jv,
        		       const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyAdjointJacobian_zp);
    #endif
    Teuchos::Array<Real> val(1,0);
    const size_t size = static_cast<size_t>(Jv->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      vecJzp_->subView(col)->dot(*v, val.view(0,1));
      (*Jv)[i] = val[0];
    }
  }
  
  void applyHessian_uo_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_uo_uo);
    #endif
    if (!isHuo_uoNotImplemented_) {
      if ( isHuo_uoZero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matHuo_uo_->apply(*v,*Hv);
      }
    }
    else {
      throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_uo_uo: Not Implemented!");
    }
  }
  
  void applyHessian_uo_un(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_uo_un);
    #endif
    if (!isHuo_unNotImplemented_) {
      if ( isHuo_unZero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matHuo_un_->apply(*v,*Hv);
      }
    }
    else {
      throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_uo_un: Not Implemented!");
    }
  }
  
  void applyHessian_un_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_un_uo);
    #endif
    if (!isHun_uoNotImplemented_) {
      if ( isHun_uoZero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matHun_uo_->apply(*v,*Hv);
      }
    }
    else {
      throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_un_uo: Not Implemented!");
    }
  }
  
  void applyHessian_un_un(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_un_un);
    #endif
    if (!isHun_unNotImplemented_) {
      if ( isHun_unZero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matHun_un_->apply(*v,*Hv);
      }
    }
    else {
      throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_un_un: Not Implemented!");
    }
  }
  
  void applyHessian_uo_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_uo_zf);
    #endif
    if ( Hv != ROL::nullPtr ) {
      if ( !isHuo_zfNotImplemented_ ) {
        if ( isHuo_zfZero_ ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        else {
          matHuo_zf_->apply(*v,*Hv);
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_uo_zf: Not Implemented!");
      }
    }
  }
  
  void applyHessian_un_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_un_zf);
    #endif
    if ( Hv != ROL::nullPtr ) {
      if ( !isHun_zfNotImplemented_ ) {
        if ( isHun_zfZero_ ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        else {
          matHun_zf_->apply(*v,*Hv);
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_un_zf: Not Implemented!");
      }
    }
  }
  
  void applyHessian_uo_zp(const ROL::Ptr<std::vector<Real>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_uo_zp);
    #endif
    if ( Hv != ROL::nullPtr ) {
      if ( !isHuo_zpNotImplemented_ ) {
        const size_t size = static_cast<size_t>(Hv->size());
        if ( isHuo_zpZero_ ) {
          Hv->assign(size,static_cast<Real>(0));
        }
        else {
          Teuchos::Array<Real> val(1,0);
          for (size_t i = 0; i < size; ++i) {
            Teuchos::ArrayView<const size_t> col(&i,1);
            vecHuo_zp_->subView(col)->dot(*v, val.view(0,1));
            (*Hv)[i] += val[0];
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_uo_zp: Not Implemented!");
      }
    }
  }
  
  void applyHessian_un_zp(const ROL::Ptr<std::vector<Real>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_un_zp);
    #endif
    if ( Hv != ROL::nullPtr ) {
      if ( !isHun_zpNotImplemented_ ) {
        const size_t size = static_cast<size_t>(Hv->size());
        if ( isHun_zpZero_ ) {
          Hv->assign(size,static_cast<Real>(0));
        }
        else {
          Teuchos::Array<Real> val(1,0);
          for (size_t i = 0; i < size; ++i) {
            Teuchos::ArrayView<const size_t> col(&i,1);
            vecHun_zp_->subView(col)->dot(*v, val.view(0,1));
            (*Hv)[i] += val[0];
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_un_zp: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zf_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zf_uo);
    #endif
    if ( v != ROL::nullPtr ) {
      if ( !isHzf_uoNotImplemented_ ) {
        if ( isHzf_uoZero_ ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        else {
          matHzf_uo_->apply(*v,*Hv);
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zf_uo: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zf_un(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zf_un);
    #endif
    if ( v != ROL::nullPtr ) {
      if ( !isHzf_unNotImplemented_ ) {
        if ( isHzf_unZero_ ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        else {
          matHzf_un_->apply(*v,*Hv);
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zf_un: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zp_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const std::vector<Real>> &v,
        		  const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zp_uo);
    #endif
    if ( v != ROL::nullPtr ) {
      if ( !isHzp_uoNotImplemented_ ) {
        const size_t size = static_cast<size_t>(v->size());
        if ( zeroOut || (isHzp_uoZero_ && !zeroOut) ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        if ( !isHzp_uoZero_ ) {
          for (size_t i = 0; i < size; ++i) {
            Teuchos::ArrayView<const size_t> col(&i,1);
            Hv->update((*v)[i],*(vecHzp_uo_->subView(col)),static_cast<Real>(1));
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zp_uo: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zp_un(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const std::vector<Real>> &v,
        		  const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zp_un);
    #endif
    if ( v != ROL::nullPtr ) {
      if ( !isHzp_unNotImplemented_ ) {
        const size_t size = static_cast<size_t>(v->size());
        if ( zeroOut || (isHzp_unZero_ && !zeroOut) ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        if ( !isHzp_unZero_ ) {
          for (size_t i = 0; i < size; ++i) {
            Teuchos::ArrayView<const size_t> col(&i,1);
            Hv->update((*v)[i],*(vecHzp_un_->subView(col)),static_cast<Real>(1));
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zp_un: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zf_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zf_zf);
    #endif
    if ( v != ROL::nullPtr ) {
      if ( !isHzf_zfNotImplemented_ ) {
        if ( isHzf_zfZero_ ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        else {
          matHzf_zf_->apply(*v,*Hv);
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zf_zf: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zf_zp(const ROL::Ptr<std::vector<Real>> &Hv,
        		  const ROL::Ptr<const Tpetra::MultiVector<>> &v,
        		  const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zf_zp);
    #endif
    if ( Hv != ROL::nullPtr && v != ROL::nullPtr ) {
      if ( !isHzf_zpNotImplemented_ ) {
        const size_t size = static_cast<size_t>(Hv->size());
        if ( zeroOut || (isHzf_zpZero_ && !zeroOut) ) {
          Hv->assign(size,static_cast<Real>(0));
        }
        if ( !isHzf_zpZero_ ) {
          Teuchos::Array<Real> val(1,0);
          for (size_t i = 0; i < size; ++i) {
            Teuchos::ArrayView<const size_t> col(&i,1);
            vecHzf_zp_->subView(col)->dot(*v, val.view(0,1));
            (*Hv)[i] += val[0];
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zf_zp: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zp_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Hv,
        		  const ROL::Ptr<const std::vector<Real>> &v,
        		  const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zp_zf);
    #endif
    if ( Hv != ROL::nullPtr && v != ROL::nullPtr ) {
      if ( !isHzp_zfNotImplemented_ ) {
        const size_t size = static_cast<size_t>(v->size());
        if ( zeroOut || (isHzp_zfZero_ && !zeroOut) ) {
          Hv->putScalar(static_cast<Real>(0));
        }
        if ( !isHzp_zfZero_ ) {
          for (size_t i = 0; i < size; ++i) {
            Teuchos::ArrayView<const size_t> col(&i,1);
            Hv->update((*v)[i],*(vecHzp_zf_->subView(col)),static_cast<Real>(1));
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zp_zf: Not Implemented!");
      }
    }
  }
  
  void applyHessian_zp_zp(const ROL::Ptr<std::vector<Real>> &Hv,
        		  const ROL::Ptr<const std::vector<Real>> &v,
        		  const bool zeroOut = true ) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::DynConstraintApplyHessian_zp_zp);
    #endif
    if ( Hv != ROL::nullPtr ) {
      if ( !isHzp_zpNotImplemented_ ) {
        const size_t size = static_cast<size_t>(Hv->size());
        if ( zeroOut || (isHzp_zpZero_ && !zeroOut) ) {
          Hv->assign(size,static_cast<Real>(0));
        }
        if ( !isHzp_zpZero_ ) {
          for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
              (*Hv)[i] += (*matHzp_zp_)[i][j]*(*v)[j];
            }
          }
        }
      }
      else {
        throw ROL::Exception::NotImplemented(">>> DynConstraint<Real>::applyHessian_zp_zp: Not Implemented!");
      }
    }
  }
  
public:
  
  DynConstraint(const ROL::Ptr<DynamicPDE<Real>> &pde,
        	const ROL::Ptr<MeshManager<Real>> &meshMgr,
        	const ROL::Ptr<const Teuchos::Comm<int>> &comm,
        	Teuchos::ParameterList &parlist,
        	std::ostream &outStream = std::cout)
    : pde_                    (   pde ),
      isResAssembled_         ( false ),
      isJuoAssembled_         ( false ),
      isJuoZero_              ( false ),
      isJuoNotImplemented_    ( false ),
      isJunAssembled_         ( false ),
      isJunZero_              ( false ),
      isJunNotImplemented_    ( false ),
      isJzfAssembled_         ( false ),
      isJzfZero_              ( false ),
      isJzfNotImplemented_    ( false ),
      isJzpAssembled_         ( false ),
      isJzpZero_              ( false ),
      isJzpNotImplemented_    ( false ),
      isHuo_uoAssembled_      ( false ),
      isHuo_uoZero_           ( false ),
      isHuo_uoNotImplemented_ ( false ),
      isHuo_zfAssembled_      ( false ),
      isHuo_zfZero_           ( false ),
      isHuo_zfNotImplemented_ ( false ),
      isHzf_uoAssembled_      ( false ),
      isHzf_uoZero_           ( false ),
      isHzf_uoNotImplemented_ ( false ),
      isHuo_zpAssembled_      ( false ),
      isHuo_zpZero_           ( false ),
      isHuo_zpNotImplemented_ ( false ),
      isHzp_uoAssembled_      ( false ),
      isHzp_uoZero_           ( false ),
      isHzp_uoNotImplemented_ ( false ),
      isHuo_unAssembled_      ( false ),
      isHuo_unZero_           ( false ),
      isHuo_unNotImplemented_ ( false ),
      isHun_uoAssembled_      ( false ),
      isHun_uoZero_           ( false ),
      isHun_uoNotImplemented_ ( false ),
      isHun_unAssembled_      ( false ),
      isHun_unZero_           ( false ),
      isHun_unNotImplemented_ ( false ),
      isHun_zfAssembled_      ( false ),
      isHun_zfZero_           ( false ),
      isHun_zfNotImplemented_ ( false ),
      isHzf_unAssembled_      ( false ),
      isHzf_unZero_           ( false ),
      isHzf_unNotImplemented_ ( false ),
      isHun_zpAssembled_      ( false ),
      isHun_zpZero_           ( false ),
      isHun_zpNotImplemented_ ( false ),
      isHzp_unAssembled_      ( false ),
      isHzp_unZero_           ( false ),
      isHzp_unNotImplemented_ ( false ),
      isHzf_zfAssembled_      ( false ),
      isHzf_zfZero_           ( false ),
      isHzf_zfNotImplemented_ ( false ),
      isHzp_zpAssembled_      ( false ),
      isHzp_zpZero_           ( false ),
      isHzp_zpNotImplemented_ ( false ),
      isHzf_zpAssembled_      ( false ),
      isHzf_zpZero_           ( false ),
      isHzf_zpNotImplemented_ ( false ),
      isHzp_zfAssembled_      ( false ),
      isHzp_zfZero_           ( false ),
      isHzp_zfNotImplemented_ ( false ) {
    // Construct assembler.
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    // Construct solver.
    solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
    // Initialize state and constraint vectors.
    cvec_ = assembler_->createResidualVector();
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<DynamicPDE<Real>> getPDE(void) const {
    return pde_;
  }

  void update_uo(const ROL::Vector<Real>    &uo,
                 const ROL::TimeStamp<Real> &ts) {
    isResAssembled_    = false;
    isJuoAssembled_    = (isJuoZero_    ? isJuoAssembled_    : false);
    isJunAssembled_    = (isJunZero_    ? isJunAssembled_    : false);
    isJzfAssembled_    = (isJzfZero_    ? isJzfAssembled_    : false);
    isJzpAssembled_    = (isJzpZero_    ? isJzpAssembled_    : false);
    isHuo_uoAssembled_ = (isHuo_uoZero_ ? isHuo_uoAssembled_ : false);
    isHuo_zfAssembled_ = (isHuo_zfZero_ ? isHuo_zfAssembled_ : false);
    isHuo_zpAssembled_ = (isHuo_zpZero_ ? isHuo_zpAssembled_ : false);
    isHzf_uoAssembled_ = (isHzf_uoZero_ ? isHzf_uoAssembled_ : false);
    isHzp_uoAssembled_ = (isHzp_uoZero_ ? isHzp_uoAssembled_ : false);
    isHun_unAssembled_ = (isHun_unZero_ ? isHun_unAssembled_ : false);
    isHuo_unAssembled_ = (isHuo_unZero_ ? isHuo_unAssembled_ : false);
    isHun_uoAssembled_ = (isHun_uoZero_ ? isHun_uoAssembled_ : false);
    isHun_zfAssembled_ = (isHun_zfZero_ ? isHun_zfAssembled_ : false);
    isHun_zpAssembled_ = (isHun_zpZero_ ? isHun_zpAssembled_ : false);
    isHzf_unAssembled_ = (isHzf_unZero_ ? isHzf_unAssembled_ : false);
    isHzp_unAssembled_ = (isHzp_unZero_ ? isHzp_unAssembled_ : false);
    isHzf_zfAssembled_ = (isHzf_zfZero_ ? isHzf_zfAssembled_ : false);
    isHzf_zpAssembled_ = (isHzf_zpZero_ ? isHzf_zpAssembled_ : false);
    isHzp_zfAssembled_ = (isHzp_zfZero_ ? isHzp_zfAssembled_ : false);
    isHzp_zpAssembled_ = (isHzp_zpZero_ ? isHzp_zpAssembled_ : false);
  }

  void update_un(const ROL::Vector<Real>    &un,
                 const ROL::TimeStamp<Real> &ts) {
    isResAssembled_    = false;
    isJuoAssembled_    = (isJuoZero_    ? isJuoAssembled_    : false);
    isJunAssembled_    = (isJunZero_    ? isJunAssembled_    : false);
    isJzfAssembled_    = (isJzfZero_    ? isJzfAssembled_    : false);
    isJzpAssembled_    = (isJzpZero_    ? isJzpAssembled_    : false);
    isHuo_uoAssembled_ = (isHuo_uoZero_ ? isHuo_uoAssembled_ : false);
    isHuo_zfAssembled_ = (isHuo_zfZero_ ? isHuo_zfAssembled_ : false);
    isHuo_zpAssembled_ = (isHuo_zpZero_ ? isHuo_zpAssembled_ : false);
    isHzf_uoAssembled_ = (isHzf_uoZero_ ? isHzf_uoAssembled_ : false);
    isHzp_uoAssembled_ = (isHzp_uoZero_ ? isHzp_uoAssembled_ : false);
    isHun_unAssembled_ = (isHun_unZero_ ? isHun_unAssembled_ : false);
    isHuo_unAssembled_ = (isHuo_unZero_ ? isHuo_unAssembled_ : false);
    isHun_uoAssembled_ = (isHun_uoZero_ ? isHun_uoAssembled_ : false);
    isHun_zfAssembled_ = (isHun_zfZero_ ? isHun_zfAssembled_ : false);
    isHun_zpAssembled_ = (isHun_zpZero_ ? isHun_zpAssembled_ : false);
    isHzf_unAssembled_ = (isHzf_unZero_ ? isHzf_unAssembled_ : false);
    isHzp_unAssembled_ = (isHzp_unZero_ ? isHzp_unAssembled_ : false);
    isHzf_zfAssembled_ = (isHzf_zfZero_ ? isHzf_zfAssembled_ : false);
    isHzf_zpAssembled_ = (isHzf_zpZero_ ? isHzf_zpAssembled_ : false);
    isHzp_zfAssembled_ = (isHzp_zfZero_ ? isHzp_zfAssembled_ : false);
    isHzp_zpAssembled_ = (isHzp_zpZero_ ? isHzp_zpAssembled_ : false);
  }

  void update_z(const ROL::Vector<Real>    &z,
                const ROL::TimeStamp<Real> &ts) {
    isResAssembled_    = false;
    isJuoAssembled_    = (isJuoZero_    ? isJuoAssembled_    : false);
    isJunAssembled_    = (isJunZero_    ? isJunAssembled_    : false);
    isJzfAssembled_    = (isJzfZero_    ? isJzfAssembled_    : false);
    isJzpAssembled_    = (isJzpZero_    ? isJzpAssembled_    : false);
    isHuo_uoAssembled_ = (isHuo_uoZero_ ? isHuo_uoAssembled_ : false);
    isHuo_zfAssembled_ = (isHuo_zfZero_ ? isHuo_zfAssembled_ : false);
    isHuo_zpAssembled_ = (isHuo_zpZero_ ? isHuo_zpAssembled_ : false);
    isHzf_uoAssembled_ = (isHzf_uoZero_ ? isHzf_uoAssembled_ : false);
    isHzp_uoAssembled_ = (isHzp_uoZero_ ? isHzp_uoAssembled_ : false);
    isHun_unAssembled_ = (isHun_unZero_ ? isHun_unAssembled_ : false);
    isHuo_unAssembled_ = (isHuo_unZero_ ? isHuo_unAssembled_ : false);
    isHun_uoAssembled_ = (isHun_uoZero_ ? isHun_uoAssembled_ : false);
    isHun_zfAssembled_ = (isHun_zfZero_ ? isHun_zfAssembled_ : false);
    isHun_zpAssembled_ = (isHun_zpZero_ ? isHun_zpAssembled_ : false);
    isHzf_unAssembled_ = (isHzf_unZero_ ? isHzf_unAssembled_ : false);
    isHzp_unAssembled_ = (isHzp_unZero_ ? isHzp_unAssembled_ : false);
    isHzf_zfAssembled_ = (isHzf_zfZero_ ? isHzf_zfAssembled_ : false);
    isHzf_zpAssembled_ = (isHzf_zpZero_ ? isHzf_zpAssembled_ : false);
    isHzp_zfAssembled_ = (isHzp_zfZero_ ? isHzp_zfAssembled_ : false);
    isHzp_zpAssembled_ = (isHzp_zpZero_ ? isHzp_zpAssembled_ : false);
  }

  void update(const ROL::Vector<Real>    &uo,
              const ROL::Vector<Real>    &un,
              const ROL::Vector<Real>    &z,
              const ROL::TimeStamp<Real> &ts) {
    update_uo(uo,ts);
    update_un(un,ts);
    update_z(z,ts);
  }

  void value(ROL::Vector<Real>    &c,
       const ROL::Vector<Real>    &uo,
       const ROL::Vector<Real>    &un,
       const ROL::Vector<Real>    &z,
       const ROL::TimeStamp<Real> &ts) const {
    const Real zero(0), one(1);
    ROL::Ptr<Tpetra::MultiVector<>> cf = getField(c);
    assembleResidual(uo,un,z,ts);
    cf->update(one, *vecR_, zero);
  }

  void solve(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &uo,
             ROL::Vector<Real> &un,
       const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) {
    ROL::DynamicConstraint<Real>::solve(c,uo,un,z,ts);
  }

  void applyJacobian_uo(ROL::Vector<Real>    &jv,
                  const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real>    &uo,
                  const ROL::Vector<Real>    &un,
                  const ROL::Vector<Real>    &z,
                  const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleJuo(uo,un,z,ts);
    applyJacobian_uo(jvf,vf,false);
  }

  void applyJacobian_un(ROL::Vector<Real>    &jv,
                  const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real>    &uo,
                  const ROL::Vector<Real>    &un,
                  const ROL::Vector<Real>    &z,
                  const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleJun(uo,un,z,ts);
    applyJacobian_un(jvf,vf,false);
  }


  void applyJacobian_z(ROL::Vector<Real>    &jv,
                 const ROL::Vector<Real>    &v,
                 const ROL::Vector<Real>    &uo,
                 const ROL::Vector<Real>    &un,
                 const ROL::Vector<Real>    &z,
                 const ROL::TimeStamp<Real> &ts) const {
    jv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);
    if (vf != ROL::nullPtr) {
      assembleJzf(uo,un,z,ts);
      applyJacobian_zf(jvf,vf,false);
    }
    if (vp != ROL::nullPtr) {
      assembleJzp(uo,un,z,ts);
      applyJacobian_zp(jvf,vp);
    }
  }


  void applyAdjointJacobian_uo(ROL::Vector<Real>    &ajv,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>     ajvf = getField(ajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleJuo(uo,un,z,ts);
    applyJacobian_uo(ajvf,vf,true);
  }


  void applyAdjointJacobian_un(ROL::Vector<Real>    &ajv,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>     ajvf = getField(ajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleJun(uo,un,z,ts);
    applyJacobian_un(ajvf,vf,true);
  }


  void applyAdjointJacobian_z(ROL::Vector<Real>    &ajv,
                        const ROL::Vector<Real>    &v,
                        const ROL::Vector<Real>    &uo,
                        const ROL::Vector<Real>    &un,
                        const ROL::Vector<Real>    &z,
                        const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>     ajvf = getField(ajv);
    ROL::Ptr<std::vector<Real>>         ajvp = getParameter(ajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    if (ajvf != ROL::nullPtr) {
      assembleJzf(uo,un,z,ts);
      applyJacobian_zf(ajvf,vf,true);
    }
    if (ajvp != ROL::nullPtr) {
      assembleJzp(uo,un,z,ts);
      applyAdjointJacobian_zp(ajvp,vf);
    }
  }


  void applyInverseJacobian_un(ROL::Vector<Real>    &ijv,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>     ijvf = getField(ijv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleJun(uo,un,z,ts);
    solveForward(ijvf,vf);
  }


  void applyInverseAdjointJacobian_un(ROL::Vector<Real>    &iajv,
                                const ROL::Vector<Real>    &v,
                                const ROL::Vector<Real>    &uo,
                                const ROL::Vector<Real>    &un,
                                const ROL::Vector<Real>    &z,
                                const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>    iajvf = getField(iajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleJun(uo,un,z,ts);
    solveAdjoint(iajvf,vf);
  }


  void applyAdjointHessian_uo_uo(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleHuo_uo(w,uo,un,z,ts);
    applyHessian_uo_uo(ahwvf,vf);
  }


  void applyAdjointHessian_uo_un(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    assembleHuo_un(w,uo,un,z,ts);
    applyHessian_uo_un(ahwvf,vf);
  }


  void applyAdjointHessian_uo_z(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>        ahwvp = getParameter(ahwv);
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    if (ahwvf != ROL::nullPtr) {
      assembleHuo_zf(w,uo,un,z,ts);
      applyHessian_uo_zf(ahwvf,vf);
    }
    if (ahwvp != ROL::nullPtr) {
      assembleHuo_zp(w,uo,un,z,ts);
      applyHessian_uo_zp(ahwvp,vf);
    }
  }


  void applyAdjointHessian_un_uo(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    assembleHun_uo(w,uo,un,z,ts);
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    applyHessian_un_uo(ahwvf,vf);
  }


  void applyAdjointHessian_un_un(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    assembleHun_un(w,uo,un,z,ts);
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    applyHessian_un_un(ahwvf,vf);
  }


  void applyAdjointHessian_un_z(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<std::vector<Real>>        ahwvp = getParameter(ahwv);
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    if (ahwvf != ROL::nullPtr) {
      assembleHun_zf(w,uo,un,z,ts);
      applyHessian_un_zf(ahwvf,vf);
    }
    if (ahwvp != ROL::nullPtr) {
      assembleHun_zp(w,uo,un,z,ts);
      applyHessian_un_zp(ahwvp,vf);
    }
  }


  void applyAdjointHessian_z_uo(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);
    if (vf != ROL::nullPtr) {
      assembleHzf_uo(w,uo,un,z,ts);
      applyHessian_zf_uo(ahwvf,vf);
    }
    if (vp != ROL::nullPtr) {
      bool zeroOut = (vf == ROL::nullPtr);
      assembleHzp_uo(w,uo,un,z,ts);
      applyHessian_zp_uo(ahwvf,vp,zeroOut);
    }
  }


  void applyAdjointHessian_z_un(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);
    if (vf != ROL::nullPtr) {
      assembleHzf_un(w,uo,un,z,ts);
      applyHessian_zf_un(ahwvf,vf);
    }
    if (vp != ROL::nullPtr) {
      bool zeroOut = (vf == ROL::nullPtr);
      assembleHzp_un(w,uo,un,z,ts);
      applyHessian_zp_un(ahwvf,vp,zeroOut);
    }
  }


  void applyAdjointHessian_z_z(ROL::Vector<Real>    &ahwv,
                         const ROL::Vector<Real>    &w,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<std::vector<Real>>        ahwvp = getParameter(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);

    bool zeroOut = (vf == ROL::nullPtr);
    if (vf != ROL::nullPtr) {
      assembleHzf_zf(w,uo,un,z,ts);
      applyHessian_zf_zf(ahwvf,vf);
    }
    if (vp != ROL::nullPtr) {
      assembleHzp_zp(w,uo,un,z,ts);
      applyHessian_zp_zp(ahwvp,vp,zeroOut);
    }
    if (vf != ROL::nullPtr && vp != ROL::nullPtr) {
      assembleHzf_zp(w,uo,un,z,ts);
      applyHessian_zf_zp(ahwvp,vf,zeroOut);
      assembleHzp_zf(w,uo,un,z,ts);
      applyHessian_zp_zf(ahwvf,vp);
    }
  }

  /***************************************************************************/
  /* Output routines.                                                        */
  /***************************************************************************/
  void printMeshData(std::ostream &outStream) const {
    assembler_->printMeshData(outStream);
  }

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> matWriter;
    if (matJuo_ != ROL::nullPtr) {
      matWriter.writeSparseFile("jacobian_uo.txt", matJuo_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian_uo.txt");
      emptyfile.close();
    }
    if (matJun_ != ROL::nullPtr) {
      matWriter.writeSparseFile("jacobian_un.txt", matJun_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian_un.txt");
      emptyfile.close();
    }
    if (matJzf_ != ROL::nullPtr) {
      matWriter.writeSparseFile("jacobian_zf.txt", matJzf_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian_zf.txt");
      emptyfile.close();
    }
    if (vecR_ != ROL::nullPtr) {
      matWriter.writeDenseFile("residual.txt", vecR_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("residual.txt");
      emptyfile.close();
    }
  }

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<>> &vec,
                          const std::string &filename) const {
    assembler_->outputTpetraVector(vec, filename);
  }

  void inputTpetraVector(ROL::Ptr<Tpetra::MultiVector<>> &vec,
                         const std::string &filename) const {
    assembler_->inputTpetraVector(vec, filename);
  }
  /***************************************************************************/
  /* End of output routines.                                                 */
  /***************************************************************************/

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<> > xp;
    try {
      xp = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::TpetraMultiVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
        if (xvec == ROL::nullPtr) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<> > xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::TpetraMultiVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getField();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::StdVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::StdVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }
};

#endif
