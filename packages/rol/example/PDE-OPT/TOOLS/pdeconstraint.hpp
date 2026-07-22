// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  pdeconstraint.hpp
    \brief Defines the SimOpt constraint for the PDE-OPT examples.
*/

#ifndef ROL_PDEOPT_PDECONSTRAINT_H
#define ROL_PDEOPT_PDECONSTRAINT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "pde.hpp"
#include "assembler.hpp"
#include "solver.hpp"
#include "pdevector.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> PDEConstraintSolverConstruct_Jacobian1        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Construction Time for Jacobian1");
    ROL::Ptr<Teuchos::Time> PDEConstraintSolverConstruct_AdjointJacobian1 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Construction Time for Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> PDEConstraintSolverSolve_Jacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Solution Time for Jacobian1");
    ROL::Ptr<Teuchos::Time> PDEConstraintSolverSolve_AdjointJacobian1     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Solution Time for Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyJacobian1                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Jacobian1");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyJacobian2                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Jacobian2");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyJacobian3                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Jacobian3");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyAdjointJacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyAdjointJacobian2            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Adjoint Jacobian2");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyAdjointJacobian3            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Adjoint Jacobian3");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian11                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian11");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian12                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian12");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian13                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian13");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian21                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian21");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian22                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian22");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian23                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian23");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian31                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian31");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian32                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian32");
    ROL::Ptr<Teuchos::Time> PDEConstraintApplyHessian33                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian33");
  }
}
#endif


template<class Real>
class PDE_Constraint : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<PDE<Real> > pde_;
  ROL::Ptr<Assembler<Real> > assembler_;
  ROL::Ptr<Solver<Real> > solver_;

  ROL::Ptr<Tpetra::MultiVector<> >           vecR_;
  ROL::Ptr<Tpetra::CrsMatrix<> >             matJ1_;
  ROL::Ptr<Tpetra::CrsMatrix<> >             matJ2_;
  ROL::Ptr<Tpetra::MultiVector<> >           vecJ3_;
  ROL::Ptr<Tpetra::CrsMatrix<> >             matH11_;
  ROL::Ptr<Tpetra::CrsMatrix<> >             matH12_;
  ROL::Ptr<Tpetra::MultiVector<> >           vecH13_;
  ROL::Ptr<Tpetra::CrsMatrix<> >             matH21_;
  ROL::Ptr<Tpetra::CrsMatrix<> >             matH22_;
  ROL::Ptr<Tpetra::MultiVector<> >           vecH23_;
  ROL::Ptr<Tpetra::MultiVector<> >           vecH31_;
  ROL::Ptr<Tpetra::MultiVector<> >           vecH32_;
  ROL::Ptr<std::vector<std::vector<Real> > > matH33_;

  bool computeJ1_,  computeJ2_,  computeJ3_;
  bool setSolver_;
  bool isJ1zero_, isJ1notImplemented_;
  bool isJ2zero_, isJ2notImplemented_;
  bool isJ3zero_, isJ3notImplemented_;
  bool isH11zero_, isH11notImplemented_;
  bool isH12zero_, isH12notImplemented_;
  bool isH13zero_, isH13notImplemented_;
  bool isH21zero_, isH21notImplemented_;
  bool isH22zero_, isH22notImplemented_;
  bool isH23zero_, isH23notImplemented_;
  bool isH31zero_, isH31notImplemented_;
  bool isH32zero_, isH32notImplemented_;
  bool isH33zero_, isH33notImplemented_;

  // Assemble the PDE Jacobian.
  void assembleJ1(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isJ1zero_ && !isJ1notImplemented_ ) {
      try {
        if (computeJ1_) {
          ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
          ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
          ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

          assembler_->assemblePDEJacobian1(matJ1_,pde_,uf,zf,zp);
          computeJ1_ = false;
          setSolver_ = true;
        }
      }
      catch ( Exception::Zero & ez ) {
        isJ1zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isJ1notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE Jacobian.
  void assembleJ2(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isJ2zero_ && !isJ2notImplemented_ ) {
      try {
        if (computeJ2_) {
          ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
          ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
          ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

          assembler_->assemblePDEJacobian2(matJ2_,pde_,uf,zf,zp);
          computeJ2_ = false;
        }
      }
      catch ( Exception::Zero & ez ) {
        isJ2zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isJ2notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE Jacobian.
  void assembleJ3(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isJ3zero_ && !isJ3notImplemented_ ) {
      try {
        if (computeJ3_) {
          ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
          ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
          ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

          assembler_->assemblePDEJacobian3(vecJ3_,pde_,uf,zf,zp);
          computeJ3_ = false;
        }
      }
      catch ( Exception::Zero & ez ) {
        isJ3zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isJ3notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH11(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH11zero_ && !isH11notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian11(matH11_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH11zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH11notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH12(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH12zero_ && !isH12notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian12(matH12_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH12zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH12notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH13(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH13zero_ && !isH13notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian13(vecH13_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH13zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH13notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH21(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH21zero_ && !isH21notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian21(matH21_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH21zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH21notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH31(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH31zero_ && !isH31notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian31(vecH31_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH31zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH31notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH22(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH22zero_ && !isH22notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian22(matH22_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH22zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH22notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH23(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH23zero_ && !isH23notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian23(vecH23_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH23zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH23notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH32(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH32zero_ && !isH32notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian32(vecH32_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH32zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH32notImplemented_ = true;
      }
    }
  }

  // Assemble the PDE adjoint Hessian.
  void assembleH33(const ROL::Vector<Real> &w, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isH33zero_ && !isH33notImplemented_ ) {
      try {
        ROL::Ptr<const Tpetra::MultiVector<> > wf = getConstField(w);
        ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
        ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
        ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian33(matH33_,pde_,wf,uf,zf,zp);
      }
      catch ( Exception::Zero & ez ) {
        isH33zero_ = true;
      }
      catch ( Exception::NotImplemented & eni ) {
        isH33notImplemented_ = true;
      }
    }
  }

  // Application routines
  void applyJacobian1(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian1);
    #endif
    if (!isJ1notImplemented_) {
      if (isJ1zero_) {
        Jv->putScalar(static_cast<Real>(0));
      }
      else {
        matJ1_->apply(*v,*Jv);
      }
    }
  }

  void applyAdjointJacobian1(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyAdjointJacobian1);
    #endif
    if ( !isJ1notImplemented_ ) {
      if ( isJ1zero_ ) {
        Jv->putScalar(static_cast<Real>(0));
      }
      else {
        matJ1_->apply(*v,*Jv,Teuchos::TRANS);
      }
    }
  }

  void applyJacobian2(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian2);
    #endif
    if (v != ROL::nullPtr && !isJ2notImplemented_) {
      if (isJ2zero_) {
        Jv->putScalar(static_cast<Real>(0));
      }
      else {
        matJ2_->apply(*v,*Jv);
      }
    }
  }

  void applyAdjointJacobian2(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyAdjointJacobian2);
    #endif
    if ( Jv != ROL::nullPtr && !isJ2notImplemented_ ) {
      if ( isJ2zero_ ) {
        Jv->putScalar(static_cast<Real>(0));
      }
      else {
        matJ2_->apply(*v,*Jv,Teuchos::TRANS);
      }
    }
  }

  void applyJacobian3(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                      const ROL::Ptr<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian3);
    #endif
    if ( v != ROL::nullPtr && !isJ3notImplemented_ ) {
      const size_t size = static_cast<size_t>(v->size());
      if ( zeroOut ) { // || (isJ3zero_ && !zeroOut) ) {
        Jv->putScalar(static_cast<Real>(0));
      }
      if ( !isJ3zero_ ) {
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          Jv->update((*v)[i],*(vecJ3_->subView(col)),static_cast<Real>(1));
        }
      }
    }
  }

  void applyAdjointJacobian3(const ROL::Ptr<std::vector<Real> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyAdjointJacobian3);
    #endif
    if ( Jv != ROL::nullPtr && !isJ3notImplemented_ ) {
      const size_t size = static_cast<size_t>(Jv->size());
      if ( isJ3zero_ ) {
        Jv->assign(size,static_cast<Real>(0));
      }
      else {
        Teuchos::Array<Real> val(1,0);
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          vecJ3_->subView(col)->dot(*v, val.view(0,1));
          (*Jv)[i] = val[0];
        }
      }
    }
  }

  void applyHessian11(const ROL::Ptr<Tpetra::MultiVector<> > &Hv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian11);
    #endif
    if (!isH11notImplemented_) {
      if ( isH11zero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matH11_->apply(*v,*Hv);
      }
    }
  }

  void applyHessian12(const ROL::Ptr<Tpetra::MultiVector<> > &Hv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian12);
    #endif
    if ( Hv != ROL::nullPtr && !isH12notImplemented_ ) {
      if ( isH12zero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matH12_->apply(*v,*Hv);
      }
    }
  }

  void applyHessian13(const ROL::Ptr<std::vector<Real> > &Hv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian13);
    #endif
    if ( Hv != ROL::nullPtr && !isH13notImplemented_ ) {
      const size_t size = static_cast<size_t>(Hv->size());
      if ( isH13zero_ ) {
        Hv->assign(size,static_cast<Real>(0));
      }
      else {
        Teuchos::Array<Real> val(1,0);
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          vecH13_->subView(col)->dot(*v, val.view(0,1));
          (*Hv)[i] += val[0];
        }
      }
    }
  }

  void applyHessian21(const ROL::Ptr<Tpetra::MultiVector<> > &Hv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian21);
    #endif
    if ( v != ROL::nullPtr && !isH21notImplemented_ ) {
      if ( isH21zero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matH21_->apply(*v,*Hv);
      }
    }
  }

  void applyHessian31(const ROL::Ptr<Tpetra::MultiVector<> > &Hv,
                      const ROL::Ptr<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian31);
    #endif
    if ( v != ROL::nullPtr && !isH31notImplemented_ ) {
      const size_t size = static_cast<size_t>(v->size());
      if ( zeroOut ) { //|| (isH31zero_ && !zeroOut) ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      if ( !isH31zero_ ) {
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          Hv->update((*v)[i],*(vecH31_->subView(col)),static_cast<Real>(1));
        }
      }
    }
  }

  void applyHessian22(const ROL::Ptr<Tpetra::MultiVector<> > &Hv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian22);
    #endif
    if ( v != ROL::nullPtr && !isH22notImplemented_ ) {
      if ( isH22zero_ ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      else {
        matH22_->apply(*v,*Hv);
      }
    }
  }

  void applyHessian23(const ROL::Ptr<std::vector<Real> > &Hv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian23);
    #endif
    if ( Hv != ROL::nullPtr && v != ROL::nullPtr && !isH23notImplemented_ ) {
      const size_t size = static_cast<size_t>(Hv->size());
      if ( zeroOut ) { // || (isH23zero_ && !zeroOut) ) {
        Hv->assign(size,static_cast<Real>(0));
      }
      if ( !isH23zero_ ) {
        Teuchos::Array<Real> val(1,0);
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          vecH23_->subView(col)->dot(*v, val.view(0,1));
          (*Hv)[i] += val[0];
        }
      }
    }
  }

  void applyHessian32(const ROL::Ptr<Tpetra::MultiVector<> > &Hv,
                      const ROL::Ptr<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian32);
    #endif
    if ( Hv != ROL::nullPtr && v != ROL::nullPtr && !isH32notImplemented_ ) {
      const size_t size = static_cast<size_t>(v->size());
      if ( zeroOut ) { // || (isH32zero_ && !zeroOut) ) {
        Hv->putScalar(static_cast<Real>(0));
      }
      if ( !isH32zero_ ) {
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          Hv->update((*v)[i],*(vecH32_->subView(col)),static_cast<Real>(1));
        }
      }
    }
  }

  void applyHessian33(const ROL::Ptr<std::vector<Real> > &Hv,
                      const ROL::Ptr<const std::vector<Real> > &v,
                      const bool zeroOut = true ) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian33);
    #endif
    if ( Hv != ROL::nullPtr && !isH33notImplemented_ ) {
      const size_t size = static_cast<size_t>(Hv->size());
      if ( zeroOut ) { // || (isH33zero_ && !zeroOut) ) {
        Hv->assign(size,static_cast<Real>(0));
      }
      if ( !isH33zero_ ) {
        for (size_t i = 0; i < size; ++i) {
          for (size_t j = 0; j < size; ++j) {
            (*Hv)[i] += (*matH33_)[i][j]*(*v)[j];
          }
        }
      }
    }
  }

  // Set the Jacobian matrix in the solver object.
  // Assumes assembleJ1 has already been called.
  void setSolver(void) {
    solver_->setA(matJ1_);
    setSolver_= false;
  }

  // Solve using the Jacobian.
  // Assumes assembleJ1 has already been called.
  void solveForward(ROL::Ptr<Tpetra::MultiVector<> > &x,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &b) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintSolverConstruct_Jacobian1);
      #endif
      setSolver();
    }

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintSolverSolve_Jacobian1);
    #endif
    solver_->solve(x,b,false);
  }

  // Solve using the adjoint Jacobian.
  // Assumes assembleJ1 has already been called.
  void solveAdjoint(ROL::Ptr<Tpetra::MultiVector<> > &x,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &b) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintSolverConstruct_AdjointJacobian1);
      #endif
      setSolver();
    }

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintSolverSolve_AdjointJacobian1);
    #endif
    solver_->solve(x,b,true);
  }

public:
  PDE_Constraint(const ROL::Ptr<PDE<Real> > &pde,
                 const ROL::Ptr<MeshManager<Real> > &meshMgr,
                 const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                 Teuchos::ParameterList &parlist,
                 std::ostream &outStream = std::cout)
    : ROL::Constraint_SimOpt<Real>(),
      pde_(pde),
      computeJ1_(true),  computeJ2_(true),  computeJ3_(true),
      setSolver_(true),
      isJ1zero_(false),  isJ1notImplemented_(false),
      isJ2zero_(false),  isJ2notImplemented_(false),
      isJ3zero_(false),  isJ3notImplemented_(false),
      isH11zero_(false),  isH11notImplemented_(false),
      isH12zero_(false),  isH12notImplemented_(false),
      isH13zero_(false),  isH13notImplemented_(false),
      isH21zero_(false),  isH21notImplemented_(false),
      isH22zero_(false),  isH22notImplemented_(false),
      isH23zero_(false),  isH23notImplemented_(false),
      isH31zero_(false),  isH31notImplemented_(false),
      isH32zero_(false),  isH32notImplemented_(false),
      isH33zero_(false),  isH33notImplemented_(false) {
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),pde_->getFields2(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
  }

  PDE_Constraint(const ROL::Ptr<PDE<Real> >       &pde,
                 const ROL::Ptr<Assembler<Real> > &assembler,
                 const ROL::Ptr<Solver<Real> >    &solver)
    : ROL::Constraint_SimOpt<Real>(),
      pde_(pde), assembler_(assembler), solver_(solver),
      computeJ1_(true),  computeJ2_(true),  computeJ3_(true),
      setSolver_(true),
      isJ1zero_(false),  isJ1notImplemented_(false),
      isJ2zero_(false),  isJ2notImplemented_(false),
      isJ3zero_(false),  isJ3notImplemented_(false),
      isH11zero_(false),  isH11notImplemented_(false),
      isH12zero_(false),  isH12notImplemented_(false),
      isH13zero_(false),  isH13notImplemented_(false),
      isH21zero_(false),  isH21notImplemented_(false),
      isH22zero_(false),  isH22notImplemented_(false),
      isH23zero_(false),  isH23notImplemented_(false),
      isH31zero_(false),  isH31notImplemented_(false),
      isH32zero_(false),  isH32notImplemented_(false),
      isH33zero_(false),  isH33notImplemented_(false) {
    assembler_->setCellNodes(*pde_);
  }

  void setParameter(const std::vector<Real> &param) {
    computeJ1_  = true; computeJ2_  = true; computeJ3_  = true;
    setSolver_  = true;
    ROL::Constraint_SimOpt<Real>::setParameter(param);
    pde_->setParameter(param);
  }

  const ROL::Ptr<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<PDE<Real> > getPDE(void) const {
    return pde_;
  }

  using ROL::Constraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {
    computeJ1_ = (flag ? true : computeJ1_);
    computeJ2_ = (flag ? true : computeJ2_);
    computeJ3_ = (flag ? true : computeJ3_);
  }

  using ROL::Constraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    computeJ1_ = (flag ? true : computeJ1_);
    computeJ2_ = (flag ? true : computeJ2_);
    computeJ3_ = (flag ? true : computeJ3_);
  }

  using ROL::Constraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
  }

  // Do something smarter here
  void update_1(const ROL::Vector<Real> &u, ROL::UpdateType type, int iter = -1) {
    computeJ1_ = true;
    computeJ2_ = true;
    computeJ3_ = true;
  }

  void update_2(const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    computeJ1_ = true;
    computeJ2_ = true;
    computeJ3_ = true;
  }

  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    update_1(u,type,iter);
    update_2(z,type,iter);
  }

  void solve_update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    update(u,z,type,iter);
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Constraint_SimOpt<Real>::solve(c,u,z,tol);
  }

  using ROL::Constraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >       cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assembler_->assemblePDEResidual(vecR_,pde_,uf,zf,zp);
    cf->scale(static_cast<Real>(1),*vecR_);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::Constraint_SimOpt<Real>::applyJacobian_1(jv,v,u,z,tol);
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> >      jvf = getField(jv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      applyJacobian1(jvf,vf);
    }
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    assembleJ2(u,z);
    assembleJ3(u,z);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
    bool useFD2 = (isJ2notImplemented_ && vf != ROL::nullPtr);
    bool useFD3 = (isJ3notImplemented_ && vp != ROL::nullPtr);
    jv.zero();
    if (useFD2 || useFD3) {
      ROL::Constraint_SimOpt<Real>::applyJacobian_2(jv,v,u,z,tol);
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> > jvf = getField(jv);
      bool zeroOut = true;
      if (!useFD2) {
        applyJacobian2(jvf,vf);
        zeroOut = (vf == ROL::nullPtr);
      }
      if (!useFD3) {
        applyJacobian3(jvf,vp,zeroOut);
      }
    }
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::Constraint_SimOpt<Real>::applyAdjointJacobian_1(ajv,v,u,z,tol);
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> >     ajvf = getField(ajv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      applyAdjointJacobian1(ajvf,vf);
    }
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleJ2(u,z);
    assembleJ3(u,z);
    ROL::Ptr<Tpetra::MultiVector<> > ajvf = getField(ajv);
    ROL::Ptr<std::vector<Real> >     ajvp = getParameter(ajv);
    bool useFD2 = (isJ2notImplemented_ && ajvf != ROL::nullPtr);
    bool useFD3 = (isJ3notImplemented_ && ajvp != ROL::nullPtr);
    if (useFD2 || useFD3) {
      ROL::Constraint_SimOpt<Real>::applyAdjointJacobian_2(ajv,v,u,z,tol);
    }
    else {
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      if (!useFD2) {
        applyAdjointJacobian2(ajvf,vf);
      }
      if (!useFD3) {
        applyAdjointJacobian3(ajvp,vf);
      }
    }
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleH11(w,u,z);
    if (isH11notImplemented_) {
      ROL::Constraint_SimOpt<Real>::applyAdjointHessian_11(ahwv,w,v,u,z,tol);
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      applyHessian11(ahwvf,vf);
    }
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleH12(w,u,z);
    assembleH13(w,u,z);
    ROL::Ptr<std::vector<Real> >     ahwvp = getParameter(ahwv);
    ROL::Ptr<Tpetra::MultiVector<> > ahwvf = getField(ahwv);
    bool useFD2 = (isH12notImplemented_ && ahwvf != ROL::nullPtr);
    bool useFD3 = (isH13notImplemented_ && ahwvp != ROL::nullPtr);
    if (useFD2 || useFD3) {
      ROL::Constraint_SimOpt<Real>::applyAdjointHessian_12(ahwv,w,v,u,z,tol);
    }
    else {
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      if (!useFD2) {
        applyHessian12(ahwvf,vf);
      }
      if (!useFD3) {
        applyHessian13(ahwvp,vf);
      }
    }
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleH21(w,u,z);
    assembleH31(w,u,z);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
    bool useFD2 = (isH21notImplemented_ && vf != ROL::nullPtr);
    bool useFD3 = (isH31notImplemented_ && vp != ROL::nullPtr);
    ahwv.zero();
    if (useFD2 || useFD3) {
      ROL::Constraint_SimOpt<Real>::applyAdjointHessian_21(ahwv,w,v,u,z,tol);
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> > ahwvf = getField(ahwv);
      bool zeroOut = true;
      if (!useFD2) {
        applyHessian21(ahwvf,vf);
        zeroOut = (vf == ROL::nullPtr);
      }
      if (!useFD3) {
        applyHessian31(ahwvf,vp,zeroOut);
      }
    }
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleH22(w,u,z);
    assembleH23(w,u,z);
    assembleH32(w,u,z);
    assembleH33(w,u,z);
    ahwv.zero();

    ROL::Ptr<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
    ROL::Ptr<std::vector<Real> >        ahwvp = getParameter(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
    bool useFD22 = (isH22notImplemented_ && vf != ROL::nullPtr);
    bool useFD23 = (isH23notImplemented_ && vp != ROL::nullPtr);
    bool useFD32 = (isH32notImplemented_ && vp != ROL::nullPtr);
    bool useFD33 = (isH33notImplemented_ && vp != ROL::nullPtr);
    if (useFD22 || useFD23 || useFD32 || useFD33) {
      ROL::Constraint_SimOpt<Real>::applyAdjointHessian_22(ahwv,w,v,u,z,tol);
    }
    else {
      bool zeroOut = true;
      if (!useFD22) {
        applyHessian22(ahwvf,vf);
        zeroOut = (vf == ROL::nullPtr);
      }
      if (!useFD23) {
        applyHessian23(ahwvp,vf,zeroOut);
      }
      if (!useFD32) {
        applyHessian32(ahwvf,vp,zeroOut);
        zeroOut = (vf == ROL::nullPtr);
      }
      if (!useFD33) {
        applyHessian33(ahwvp,vp,zeroOut);
      }
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::Constraint_SimOpt<Real>::applyInverseJacobian_1(ijv,v,u,z,tol);
    }
    else {
      if (isJ1zero_) {
        throw Exception::Zero(">>> (PDE_Constraint::applyInverseJacobian_1): Jacobian is zero!");
      }
      else {
        ROL::Ptr<Tpetra::MultiVector<> >     ijvf = getField(ijv);
        ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
        solveForward(ijvf,vf);
      }
    }
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,
                               const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::Constraint_SimOpt<Real>::applyInverseAdjointJacobian_1(iajv,v,u,z,tol);
    }
    else {
      if (isJ1zero_) {
        throw Exception::Zero(">>> (PDE_Constraint::applyAdjointInverseJacobian_1): Jacobian is zero!");
      }
      else {
        ROL::Ptr<Tpetra::MultiVector<> >    iajvf = getField(iajv);
        ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
        solveAdjoint(iajvf,vf);
      }
    }
  }


  /***************************************************************************/
  /* Output routines.                                                        */
  /***************************************************************************/
  void printMeshData(std::ostream &outStream) const {
    assembler_->printMeshData(outStream);
  }

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > matWriter;
    if (matJ1_ != ROL::nullPtr) {
      matWriter.writeSparseFile("jacobian1.txt", matJ1_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian1.txt");
      emptyfile.close();
    }
    if (matJ2_ != ROL::nullPtr) {
      matWriter.writeSparseFile("jacobian2.txt", matJ2_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian2.txt");
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

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<> > &vec,
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
