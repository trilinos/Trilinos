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

/*! \file  pdeconstraint.hpp
    \brief Defines the SimOpt constraint for the PDE-OPT examples.
*/

#ifndef ROL_PDEOPT_PDECONSTRAINT_H
#define ROL_PDEOPT_PDECONSTRAINT_H

#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"
#include "pde.hpp"
#include "assembler.hpp"
#include "solver.hpp"
#include "pdevector.hpp"


//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    Teuchos::RCP<Teuchos::Time> PDEConstraintSolverConstruct_Jacobian1        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Construction Time for Jacobian1");
    Teuchos::RCP<Teuchos::Time> PDEConstraintSolverConstruct_AdjointJacobian1 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Construction Time for Adjoint Jacobian1");
    Teuchos::RCP<Teuchos::Time> PDEConstraintSolverSolve_Jacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Solution Time for Jacobian1");
    Teuchos::RCP<Teuchos::Time> PDEConstraintSolverSolve_AdjointJacobian1     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Solver Solution Time for Adjoint Jacobian1");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyJacobian1                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Jacobian1");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyJacobian2                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Jacobian2");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyJacobian3                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Jacobian3");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyAdjointJacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Adjoint Jacobian1");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyAdjointJacobian2            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Adjoint Jacobian2");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyAdjointJacobian3            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Adjoint Jacobian3");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian11                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian11");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian12                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian12");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian13                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian13");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian21                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian21");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian22                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian22");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian23                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian23");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian31                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian31");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian32                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian32");
    Teuchos::RCP<Teuchos::Time> PDEConstraintApplyHessian33                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Constraint Apply Hessian33");
  }
}
#endif


template<class Real>
class PDE_Constraint : public ROL::ParametrizedEqualityConstraint_SimOpt<Real> {
private:
  const Teuchos::RCP<PDE<Real> > pde_;
  Teuchos::RCP<Assembler<Real> > assembler_;
  Teuchos::RCP<Solver<Real> > solver_;

  Teuchos::RCP<Tpetra::MultiVector<> >           vecR_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             matJ1_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             matJ2_;
  Teuchos::RCP<Tpetra::MultiVector<> >           vecJ3_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             matH11_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             matH12_;
  Teuchos::RCP<Tpetra::MultiVector<> >           vecH13_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             matH21_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             matH22_;
  Teuchos::RCP<Tpetra::MultiVector<> >           vecH23_;
  Teuchos::RCP<std::vector<std::vector<Real> > > matH33_;

  bool computeJ1_,  computeJ2_,  computeJ3_;
  bool computeH11_, computeH12_, computeH13_;
  bool computeH21_, computeH22_, computeH23_;
  bool computeH31_, computeH32_, computeH33_;
  bool setSolver_;
  bool isJ1zero_, isJ1notImplemented_;
  bool isJ2zero_, isJ2notImplemented_;
  bool isJ3zero_, isJ3notImplemented_;

  // Assemble the PDE Jacobian.
  void assembleJ1(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    if ( !isJ1zero_ && !isJ1notImplemented_ ) {
      try {
        if (computeJ1_) {
          Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
          Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
          Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

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
          Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
          Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
          Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

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
          Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
          Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
          Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

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

  void applyJacobian1(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                      const Teuchos::RCP<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian1);
    #endif
    matJ1_->apply(*v,*Jv);
  }

  void applyAdjointJacobian1(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                             const Teuchos::RCP<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyAdjointJacobian1);
    #endif
    matJ1_->apply(*v,*Jv,Teuchos::TRANS);
  }

  void applyJacobian2(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                      const Teuchos::RCP<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian2);
    #endif
    matJ2_->apply(*v,*Jv);
  }

  void applyAdjointJacobian2(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                             const Teuchos::RCP<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyAdjointJacobian2);
    #endif
    matJ2_->apply(*v,*Jv,Teuchos::TRANS);
  }

  // Application routines
  void applyJacobian3(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                      const Teuchos::RCP<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian3);
    #endif
    if ( zeroOut ) {
      Jv->putScalar(static_cast<Real>(0));
    }
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Jv->update((*v)[i],*(vecJ3_->subView(col)),static_cast<Real>(1));
    }
  }

  void applyAdjointJacobian3(const Teuchos::RCP<std::vector<Real> > &Jv,
                             const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyJacobian3);
    #endif
    Teuchos::Array<Real> val(1,0);
    const size_t size = static_cast<size_t>(Jv->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      vecJ3_->subView(col)->dot(*v, val.view(0,1));
      (*Jv)[i] = val[0];
    }
  }

  void applyHessian13(const Teuchos::RCP<std::vector<Real> > &Hv,
                      const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian13);
    #endif
    const size_t size = static_cast<size_t>(Hv->size());
    if ( zeroOut ) {
      Hv->assign(size,static_cast<Real>(0));
    }
    Teuchos::Array<Real> val(1,0);
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      vecH13_->subView(col)->dot(*v, val.view(0,1));
      (*Hv)[i] += val[0];
    }
  }

  void applyHessian31(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                      const Teuchos::RCP<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian31);
    #endif
    if ( zeroOut ) {
      Hv->putScalar(static_cast<Real>(0));
    }
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Hv->update((*v)[i],*(vecH13_->subView(col)),static_cast<Real>(1));
    }
  }

  void applyHessian23(const Teuchos::RCP<std::vector<Real> > &Hv,
                      const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian23);
    #endif
    const size_t size = static_cast<size_t>(Hv->size());
    if ( zeroOut ) {
      Hv->assign(size,static_cast<Real>(0));
    }
    Teuchos::Array<Real> val(1,0);
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      vecH23_->subView(col)->dot(*v, val.view(0,1));
      (*Hv)[i] += val[0];
    }
  }

  void applyHessian32(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                      const Teuchos::RCP<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian32);
    #endif
    if ( zeroOut ) {
      Hv->putScalar(static_cast<Real>(0));
    }
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Hv->update((*v)[i],*(vecH23_->subView(col)),static_cast<Real>(1));
    }
  }

  void applyHessian33(const Teuchos::RCP<std::vector<Real> > &Hv,
                      const Teuchos::RCP<const std::vector<Real> > &v,
                      const bool zeroOut = true ) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian33);
    #endif
    const int size = Hv->size();
    if ( zeroOut ) {
      Hv->assign(size,static_cast<Real>(0));
    }
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        (*Hv)[i] += (*matH33_)[i][j]*(*v)[j];
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
  void solveForward(Teuchos::RCP<Tpetra::MultiVector<> > &x,
                    const Teuchos::RCP<const Tpetra::MultiVector<> > &b) {
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
  void solveAdjoint(Teuchos::RCP<Tpetra::MultiVector<> > &x,
                    const Teuchos::RCP<const Tpetra::MultiVector<> > &b) {
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
  PDE_Constraint(const Teuchos::RCP<PDE<Real> > &pde,
                 const Teuchos::RCP<MeshManager<Real> > &meshMgr,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 Teuchos::ParameterList &parlist,
                 std::ostream &outStream = std::cout)
    : ROL::ParametrizedEqualityConstraint_SimOpt<Real>(),
      pde_(pde),
      computeJ1_(true),  computeJ2_(true),  computeJ3_(true),
      computeH11_(true), computeH12_(true), computeH13_(true),
      computeH21_(true), computeH22_(true), computeH23_(true),
      computeH31_(true), computeH32_(true), computeH33_(true),
      setSolver_(true),
      isJ1zero_(false),  isJ1notImplemented_(false),
      isJ2zero_(false),  isJ2notImplemented_(false),
      isJ3zero_(false),  isJ3notImplemented_(false) {
    assembler_ = Teuchos::rcp(new Assembler<Real>(pde_->getFields(),meshMgr,comm,parlist,outStream));
    assembler_->setCellNodes(*pde_);
    solver_ = Teuchos::rcp(new Solver<Real>(parlist.sublist("Solver")));
  }

  PDE_Constraint(const Teuchos::RCP<PDE<Real> >       &pde,
                 const Teuchos::RCP<Assembler<Real> > &assembler,
                 const Teuchos::RCP<Solver<Real> >    &solver)
    : ROL::ParametrizedEqualityConstraint_SimOpt<Real>(),
      pde_(pde), assembler_(assembler), solver_(solver),
      computeJ1_(true),  computeJ2_(true),  computeJ3_(true),
      computeH11_(true), computeH12_(true), computeH13_(true),
      computeH21_(true), computeH22_(true), computeH23_(true),
      computeH31_(true), computeH32_(true), computeH33_(true),
      setSolver_(true),
      isJ1zero_(false),  isJ1notImplemented_(false),
      isJ2zero_(false),  isJ2notImplemented_(false),
      isJ3zero_(false),  isJ3notImplemented_(false) {
    assembler_->setCellNodes(*pde_);
  }

  void setParameter(const std::vector<Real> &param) {
    computeJ1_  = true; computeJ2_  = true; computeJ3_  = true;
    computeH11_ = true; computeH12_ = true; computeH13_ = true;
    computeH21_ = true; computeH22_ = true; computeH23_ = true;
    computeH31_ = true; computeH32_ = true; computeH33_ = true;
    setSolver_  = true;
    ROL::ParametrizedEqualityConstraint_SimOpt<Real>::setParameter(param);
    pde_->setParameter(param);
  }

  const Teuchos::RCP<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {
    computeJ1_ = (flag ? true : computeJ1_);
    computeJ2_ = (flag ? true : computeJ2_);
    computeJ3_ = (flag ? true : computeJ3_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    computeJ1_ = (flag ? true : computeJ1_);
    computeJ2_ = (flag ? true : computeJ2_);
    computeJ3_ = (flag ? true : computeJ3_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
    computeH11_ = (flag ? true : computeH11_);
    computeH12_ = (flag ? true : computeH12_);
    computeH13_ = (flag ? true : computeH13_);
    computeH21_ = (flag ? true : computeH21_);
    computeH22_ = (flag ? true : computeH22_);
    computeH23_ = (flag ? true : computeH23_);
    computeH31_ = (flag ? true : computeH31_);
    computeH32_ = (flag ? true : computeH32_);
    computeH33_ = (flag ? true : computeH33_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> >       cf = getField(c);
    Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assembler_->assemblePDEResidual(vecR_,pde_,uf,zf,zp);
    cf->scale(static_cast<Real>(1),*vecR_);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_1(jv,v,u,z,tol);
    }
    else {
      if (isJ1zero_) {
        jv.zero();
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >      jvf = getField(jv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
        applyJacobian1(jvf,vf);
      }
    }
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    int NotImplemented(0), IsZero(0);
    // Apply Jacobian of field controls to vector
    assembleJ2(u,z);
    if (isJ2notImplemented_) {
      NotImplemented++;
    }
    else {
      if (isJ2zero_) {
        IsZero++;
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >      jvf = getField(jv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
        applyJacobian2(jvf,vf);
      }
    }
    // Apply Jacobian of parametric controls to vector
    assembleJ3(u,z);
    if (isJ3notImplemented_) {
      NotImplemented++;
    }
    else {
      if (isJ3zero_) {
        IsZero++;
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >  jvf = getField(jv);
        Teuchos::RCP<const std::vector<Real> > vp = getConstParameter(v);
        applyJacobian3(jvf,vp);
      }
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      jv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented == 2 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_2(jv,v,u,z,tol);
    }
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointJacobian_1(ajv,v,u,z,tol);
    }
    else {
      if (isJ1zero_) {
        ajv.zero();
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >     ajvf = getField(ajv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
        applyAdjointJacobian1(ajvf,vf);
      }
    }
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    int NotImplemented(0), IsZero(0);
    // Apply Jacobian of field controls to vector
    assembleJ2(u,z);
    if (isJ2notImplemented_) {
      NotImplemented++;
    }
    else {
      if (isJ2zero_) {
        IsZero++;
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >     ajvf = getField(ajv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
        applyAdjointJacobian2(ajvf,vf);
      }
    }
    // Apply Jacobian of parametric controls to vector
    assembleJ3(u,z);
    if (isJ3notImplemented_) {
      NotImplemented++;
    }
    else {
      if (isJ3zero_) {
        IsZero++;
      }
      else {
        Teuchos::RCP<std::vector<Real> >         ajvp = getParameter(ajv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
        applyAdjointJacobian3(ajvp,vf);
      }
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      ajv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented == 2 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointJacobian_2(ajv,v,u,z,tol);
    }
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
    try {
      if (computeH11_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian11(matH11_,pde_,wf,uf,zf,zp);
        computeH11_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH11_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_11(ahwv,w,v,u,z,tol);
    }
    if ( computeH11_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian11);
      #endif
      matH11_->apply(*vf,*ahwvf);
    }
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    // Apply Jacobian of field controls to vector
    ahwv.zero();
    int NotImplemented(0), IsZero(0);
    try {
      if (computeH12_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian12(matH12_,pde_,wf,uf,zf,zp);
        computeH12_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeH12_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeH12_ = true;
      NotImplemented++;
    }
    if ( !computeH12_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian12);
      #endif
      matH12_->apply(*vf,*ahwvf);
    }
    // Apply Jacobian of parametric controls to vector
    try {
      if (computeH13_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian13(vecH13_,pde_,wf,uf,zf,zp);
        computeH13_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeH13_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeH13_ = true;
      NotImplemented++;
    }
    if ( !computeH13_ ) {
      Teuchos::RCP<std::vector<Real> >        ahwvp = getParameter(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      applyHessian13(ahwvp,vf,true);
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      ahwv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented == 2 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_12(ahwv,w,v,u,z,tol);
    }
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    // Apply Jacobian of field controls to vector
    ahwv.zero();
    int NotImplemented(0), IsZero(0);
    try {
      if (computeH21_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian21(matH21_,pde_,wf,uf,zf,zp);
        computeH21_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      ahwv.zero();
      computeH21_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      ahwv.zero();
      computeH21_ = true;
      NotImplemented++;
    }
    if ( !computeH21_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian21);
      #endif
      matH21_->apply(*vf,*ahwvf);
    }
    // Apply Jacobian of parametric controls to vector
    try {
      if (computeH31_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian31(vecH13_,pde_,wf,uf,zf,zp);
        computeH31_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeH31_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeH31_ = true;
      NotImplemented++;
    }
    if ( !computeH31_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      applyHessian31(ahwvf,vp,false);
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      ahwv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented == 2 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_21(ahwv,w,v,u,z,tol);
    }
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    // Apply Hessian wrt field controls to field control vector
    ahwv.zero();
    int NotImplemented(0), IsZero(0);
    try {
      if (computeH22_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian22(matH22_,pde_,wf,uf,zf,zp);
        computeH22_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      ahwv.zero();
      computeH22_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      ahwv.zero();
      computeH22_ = true;
      NotImplemented++;
    }
    if ( !computeH22_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEConstraintApplyHessian22);
      #endif
      matH22_->apply(*vf,*ahwvf);
    }
    // Apply Hessian wrt field controls to parametric control vector
    try {
      if (computeH23_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian23(vecH23_,pde_,wf,uf,zf,zp);
        computeH23_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      Teuchos::RCP<std::vector<Real> > ahwvp = getParameter(ahwv);
      if ( ahwvp != Teuchos::null ) {
        const int size = ahwvp->size();
        ahwvp->assign(size,static_cast<Real>(0));
      }
      computeH23_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      Teuchos::RCP<std::vector<Real> > ahwvp = getParameter(ahwv);
      if ( ahwvp != Teuchos::null ) {
        const int size = ahwvp->size();
        ahwvp->assign(size,static_cast<Real>(0));
      }
      computeH23_ = true;
      NotImplemented++;
    }
    if ( !computeH23_ ) {
      Teuchos::RCP<std::vector<Real> >        ahwvp = getParameter(ahwv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      applyHessian23(ahwvp,vf,true);
    }
    // Apply Hessian wrt parametric controls to field control vector
    try {
      if (computeH32_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian32(vecH23_,pde_,wf,uf,zf,zp);
        computeH32_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeH32_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeH32_ = true;
      NotImplemented++;
    }
    if ( !computeH32_ ) {
      Teuchos::RCP<Tpetra::MultiVector<> >    ahwvf = getField(ahwv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      applyHessian32(ahwvf,vp,false);
    }
    // Apply Hessian wrt parametric controls to field control vector
    try {
      if (computeH33_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wf = getConstField(w);
        Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
        Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
        Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

        assembler_->assemblePDEHessian33(matH33_,pde_,wf,uf,zf,zp);
        computeH33_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeH33_ = true;
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      computeH33_ = true;
      NotImplemented++;
    }
    if ( !computeH33_ ) {
      Teuchos::RCP<std::vector<Real> >        ahwvp = getParameter(ahwv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      applyHessian33(ahwvp,vp,false);
    }
    // Zero Jacobian if all routines return Exception::Zero
    if ( IsZero > 0 && (IsZero + NotImplemented == 4) ) {
      ahwv.zero();
    }
    // Default to finite differences if all routines return Exception::NotImplemented
    if ( NotImplemented == 4 ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_22(ahwv,w,v,u,z,tol);
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    assembleJ1(u,z);
    if (isJ1notImplemented_) {
      ROL::EqualityConstraint_SimOpt<Real>::applyInverseJacobian_1(ijv,v,u,z,tol);
    }
    else {
      if (isJ1zero_) {
        throw Exception::Zero(">>> (PDE_Constraint::applyInverseJacobian_1): Jacobian is zero!");
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >     ijvf = getField(ijv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
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
      ROL::EqualityConstraint_SimOpt<Real>::applyInverseAdjointJacobian_1(iajv,v,u,z,tol);
    }
    else {
      if (isJ1zero_) {
        throw Exception::Zero(">>> (PDE_Constraint::applyAdjointInverseJacobian_1): Jacobian is zero!");
      }
      else {
        Teuchos::RCP<Tpetra::MultiVector<> >    iajvf = getField(iajv);
        Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
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
    if (matJ1_ != Teuchos::null) {
      matWriter.writeSparseFile("jacobian1.txt", matJ1_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian1.txt");
      emptyfile.close();
    }
    if (matJ2_ != Teuchos::null) {
      matWriter.writeSparseFile("jacobian2.txt", matJ2_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian2.txt");
      emptyfile.close();
    }
    if (vecR_ != Teuchos::null) {
      matWriter.writeDenseFile("residual.txt", vecR_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("residual.txt");
      emptyfile.close();
    }
  }

  void outputTpetraVector(const Teuchos::RCP<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }
  /***************************************************************************/
  /* End of output routines.                                                 */
  /***************************************************************************/

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<const ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getField();
      if (xvec == Teuchos::null) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    Teuchos::RCP<Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getField();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const std::vector<Real> > xp;
    try {
      Teuchos::RCP<const ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }

  Teuchos::RCP<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    Teuchos::RCP<std::vector<Real> > xp;
    try {
      Teuchos::RCP<ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }
};

#endif
