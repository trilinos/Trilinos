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

#ifndef ROL_PDEOPT_LINDYNCONSTRAINT_H
#define ROL_PDEOPT_LINDYNCONSTRAINT_H

#include "ROL_DynamicConstraint.hpp"
#include "assembler.hpp"
#include "solver.hpp"
#include "pdevector.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> LinDynConstraintSolverConstruct_Jacobian_un        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Solver Construction Time for Jacobian un");
    ROL::Ptr<Teuchos::Time> LinDynConstraintSolverSolve_Jacobian_un            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Solver Solution Time for Jacobian un");
    ROL::Ptr<Teuchos::Time> LinDynConstraintSolverSolve_AdjointJacobian_un     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Solver Solution Time for Adjoint Jacobian un");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyJacobian_uo                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Jacobian uo");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyJacobian_un                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Jacobian un");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyJacobian_zf                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Jacobian zf");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyJacobian_zp                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Jacobian zp");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyAdjointJacobian_uo            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Adjoint Jacobian uo");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyAdjointJacobian_un            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Adjoint Jacobian un");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyAdjointJacobian_zf            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Adjoint Jacobian zf");
    ROL::Ptr<Teuchos::Time> LinDynConstraintApplyAdjointJacobian_zp            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LinDyn Constraint Apply Adjoint Jacobian zp");
  }
}
#endif


template<class Real>
class LinDynConstraint : public ROL::DynamicConstraint<Real> {
private:
  const ROL::Ptr<DynamicPDE<Real>>    pde_;
  ROL::Ptr<Solver<Real>>           solver_;
  ROL::Ptr<Assembler<Real>>     assembler_;

  ROL::Ptr<Tpetra::MultiVector<>> uvec_;

  mutable ROL::Ptr<Tpetra::MultiVector<>> zvec_;
  mutable ROL::Ptr<std::vector<Real>>     zpar_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matJuo_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matJun_;
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecR_;
  mutable ROL::Ptr<Tpetra::CrsMatrix<>>   matJzf_;
  mutable ROL::Ptr<Tpetra::MultiVector<>> vecJzp_;
  mutable ROL::Ptr<Tpetra::MultiVector<>> cvec_;

  const bool isLTI_;
  mutable bool isAssembled_;
  mutable bool initialize_;

  void initialize(const ROL::Vector<Real> &z) const {
    if (initialize_) {
      // Initialize control vectors
      ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
      ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);
      if (zf != ROL::nullPtr) {
        zvec_ = assembler_->createControlVector();
        zvec_->putScalar(static_cast<Real>(0));
      }
      if (zp != ROL::nullPtr) {
        zpar_ = ROL::makePtr<std::vector<Real>>(zp->size(),static_cast<Real>(0));
      }
      initialize_ = false;
    }
  }
  

  void assemble(const ROL::Vector<Real> &z, const ROL::TimeStamp<Real> &ts) const {
    initialize(z);
    if (!isAssembled_) {
      // Assemble uold Jacobian.
      assembler_->assembleDynPDEJacobian_uo(matJuo_,pde_,ts,uvec_,uvec_,zvec_,zpar_);
      // Assemble unew Jacobian and initialize linear solver.
      assembler_->assembleDynPDEJacobian_un(matJun_,pde_,ts,uvec_,uvec_,zvec_,zpar_);
      // Assemble old affine term.
      assembler_->assembleDynPDEResidual(vecR_,pde_,ts,uvec_,uvec_,zvec_,zpar_);
      // Assemble old control Jacobian.
      if (zvec_ != ROL::nullPtr) {
        assembler_->assembleDynPDEJacobian_zf(matJzf_,pde_,ts,uvec_,uvec_,zvec_,zpar_);
      }
      if (zpar_ != ROL::nullPtr) {
        assembler_->assembleDynPDEJacobian_zp(vecJzp_,pde_,ts,uvec_,uvec_,zvec_,zpar_);
      }
      // Set solver with unew Jacobian
      setSolver();
      isAssembled_ = true;
    }
  }

  void setSolver(void) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintSolverConstruct_Jacobian_un);
    #endif
    solver_->setA(matJun_);
  }

  void solveForward(ROL::Ptr<Tpetra::MultiVector<>> &x,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &b) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintSolverSolve_Jacobian_un);
    #endif
    solver_->solve(x,b,false);
  }

  void solveAdjoint(ROL::Ptr<Tpetra::MultiVector<>> &x,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &b) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintSolverSolve_AdjointJacobian_un);
    #endif
    solver_->solve(x,b,true);
  }

  void applyJacobian_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                        const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyJacobian_uo);
      #endif
      matJuo_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyAdjointJacobian_uo);
      #endif
      matJuo_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }

  void applyJacobian_un(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                        const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyJacobian_un);
      #endif
      matJun_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyAdjointJacobian_un);
      #endif
      matJun_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }

  void applyJacobian_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                        const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyJacobian_zf);
      #endif
      matJzf_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyAdjointJacobian_zf);
      #endif
      matJzf_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }

  void applyJacobian_zp(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const std::vector<Real>> &v,
                        const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyJacobian_zp);
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

  void applyAdjointJacobian_zp(const ROL::Ptr<std::vector<Real> > &Jv,
                               const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinDynConstraintApplyAdjointJacobian_zp);
    #endif
    Teuchos::Array<Real> val(1,0);
    const size_t size = static_cast<size_t>(Jv->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      vecJzp_->subView(col)->dot(*v, val.view(0,1));
      (*Jv)[i] = val[0];
    }
  }

public:

  LinDynConstraint(const ROL::Ptr<DynamicPDE<Real>> &pde,
                   const ROL::Ptr<MeshManager<Real>> &meshMgr,
                   const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                   Teuchos::ParameterList &parlist,
                   bool isLTI = false,
                   std::ostream &outStream = std::cout)
    : pde_         (   pde ),
      isLTI_       ( isLTI ),
      isAssembled_ ( false ),
      initialize_  (  true ) {
    // Construct assembler.
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    // Construct solver.
    solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
    // Initialize state and constraint vectors.
    cvec_ = assembler_->createResidualVector();
    uvec_ = assembler_->createStateVector();
    uvec_->putScalar(static_cast<Real>(0));
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<DynamicPDE<Real>> getPDE(void) const {
    return pde_;
  }

  void update(const ROL::Vector<Real>    &uo,
              const ROL::Vector<Real>    &un,
              const ROL::Vector<Real>    &z,
              const ROL::TimeStamp<Real> &ts) {
    isAssembled_ = (!isAssembled_ ? isAssembled_ : isLTI_);
    assemble(z,ts);
  }

  void value(ROL::Vector<Real>    &c,
       const ROL::Vector<Real>    &uo,
       const ROL::Vector<Real>    &un,
       const ROL::Vector<Real>    &z,
       const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>        cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
    ROL::Ptr<const Tpetra::MultiVector<>> unf = getConstField(un);
    ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);

    const Real one(1);
    assemble(z,ts);
    c.zero();
    cvec_->putScalar(static_cast<Real>(0));

    // Apply old state contribution
    applyJacobian_uo(cf, uof, false);
    // Add affine/control terms
    cf->update(one, *vecR_, one);
    if (zf != ROL::nullPtr) {
      applyJacobian_zf(cvec_, zf, false);
      cf->update(one, *cvec_, one);
    }
    if (zp != ROL::nullPtr) {
      applyJacobian_zp(cvec_, zp, true);
      cf->update(one, *cvec_, one);
    }
    // Add new state contribution
    applyJacobian_un(cvec_, unf, false);
    cf->update(one, *cvec_, one);
  }

  void solve(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &uo,
             ROL::Vector<Real> &un,
       const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) {
    ROL::Ptr<Tpetra::MultiVector<>>        cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
    ROL::Ptr<Tpetra::MultiVector<>>       unf = getField(un);
    ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);

    const Real one(1);
    assemble(z,ts);
    c.zero();
    cvec_->putScalar(static_cast<Real>(0));

    // Apply old state contribution
    applyJacobian_uo(cf, uof, false);
    // Add affine/control terms
    cf->update(one, *vecR_, one);
    if (zf != ROL::nullPtr) {
      applyJacobian_zf(cvec_, zf, false);
      cf->update(one, *cvec_, one);
    }
    if (zp != ROL::nullPtr) {
      applyJacobian_zp(cvec_, zp, true);
      cf->update(one, *cvec_, one);
    }
    // Apply inverse of new state jacobian
    solveForward(unf, cf);
    unf->scale(static_cast<Real>(-1));
    // Compute residual
    applyJacobian_un(cvec_, unf, false);
    cf->update(one, *cvec_, one);
  }

  void applyJacobian_uo(ROL::Vector<Real>    &jv,
                  const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real>    &uo,
                  const ROL::Vector<Real>    &un,
                  const ROL::Vector<Real>    &z,
                  const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);

    assemble(z,ts);
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

    assemble(z,ts);
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

    const Real one(1);
    assemble(z,ts);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    if (vf != ROL::nullPtr) {
      applyJacobian_zf(cvec_,vf,false);
      jvf->update(one,*cvec_,one);
    }
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);
    bool zeroOut = (vf == ROL::nullPtr);
    if (vp != ROL::nullPtr) {
      applyJacobian_zp(cvec_,vp,zeroOut);
      jvf->update(one,*cvec_,one);
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

    assemble(z,ts);
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

    assemble(z,ts);
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
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    assemble(z,ts);
    if (zf != ROL::nullPtr) {
      applyJacobian_zf(ajvf,vf,true);
    }
    if (zp != ROL::nullPtr) {
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

    assemble(z,ts);
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

    assemble(z,ts);
    solveAdjoint(iajvf,vf);
  }


  void applyAdjointHessian_uo_uo(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_uo_un(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_uo_z(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_un_uo(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_un_un(ROL::Vector<Real>    &ahwv,
                           const ROL::Vector<Real>    &w,
                           const ROL::Vector<Real>    &v,
                           const ROL::Vector<Real>    &uo,
                           const ROL::Vector<Real>    &un,
                           const ROL::Vector<Real>    &z,
                           const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_un_z(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_z_uo(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_z_un(ROL::Vector<Real>    &ahwv,
                          const ROL::Vector<Real>    &w,
                          const ROL::Vector<Real>    &v,
                          const ROL::Vector<Real>    &uo,
                          const ROL::Vector<Real>    &un,
                          const ROL::Vector<Real>    &z,
                          const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
  }


  void applyAdjointHessian_z_z(ROL::Vector<Real>    &ahwv,
                         const ROL::Vector<Real>    &w,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ahwv.zero();
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
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> vecWriter;
    vecWriter.writeDenseFile(filename, vec);
    std::string mapfile = "map_" + filename;
    vecWriter.writeMapFile(mapfile, *(vec->getMap()));
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
