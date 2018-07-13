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

/*! \file  constraint.hpp
    \brief Defines the SimOpt constraint for the 'poisson' example.
*/

#ifndef ROL_PDEOPT_LTICONSTRAINT_H
#define ROL_PDEOPT_LTICONSTRAINT_H

#include "ROL_DynamicConstraint.hpp"
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
    ROL::Ptr<Teuchos::Time> LTIConstraintSolverConstruct_Jacobian_un        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Solver Construction Time for Jacobian un");
    ROL::Ptr<Teuchos::Time> LTIConstraintSolverSolve_Jacobian_un            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Solver Solution Time for Jacobian un");
    ROL::Ptr<Teuchos::Time> LTIConstraintSolverSolve_AdjointJacobian_un     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Solver Solution Time for Adjoint Jacobian un");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyJacobian_uo                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Jacobian uo");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyJacobian_un                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Jacobian un");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyJacobian_zf                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Jacobian zf");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyJacobian_zp                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Jacobian zp");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyAdjointJacobian_uo            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Adjoint Jacobian uo");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyAdjointJacobian_un            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Adjoint Jacobian un");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyAdjointJacobian_zf            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Adjoint Jacobian zf");
    ROL::Ptr<Teuchos::Time> LTIConstraintApplyAdjointJacobian_zp            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: LTI Constraint Apply Adjoint Jacobian zp");
  }
}
#endif


template<class Real>
class LTI_Constraint : public ROL::DynamicConstraint<Real> {
private:
  const ROL::Ptr<PDE<Real>>       pde_;
  const ROL::Ptr<PDE<Real>>      mass_;
  ROL::Ptr<Solver<Real>>       solver_;
  ROL::Ptr<Assembler<Real>> assembler_;

  Real theta_, dt_;

  ROL::Ptr<Tpetra::MultiVector<>> uvec_;
  ROL::Ptr<Tpetra::MultiVector<>> zvec_;
  ROL::Ptr<std::vector<Real>>     zpar_;

  ROL::Ptr<Tpetra::MultiVector<>>   vecR_;
  ROL::Ptr<Tpetra::CrsMatrix<>>     matM_;
  ROL::Ptr<Tpetra::CrsMatrix<>>     matA_;
  ROL::Ptr<Tpetra::CrsMatrix<>>   matJuo_;
  ROL::Ptr<Tpetra::CrsMatrix<>>   matJun_;
  ROL::Ptr<Tpetra::CrsMatrix<>>   matJzf_;
  ROL::Ptr<Tpetra::MultiVector<>> vecJzp_;

  bool initZvec_, initZpar_;
  bool assembleRHS_, assembleM_, assembleA_;
  bool assembleJuo_, assembleJun_, assembleJzf_;
  bool assembleJzp_, setSolver_;

  mutable ROL::Ptr<Tpetra::MultiVector<>> cvec_;

  void assemble(const ROL::Vector<Real> &z) {
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    // Initialize field component of z.
    if (initZvec_ && zf != ROL::nullPtr) {
      zvec_ = assembler_->createControlVector();
      zvec_->putScalar(static_cast<Real>(0));
    }
    initZvec_ = false;
    // Initialize parameter component of z.
    if (initZpar_ && zp != ROL::nullPtr) {
      zpar_ = ROL::makePtr<std::vector<Real>>(zp->size(),static_cast<Real>(0));
    }
    initZpar_ = false;
    // Assemble affine term.
    if (assembleRHS_) {
      assembler_->assemblePDEResidual(vecR_,pde_,uvec_,zvec_,zpar_);
    }
    assembleRHS_ = false;
    // Assemble mass matrix.
    if (assembleM_) {
      assembler_->assemblePDEJacobian1(matM_,mass_,uvec_,zvec_,zpar_);
    }
    assembleM_ = false;
    // Assemble stiffness matrix.
    if (assembleA_) {
      assembler_->assemblePDEJacobian1(matA_,pde_,uvec_,zvec_,zpar_);
    }
    assembleA_ = false;
    // Assemble jacobian_uo.
    if (assembleJuo_) {
      const Real one(1);
      matJuo_ = ROL::dynamicPtrCast<Tpetra::CrsMatrix<>>(matM_->add(-dt_*(one-theta_),*matA_,one,matM_->getDomainMap(),matM_->getRangeMap(),Teuchos::null));
    }
    assembleJuo_ = false;
    // Assemble jacobian_un.
    if (assembleJun_) {
      const Real one(1);
      matJun_ = ROL::dynamicPtrCast<Tpetra::CrsMatrix<>>(matM_->add(dt_*theta_,*matA_,one,matM_->getDomainMap(),matM_->getRangeMap(),Teuchos::null));
    }
    assembleJun_ = false;
    // Assemble jacobian_zf.
    if (assembleJzf_ && zf != ROL::nullPtr) {
      assembler_->assemblePDEJacobian2(matJzf_,pde_,uvec_,zvec_,zpar_);
    }
    assembleJzf_ = false;
    // Assemble jacobian_3.
    if (assembleJzp_ && zp != ROL::nullPtr) {
      assembler_->assemblePDEJacobian3(vecJzp_,pde_,uvec_,zvec_,zpar_);
    }
    assembleJzp_ = false;
  }

  void setSolver(void) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintSolverConstruct_Jacobian_un);
      #endif
      solver_->setA(matJun_);
    }
    setSolver_ = false;
  }

  void solveForward(ROL::Ptr<Tpetra::MultiVector<>> &x,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &b) const {

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintSolverSolve_Jacobian_un);
    #endif
    solver_->solve(x,b,false);
  }

  void solveAdjoint(ROL::Ptr<Tpetra::MultiVector<>> &x,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &b) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintSolverSolve_AdjointJacobian_un);
    #endif
    solver_->solve(x,b,true);
  }

  void applyJacobian_uo(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                        const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyJacobian_uo);
      #endif
      matJuo_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyAdjointJacobian_uo);
      #endif
      matJuo_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }

  void applyJacobian_un(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                        const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyJacobian_un);
      #endif
      matJun_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyAdjointJacobian_un);
      #endif
      matJun_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }

  void applyJacobian_zf(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                        const bool trans = false) const {
    if (!trans) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyJacobian_zf);
      #endif
      matJzf_->apply(*v,*Jv,Teuchos::NO_TRANS);
    }
    else {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyAdjointJacobian_zf);
      #endif
      matJzf_->apply(*v,*Jv,Teuchos::TRANS);
    }
  }

  void applyJacobian_zp(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                        const ROL::Ptr<const std::vector<Real>> &v,
                        const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyJacobian_zp);
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
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LTIConstraintApplyAdjointJacobian_zp);
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

  LTI_Constraint(const ROL::Ptr<PDE<Real>> &pde,
                 const ROL::Ptr<PDE<Real>> &mass,
                 const ROL::Ptr<MeshManager<Real>> &meshMgr,
                 const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                 const ROL::Vector<Real> &z,
                 Teuchos::ParameterList &parlist,
                 std::ostream &outStream = std::cout)
    : pde_         (  pde ),
      mass_        ( mass ),
      initZvec_    ( true ),
      initZpar_    ( true ),
      assembleRHS_ ( true ),
      assembleM_   ( true ),
      assembleA_   ( true ), 
      assembleJuo_ ( true ),
      assembleJun_ ( true ),
      assembleJzf_ ( true ),
      assembleJzp_ ( true ),
      setSolver_   ( true ) {
    // Get time discretization parameters
    Real  T = parlist.sublist("Time Discretization").get("End Time",             1.0);
    int  nt = parlist.sublist("Time Discretization").get("Number of Time Steps", 100);
    theta_  = parlist.sublist("Time Discretization").get("Theta",                1.0);
    dt_     = T/static_cast<Real>(nt-1);
    // Construct assembler.
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    assembler_->setCellNodes(*mass_);
    // Construct solver.
    solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
    // Initialize zero vectors.
    cvec_ = assembler_->createResidualVector();
    uvec_ = assembler_->createStateVector();
    uvec_->putScalar(static_cast<Real>(0));
    // Assemble matrices
    assemble(z);
    setSolver();
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<PDE<Real>> getPDE(void) const {
    return pde_;
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

    c.zero();
    cvec_->putScalar(static_cast<Real>(0));

    const Real one(1);
    // Old state contribution
    applyJacobian_uo(cf,uof,false);
    // Add load contribution
    cf->update(dt_,*vecR_,one);
    if (zf != ROL::nullPtr) {
      applyJacobian_zf(cvec_,zf,false);
      cf->update(dt_,*cvec_,one);
    }
    if (zp != ROL::nullPtr) {
      applyJacobian_zp(cvec_,zp,false);
      cf->update(dt_,*cvec_,one);
    }
    // Add new contribution
    applyJacobian_un(cvec_,unf,false);
    cf->update(one,*cvec_,-one);
  }

  void solve(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &uo,
             ROL::Vector<Real> &un,
       const ROL::Vector<Real> &z,
       const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>        cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<>> uof = getConstField(uo);
    ROL::Ptr<Tpetra::MultiVector<>>       unf = getField(un);
    ROL::Ptr<const Tpetra::MultiVector<>>  zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>      zp = getConstParameter(z);

    c.zero();
    cvec_->putScalar(static_cast<Real>(0));

    const Real one(1);
    // Old state contribution
    applyJacobian_uo(cf,uof,false);
    // Load contribution
    cf->update(dt_,*vecR_,one);
    if (zf != ROL::nullPtr) {
      applyJacobian_zf(cvec_,zf,false);
      cf->update(dt_,*cvec_,one);
    }
    if (zp != ROL::nullPtr) {
      applyJacobian_zp(cvec_,zp,false);
      cf->update(dt_,*cvec_,one);
    }
    // Apply inverse of new state jacobian
    solveForward(unf,cf);
    // Compute residual
    applyJacobian_un(cvec_,unf,false);
    cf->update(one,*cvec_,-one);
  }

  void applyJacobian_uo(ROL::Vector<Real>    &jv,
                  const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real>    &uo,
                  const ROL::Vector<Real>    &un,
                  const ROL::Vector<Real>    &z,
                  const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);

    const Real one(1);
    applyJacobian_uo(jvf,vf,false);
    jvf->scale(-one);
  }

  void applyJacobian_un(ROL::Vector<Real>    &jv,
                  const ROL::Vector<Real>    &v,
                  const ROL::Vector<Real>    &uo,
                  const ROL::Vector<Real>    &un,
                  const ROL::Vector<Real>    &z,
                  const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);

    applyJacobian_un(jvf,vf,false);
  }


  void applyJacobian_z(ROL::Vector<Real>    &jv,
                 const ROL::Vector<Real>    &v,
                 const ROL::Vector<Real>    &uo,
                 const ROL::Vector<Real>    &un,
                 const ROL::Vector<Real>    &z,
                 const ROL::TimeStamp<Real> &tol) const {
    jv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);

    const Real one(1);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    if (vf != ROL::nullPtr) {
      applyJacobian_zf(cvec_,vf,false);
      jvf->update(-dt_,*cvec_,one);
    }
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);
    bool zeroOut = (vf == ROL::nullPtr);
    if (vp != ROL::nullPtr) {
      applyJacobian_zp(cvec_,vp,zeroOut);
      jvf->update(-dt_,*cvec_,one);
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

    const Real one(1);
    applyJacobian_uo(ajvf,vf,true);
    ajvf->scale(-one);
  }


  void applyAdjointJacobian_un(ROL::Vector<Real>    &ajv,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>     ajvf = getField(ajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);

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

    if (zf != ROL::nullPtr) {
      applyJacobian_zf(ajvf,vf,true);
    }
    if (zp != ROL::nullPtr) {
      applyAdjointJacobian_zp(ajvp,vf);
    }
    ajv.scale(-dt_);
  }


  void applyInverseJacobian_un(ROL::Vector<Real>    &ijv,
                         const ROL::Vector<Real>    &v,
                         const ROL::Vector<Real>    &uo,
                         const ROL::Vector<Real>    &un,
                         const ROL::Vector<Real>    &z,
                         const ROL::TimeStamp<Real> &ts) const {
    ROL::Ptr<Tpetra::MultiVector<>>     ijvf = getField(ijv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);

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
    if (matM_ != ROL::nullPtr) {
      matWriter.writeSparseFile("mass.txt", matM_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("mass.txt");
      emptyfile.close();
    }
    if (matA_ != ROL::nullPtr) {
      matWriter.writeSparseFile("stiffness.txt", matA_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("stiffness.txt");
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
  }
  /***************************************************************************/
  /* End of output routines.                                                 */
  /***************************************************************************/

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<>> getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::Ptr<const ROL::TpetraMultiVector<Real>> xvec
        = dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
      if (xvec == ROL::nullPtr) {
        xp = ROL::nullPtr;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  ROL::Ptr<Tpetra::MultiVector<>> getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::Ptr<ROL::TpetraMultiVector<Real>> xvec
        = dynamic_cast<PDE_OptVector<Real>&>(x).getField();
      if ( xvec == ROL::nullPtr ) {
        xp = ROL::nullPtr;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  ROL::Ptr<const std::vector<Real>> getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real>> xp;
    try {
      ROL::Ptr<const ROL::StdVector<Real>> xvec
        = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
      if ( xvec == ROL::nullPtr ) {
        xp = ROL::nullPtr;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = ROL::nullPtr;
    }
    return xp;
  }

  ROL::Ptr<std::vector<Real>> getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real>> xp;
    try {
      ROL::Ptr<ROL::StdVector<Real>> xvec
        = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
      if ( xvec == ROL::nullPtr ) {
        xp = ROL::nullPtr;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = ROL::nullPtr;
    }
    return xp;
  }
};

#endif
