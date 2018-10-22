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

#ifndef ROL_PDEOPT_LINEARPDECONSTRAINT_H
#define ROL_PDEOPT_LINEARPDECONSTRAINT_H

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
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintSolverConstruct_Jacobian1        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Solver Construction Time for Jacobian1");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintSolverConstruct_AdjointJacobian1 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Solver Construction Time for Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintSolverSolve_Jacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Solver Solution Time for Jacobian1");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintSolverSolve_AdjointJacobian1     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Solver Solution Time for Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintApplyJacobian1                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Apply Jacobian1");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintApplyJacobian2                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Apply Jacobian2");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintApplyJacobian3                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Apply Jacobian3");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintApplyAdjointJacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Apply Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintApplyAdjointJacobian2            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Apply Adjoint Jacobian2");
    ROL::Ptr<Teuchos::Time> LinearPDEConstraintApplyAdjointJacobian3            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Linear PDE Constraint Apply Adjoint Jacobian3");
  }
}
#endif


template<class Real>
class Linear_PDE_Constraint : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<PDE<Real> > pde_;
  ROL::Ptr<Assembler<Real> > assembler_;
  ROL::Ptr<Solver<Real> >    solver_;
  bool initZvec_, initZpar_, isStochastic_;
  bool assembleRHS_, assembleJ1_, assembleJ2_, assembleJ3_, setSolver_;

  ROL::Ptr<Tpetra::MultiVector<> > cvec_;
  ROL::Ptr<Tpetra::MultiVector<> > uvec_;
  ROL::Ptr<Tpetra::MultiVector<> > zvec_;
  ROL::Ptr<std::vector<Real> >     zpar_;

  ROL::Ptr<Tpetra::MultiVector<> > vecR_;
  ROL::Ptr<Tpetra::CrsMatrix<> >   matJ1_;
  ROL::Ptr<Tpetra::CrsMatrix<> >   matJ2_;
  ROL::Ptr<Tpetra::MultiVector<> > vecJ3_;

  void assemble(const ROL::Ptr<const Tpetra::MultiVector<> > &zf,
                const ROL::Ptr<const std::vector<Real> > &zp) {
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
    // Assemble jacobian_1.
    if (assembleJ1_) {
      assembler_->assemblePDEJacobian1(matJ1_,pde_,uvec_,zvec_,zpar_);
    }
    assembleJ1_ = false;
    // Assemble jacobian_2.
    if (assembleJ2_ && zf != ROL::nullPtr) {
      assembler_->assemblePDEJacobian2(matJ2_,pde_,uvec_,zvec_,zpar_);
    }
    assembleJ2_ = false;
    // Assemble jacobian_3.
    if (assembleJ3_ && zp != ROL::nullPtr) {
      assembler_->assemblePDEJacobian3(vecJ3_,pde_,uvec_,zvec_,zpar_);
    }
    assembleJ3_ = false;
  }

  void setSolver(void) {
    solver_->setA(matJ1_);
    setSolver_ = false;
  }

  void solveForward(ROL::Ptr<Tpetra::MultiVector<> > &x,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &b) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintSolverConstruct_Jacobian1);
      #endif
      setSolver();
    }

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintSolverSolve_Jacobian1);
    #endif
    solver_->solve(x,b,false);
  }

  void solveAdjoint(ROL::Ptr<Tpetra::MultiVector<> > &x,
                    const ROL::Ptr<const Tpetra::MultiVector<> > &b) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintSolverConstruct_AdjointJacobian1);
      #endif
      setSolver();
    }

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintSolverSolve_AdjointJacobian1);
    #endif
    solver_->solve(x,b,true);
  }

  void applyJacobian1(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintApplyJacobian1);
    #endif
    matJ1_->apply(*v,*Jv);
  }

  void applyAdjointJacobian1(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintApplyAdjointJacobian1);
    #endif
    matJ1_->apply(*v,*Jv,Teuchos::TRANS);
  }

  void applyJacobian2(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                      const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintApplyJacobian2);
    #endif
    matJ2_->apply(*v,*Jv);
  }

  void applyAdjointJacobian2(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintApplyAdjointJacobian2);
    #endif
    matJ2_->apply(*v,*Jv,Teuchos::TRANS);
  }

  // Application routines
  void applyJacobian3(const ROL::Ptr<Tpetra::MultiVector<> > &Jv,
                      const ROL::Ptr<const std::vector<Real> > &v,
                      const bool zeroOut = true) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintApplyJacobian3);
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

  void applyAdjointJacobian3(const ROL::Ptr<std::vector<Real> > &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<> > &v) const {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::LinearPDEConstraintApplyAdjointJacobian3);
    #endif
    Teuchos::Array<Real> val(1,0);
    const size_t size = static_cast<size_t>(Jv->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      vecJ3_->subView(col)->dot(*v, val.view(0,1));
      (*Jv)[i] = val[0];
    }
  }

public:
  virtual ~Linear_PDE_Constraint() {}

  Linear_PDE_Constraint(const ROL::Ptr<PDE<Real> > &pde,
                        const ROL::Ptr<MeshManager<Real> > &meshMgr,
                        const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                        Teuchos::ParameterList &parlist,
                        std::ostream &outStream = std::cout,
                        const bool isStochastic = false)
    : ROL::Constraint_SimOpt<Real>(), pde_(pde),
      initZvec_(true), initZpar_(true), isStochastic_(isStochastic),
      assembleRHS_(true), assembleJ1_(true), assembleJ2_(true), assembleJ3_(true), setSolver_(true) {
    // Construct assembler.
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    // Construct solver.
    solver_ = ROL::makePtr<Solver<Real>>(parlist.sublist("Solver"));
    // Initialize zero vectors.
    cvec_ = assembler_->createResidualVector();
    uvec_ = assembler_->createStateVector();
    uvec_->putScalar(static_cast<Real>(0));
  }

  void setParameter(const std::vector<Real> &param) {
    if (isStochastic_) {
      ROL::Constraint_SimOpt<Real>::setParameter(param);
      pde_->setParameter(param);
      assembleRHS_ = true;
      assembleJ1_  = true;
      assembleJ2_  = true;
      assembleJ3_  = true;
      setSolver_   = true;
    }
  }

  const ROL::Ptr<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<PDE<Real> > getPDE(void) const {
    return pde_;
  }

  using ROL::Constraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {}

  using ROL::Constraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {}

  using ROL::Constraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
  }

  using ROL::Constraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >       cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);
    c.zero();
    cvec_->putScalar(static_cast<Real>(0));

    const Real one(1);
    cf->scale(one,*vecR_);
    if (zf != ROL::nullPtr) {
      applyJacobian2(cvec_,zf);
      cf->update(one,*cvec_,one);
    }
    if (zp != ROL::nullPtr) {
      applyJacobian3(cvec_,zp,false);
      cf->update(one,*cvec_,one);
    }
    applyJacobian1(cvec_,uf);
    cf->update(one,*cvec_,one);
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >       cf = getField(c);
    ROL::Ptr<Tpetra::MultiVector<> >       uf = getField(u);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);
    c.zero();
    cvec_->putScalar(static_cast<Real>(0));

    const Real one(1);
    cf->scale(one,*vecR_);
    if (zf != ROL::nullPtr) {
      applyJacobian2(cvec_,zf);
      cf->update(one,*cvec_,one);
    }
    if (zp != ROL::nullPtr) {
      applyJacobian3(cvec_,zp,false);
      cf->update(one,*cvec_,one);
    }

    solveForward(uf,cf);
    uf->scale(-one);

    applyJacobian1(cvec_,uf);
    cf->update(one,*cvec_,one);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);
    applyJacobian1(jvf,vf);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
    ROL::Ptr<Tpetra::MultiVector<> >      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    const Real one(1);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    if (vf != ROL::nullPtr) {
      applyJacobian2(cvec_,vf);
      jvf->update(one,*cvec_,one);
    }
    ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
    bool zeroOut = (vf == ROL::nullPtr);
    if (vp != ROL::nullPtr) {
      applyJacobian3(cvec_,vp,zeroOut);
      jvf->update(one,*cvec_,one);
    }
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >     ajvf = getField(ajv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);
    applyAdjointJacobian1(ajvf,vf);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >     ajvf = getField(ajv);
    ROL::Ptr<std::vector<Real> >         ajvp = getParameter(ajv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    if (zf != ROL::nullPtr) {
      applyAdjointJacobian2(ajvf,vf);
    }
    if (zp != ROL::nullPtr) {
      applyAdjointJacobian3(ajvp,vf);
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >     ijvf = getField(ijv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);
    solveForward(ijvf,vf);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,
                               const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> >    iajvf = getField(iajv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);
    solveAdjoint(iajvf,vf);
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
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
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
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
      ROL::Ptr<const ROL::TpetraMultiVector<Real> > xvec
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

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<> > xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      ROL::Ptr<ROL::TpetraMultiVector<Real> > xvec
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

  ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
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
    catch (std::exception &e) {
      xp = ROL::nullPtr;
    }
    return xp;
  }

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
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
    catch (std::exception &e) {
      xp = ROL::nullPtr;
    }
    return xp;
  }
};

#endif
