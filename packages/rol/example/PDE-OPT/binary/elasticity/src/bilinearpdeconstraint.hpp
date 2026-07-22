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

#ifndef ROL_PDEOPT_BILINEARPDECONSTRAINT_H
#define ROL_PDEOPT_BILINEARPDECONSTRAINT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "../../../TOOLS/solver.hpp"
#include "../../../TOOLS/pdevector.hpp"
#include "femdata.hpp"

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintSolverConstruct_Jacobian1        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Solver Construction Time for Jacobian1");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintSolverConstruct_AdjointJacobian1 = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Solver Construction Time for Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintSolverSolve_Jacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Solver Solution Time for Jacobian1");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintSolverSolve_AdjointJacobian1     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Solver Solution Time for Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintApplyJacobian1                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Apply Jacobian1");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintApplyJacobian2                   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Apply Jacobian2");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintApplyAdjointJacobian1            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Apply Adjoint Jacobian1");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintApplyAdjointJacobian2            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Apply Adjoint Jacobian2");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintApplyHessian                     = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Apply Hessian");
    ROL::Ptr<Teuchos::Time> BilinearPDEConstraintApplyAdjontHessian               = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Bilinear PDE Constraint Apply Adjoint Hessian");
  }
}
#endif


template<class Real>
class Bilinear_PDE_Constraint : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<FEM_Data<Real>> fem_;
  ROL::Ptr<Solver<Real>>         solver_;
  bool setSolver_, assembleJ1_, buildState_, buildAdjoint_;

  ROL::Ptr<Tpetra::MultiVector<>>                                        cvec_;
  ROL::Ptr<Tpetra::RowMatrix<>>                                          matJ1_;
  ROL::Ptr<Tpetra::MultiVector<>>                                        S0u_;
  std::vector<std::vector<std::vector<ROL::Ptr<Tpetra::MultiVector<>>>>> Su_;
  std::vector<std::vector<std::vector<ROL::Ptr<Tpetra::MultiVector<>>>>> Sl_;

  unsigned M_, N_, T_;

  void setSolver(void) {
    ROL::Ptr<Tpetra::CrsMatrix<>> tmp = ROL::staticPtrCast<Tpetra::CrsMatrix<>>(matJ1_);
    solver_->setA(tmp);
    setSolver_ = false;
  }

  void solveForward(ROL::Ptr<Tpetra::MultiVector<>> &x,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &b) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintSolverConstruct_Jacobian1);
      #endif
      setSolver();
    }

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintSolverSolve_Jacobian1);
    #endif
    solver_->solve(x,b,false);
  }

  void solveAdjoint(ROL::Ptr<Tpetra::MultiVector<>> &x,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &b) {
    if (setSolver_) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintSolverConstruct_AdjointJacobian1);
      #endif
      setSolver();
    }

    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintSolverSolve_AdjointJacobian1);
    #endif
    solver_->solve(x,b,true);
  }

  void buildSu(const ROL::Ptr<const Tpetra::MultiVector<>> &u) {
    if (buildState_) {
      fem_->applyS0(*S0u_,*u);
      for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j) {
          for (unsigned k = 0; k < T_; ++k) {
            fem_->applyS(i,j,k,*Su_[i][j][k],*u);
          }
        }
      }
      buildState_ = false;
    }
  }

  void applySuz(const ROL::Ptr<Tpetra::MultiVector<>>       &Juz,
                const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                const ROL::Ptr<const std::vector<Real>>     &z) {
    const Real zero(0), one(1);
    buildSu(u);
    Juz->update(one,*S0u_,zero);
    for (unsigned i = 0; i < M_; ++i) {
      for (unsigned j = 0; j < N_; ++j) {
        for (unsigned k = 0; k < T_; ++k) {
          Juz->update((*z)[i+M_*(j+N_*k)],*Su_[i][j][k],one);
        }
      }
    }
  }

  void buildJacobian1(const ROL::Ptr<const std::vector<Real>> &z) {
    if (assembleJ1_) {
      const Real zero(0), one(1);
      ROL::Ptr<Teuchos::ParameterList> list = ROL::makePtr<Teuchos::ParameterList>();
      matJ1_ = matJ1_->add(one,*fem_->getS0(),zero,
        fem_->getS0()->getDomainMap(),fem_->getS0()->getRangeMap(),list);
      for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j) {
          for (unsigned k = 0; k < T_; ++k) {
            matJ1_ = matJ1_->add((*z)[i+M_*(j+N_*k)],*fem_->getS(i,j,k),one,
              fem_->getS(i,j,k)->getDomainMap(),fem_->getS(i,j,k)->getRangeMap(),list);
          }
        }
      }
      assembleJ1_ = false;
      setSolver_  = true;
    }
  }

  void applyJacobian1(const ROL::Ptr<Tpetra::MultiVector<>>       &Jv,
                      const ROL::Ptr<const Tpetra::MultiVector<>> &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BilinearPDEConstraintApplyJacobian1);
    #endif
    matJ1_->apply(*v,*Jv);
  }

  void applyAdjointJacobian1(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<>> &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintApplyAdjointJacobian1);
    #endif
    matJ1_->apply(*v,*Jv,Teuchos::TRANS);
  }

  void applyJacobian2(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                      const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                      const ROL::Ptr<const std::vector<Real>> &v) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintApplyJacobian2);
    #endif
    buildSu(u);
    Jv->scale(static_cast<Real>(0));
    const Real one(1);
    for (unsigned i = 0; i < M_; ++i) {
      for (unsigned j = 0; j < N_; ++j) {
        for (unsigned k = 0; k < T_; ++k) {
          Jv->update((*v)[i+M_*(j+N_*k)],*Su_[i][j][k],one);
        }
      }
    }
  }

  void applyAdjointJacobian2(const ROL::Ptr<std::vector<Real>> &Jv,
                             const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                             const ROL::Ptr<const Tpetra::MultiVector<>> &u) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintApplyAdjointJacobian2);
    #endif
    buildSu(u);
    Teuchos::Array<Real> val(1,0);
    for (unsigned i = 0; i < M_; ++i) {
      for (unsigned j = 0; j < N_; ++j) {
        for (unsigned k = 0; k < T_; ++k) {
          Su_[i][j][k]->dot(*v, val.view(0,1));
          (*Jv)[i+M_*(j+N_*k)] = val[0];
        }
      }
    }
  }

  void buildSl(const ROL::Ptr<const Tpetra::MultiVector<>> &l) {
    if (buildAdjoint_) {
      for (unsigned i = 0; i < M_; ++i) {
        for (unsigned j = 0; j < N_; ++j) {
          for (unsigned k = 0; k < T_; ++k) {
            fem_->applyS(i,j,k,*Sl_[i][j][k],*l,true);
          }
        }
      }
      buildAdjoint_ = false;
    }
  }

  void applyHessian(const ROL::Ptr<Tpetra::MultiVector<>> &H,
                    const ROL::Ptr<const std::vector<Real>> &v,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &l) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintApplyHessian);
    #endif
    buildSl(l);
    H->scale(static_cast<Real>(0));
    const Real one(1);
    for (unsigned i = 0; i < M_; ++i) {
      for (unsigned j = 0; j < N_; ++j) {
        for (unsigned k = 0; k < T_; ++k) {
          H->update((*v)[i+M_*(j+N_*k)],*Sl_[i][j][k],one);
        }
      }
    }
  }

  void applyAdjointHessian(const ROL::Ptr<std::vector<Real>> &H,
                           const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                           const ROL::Ptr<const Tpetra::MultiVector<>> &l) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::BinlinearPDEConstraintApplyAdjointHessian);
    #endif
    buildSl(l);
    Teuchos::Array<Real> val(1,0);
    for (unsigned i = 0; i < M_; ++i) {
      for (unsigned j = 0; j < N_; ++j) {
        for (unsigned k = 0; k < T_; ++k) {
          Sl_[i][j][k]->dot(*v, val.view(0,1));
          (*H)[i+M_*(j+N_*k)] = val[0];
        }
      }
    }
  }

public:
  Bilinear_PDE_Constraint(const ROL::Ptr<FEM_Data<Real>>           &fem,
                          Teuchos::ParameterList                   &list,
                          std::ostream                             &outStream = std::cout)
    : ROL::Constraint_SimOpt<Real>(),
      fem_(fem),
      setSolver_(true),
      assembleJ1_(true),
      buildState_(true),
      buildAdjoint_(true) {
    // Construct solver.
    solver_ = ROL::makePtr<Solver<Real>>(list.sublist("Solver"));
    // Initialize zero vectors.
    cvec_ = fem_->getAssembler()->createResidualVector();

    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(list.sublist("Problem"), "Young's Modulus");
    M_ = list.sublist("Problem").get("Number of Horizontal Cells",10);
    N_ = list.sublist("Problem").get("Number of Vertical Cells",20);
    T_ = ym.size();

    matJ1_ = ROL::makePtr<Tpetra::CrsMatrix<>>(fem_->getS0()->getCrsGraph());
    S0u_   = fem_->getAssembler()->createResidualVector();
    Su_.resize(M_);
    Sl_.resize(M_);
    for (unsigned i = 0; i < M_; ++i) {
      Su_[i].resize(N_);
      Sl_[i].resize(N_);
      for (unsigned j = 0; j < N_; ++j) {
        Su_[i][j].resize(T_);
        Sl_[i][j].resize(T_);
        for (unsigned k = 0; k < T_; ++k) {
          Su_[i][j][k] = fem_->getAssembler()->createResidualVector();
          Sl_[i][j][k] = fem_->getAssembler()->createResidualVector();
        }
      }
    }
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return fem_->getAssembler();
  }

  const ROL::Ptr<PDE<Real>> getPDE(void) const {
    return fem_->getPDE();
  }

  using ROL::Constraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {
    buildState_   = true;
    buildAdjoint_ = true;
  }

  using ROL::Constraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    assembleJ1_   = true;
  }

  using ROL::Constraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    //update_1(u,flag,iter);
    //update_2(z,flag,iter);
  }

  using ROL::Constraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<>>       cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<>> uf = getConstField(u);
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    const Real one(1);
    cvec_->putScalar(static_cast<Real>(0));
    c.zero();

    cf->scale(one,*fem_->getForce());
    applySuz(cvec_,uf,zp);
    cf->update(one,*cvec_,one);
  }

  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<>>       cf = getField(c);
    ROL::Ptr<Tpetra::MultiVector<>>       uf = getField(u);
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    c.zero();
    cvec_->putScalar(static_cast<Real>(0));
    u.zero();

    const Real one(1);
    cf->scale(one,*fem_->getForce());
    buildJacobian1(zp);
    solveForward(uf,cf);
    uf->scale(-one);

    applySuz(cvec_,uf,zp);
    cf->update(one,*cvec_,one);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    buildJacobian1(zp);
    applyJacobian1(jvf,vf);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>       jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>>  uf = getConstField(u);
    ROL::Ptr<const std::vector<Real>>      vp = getConstParameter(v);

    applyJacobian2(jvf,uf,vp);
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ajv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>     ajvf = getField(ajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    buildJacobian1(zp);
    applyAdjointJacobian1(ajvf,vf);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ajv.zero();
    ROL::Ptr<std::vector<Real>>          ajvp = getParameter(ajv);
    ROL::Ptr<const Tpetra::MultiVector<>>  uf = getConstField(u);
    ROL::Ptr<const Tpetra::MultiVector<>>  vf = getConstField(v);

    applyAdjointJacobian2(ajvp,vf,uf);
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ijv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>     ijvf = getField(ijv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    buildJacobian1(zp);
    solveForward(ijvf,vf);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,
                               const ROL::Vector<Real> &z, Real &tol) {
    iajv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>    iajvf = getField(iajv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    ROL::Ptr<const std::vector<Real>>     zp = getConstParameter(z);

    buildJacobian1(zp);
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
    ROL::Ptr<std::vector<Real>>        ahwvp = getParameter(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> wf = getConstField(w);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);

    applyAdjointHessian(ahwvp,vf,wf);
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
    ROL::Ptr<Tpetra::MultiVector<>>    ahwvf = getField(ahwv);
    ROL::Ptr<const Tpetra::MultiVector<>> wf = getConstField(w);
    ROL::Ptr<const std::vector<Real>>     vp = getConstParameter(v);

    applyHessian(ahwvf,vp,wf);
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
    fem_->printMeshData(outStream);
  }

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> matWriter;
    if (matJ1_ != ROL::nullPtr) {
      matWriter.writeSparseFile("jacobian1.txt", ROL::staticPtrCast<Tpetra::CrsMatrix<>>(matJ1_));
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("jacobian1.txt");
      emptyfile.close();
    }
    if (fem_->getForce() != ROL::nullPtr) {
      matWriter.writeDenseFile("residual.txt", fem_->getForce());
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("residual.txt");
      emptyfile.close();
    }
  }

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<>> &vec,
                          const std::string &filename) const {
    fem_->outputTpetraVector(vec,filename);
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
      try {
        ROL::Ptr<const ROL::TpetraMultiVector<Real>> xvec
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

  ROL::Ptr<Tpetra::MultiVector<>> getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::TpetraMultiVector<Real>> xvec
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

  ROL::Ptr<const std::vector<Real>> getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real>> xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
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
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<std::vector<Real>> getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real>> xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
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
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }
};

#endif
