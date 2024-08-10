// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TOPOT_COMPLIANCE_ROBJ_H
#define TOPOT_COMPLIANCE_ROBJ_H

#include "ROL_Objective.hpp"
#include "ROL_VectorController.hpp"
#include "ROL_SecantFactory.hpp"
#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/assembler.hpp"
#include "../../../../TOOLS/solver.hpp"
#include "../../../../TOOLS/pdevector.hpp"

template <class Real>
class TopOptComplianceObjective : public ROL::Objective<Real> {
private:
  // Elasticity PDE
  const ROL::Ptr<PDE<Real>> pde_;
  ROL::Ptr<Assembler<Real>> assembler_;
  ROL::Ptr<Solver<Real>> solver_, solver_cache_, solver_tmp_;
  ROL::Ptr<Tpetra::CrsMatrix<>> matJ1_, matJ1_cache_, matJ1_tmp_;
  ROL::Ptr<Tpetra::CrsMatrix<>> matJ2_, matJ2_cache_, matJ2_tmp_;
  ROL::Ptr<Tpetra::CrsMatrix<>> matH22_;
  bool assembleJ1_, assembleJ2_, assembleH22_, assembleF_;
  ROL::ParameterList solverList_;
  Real tol_, maxtol_;
  bool useInexact_;

  // Output information
  const std::string uname_, dname_;

  // Vector Storage
  ROL::Ptr<Tpetra::MultiVector<>> F_data_, state_sens_data_;
  ROL::Ptr<Tpetra::MultiVector<>> dstat_data_, dctrl_data_;
  ROL::Ptr<ROL::Vector<Real>> state_, dctrl_;
  ROL::Ptr<ROL::VectorController<Real>> stateStore_;

  bool storage_;
  Real cmpScale_;
  bool nuke_;
  int printFreq_;

  unsigned nstat_, nsens_;
  unsigned nupda_, nfval_, ngrad_, nhess_, nprec_;
  unsigned napJ1_, nasJ1_, napJ2_, nasJ2_, napH2_, nasH2_, nasLo_;

  // Hessian Approximation Data
  std::string hessAppr_;
  Real H2scale_;
  ROL::Ptr<ROL::Secant<Real>> secant_;
  ROL::Ptr<Tpetra::MultiVector<>> s_data_;
  ROL::Ptr<ROL::Vector<Real>> g_new_, g_old_, s_, g_;
  bool computeGrad_;

public:
  TopOptComplianceObjective(
      const ROL::Ptr<PDE<Real>>                    &pde,
      const ROL::Ptr<MeshManager<Real>>            &mesh,
      const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
      ROL::ParameterList                           &list,
      std::ostream                                 &stream = std::cout,
      std::string                                   uname  = "state",
      std::string                                   dname  = "density")
    : pde_(pde), assembleJ1_(true), assembleJ2_(true),
      assembleH22_(true), assembleF_(true), tol_(1e-2), maxtol_(1e-2),
      uname_(uname), dname_(dname),
      nstat_(0), nsens_(0),
      nupda_(0), nfval_(0), ngrad_(0), nhess_(0), nprec_(0),
      napJ1_(0), nasJ1_(0), napJ2_(0), nasJ2_(0), napH2_(0), nasH2_(0), nasLo_(0), H2scale_(1) {
    // Elasticity PDE
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),
                                               pde_->getFields2(),
                                               mesh,comm,list,stream);
    assembler_->setCellNodes(*pde_);
    solver_       = ROL::makePtr<Solver<Real>>(list.sublist("Solver"));
    solver_cache_ = ROL::makePtr<Solver<Real>>(list.sublist("Solver"));
    solverList_   = list.sublist("Solver");
    useInexact_   = list.sublist("Problem").get("Use Inexact Linear Solves",false);

    // Vector storage
    F_data_           = assembler_->createResidualVector();
    state_sens_data_  = assembler_->createStateVector();
    dstat_data_       = assembler_->createResidualVector();
    dctrl_data_       = assembler_->createControlVector();
    state_            = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_->createStateVector(),pde_,assembler_,list);
    dctrl_            = ROL::makePtr<PDE_DualOptVector<Real>>(dctrl_data_,pde_,assembler_,list);
    stateStore_       = ROL::makePtr<ROL::VectorController<Real>>();

    storage_   = list.sublist("Problem").get("Use state storage", true);
    cmpScale_  = list.sublist("Problem").get("Compliance Scaling", 1.0);
    hessAppr_  = list.sublist("Problem").get("Hessian Approximation", "None");
    nuke_      = list.sublist("Problem").get("Use Basic Update",false);
    printFreq_ = list.sublist("Problem").get("Output Frequency",0);

    computeGrad_ = true;
    if (hessAppr_ == "Secant") {
      s_data_ = assembler_->createControlVector();
      s_      = ROL::makePtr<PDE_PrimalOptVector<Real>>(s_data_,pde_,assembler_,list);
      g_      = ROL::makePtr<PDE_DualOptVector<Real>>(assembler_->createControlVector(),pde_,assembler_,list);
      g_old_  = ROL::makePtr<PDE_DualOptVector<Real>>(assembler_->createControlVector(),pde_,assembler_,list);
      g_new_  = ROL::makePtr<PDE_DualOptVector<Real>>(assembler_->createControlVector(),pde_,assembler_,list);
      secant_ = ROL::SecantFactory<Real>(list);
    }

    initialize();
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  void update(const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    nupda_++;
    stateStore_->objectiveUpdate(type);
    if (nuke_) {
      update_temp(z,iter);
    }
    else {
      switch (type) {
        case ROL::UpdateType::Initial: update_initial(z,iter); break;
        case ROL::UpdateType::Accept:  update_accept(z,iter);  break;
        case ROL::UpdateType::Revert:  update_revert(z,iter);  break;
        case ROL::UpdateType::Trial:   update_trial(z,iter);   break;
        case ROL::UpdateType::Temp:    update_temp(z,iter);    break;
      }
    }
    // Print
    if (printFreq_ > 0 && iter%printFreq_ == 0) {
      std::stringstream ufile, zfile;
      ufile  << uname_ << "_" << iter << ".txt";
      zfile  << dname_ << "_" << iter << ".txt";
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<const Tpetra::MultiVector<>> u_data = getConstField(*state_);
      assembler_->outputTpetraVector(u_data, ufile.str());
      assembler_->outputTpetraVector(z_data, zfile.str());
    }
  }

  void update( const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    nupda_++;
    stateStore_->objectiveUpdate(true);
    if (nuke_) {
      update_temp(z,iter);
    }
    else {
      // Trial:    flag = false, iter = -1
      // Check:    flag = true,  iter = -1
      // Reject:   flag = false, iter > -1
      // Accept:   flag = true,  iter > -1
      if (flag) {
        if (iter > -1) {
          update_accept(z,iter);
        }
        else {
          update_temp(z,iter);
        }
      }
      else {
        if (iter > -1) {
          update_revert(z,iter);
        }
        else {
          update_trial(z,iter);
        }
      }
    }
  }

  Real normalize(const ROL::Vector<Real> &z, Real &tol) {
    update(z,ROL::UpdateType::Temp);
    Real val = value(z,tol);
    cmpScale_ /= val;
    return cmpScale_;
  }

  void printToFile(const ROL::Vector<Real> &z, std::ostream &stream = std::cout,
                   const std::string ufile = "state.txt", const std::string dfile = "density.txt") {
    Real tol(1e-2*std::sqrt(ROL::ROL_EPSILON<Real>()));
    update(z,ROL::UpdateType::Temp);
    ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
    ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
    solve_state_equation(u_data, z_data, tol);
    assembler_->outputTpetraVector(u_data, ufile);
    assembler_->outputTpetraVector(z_data, dfile);
    assembler_->printMeshData(stream);
  }

  void solveState(ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    Real tol(1e-2*std::sqrt(ROL::ROL_EPSILON<Real>()));
    update(z,ROL::UpdateType::Temp);
    ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
    ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(u);
    solve_state_equation(u_data, z_data, tol);
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    nfval_++;
    ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
    ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
    // Solve state equation
    solve_state_equation(*state_,z,tol);
    // Compute compliance
    Teuchos::Array<Real> cmp(1,0);
    F_data_->dot(*u_data, cmp.view(0,1));
    return cmpScale_*cmp[0];
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    ngrad_++;
    if (hessAppr_ == "Secant") {
      // Will not work for stochastic problems
      if (computeGrad_) {
        ROL::Ptr<Tpetra::MultiVector<>>       g_data = getField(g);
        ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
        ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
        // Solve state equation
        solve_state_equation(*state_,z,tol);
        // Apply adjoint density jacobian
        assembleJ2(u_data, z_data);
        applyJacobian2(g_data, u_data, true);
        g_data->scale(-cmpScale_);

        g_->set(g);
        computeGrad_ = false;
      }
      else {
        g.set(*g_);
      }
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<>>       g_data = getField(g);
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
      // Solve state equation
      solve_state_equation(*state_,z,tol);
      // Apply adjoint density jacobian
      assembleJ2(u_data, z_data);
      applyJacobian2(g_data, u_data, true);
      g_data->scale(-cmpScale_);
    }
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nhess_++;
    if (hessAppr_ == "Zero") {
      hv.zero();
    }
    else if (hessAppr_ == "Positive Riesz") {
      hv.set(v.dual());
    }
    else if (hessAppr_ == "Negative Riesz") {
      hv.set(v.dual());
      hv.scale(static_cast<Real>(-1));
    }
    else if (hessAppr_ == "Negative") {
      ROL::Ptr<Tpetra::MultiVector<>>       h_data = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<>> v_data = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
      // Solve state equation
      solve_state_equation(*state_,z,tol_);
      // Apply density hessian to direction
      assembleH22(u_data, z_data);
      applyHessian22(h_data, v_data);
      h_data->scale(-cmpScale_);
    }
    else if (hessAppr_ == "Positive") {
      ROL::Ptr<Tpetra::MultiVector<>>       h_data = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<>> v_data = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
      // Solve state equation
      solve_state_equation(*state_,z,tol_);
      // Solve state sensitivity equation
      assembleJ2(u_data, z_data);
      solve_state_sensitivity(state_sens_data_, v_data, tol_);
      // Apply adjoint density jacobian to state sensitivity
      applyJacobian2(h_data, state_sens_data_, true);
      h_data->scale(static_cast<Real>(2)*cmpScale_);
    }
    else if (hessAppr_ == "No Solve") {
      ROL::Ptr<Tpetra::MultiVector<>>       h_data = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<>> v_data = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
      // Solve state equation
      solve_state_equation(*state_,z,tol_);
      // Apply density hessian to direction
      assembleH22(u_data, z_data);
      applyHessian22(h_data, v_data);
      // Approximate positive term
      assembleJ2(u_data, z_data);
      applyJacobian2(dstat_data_, v_data, false);
      applyJacobian2(dctrl_data_, dstat_data_, true);
      // Approximation of inv(K) by aI where a solves min||aF-U||
      // or b=1/a solves min||F-bU||
      Teuchos::Array<Real> UdotF(1,0), UdotU(1,0);
      u_data->dot(*F_data_,UdotF.view(0,1));
      F_data_->dot(*F_data_,UdotU.view(0,1));
      H2scale_ = UdotF[0]/UdotU[0];
      //u_data->dot(*u_data,UdotU.view(0,1));
      //H2scale_ = UdotU[0]/UdotF[0];
      h_data->update(static_cast<Real>(2)*cmpScale_*H2scale_,*dctrl_data_,-cmpScale_);
    }
    else if (hessAppr_ == "Secant") {
      // Will not work for stochastic problems
      ROL::Ptr<Tpetra::MultiVector<>>       h_data = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<>> v_data = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
      // Solve state equation
      solve_state_equation(*state_,z,tol_);
      // Apply density hessian to direction
      assembleH22(u_data, z_data);
      applyHessian22(h_data, v_data);
      // Apply secant approximation of positive term
      secant_->applyB(*dctrl_, v);
      h_data->update(static_cast<Real>(2)*cmpScale_,*dctrl_data_,-cmpScale_);
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<>>       h_data = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<>> v_data = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(*state_);
      // Solve state equation
      solve_state_equation(*state_,z,tol_);
      // Solve state sensitivity equation
      assembleJ2(u_data, z_data);
      solve_state_sensitivity(state_sens_data_, v_data, tol_);
      // Apply density hessian to direction
      assembleH22(u_data, z_data);
      applyHessian22(h_data, v_data);
      // Apply adjoint density jacobian to state sensitivity
      applyJacobian2(dctrl_data_, state_sens_data_, true);
      h_data->update(static_cast<Real>(2)*cmpScale_,*dctrl_data_,-cmpScale_);
    }
  }

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( ROL::Vector<Real> &Pv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nprec_++;
    Pv.set(v.dual());
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  Objective_Compliance::summarize" << std::endl;
    stream << "    Number of calls to update:          " << nupda_ << std::endl;
    stream << "    Number of calls to value:           " << nfval_ << std::endl;
    stream << "    Number of calls to gradient:        " << ngrad_ << std::endl;
    stream << "    Number of calls to hessVec:         " << nhess_ << std::endl;
    stream << "    Number of calls to precond:         " << nprec_ << std::endl;
    stream << "    Number of state solves:             " << nstat_ << std::endl;
    stream << "    Number of state sensitivity solves: " << nsens_ << std::endl;
    stream << "    Number of load assemblies:          " << nasLo_ << std::endl;
    stream << "    Number of Jacobian_1 assemblies:    " << nasJ1_ << std::endl;
    stream << "    Number of Jacobian_1 sovles:        " << napJ1_ << std::endl;
    stream << "    Number of Jacobian_2 assemblies:    " << nasJ2_ << std::endl;
    stream << "    Number of Jacobian_2 applies:       " << napJ2_ << std::endl;
    stream << "    Number of Hessian_22 assemblies:    " << nasH2_ << std::endl;
    stream << "    Number of Hessian_22 applies:       " << napH2_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
  }

  void reset(void) {
    nupda_ = 0; nfval_ = 0; ngrad_ = 0; nhess_ = 0; nprec_ = 0; nstat_ = 0; nsens_ = 0;
    nasLo_ = 0; nasJ1_ = 0; napJ1_ = 0; nasJ2_ = 0; napJ2_ = 0; nasH2_ = 0; napH2_ = 0;
  }

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    pde_->setParameter(param);
    assembleF_ = true;
  }

private:

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

  void initialize(void) {
    ROL::Ptr<Tpetra::MultiVector<>> u_data = getField(*state_);
    u_data->putScalar(0.0);
    dctrl_data_->putScalar(1.0);
    assembleJ1_     = true;
    assembleJ2_     = true;
    assembleJ1(u_data, dctrl_data_,false);
    assembleJ2(u_data, dctrl_data_);
    matJ1_tmp_      = matJ1_;
    matJ2_tmp_      = matJ2_;
    matJ1_          = matJ1_cache_;
    matJ2_          = matJ2_cache_;
    assembleJ1_     = true;
    assembleJ2_     = true;
    assembleJ1(u_data, dctrl_data_,false);
    assembleJ2(u_data, dctrl_data_);
    matJ1_cache_    = matJ1_;
    matJ2_cache_    = matJ2_;
    solver_tmp_     = solver_;
    solver_         = solver_cache_;
    assembleJ1_     = true;
    assembleJ2_     = true;
    // Upon exit, matJ1, matJ2 and solver are set to their cached versions
  }

  // Assemble load vector
  void assembleF(void) {
    if (assembleF_) {
      nasLo_++;
      ROL::Ptr<Tpetra::MultiVector<>> u0 = assembler_->createStateVector();
      ROL::Ptr<Tpetra::MultiVector<>> z1 = assembler_->createControlVector();
      u0->putScalar(0.0);
      z1->putScalar(1.0);
      assembler_->assemblePDEResidual(F_data_,pde_,u0,z1);
      F_data_->scale(static_cast<Real>(-1));
      assembleF_ = false;
    }
  }

  void assembleJ2(const ROL::Ptr<const Tpetra::MultiVector<>> &u, const ROL::Ptr<const Tpetra::MultiVector<>> &z) {
    if (assembleJ2_) {
      nasJ2_++;
      assembler_->assemblePDEJacobian2(matJ2_,pde_,u,z);
      assembleJ2_ = false;
    }
  }

  void applyJacobian2(ROL::Ptr<Tpetra::MultiVector<>> &Jv, const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                      const bool transpose) {
    napJ2_++;
    if (transpose) {
      matJ2_->apply(*v,*Jv,Teuchos::TRANS);
    }
    else {
      matJ2_->apply(*v,*Jv);
    }
  }

  void assembleH22(const ROL::Ptr<const Tpetra::MultiVector<>> &u, const ROL::Ptr<const Tpetra::MultiVector<>> &z) {
    if (assembleH22_) {
      nasH2_++;
      assembler_->assemblePDEHessian22(matH22_,pde_,u,u,z);
      assembleH22_ = false;
    }
  }

  void applyHessian22(ROL::Ptr<Tpetra::MultiVector<>> &Hv, const ROL::Ptr<const Tpetra::MultiVector<>> &v) {
    matH22_->apply(*v,*Hv);
    napH2_++;
  }

  void assembleJ1(const ROL::Ptr<const Tpetra::MultiVector<>> &u, const ROL::Ptr<const Tpetra::MultiVector<>> &z, bool setSolve = true) {
    if (assembleJ1_) {
      nasJ1_++;
      assembler_->assemblePDEJacobian1(matJ1_,pde_,u,z);
      if (setSolve) {
        solver_->setA(matJ1_);
      }
      assembleJ1_ = false;
    }
  }

  void solve(ROL::Ptr<Tpetra::MultiVector<>> &x, const ROL::Ptr<const Tpetra::MultiVector<>> &b, Real &tol) {
    napJ1_++;
    if (useInexact_) {
      solverList_.sublist("Belos").set("Convergence Tolerance",tol);
      solver_->setParameters(solverList_);
    }
    solver_->solve(x,b,false);
  }

  void solve_state_equation(ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    bool isComputed = false;
    Real mytol = std::min(maxtol_,tol);
    if (storage_ && tol_ <= mytol) {
      isComputed = stateStore_->get(u, ROL::Objective<Real>::getParameter());
    }
    if (!isComputed || !storage_ || tol_ > mytol) {
      ROL::Ptr<const Tpetra::MultiVector<>> z_data = getConstField(z);
      ROL::Ptr<Tpetra::MultiVector<>>       u_data = getField(u);
      solve_state_equation(u_data,z_data,mytol);
      tol_ = mytol;
      if (storage_) {
        stateStore_->set(u, ROL::Objective<Real>::getParameter());
      }
    }
  }

  void solve_state_equation(ROL::Ptr<Tpetra::MultiVector<>> &u, const ROL::Ptr<const Tpetra::MultiVector<>> &z, Real &tol) {
    nstat_++;
    assembleF(); assembleJ1(u,z,true);
    solve(u,F_data_,tol);
  }

  void solve_state_sensitivity(ROL::Ptr<Tpetra::MultiVector<>> &s, const ROL::Ptr<const Tpetra::MultiVector<>> &v, Real &tol) {
    // This assumes that assembleJ1 and assemble J2 were called
    nsens_++;
    applyJacobian2(dstat_data_, v, false);
    solve(s, dstat_data_, tol);
  }

  void update_initial( const ROL::Vector<Real> &z, int iter ) {
    assembleJ1_   = true;
    assembleJ2_   = true;
    assembleH22_  = true;

    // For secant approximation
    if (hessAppr_ == "Secant") {
      computeGrad_ = true;
      const Real half(0.5), one(1);
      Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      gradient(*g_old_,z,tol);       // g_old = J'(z)
      g_old_->scale(half/cmpScale_); // g_old = J'(z)/2c
      s_->set(z); s_->scale(-one);   // s     = -z
    }
  }

  void update_accept( const ROL::Vector<Real> &z, int iter ) {
    // Push trial state storage into base state storage.
    solver_tmp_       = solver_cache_;
    solver_cache_     = solver_;
    matJ1_tmp_        = matJ1_cache_;
    matJ1_cache_      = matJ1_;
    matJ2_tmp_        = matJ2_cache_;
    matJ2_cache_      = matJ2_;
    assembleJ1_       = false;
    // If assembleJ2 was not called, then assembleJ2_ = true
    assembleH22_      = true;

    // For secant update
    // Secant Equation: B s = 1/2c(J'(z) - (J'(z_old) + H s))
    if (hessAppr_ == "Secant") {
      computeGrad_ = true;
      const Real half(0.5), one(1);
      Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      s_->plus(z);                         // s     = z - z_old
      applyHessian22(dctrl_data_,s_data_); // dctrl = -H s
      g_old_->axpy(-half,*dctrl_);         // g_old = (J'(z_old) + cH s)/2c
      gradient(*g_new_,z,tol);             // gnew  = J'(z)
      g_new_->scale(half/cmpScale_);       // gnew  = J'(z)/2c
      secant_->updateStorage(z,*g_new_,*g_old_,*s_,s_->norm(),iter);
      s_->set(z); s_->scale(-one);         // s     = -z
      g_old_->set(*g_new_);                // g_old = J'(z)/2c
    }
  }

  void update_temp( const ROL::Vector<Real> &z, int iter ) {
    assembleJ1_     = true;
    assembleJ2_     = true;
    assembleH22_    = true;
    if (hessAppr_ == "Secant") {
      computeGrad_    = true;
    }
  }

  void update_trial( const ROL::Vector<Real> &z, int iter ) {
    solver_       = solver_tmp_;
    matJ1_        = matJ1_tmp_;
    matJ2_        = matJ2_tmp_;
    assembleJ1_   = true;
    assembleJ2_   = true;
    assembleH22_  = false;
  }

  void update_revert( const ROL::Vector<Real> &z, int iter ) {
    solver_       = solver_cache_;
    matJ1_        = matJ1_cache_;
    matJ2_        = matJ2_cache_;
    assembleJ1_   = false;
    // If assembleJ2 was not called, then assembleJ2_ = true
    assembleH22_  = false;
    // If gradient was not called, then computeGrad_ = true
  }
}; // class TopOptComplianceObjective

#endif
