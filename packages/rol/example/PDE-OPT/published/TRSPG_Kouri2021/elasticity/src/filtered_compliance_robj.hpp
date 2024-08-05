// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TOPOPT_FILTERED_COMPLIANCE_ROBJ_H
#define TOPOPT_FILTERED_COMPLIANCE_ROBJ_H

#include "ROL_Objective.hpp"
#include "../../../../TOOLS/pde.hpp"
#include "../../../../TOOLS/assembler.hpp"
#include "../../../../TOOLS/solver.hpp"
#include "../../../../TOOLS/pdevector.hpp"
#include "compliance_robj.hpp"

template <class Real>
class TopOptFilteredComplianceObjective : public ROL::Objective<Real> {
private:
  // Compliance objective function
  ROL::Ptr<TopOptComplianceObjective<Real>> comp_;

  // FIlter PDE
  const ROL::Ptr<PDE<Real>> pde_filter_;
  ROL::Ptr<Assembler<Real>> assembler_filter_;
  ROL::Ptr<Solver<Real>> solver_filter_;
  ROL::Ptr<Tpetra::CrsMatrix<>> matF1_, matF2_;

  // Vector Storage
  ROL::Ptr<Tpetra::MultiVector<>> f_state_, f_res_;
  ROL::Ptr<ROL::Vector<Real>> Fz_, Fzcache_, Fv_, dctrl_;

  bool nuke_;
  int printFreq_;

  unsigned nupda_, nfval_, ngrad_, nhess_, nprec_;
  unsigned nasFi_, napFi_;

public:
  TopOptFilteredComplianceObjective(
      const ROL::Ptr<PDE<Real>>                    &pde,
      const ROL::Ptr<PDE<Real>>                    &filter,
      const ROL::Ptr<MeshManager<Real>>            &mesh,
      const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
      ROL::ParameterList                           &list,
      std::ostream                                 &stream = std::cout)
    : pde_filter_(filter),
      nupda_(0), nfval_(0), ngrad_(0), nhess_(0), nprec_(0),
      nasFi_(0), napFi_(0) {
    // Compliance objective function
    comp_ = ROL::makePtr<TopOptComplianceObjective<Real>>(pde,mesh,comm,list,stream,"state","filtered_density");

    // Filter PDE
    assembler_filter_ = ROL::makePtr<Assembler<Real>>(pde_filter_->getFields(),
                                                      pde_filter_->getFields2(),
                                                      mesh,comm,list,stream);
    assembler_filter_->setCellNodes(*pde_filter_);
    ROL::ParameterList solverList = list.sublist("Solver");
    solverList.set("Preconditioner","MueLu");
    solverList.sublist("MueLu").set("number of equations",1);
    solverList.sublist("MueLu").set("problem: symmetric",true);
    solverList.sublist("Belos").set("Convergence Tolerance",1e-10);
    solver_filter_    = ROL::makePtr<Solver<Real>>(solverList);

    // Vector storage
    f_state_ = assembler_filter_->createStateVector();
    f_res_   = assembler_filter_->createResidualVector();
    Fz_      = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_filter_->createStateVector(),pde_filter_,assembler_filter_,list);
    Fzcache_ = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_filter_->createStateVector(),pde_filter_,assembler_filter_,list);
    Fv_      = ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_filter_->createStateVector(),pde_filter_,assembler_filter_,list);
    dctrl_   = ROL::makePtr<PDE_DualSimVector<Real>>(assembler_filter_->createStateVector(),pde_filter_,assembler_filter_,list);

    nuke_      = list.sublist("Problem").get("Use Basic Update",false);
    printFreq_ = list.sublist("Problem").get("Output Frequency",0);

    assembleFilter();
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return comp_->getAssembler();
  }

  const ROL::Ptr<Assembler<Real>> getFilterAssembler(void) const {
    return assembler_filter_;
  }

  void update(const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    nupda_++;
    if (nuke_) {
      update_temp(z,iter);
    }
    else {
      switch (type) {
        case ROL::UpdateType::Initial: update_initial(z,iter); break;
        case ROL::UpdateType::Accept:  update_accept(z,iter);  break;
        case ROL::UpdateType::Revert:  update_revert(z,iter);  break;
        case ROL::UpdateType::Trial:   update_trial(z,iter);   break;
        case ROL::UpdateType::Temp:    update_temp(z,iter);   break;
      }
    }
    comp_->update(*Fz_,type,iter);
    // Print
    if (printFreq_ > 0 && iter%printFreq_ == 0) {
      std::stringstream dfile;
      dfile  << "density_"          << iter << ".txt";
      ROL::Ptr<const Tpetra::MultiVector<>> d_data = getConstField(z);
      assembler_filter_->outputTpetraVector(d_data, dfile.str());
    }
  }

  void update( const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    nupda_++;
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
    comp_->update(*Fz_,flag,iter);
  }

  Real normalize(const ROL::Vector<Real> &z, Real &tol) {
    applyFilter(*Fz_,z,false);
    return comp_->normalize(*Fz_, tol);
  }

  void printToFile(const ROL::Vector<Real> &z, std::ostream &stream = std::cout,
                   std::string ufile = "state.txt", std::string dfile = "density.txt",
                   std::string ffile = "filtered_density.txt") {
    update(z,ROL::UpdateType::Temp);
    ROL::Ptr<const Tpetra::MultiVector<>> dens_data = getConstField(z);
    assembler_filter_->outputTpetraVector(dens_data,  dfile);
    comp_->printToFile(*Fz_, stream, ufile, ffile);
  }

  void solveState(ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    update(z,ROL::UpdateType::Temp);
    comp_->solveState(u,*Fz_);
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    nfval_++;
    return comp_->value(*Fz_,tol);
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    ngrad_++;
    comp_->gradient(*dctrl_, *Fz_, tol);
    applyFilter(g, *dctrl_, true);
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nhess_++;
    applyFilter(*Fv_, v, false);
    comp_->hessVec(*dctrl_, *Fv_, *Fz_, tol);
    applyFilter(hv, *dctrl_, true);
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
    stream << "  Filtered_Objective_Compliance::summarize" << std::endl;
    stream << "    Number of calls to update:          " << nupda_ << std::endl;
    stream << "    Number of calls to value:           " << nfval_ << std::endl;
    stream << "    Number of calls to gradient:        " << ngrad_ << std::endl;
    stream << "    Number of calls to hessVec:         " << nhess_ << std::endl;
    stream << "    Number of calls to precond:         " << nprec_ << std::endl;
    stream << "    Number of filter assemblies:        " << nasFi_ << std::endl;
    stream << "    Number of filter applies:           " << napFi_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
    comp_->summarize(stream);
  }

  void reset(void) {
    nupda_ = 0; nfval_ = 0; ngrad_ = 0; nhess_ = 0; nprec_ = 0; nasFi_ = 0; napFi_ = 0;
    comp_->reset();
  }

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    comp_->setParameter(param);
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

  void assembleFilter(void) {
    nasFi_++;
    ROL::Ptr<Tpetra::MultiVector<>> f_ctrl = assembler_filter_->createControlVector();
    assembler_filter_->assemblePDEJacobian1(matF1_,pde_filter_,f_state_,f_ctrl);
    assembler_filter_->assemblePDEJacobian2(matF2_,pde_filter_,f_state_,f_ctrl);
    solver_filter_->setA(matF1_);
  }

  void applyFilter(ROL::Vector<Real> &Fx, const ROL::Vector<Real> &x, const bool transpose) {
    napFi_++;
    ROL::Ptr<Tpetra::MultiVector<>>      Fx_data = getField(Fx);
    ROL::Ptr<const Tpetra::MultiVector<>> x_data = getConstField(x);
    if (transpose) {
      solver_filter_->solve(f_state_,x_data,false);
      matF2_->apply(*f_state_,*Fx_data,Teuchos::TRANS);
    }
    else {
      matF2_->apply(*x_data,*f_res_);
      solver_filter_->solve(Fx_data,f_res_,false);
    }
    Fx.scale(static_cast<Real>(-1));
  }

  void update_initial( const ROL::Vector<Real> &z, int iter ) {
    applyFilter(*Fz_,z,false);
    Fzcache_->set(*Fz_);
  }

  void update_accept( const ROL::Vector<Real> &z, int iter ) {
    Fzcache_->set(*Fz_);
  }

  void update_temp( const ROL::Vector<Real> &z, int iter ) {
    applyFilter(*Fz_,z,false);
  }

  void update_trial( const ROL::Vector<Real> &z, int iter ) {
    applyFilter(*Fz_,z,false);
  }

  void update_revert( const ROL::Vector<Real> &z, int iter ) {
    Fz_->set(*Fzcache_);
  }
}; // class TopOptComplianceObjective

#endif
