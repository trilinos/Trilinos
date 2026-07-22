// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef FILTERED_COMPLIANCE_ROBJ_H
#define FILTERED_COMPLIANCE_ROBJ_H

#include "ROL_Objective.hpp"
#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/assembler.hpp"
#include "../../../TOOLS/solver.hpp"
#include "../../../TOOLS/pdevector.hpp"
#include "compliance_robj.hpp"
#include "filter.hpp"

template <class Real>
class Filtered_Compliance_Objective : public ROL::Objective<Real> {
private:
  // Compliance objective function
  ROL::Ptr<Compliance_Objective<Real>> comp_;

  // Filter PDE
  ROL::Ptr<Filter<Real>> filter_;

  // Vector Storage
  ROL::Ptr<ROL::Vector<Real>> Fz_, Fzcache_, Fv_, dctrl_;

  bool nuke_;
  int printFreq_;

  unsigned nupda_, nfval_, ngrad_, nhess_, nprec_;

public:
  Filtered_Compliance_Objective(
      const ROL::Ptr<PDE<Real>>                    &pde,
      const ROL::Ptr<PDE<Real>>                    &filter,
      const ROL::Ptr<MeshManager<Real>>            &mesh,
      const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
      ROL::ParameterList                           &list,
      std::ostream                                 &stream = std::cout)
    : Filtered_Compliance_Objective(ROL::makePtr<Filter<Real>>(pde,mesh,comm,list,stream),pde,mesh,comm,list,stream) {}

  Filtered_Compliance_Objective(
      const ROL::Ptr<Filter<Real>>                 &filter,
      const ROL::Ptr<PDE<Real>>                    &pde,
      const ROL::Ptr<MeshManager<Real>>            &mesh,
      const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
      ROL::ParameterList                           &list,
      std::ostream                                 &stream = std::cout)
    : filter_(filter), nupda_(0), nfval_(0), ngrad_(0), nhess_(0), nprec_(0) {
    // Compliance objective function
    comp_ = ROL::makePtr<Compliance_Objective<Real>>(pde,mesh,comm,list,stream,"state","filtered_density");

    // Vector storage
    Fz_      = filter_->createPrimalFilterVector();
    Fzcache_ = filter_->createPrimalFilterVector();
    Fv_      = filter_->createPrimalFilterVector();
    dctrl_   = filter_->createDualFilterVector();

    nuke_      = list.sublist("Problem").get("Use Basic Update",false);
    printFreq_ = list.sublist("Problem").get("Output Frequency",0);
  }

  Filtered_Compliance_Objective(
      const ROL::Ptr<Filter<Real>>                 &filter,
      const ROL::Ptr<PDE<Real>>                    &pde,
      const ROL::Ptr<Assembler<Real>>              &assembler,
      ROL::ParameterList                           &list)
    : filter_(filter), nupda_(0), nfval_(0), ngrad_(0), nhess_(0), nprec_(0) {
    // Compliance objective function
    comp_ = ROL::makePtr<Compliance_Objective<Real>>(pde,assembler,list,"state","filtered_density");

    // Vector storage
    Fz_      = filter_->createPrimalFilterVector();
    Fzcache_ = filter_->createPrimalFilterVector();
    Fv_      = filter_->createPrimalFilterVector();
    dctrl_   = filter_->createDualFilterVector();

    nuke_      = list.sublist("Problem").get("Use Basic Update",false);
    printFreq_ = list.sublist("Problem").get("Output Frequency",0);
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return comp_->getAssembler();
  }

  const ROL::Ptr<Assembler<Real>> getFilterAssembler(void) const {
    return filter_->getAssembler();
  }

  ROL::Ptr<ROL::Vector<Real>> applyFilter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<ROL::Vector<Real>> Fx = Fz_->clone();
    filter_->apply(*Fx,x,false);
    return Fx;
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
      filter_->getAssembler()->outputTpetraVector(d_data, dfile.str());
    }
  }

  Real normalize(const ROL::Vector<Real> &z, Real &tol) {
    filter_->apply(*Fz_,z,false);
    return comp_->normalize(*Fz_, tol);
  }

  void printToFile(const ROL::Vector<Real> &z, std::ostream &stream = std::cout,
                   std::string ufile = "state.txt", std::string dfile = "density.txt",
                   std::string ffile = "filtered_density.txt") {
    update(z,ROL::UpdateType::Temp);
    ROL::Ptr<const Tpetra::MultiVector<>> dens_data = getConstField(z);
    filter_->getAssembler()->outputTpetraVector(dens_data,  dfile);
    comp_->printToFile(*Fz_, stream, ufile, ffile);
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    nfval_++;
    return comp_->value(*Fz_,tol);
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    ngrad_++;
    comp_->gradient(*dctrl_, *Fz_, tol);
    filter_->apply(g, *dctrl_, true);
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nhess_++;
    filter_->apply(*Fv_, v, false);
    comp_->hessVec(*dctrl_, *Fv_, *Fz_, tol);
    filter_->apply(hv, *dctrl_, true);
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
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
    filter_->summarize(stream);
    comp_->summarize(stream);
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

  void update_initial( const ROL::Vector<Real> &z, int iter ) {
    filter_->apply(*Fz_,z,false);
    Fzcache_->set(*Fz_);
  }

  void update_accept( const ROL::Vector<Real> &z, int iter ) {
    Fzcache_->set(*Fz_);
  }

  void update_temp( const ROL::Vector<Real> &z, int iter ) {
    filter_->apply(*Fz_,z,false);
  }

  void update_trial( const ROL::Vector<Real> &z, int iter ) {
    filter_->apply(*Fz_,z,false);
  }

  void update_revert( const ROL::Vector<Real> &z, int iter ) {
    Fz_->set(*Fzcache_);
  }
}; // class Compliance_Objective

#endif
