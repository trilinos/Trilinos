// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BINARY_ELASTICITY_FILTER_H
#define BINARY_ELASTICITY_FILTER_H

#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/assembler.hpp"
#include "../../../TOOLS/solver.hpp"
#include "../../../TOOLS/pdevector.hpp"

template <class Real>
class Filter {
private:
  // FIlter PDE
  const ROL::Ptr<PDE<Real>> pde_filter_;
  ROL::Ptr<Assembler<Real>> assembler_filter_;
  ROL::Ptr<Solver<Real>> solver_filter_;
  ROL::Ptr<Tpetra::CrsMatrix<>> matF1_, matF2_;
  ROL::ParameterList list_;

  // Vector Storage
  ROL::Ptr<Tpetra::MultiVector<>> f_state_, f_res_;

  unsigned nasFi_, napFi_;

public:
  Filter(const ROL::Ptr<PDE<Real>>                    &pde,
         const ROL::Ptr<MeshManager<Real>>            &mesh,
         const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
         ROL::ParameterList                           &list,
         std::ostream                                 &stream = std::cout)
    : pde_filter_(pde), list_(list), nasFi_(0), napFi_(0) {
    // Filter PDE
    assembler_filter_ = ROL::makePtr<Assembler<Real>>(pde_filter_->getFields(),
                                                      pde_filter_->getFields2(),
                                                      mesh,comm,list_,stream);
    assembler_filter_->setCellNodes(*pde_filter_);
    solver_filter_ = ROL::makePtr<Solver<Real>>(list.sublist("Solver"));
    // Vector storage
    f_state_ = assembler_filter_->createStateVector();
    f_res_   = assembler_filter_->createResidualVector();

    assemble();
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_filter_;
  }

  void apply(ROL::Vector<Real> &Fx, const ROL::Vector<Real> &x, const bool transpose) {
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

  ROL::Ptr<ROL::Vector<Real>> createPrimalFilterVector(void) {
    return ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_filter_->createStateVector(),pde_filter_,assembler_filter_,list_);
  }

  ROL::Ptr<ROL::Vector<Real>> createDualFilterVector(void) {
    return ROL::makePtr<PDE_DualSimVector<Real>>(assembler_filter_->createStateVector(),pde_filter_,assembler_filter_,list_);
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  Filter::summarize" << std::endl;
    stream << "    Number of filter assemblies:        " << nasFi_ << std::endl;
    stream << "    Number of filter applies:           " << napFi_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
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

  void assemble(void) {
    ROL::Ptr<Tpetra::MultiVector<>> f_ctrl = assembler_filter_->createControlVector();
    assembler_filter_->assemblePDEJacobian1(matF1_,pde_filter_,f_state_,f_ctrl);
    assembler_filter_->assemblePDEJacobian2(matF2_,pde_filter_,f_state_,f_ctrl);
    solver_filter_->setA(matF1_);
  }
}; // class Filter

#endif
