// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TOPOPT_VOLUME_OBJ_H
#define TOPOPT_VOLUME_OBJ_H

#include "ROL_Objective.hpp"
#include "../../TOOLS/qoi.hpp"
#include "../../TOOLS/assembler.hpp"
#include "../../TOOLS/pdevector.hpp"

template <class Real>
class TopOptVolumeObjective : public ROL::Objective<Real> {
private:
  const ROL::Ptr<QoI<Real>> qoi_;
  const ROL::Ptr<Assembler<Real>> assembler_;
  bool assemble_;
  ROL::Ptr<Tpetra::MultiVector<>> controlZero_;
  ROL::Ptr<ROL::Vector<Real>> V_, Vdual_;

  unsigned nupda_;
  unsigned nfval_;
  unsigned ngrad_;
  unsigned nhess_;
  unsigned nasse_;

public:
  TopOptVolumeObjective(
      const ROL::Ptr<QoI<Real>>                    &qoi,
      const ROL::Ptr<Assembler<Real>>              &assembler,
      const ROL::Ptr<const ROL::Vector<Real>>      &z)
    : qoi_(qoi), assembler_(assembler), assemble_(true),
      nupda_(0), nfval_(0), ngrad_(0), nhess_(0), nasse_(0) {
    controlZero_ = assembler_->createControlVector(); controlZero_->putScalar(0.0);
    V_     = z->dual().clone();
    Vdual_ = z->clone();
    assemble();
  }

  void update( const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    nupda_++;
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    nfval_++;
    return Vdual_->dot(z);
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    ngrad_++;
    g.set(*V_);
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nhess_++;
    hv.zero();
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  Volume_Objective::summarize" << std::endl;
    stream << "    Number of calls to update:   " << nupda_ << std::endl;
    stream << "    Number of calls to value:    " << nfval_ << std::endl;
    stream << "    Number of calls to gradient: " << ngrad_ << std::endl;
    stream << "    Number of calls to hessVec:  " << nhess_ << std::endl;
    stream << "    Number of assemblies:        " << nasse_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
  }

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    qoi_->setParameter(param);
  }

private:

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
    nasse_++;
    if (assemble_) {
      ROL::Ptr<Tpetra::MultiVector<>> V = getField(*V_);
      assembler_->assembleQoIGradient2(V,qoi_,ROL::nullPtr,controlZero_);
      Vdual_->set(V_->dual());
      assemble_ = false;
    }
  }
}; // class TopOptVolumeObjective

#endif
