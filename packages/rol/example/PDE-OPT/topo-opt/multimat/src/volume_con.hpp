// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef VOLUME_CON_MULTIMAT_H
#define VOLUME_CON_MULTIMAT_H

#include "ROL_Constraint.hpp"
#include "ROL_SingletonVector.hpp"
#include "../../../TOOLS/qoi.hpp"
#include "../../../TOOLS/assembler.hpp"
#include "../../../TOOLS/pdevector.hpp"

template <class Real>
class MultiMat_Volume_Constraint : public ROL::Constraint<Real> {
private:
  const ROL::Ptr<QoI<Real>> qoi_;
  const ROL::Ptr<Assembler<Real>> assembler_;
  bool assemble_;
  ROL::Ptr<ROL::Vector<Real>> V_, Vdual_;
  Real V0_;

  unsigned nupda_;
  unsigned nfval_;
  unsigned njaco_;
  unsigned najac_;
  unsigned nhess_;
  unsigned nasse_;

public:
  MultiMat_Volume_Constraint(
      const ROL::Ptr<QoI<Real>>                    &qoi,
      const ROL::Ptr<Assembler<Real>>              &assembler,
      const ROL::Ptr<const ROL::Vector<Real>>      &z)
    : qoi_(qoi), assembler_(assembler), assemble_(true),
      nupda_(0), nfval_(0), njaco_(0), najac_(0), nhess_(0), nasse_(0) {
    V_     = z->dual().clone();
    Vdual_ = z->clone();
    assemble();
  }

  void update( const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    nupda_++;
  }

  void value( ROL::Vector<Real> &c, const ROL::Vector<Real> &z, Real &tol ) {
    nfval_++;
    dynamic_cast<ROL::SingletonVector<Real>&>(c).setValue(V0_ + Vdual_->dot(z));
  }

  void applyJacobian(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &z, Real &tol) {
    njaco_++;
    dynamic_cast<ROL::SingletonVector<Real>&>(jv).setValue(Vdual_->dot(v));
  }

  void applyAdjointJacobian( ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    najac_++;
    Real vval = dynamic_cast<const ROL::SingletonVector<Real>&>(v).getValue();
    jv.set(*V_);
    jv.scale(vval);
  }

  void applyAdjointHessian( ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nhess_++;
    ahwv.zero();
  }

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > matWriter;
    if (V_ != ROL::nullPtr) {
      ROL::Ptr<Tpetra::MultiVector<>> V = getField(*V_);
      matWriter.writeDenseFile("volume_jacobian.txt", V);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("volume_jacobian.txt");
      emptyfile.close();
    }
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  Volume_Constraint::summarize" << std::endl;
    stream << "    Number of calls to update:               " << nupda_ << std::endl;
    stream << "    Number of calls to value:                " << nfval_ << std::endl;
    stream << "    Number of calls to applyJacobian:        " << njaco_ << std::endl;
    stream << "    Number of calls to applyAdjointJacobian: " << najac_ << std::endl;
    stream << "    Number of calls to applyAdjointHessian:  " << nhess_ << std::endl;
    stream << "    Number of assemblies:                    " << nasse_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
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
      ROL::Ptr<Tpetra::MultiVector<>> z0 = assembler_->createControlVector();
      z0->putScalar(0.0);
      V0_ = assembler_->assembleQoIValue(qoi_,ROL::nullPtr,z0);
      assembler_->assembleQoIGradient2(V,qoi_,ROL::nullPtr,z0);
      Vdual_->set(V_->dual());
      assemble_ = false;
      outputTpetraData();
    }
  }
}; // class Volume_Constraint

#endif
