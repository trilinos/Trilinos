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

/*! \file  con_selection.hpp
    \brief Defines the selectiont constraint.
*/

#ifndef ROL_PDEOPT_SELECTION_H
#define ROL_PDEOPT_SELECTION_H

#include "ROL_Constraint.hpp"
#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/assembler.hpp"
#include "../../../TOOLS/pdevector.hpp"

template<class Real>
class MultiMat_Selection_Constraint : public ROL::Constraint<Real> {
private:
  const ROL::Ptr<PDE<Real>> pde_;
  ROL::Ptr<Assembler<Real>> assembler_;
  ROL::Ptr<Tpetra::MultiVector<>> vecR_;
  ROL::Ptr<Tpetra::CrsMatrix<>>   matJ_;

  unsigned nupda_;
  unsigned nfval_;
  unsigned njaco_;
  unsigned najac_;
  unsigned nhess_;
  unsigned nasse_;

public:
  MultiMat_Selection_Constraint(const ROL::Ptr<PDE<Real>> &pde,
                                const ROL::Ptr<MeshManager<Real>> &meshMgr,
                                const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                                ROL::ParameterList &parlist,
                                std::ostream &outStream = std::cout)
    : pde_(pde),
      nupda_(0), nfval_(0), njaco_(0), najac_(0), nhess_(0), nasse_(0) {
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),pde_->getFields2(),meshMgr,comm,parlist,outStream);
    assembler_->setCellNodes(*pde_);
    assemble();
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<PDE<Real>> getPDE(void) const {
    return pde_;
  }

  void update( const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    nupda_++;
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &z, Real &tol) {
    nfval_++;
    ROL::Ptr<Tpetra::MultiVector<>>       cf = getField(c);
    ROL::Ptr<const Tpetra::MultiVector<>> zf = getConstField(z);
    buildResidual(cf,zf);
  }

  void applyJacobian(ROL::Vector<Real> &jv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &z, Real &tol) {
    njaco_++;
    ROL::Ptr<Tpetra::MultiVector<>>      jvf = getField(jv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
    applyJacobian(jvf,vf);
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                      const ROL::Vector<Real> &v,
                      const ROL::Vector<Real> &z, Real &tol) {
    najac_++;
    ROL::Ptr<Tpetra::MultiVector<> >     ajvf = getField(ajv);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    applyJacobian(ajvf,vf,true);
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                     const ROL::Vector<Real> &w,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &z, Real &tol) {
    nhess_++;
    ahwv.zero();
  }

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > matWriter;
    if (matJ_ != ROL::nullPtr) {
      matWriter.writeSparseFile("selection_jacobian.txt", matJ_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("selection_jacobian.txt");
      emptyfile.close();
    }
    if (vecR_ != ROL::nullPtr) {
      matWriter.writeDenseFile("selection_rhs.txt", vecR_);
    }
    else {
      std::ofstream emptyfile;
      emptyfile.open("selection_rhs.txt");
      emptyfile.close();
    }
  }

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    assembler_->outputTpetraVector(vec, filename);
  }

  void inputTpetraVector(ROL::Ptr<Tpetra::MultiVector<>> &vec,
                         const std::string &filename) const {
    assembler_->inputTpetraVector(vec, filename);
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  Selection_Constraint::summarize" << std::endl;
    stream << "    Number of calls to update:               " << nupda_ << std::endl;
    stream << "    Number of calls to value:                " << nfval_ << std::endl;
    stream << "    Number of calls to applyJacobian:        " << njaco_ << std::endl;
    stream << "    Number of calls to applyAdjointJacobian: " << najac_ << std::endl;
    stream << "    Number of calls to applyAdjointHessian:  " << nhess_ << std::endl;
    stream << "    Number of assemblies:                    " << nasse_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<>> getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<>> xp;
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

  ROL::Ptr<Tpetra::MultiVector<>> getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<>> xp;
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

  void assemble(void) {
    nasse_++;
    ROL::Ptr<Tpetra::MultiVector<>> u = assembler_->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> z = assembler_->createControlVector();
    u->putScalar(static_cast<Real>(0)); z->putScalar(static_cast<Real>(0));
    assembler_->assemblePDEResidual(vecR_, pde_, u, z);
    assembler_->assemblePDEJacobian2(matJ_, pde_, u, z);
    outputTpetraData();
  }

  // Application routines
  void applyJacobian(const ROL::Ptr<Tpetra::MultiVector<>> &Jv,
                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                     bool transpose = false) const {
    if (transpose) {
      matJ_->apply(*v,*Jv,Teuchos::TRANS);
    }
    else {
      matJ_->apply(*v,*Jv);
    }
  }

  void buildResidual(const ROL::Ptr<Tpetra::MultiVector<>> &r,
                     const ROL::Ptr<const Tpetra::MultiVector<>> &z) const {
    applyJacobian(r,z);
    r->update(static_cast<Real>(1),*vecR_,static_cast<Real>(1));
  }
};

#endif
