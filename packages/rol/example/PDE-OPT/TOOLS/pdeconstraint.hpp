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

#ifndef ROL_PDEOPT_PDECONSTRAINT_H
#define ROL_PDEOPT_PDECONSTRAINT_H

#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "pde.hpp"
#include "assembler.hpp"

template<class Real>
class PDE_Constraint : public ROL::ParametrizedEqualityConstraint_SimOpt<Real> {
private:
  const Teuchos::RCP<PDE<Real> > pde_;
  Teuchos::RCP<Assembler<Real> > assembler_;

  bool computeJ1_, computeJ2_;
  bool computeH11_, computeH12_, computeH21_, computeH22_;

public:
  PDE_Constraint(const Teuchos::RCP<PDE<Real> > &pde,
                 const Teuchos::RCP<MeshManager<Real> > &meshMgr,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 Teuchos::ParameterList &parlist,
                 std::ostream &outStream = std::cout)
    : pde_(pde), computeJ1_(true), computeJ2_(true),
      computeH11_(true), computeH12_(true), computeH21_(true), computeH22_(true) {
    assembler_ = Teuchos::rcp(new Assembler<Real>(pde_->getFields(),meshMgr,comm,parlist,outStream));
    assembler_->setCellNodes(*pde_);
  }

  PDE_Constraint(const Teuchos::RCP<PDE<Real> > &pde,
                 const Teuchos::RCP<Assembler<Real> > &assembler)
    : pde_(pde), assembler_(assembler), computeJ1_(true), computeJ2_(true),
      computeH11_(true), computeH12_(true), computeH21_(true), computeH22_(true) {
    assembler_->setCellNodes(*pde_);
  }

  void setParameter(const std::vector<Real> &param) {
    pde_->setParameter(param);
    ROL::ParametrizedEqualityConstraint_SimOpt<Real>::setParameter(param);
  }

  const Teuchos::RCP<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {
    computeJ1_ = (flag ? true : computeJ1_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    computeJ2_ = (flag ? true : computeJ2_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
    computeH11_ = (flag ? true : computeH11_);
    computeH12_ = (flag ? true : computeH12_);
    computeH21_ = (flag ? true : computeH21_);
    computeH22_ = (flag ? true : computeH22_);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > cp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(c)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

    Real one(1);
    assembler_->assemblePDEResidual(pde_,up,zp);
    cp->scale(one,*(assembler_->getPDEResidual()));
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEJacobian1(pde_,up,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
     computeJ1_ = true;
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_1(jv,v,u,z,tol);
    }
    if ( computeJ1_ ) {
      jv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > jvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEJacobian1(jvp,vp,false);
    }
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ2_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEJacobian2(pde_,up,zp);
        computeJ2_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
     computeJ2_ = true;
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyJacobian_2(jv,v,u,z,tol);
    }
    if ( computeJ2_ ) {
      jv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > jvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEJacobian2(jvp,vp,false);
    }
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEJacobian1(pde_,up,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
     computeJ1_ = true;
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointJacobian_1(ajv,v,u,z,tol);
    }
    if ( computeJ1_ ) {
      ajv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > ajvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEJacobian1(ajvp,vp,true);
    }
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ2_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEJacobian2(pde_,up,zp);
        computeJ2_ = false;
      }
    }
    catch ( Exception::Zero & ez ) {
      computeJ2_ = true;
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointJacobian_2(ajv,v,u,z,tol);
    }
    if ( computeJ2_ ) {
      ajv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > ajvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEJacobian2(ajvp,vp,true);
    }
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeH11_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian11(pde_,wp,up,zp);
        computeH11_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH11_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_11(ahwv,w,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (PDE_Constraint::applyAdjointHessian_11): Hessian not implemented.");
    }
    if ( computeH11_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEHessian11(ahwvp,vp);
    }
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      if ( computeH12_ ) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian12(pde_,wp,up,zp);
        computeH12_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH12_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_12(ahwv,w,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (PDE_Constraint::applyAdjointHessian_12): Hessian not implemented.");
    }
    if ( computeH12_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEHessian12(ahwvp,vp);
    }
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeH21_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian21(pde_,wp,up,zp);
        computeH21_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH21_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_21(ahwv,w,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (PDE_Constraint::applyAdjointHessian_21): Hessian not implemented.");
    }
    if ( computeH21_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEHessian21(ahwvp,vp);
    }
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeH22_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian22(pde_,wp,up,zp);
        computeH22_ = false;
      }
    }
    catch (Exception::Zero &ez) {
      computeH22_ = true;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::EqualityConstraint_SimOpt<Real>::applyAdjointHessian_22(ahwv,w,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (PDE_Constraint::applyAdjointHessian_22): Hessian not implemented.");
    }
    if ( computeH22_ ) {
      ahwv.zero();
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      assembler_->applyPDEHessian22(ahwvp,vp);
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEJacobian1(pde_,up,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyInverseJacobian_1(ijv,v,u,z,tol);
    }
    Teuchos::RCP<Tpetra::MultiVector<> > ijvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ijv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    assembler_->applyInverseJacobian1(ijvp,vp,false);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                                     const ROL::Vector<Real> &z, Real &tol) {
    try {
      if (computeJ1_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEJacobian1(pde_,up,zp);
        computeJ1_ = false;
      }
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::EqualityConstraint_SimOpt<Real>::applyInverseAdjointJacobian_1(iajv,v,u,z,tol);
    }
    Teuchos::RCP<Tpetra::MultiVector<> > iajvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(iajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    assembler_->applyInverseJacobian1(iajvp,vp,true);
  }

};

#endif
