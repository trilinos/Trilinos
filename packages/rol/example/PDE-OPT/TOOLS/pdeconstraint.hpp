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

  void setParameter(const std::vector<Real> &param) {
    pde_->setParameter(param);
    ROL::ParametrizedEqualityConstraint_SimOpt<Real>::setParameter(param);
  }

  Teuchos::RCP<Assembler<Real> > getAssembler(void) const {
    return assembler_;
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_1;
  void update_1(const ROL::Vector<Real> &u, bool flag = true, int iter = -1) {
    if (flag && iter > -1) {
      computeJ1_ = true;
    }
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update_2;
  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    if (flag && iter > -1) {
      computeJ2_ = true;
    }
  }

  using ROL::EqualityConstraint_SimOpt<Real>::update;
  void update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    update_1(u,flag,iter);
    update_2(z,flag,iter);
    if (flag && iter > -1) {
      computeH11_ = true;
      computeH12_ = true;
      computeH21_ = true;
      computeH22_ = true;
    }
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
    assembler_->assemblePDEResidual(*up,*zp,*pde_);
    cp->scale(one,*(assembler_->getPDEResidual()));
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if (computeJ1_) {
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

      assembler_->assemblePDEJacobian1(*up,*zp,*pde_);
      computeJ1_ = false;
    }
    assembler_->getPDEJacobian1(false)->apply(*vp,*jvp);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if (computeJ2_) {
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

      assembler_->assemblePDEJacobian2(*up,*zp,*pde_);
      computeJ2_ = false;
    }
    assembler_->getPDEJacobian2(false)->apply(*vp,*jvp);
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ajvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if (computeJ1_) {
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

      assembler_->assemblePDEJacobian1(*up,*zp,*pde_);
      computeJ1_ = false;
    }
    assembler_->getPDEJacobian1(true)->apply(*vp,*ajvp);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ajvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if (computeJ2_) {
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

      assembler_->assemblePDEJacobian2(*up,*zp,*pde_);
      computeJ2_ = false;
    }
    assembler_->getPDEJacobian2(true)->apply(*vp,*ajvp);
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      if (computeH11_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian11(*up,*zp,*wp,*pde_);
        computeH11_ = false;
      }
      assembler_->getPDEHessian11()->apply(*vp,*ahwvp);
    }
    catch (Exception::Zero &zero) {
      ahwv.zero();
    }
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      if (computeH12_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian12(*up,*zp,*wp,*pde_);
        computeH12_ = false;
      }
      assembler_->getPDEHessian12()->apply(*vp,*ahwvp);
    }
    catch (Exception::Zero &zero) {
      ahwv.zero();
    }
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      if (computeH21_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian21(*up,*zp,*wp,*pde_);
        computeH21_ = false;
      }
      assembler_->getPDEHessian21()->apply(*vp,*ahwvp);
    }
    catch (Exception::Zero &zero) {
      ahwv.zero();
    }
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      if (computeH22_) {
        Teuchos::RCP<const Tpetra::MultiVector<> > wp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > up =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
        Teuchos::RCP<const Tpetra::MultiVector<> > zp =
          (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

        assembler_->assemblePDEHessian22(*up,*zp,*wp,*pde_);
        computeH22_ = false;
      }
      assembler_->getPDEHessian22()->apply(*vp,*ahwvp);
    }
    catch (Exception::Zero &zero) {
      ahwv.zero();
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ijvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ijv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if (computeJ1_) {
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

      assembler_->assemblePDEJacobian1(*up,*zp,*pde_);
      computeJ1_ = false;
    }
    assembler_->linearPDEsolve(ijvp,vp,false);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                                     const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > iajvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(iajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if (computeJ1_) {
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

      assembler_->assemblePDEJacobian1(*up,*zp,*pde_);
      computeJ1_ = false;
    }
    assembler_->linearPDEsolve(iajvp,vp,true);
  }

};

#endif
