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

#ifndef PDE_INTEGRALCONSTRAINT_HPP
#define PDE_INTEGRALCONSTRAINT_HPP

#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"
#include "ROL_StdVector.hpp"
#include "qoi.hpp"
#include "assembler.hpp"

template<class Real>
class IntegralConstraint : public ROL::ParametrizedEqualityConstraint_SimOpt<Real> {
private:
  const Teuchos::RCP<QoI<Real> > qoi_;
  const Teuchos::RCP<Assembler<Real> > assembler_;
public:
  IntegralConstraint(const Teuchos::RCP<QoI<Real> > &qoi,
                     const Teuchos::RCP<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    qoi_->setParameter(param);
    ROL::ParametrizedEqualityConstraint_SimOpt<Real>::setParameter(param);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<std::vector<Real> > cp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(c)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    (*cp)[0] = assembler_->assembleQoIValue(*up,*zp,*qoi_);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<std::vector<Real> > jvp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    Teuchos::Array<Real> jdotv(1,0);
    assembler_->assembleQoIGradient1(*up,*zp,*qoi_);
    assembler_->getQoIGradient1()->dot(*vp,jdotv.view(0,1));
    (*jvp)[0] = jdotv[0];
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<std::vector<Real> > jvp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    Teuchos::Array<Real> jdotv(1,0);
    assembler_->assembleQoIGradient2(*up,*zp,*qoi_);
    assembler_->getQoIGradient2()->dot(*vp,jdotv.view(0,1));
    (*jvp)[0] = jdotv[0];
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    assembler_->assembleQoIGradient1(*up,*zp,*qoi_);
    jvp->scale((*vp)[0],*(assembler_->getQoIGradient1()));
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    assembler_->assembleQoIGradient2(*up,*zp,*qoi_);
    jvp->scale((*vp)[0],*(assembler_->getQoIGradient2()));
  }

  void applyAdjointHessian_11( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const std::vector<Real> > wp =
        (Teuchos::dyn_cast<const ROL::StdVector<Real> >(w)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec11(*vp,*up,*zp,*qoi_);
      ahwvp->scale((*wp)[0],*(assembler_->getQoIHessVec11()));
    }
    catch (Exception::Zero &ez) {
      ahwv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (IntegralConstraint::hessVec_11): Hessian not implemented.");
    }
  }

  void applyAdjointHessian_12( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const std::vector<Real> > wp =
        (Teuchos::dyn_cast<const ROL::StdVector<Real> >(w)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec12(*vp,*up,*zp,*qoi_);
      ahwvp->scale((*wp)[0],*(assembler_->getQoIHessVec12()));
    }
    catch (Exception::Zero &ez) {
      ahwv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (IntegralConstraint::hessVec_12): Hessian not implemented.");
    }
  }

  void applyAdjointHessian_21( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const std::vector<Real> > wp =
        (Teuchos::dyn_cast<const ROL::StdVector<Real> >(w)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec21(*vp,*up,*zp,*qoi_);
      ahwvp->scale((*wp)[0],*(assembler_->getQoIHessVec21()));
    }
    catch (Exception::Zero &ez) {
      ahwv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (IntegralConstraint::hessVec_21): Hessian not implemented.");
    }
  }

  void applyAdjointHessian_22( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > ahwvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
      Teuchos::RCP<const std::vector<Real> > wp =
        (Teuchos::dyn_cast<const ROL::StdVector<Real> >(w)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec22(*vp,*up,*zp,*qoi_);
      ahwvp->scale((*wp)[0],*(assembler_->getQoIHessVec22()));
    }
    catch (Exception::Zero &ez) {
      ahwv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (IntegralConstraint::hessVec_22): Hessian not implemented.");
    }
  }
};

#endif
