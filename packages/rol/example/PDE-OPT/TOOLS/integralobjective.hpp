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

#ifndef PDE_INTEGRALOBJECTIVE_HPP
#define PDE_INTEGRALOBJECTIVE_HPP

#include "ROL_ParametrizedObjective_SimOpt.hpp"
#include "qoi.hpp"
#include "assembler.hpp"

template<class Real>
class IntegralObjective : public ROL::ParametrizedObjective_SimOpt<Real> {
private:
  const Teuchos::RCP<QoI<Real> > qoi_;
  const Teuchos::RCP<Assembler<Real> > assembler_;
public:
  IntegralObjective(const Teuchos::RCP<QoI<Real> > &qoi,
                    const Teuchos::RCP<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    qoi_->setParameter(param);
    ROL::ParametrizedObjective_SimOpt<Real>::setParameter(param);
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    return assembler_->assembleQoIValue(qoi_,up,zp);
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > gp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIGradient1(qoi_,up,zp);
      gp->scale(static_cast<Real>(1),*(assembler_->getQoIGradient1()));
    }
    catch ( Exception::Zero & ez ) {
      g.zero();
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::ParametrizedObjective_SimOpt<Real>::gradient_1(g,u,z,tol);
    }
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > gp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIGradient2(qoi_,up,zp);
      gp->scale(static_cast<Real>(1),*(assembler_->getQoIGradient2()));
    }
    catch ( Exception::Zero & ez ) {
      g.zero();
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::ParametrizedObjective_SimOpt<Real>::gradient_2(g,u,z,tol);
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > hvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec11(qoi_,vp,up,zp);
      hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec11()));
    }
    catch (Exception::Zero &ez) {
      hv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      ROL::ParametrizedObjective_SimOpt<Real>::hessVec_11(hv,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (IntegratedObjective::hessVec_11): Hessian not implemented.");
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > hvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec12(qoi_,vp,up,zp);
      hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec12()));
    }
    catch (Exception::Zero &ez) {
      hv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      ROL::ParametrizedObjective_SimOpt<Real>::hessVec_12(hv,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (IntegratedObjective::hessVec_12): Hessian not implemented.");
    }
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > hvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec21(qoi_,vp,up,zp);
      hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec21()));
    }
    catch (Exception::Zero &ez) {
      hv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      ROL::ParametrizedObjective_SimOpt<Real>::hessVec_21(hv,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (IntegratedObjective::hessVec_21): Hessian not implemented.");
    }
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    try {
      Teuchos::RCP<Tpetra::MultiVector<> > hvp =
        (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > zp =
        (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
      assembler_->assembleQoIHessVec22(qoi_,vp,up,zp);
      hvp->scale(static_cast<Real>(1),*(assembler_->getQoIHessVec22()));
    }
    catch (Exception::Zero &ez) {
      hv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      ROL::ParametrizedObjective_SimOpt<Real>::hessVec_22(hv,v,u,z,tol);
      //throw Exception::NotImplemented(">>> (IntegratedObjective::hessVec_22): Hessian not implemented.");
    }
  }
};

#endif
