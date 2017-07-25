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

#include "ROL_Objective_SimOpt.hpp"
#include "qoi.hpp"
#include "assembler.hpp"
#include "pdevector.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

template<class Real>
class IntegralObjective : public ROL::Objective_SimOpt<Real> {
private:
  const Teuchos::RCP<QoI<Real> > qoi_;
  const Teuchos::RCP<Assembler<Real> > assembler_;

  Teuchos::RCP<Tpetra::MultiVector<> > vecG1_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecG2_;
  Teuchos::RCP<std::vector<Real> >     vecG3_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH11_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH12_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH13_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH21_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH22_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH23_;
  Teuchos::RCP<std::vector<Real> >     vecH31_;
  Teuchos::RCP<std::vector<Real> >     vecH32_;
  Teuchos::RCP<std::vector<Real> >     vecH33_;

public:
  IntegralObjective(const Teuchos::RCP<QoI<Real> > &qoi,
                    const Teuchos::RCP<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective_SimOpt<Real>::setParameter(param);
    qoi_->setParameter(param);
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
    return assembler_->assembleQoIValue(qoi_,uf,zf,zp);
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
    try {
      Teuchos::RCP<Tpetra::MultiVector<> >       gf = getField(g);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIGradient1(vecG1_,qoi_,uf,zf,zp);
      gf->scale(static_cast<Real>(1),*vecG1_);
    }
    catch ( Exception::Zero & ez ) {
      g.zero();
    }
    catch ( Exception::NotImplemented & eni ) {
      ROL::Objective_SimOpt<Real>::gradient_1(g,u,z,tol);
    }
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    g.zero();
    // Compute control field gradient
    try {
      // Get state and control vectors
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      Teuchos::RCP<Tpetra::MultiVector<> >       gf = getField(g);
      assembler_->assembleQoIGradient2(vecG2_,qoi_,uf,zf,zp);
      gf->scale(static_cast<Real>(1),*vecG2_);
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Compute control parameter gradient
    try {
      // Get state and control vectors
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      Teuchos::RCP<std::vector<Real> >           gp = getParameter(g);
      assembler_->assembleQoIGradient3(vecG3_,qoi_,uf,zf,zp);
      gp->assign(vecG3_->begin(),vecG3_->end());
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Zero gradient
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      g.zero();
    }
    // Not Implemented
    if ( NotImplemented == 2 ) {
      ROL::Objective_SimOpt<Real>::gradient_2(g,u,z,tol);
    }
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
    try {
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec11(vecH11_,qoi_,vf,uf,zf,zp);
      hvf->scale(static_cast<Real>(1),*vecH11_);
    }
    catch (Exception::Zero &ez) {
      hv.zero();
    }
    catch (Exception::NotImplemented &eni) {
      ROL::Objective_SimOpt<Real>::hessVec_11(hv,v,u,z,tol);
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    hv.zero();
    // Compute state field/control field hessvec
    try {
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec12(vecH12_,qoi_,vf,uf,zf,zp);
      hvf->scale(static_cast<Real>(1),*vecH12_);
    }
    catch (Exception::Zero &ez) {
      hv.zero();
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      hv.zero();
      NotImplemented++;
    }
    // Compute state field/control parameter hessvec
    try {
      // Get state and control vectors
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec13(vecH13_,qoi_,vp,uf,zf,zp);
      hvf->update(static_cast<Real>(1),*vecH13_,static_cast<Real>(1));
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Zero hessvec
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      hv.zero();
    }
    // Not Implemented
    if ( NotImplemented == 2 ) {
      ROL::Objective_SimOpt<Real>::hessVec_12(hv,v,u,z,tol);
    }
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    hv.zero();
    // Compute control field/state field hessvec
    try {
      // Get state and control vectors
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec21(vecH21_,qoi_,vf,uf,zf,zp);
      hvf->scale(static_cast<Real>(1),*vecH21_);
    }
    catch ( Exception::Zero & ez ) {
      hv.zero();
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      hv.zero();
      NotImplemented++;
    }
    // Compute control parameter/state field hessvec
    try {
      // Get state and control vectors
      Teuchos::RCP<std::vector<Real> >          hvp = getParameter(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec31(vecH31_,qoi_,vf,uf,zf,zp);
      hvp->assign(vecH31_->begin(),vecH31_->end());
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Zero hessvec
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      hv.zero();
    }
    // Not Implemented
    if ( NotImplemented == 2 ) {
      ROL::Objective_SimOpt<Real>::hessVec_21(hv,v,u,z,tol);
    }
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    hv.zero();
    // Compute control field/field hessvec
    try {
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec22(vecH22_,qoi_,vf,uf,zf,zp);
      hvf->scale(static_cast<Real>(1),*vecH22_);
    }
    catch (Exception::Zero &ez) {
      hv.zero();
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      hv.zero();
      NotImplemented++;
    }
    // Compute control field/parameter hessvec
    try {
      // Get state and control vectors
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec23(vecH23_,qoi_,vp,uf,zf,zp);
      hvf->update(static_cast<Real>(1),*vecH23_,static_cast<Real>(1));
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Compute control parameter/field hessvec
    try {
      Teuchos::RCP<std::vector<Real> >          hvp = getParameter(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec32(vecH32_,qoi_,vf,uf,zf,zp);
      hvp->assign(vecH32_->begin(),vecH32_->end());
    }
    catch (Exception::Zero &ez) {
      Teuchos::RCP<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != Teuchos::null ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      Teuchos::RCP<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != Teuchos::null ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      NotImplemented++;
    }
    // Compute control parameter/parameter hessvec
    try {
      Teuchos::RCP<std::vector<Real> >          hvp = getParameter(hv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec33(vecH33_,qoi_,vp,uf,zf,zp);
      const int size = hvp->size();
      for (int i = 0; i < size; ++i) {
        (*hvp)[i] += (*vecH33_)[i];
      }
    }
    catch (Exception::Zero &ez) {
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      NotImplemented++;
    }

    // Zero hessvec
    if ( IsZero > 0 && (IsZero + NotImplemented == 4) ) {
      hv.zero();
    }
    // Not Implemented
    if ( NotImplemented == 4 ) {
      ROL::Objective_SimOpt<Real>::hessVec_22(hv,v,u,z,tol);
    }
  }

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<const ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getField();
      if (xvec == Teuchos::null) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    Teuchos::RCP<Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getField();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const std::vector<Real> > xp;
    try {
      Teuchos::RCP<const ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }

  Teuchos::RCP<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    Teuchos::RCP<std::vector<Real> > xp;
    try {
      Teuchos::RCP<ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }
};


template<class Real>
class IntegralOptObjective : public ROL::Objective<Real> {
private:
  const Teuchos::RCP<QoI<Real> > qoi_;
  const Teuchos::RCP<Assembler<Real> > assembler_;

  Teuchos::RCP<Tpetra::MultiVector<> > vecG2_;
  Teuchos::RCP<std::vector<Real> >     vecG3_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH22_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecH23_;
  Teuchos::RCP<std::vector<Real> >     vecH32_;
  Teuchos::RCP<std::vector<Real> >     vecH33_;

public:
  IntegralOptObjective(const Teuchos::RCP<QoI<Real> > &qoi,
                       const Teuchos::RCP<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    qoi_->setParameter(param);
  }

  Real value(const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
    return assembler_->assembleQoIValue(qoi_,Teuchos::null,zf,zp);
  }

  void gradient(ROL::Vector<Real> &g,
                const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    g.zero();
    // Compute control field gradient
    try {
      // Get state and control vectors
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      Teuchos::RCP<Tpetra::MultiVector<> >       gf = getField(g);
      assembler_->assembleQoIGradient2(vecG2_,qoi_,Teuchos::null,zf,zp);
      gf->scale(static_cast<Real>(1),*vecG2_);
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Compute control parameter gradient
    try {
      // Get state and control vectors
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      Teuchos::RCP<std::vector<Real> >           gp = getParameter(g);
      assembler_->assembleQoIGradient3(vecG3_,qoi_,Teuchos::null,zf,zp);
      gp->assign(vecG3_->begin(),vecG3_->end());
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Zero gradient
    if ( IsZero == 2 || (IsZero == 1 && NotImplemented == 1) ) {
      g.zero();
    }
    // Not Implemented
    if ( NotImplemented == 2 ) {
      ROL::Objective<Real>::gradient(g,z,tol);
    }
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
               const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    hv.zero();
    // Compute control field/field hessvec
    try {
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec22(vecH22_,qoi_,vf,Teuchos::null,zf,zp);
      hvf->scale(static_cast<Real>(1),*vecH22_);
    }
    catch (Exception::Zero &ez) {
      hv.zero();
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      hv.zero();
      NotImplemented++;
    }
    // Compute control field/parameter hessvec
    try {
      // Get state and control vectors
      Teuchos::RCP<Tpetra::MultiVector<> >      hvf = getField(hv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec23(vecH23_,qoi_,vp,Teuchos::null,zf,zp);
      hvf->update(static_cast<Real>(1),*vecH23_,static_cast<Real>(1));
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Compute control parameter/field hessvec
    try {
      Teuchos::RCP<std::vector<Real> >          hvp = getParameter(hv);
      Teuchos::RCP<const Tpetra::MultiVector<> > vf = getConstField(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec32(vecH32_,qoi_,vf,Teuchos::null,zf,zp);
      hvp->assign(vecH32_->begin(),vecH32_->end());
    }
    catch (Exception::Zero &ez) {
      Teuchos::RCP<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != Teuchos::null ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      Teuchos::RCP<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != Teuchos::null ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      NotImplemented++;
    }
    // Compute control parameter/parameter hessvec
    try {
      Teuchos::RCP<std::vector<Real> >          hvp = getParameter(hv);
      Teuchos::RCP<const std::vector<Real> >     vp = getConstParameter(v);
      Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
      Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec33(vecH33_,qoi_,vp,Teuchos::null,zf,zp);
      const int size = hvp->size();
      for (int i = 0; i < size; ++i) {
        (*hvp)[i] += (*vecH33_)[i];
      }
    }
    catch (Exception::Zero &ez) {
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      NotImplemented++;
    }

    // Zero hessvec
    if ( IsZero > 0 && (IsZero + NotImplemented == 4) ) {
      hv.zero();
    }
    // Not Implemented
    if ( NotImplemented == 4 ) {
      ROL::Objective<Real>::hessVec(hv,v,z,tol);
    }
  }

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<const ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getField();
      if (xvec == Teuchos::null) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    Teuchos::RCP<Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getField();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const std::vector<Real> > xp;
    try {
      Teuchos::RCP<const ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }

  Teuchos::RCP<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    Teuchos::RCP<std::vector<Real> > xp;
    try {
      Teuchos::RCP<ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }
};

template<class Real>
class IntegralAffineSimObjective : public ROL::Objective_SimOpt<Real> {
private:
  const Teuchos::RCP<QoI<Real> > qoi_;
  const Teuchos::RCP<Assembler<Real> > assembler_;

  Real shift_;
  Teuchos::RCP<Tpetra::MultiVector<> > vecG1_;
  Teuchos::RCP<Tpetra::MultiVector<> > ufield_;
  Teuchos::RCP<Tpetra::MultiVector<> > zfield_;
  Teuchos::RCP<std::vector<Real> >     zparam_;

  bool isAssembled_;

  void assemble(const Teuchos::RCP<const Tpetra::MultiVector<> > &zf,
                const Teuchos::RCP<const std::vector<Real> >     &zp) {
    if ( !isAssembled_ ) {
      ufield_ = assembler_->createStateVector();
      ufield_->putScalar(static_cast<Real>(0));
      if ( zf != Teuchos::null ) {
        zfield_ = assembler_->createControlVector();
        zfield_->putScalar(static_cast<Real>(0));
      }
      if ( zp != Teuchos::null ) {
       zparam_ = Teuchos::rcp(new std::vector<Real>(zp->size(),static_cast<Real>(0)));
      }
      shift_ = assembler_->assembleQoIValue(qoi_,ufield_,zfield_,zparam_);
      assembler_->assembleQoIGradient1(vecG1_,qoi_,ufield_,zfield_,zparam_);
      isAssembled_ = true;
    }
  }

public:
  IntegralAffineSimObjective(const Teuchos::RCP<QoI<Real> > &qoi,
                             const Teuchos::RCP<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler), isAssembled_(false) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective_SimOpt<Real>::setParameter(param);
    qoi_->setParameter(param);
    isAssembled_ = false;
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > uf = getConstField(u);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    const size_t n = uf->getNumVectors();
    Teuchos::Array<Real> val(n,0);
    vecG1_->dot(*uf, val.view(0,n));
    
    return val[0] + shift_;
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::RCP<Tpetra::MultiVector<> >       gf = getField(g);
    Teuchos::RCP<const Tpetra::MultiVector<> > zf = getConstField(z);
    Teuchos::RCP<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    gf->scale(static_cast<Real>(1),*vecG1_);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

private: // Vector accessor functions

  Teuchos::RCP<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<const ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getField();
      if (xvec == Teuchos::null) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    Teuchos::RCP<Tpetra::MultiVector<> > xp;
    try {
      xp = Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(x).getVector();
    }
    catch (std::exception &e) {
      Teuchos::RCP<ROL::TpetraMultiVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getField();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    return xp;
  }

  Teuchos::RCP<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    Teuchos::RCP<const std::vector<Real> > xp;
    try {
      Teuchos::RCP<const ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }

  Teuchos::RCP<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    Teuchos::RCP<std::vector<Real> > xp;
    try {
      Teuchos::RCP<ROL::StdVector<Real> > xvec
        = Teuchos::dyn_cast<PDE_OptVector<Real> >(x).getParameter();
      if ( xvec == Teuchos::null ) {
        xp = Teuchos::null;
      }
      else {
        xp = xvec->getVector();
      }
    }
    catch (std::exception &e) {
      xp = Teuchos::null;
    }
    return xp;
  }
};

#endif
