// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const ROL::Ptr<QoI<Real> > qoi_;
  const ROL::Ptr<Assembler<Real> > assembler_;

  ROL::Ptr<Tpetra::MultiVector<> > vecG1_;
  ROL::Ptr<Tpetra::MultiVector<> > vecG2_;
  ROL::Ptr<std::vector<Real> >     vecG3_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH11_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH12_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH13_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH21_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH22_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH23_;
  ROL::Ptr<std::vector<Real> >     vecH31_;
  ROL::Ptr<std::vector<Real> >     vecH32_;
  ROL::Ptr<std::vector<Real> >     vecH33_;

public:
  IntegralObjective(const ROL::Ptr<QoI<Real> > &qoi,
                    const ROL::Ptr<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective_SimOpt<Real>::setParameter(param);
    qoi_->setParameter(param);
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
    return assembler_->assembleQoIValue(qoi_,uf,zf,zp);
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
    try {
      ROL::Ptr<Tpetra::MultiVector<> >       gf = getField(g);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      ROL::Ptr<Tpetra::MultiVector<> >       gf = getField(g);
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
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      ROL::Ptr<std::vector<Real> >           gp = getParameter(g);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<std::vector<Real> >          hvp = getParameter(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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
      ROL::Ptr<std::vector<Real> >          hvp = getParameter(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec32(vecH32_,qoi_,vf,uf,zf,zp);
      hvp->assign(vecH32_->begin(),vecH32_->end());
    }
    catch (Exception::Zero &ez) {
      ROL::Ptr<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != ROL::nullPtr ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::Ptr<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != ROL::nullPtr ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      NotImplemented++;
    }
    // Compute control parameter/parameter hessvec
    try {
      ROL::Ptr<std::vector<Real> >          hvp = getParameter(hv);
      ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
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

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::StdVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
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

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::StdVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
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
};


template<class Real>
class IntegralSimObjective : public ROL::Objective<Real> {
private:
  const ROL::Ptr<QoI<Real>> qoi_;
  const ROL::Ptr<Assembler<Real>> assembler_;

  ROL::Ptr<Tpetra::MultiVector<>> vecG1_;
  ROL::Ptr<Tpetra::MultiVector<>> vecH11_;

public:
  IntegralSimObjective(const ROL::Ptr<QoI<Real>> &qoi,
                       const ROL::Ptr<Assembler<Real>> &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    qoi_->setParameter(param);
  }

  Real value(const ROL::Vector<Real> &u, Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<>> uf = getConstField(u);
    return assembler_->assembleQoIValue(qoi_,uf,ROL::nullPtr,ROL::nullPtr);
  }

  void gradient(ROL::Vector<Real> &g,
                const ROL::Vector<Real> &u, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    g.zero();
    // Compute control field gradient
    try {
      // Get state and control vectors
      ROL::Ptr<const Tpetra::MultiVector<>> uf = getConstField(u);
      ROL::Ptr<Tpetra::MultiVector<>>       gf = getField(g);
      assembler_->assembleQoIGradient1(vecG1_,qoi_,uf,ROL::nullPtr,ROL::nullPtr);
      gf->scale(static_cast<Real>(1),*vecG1_);
    }
    catch ( Exception::Zero & ez ) {
      IsZero++;
    }
    catch ( Exception::NotImplemented & eni ) {
      NotImplemented++;
    }
    // Zero gradient
    if ( IsZero == 1 ) {
      g.zero();
    }
    // Not Implemented
    if ( NotImplemented == 1 ) {
      ROL::Objective<Real>::gradient(g,u,tol);
    }
  }

  void hessVec(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
               const ROL::Vector<Real> &u, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    hv.zero();
    // Compute control field/field hessvec
    try {
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
      assembler_->assembleQoIHessVec11(vecH11_,qoi_,vf,uf,ROL::nullPtr,ROL::nullPtr);
      hvf->scale(static_cast<Real>(1),*vecH11_);
    }
    catch (Exception::Zero &ez) {
      hv.zero();
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      hv.zero();
      NotImplemented++;
    }
    // Zero hessvec
    if ( IsZero == 1 ) {
      hv.zero();
    }
    // Not Implemented
    if ( NotImplemented == 1 ) {
      ROL::Objective<Real>::hessVec(hv,v,u,tol);
    }
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<> > xp;
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
};

template<class Real>
class IntegralOptObjective : public ROL::Objective<Real> {
private:
  const ROL::Ptr<QoI<Real> > qoi_;
  const ROL::Ptr<Assembler<Real> > assembler_;

  ROL::Ptr<Tpetra::MultiVector<> > vecG2_;
  ROL::Ptr<std::vector<Real> >     vecG3_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH22_;
  ROL::Ptr<Tpetra::MultiVector<> > vecH23_;
  ROL::Ptr<std::vector<Real> >     vecH32_;
  ROL::Ptr<std::vector<Real> >     vecH33_;

public:
  IntegralOptObjective(const ROL::Ptr<QoI<Real> > &qoi,
                       const ROL::Ptr<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    qoi_->setParameter(param);
  }

  Real value(const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
    return assembler_->assembleQoIValue(qoi_,ROL::nullPtr,zf,zp);
  }

  void gradient(ROL::Vector<Real> &g,
                const ROL::Vector<Real> &z, Real &tol ) {
    int NotImplemented(0), IsZero(0);
    g.zero();
    // Compute control field gradient
    try {
      // Get state and control vectors
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      ROL::Ptr<Tpetra::MultiVector<> >       gf = getField(g);
      assembler_->assembleQoIGradient2(vecG2_,qoi_,ROL::nullPtr,zf,zp);
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
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      ROL::Ptr<std::vector<Real> >           gp = getParameter(g);
      assembler_->assembleQoIGradient3(vecG3_,qoi_,ROL::nullPtr,zf,zp);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec22(vecH22_,qoi_,vf,ROL::nullPtr,zf,zp);
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
      ROL::Ptr<Tpetra::MultiVector<> >      hvf = getField(hv);
      ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec23(vecH23_,qoi_,vp,ROL::nullPtr,zf,zp);
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
      ROL::Ptr<std::vector<Real> >          hvp = getParameter(hv);
      ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec32(vecH32_,qoi_,vf,ROL::nullPtr,zf,zp);
      hvp->assign(vecH32_->begin(),vecH32_->end());
    }
    catch (Exception::Zero &ez) {
      ROL::Ptr<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != ROL::nullPtr ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      IsZero++;
    }
    catch (Exception::NotImplemented &eni) {
      ROL::Ptr<std::vector<Real> > hvp = getParameter(hv);
      if ( hvp != ROL::nullPtr ) {
        const int size = hvp->size();
        hvp->assign(size,static_cast<Real>(0));
      }
      NotImplemented++;
    }
    // Compute control parameter/parameter hessvec
    try {
      ROL::Ptr<std::vector<Real> >          hvp = getParameter(hv);
      ROL::Ptr<const std::vector<Real> >     vp = getConstParameter(v);
      ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
      ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);
      assembler_->assembleQoIHessVec33(vecH33_,qoi_,vp,ROL::nullPtr,zf,zp);
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

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::StdVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
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

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::StdVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
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
};

template<class Real>
class IntegralAffineSimObjective : public ROL::Objective_SimOpt<Real> {
private:
  const ROL::Ptr<QoI<Real> > qoi_;
  const ROL::Ptr<Assembler<Real> > assembler_;

  Real shift_;
  ROL::Ptr<Tpetra::MultiVector<> > vecG1_;
  ROL::Ptr<Tpetra::MultiVector<> > ufield_;
  ROL::Ptr<Tpetra::MultiVector<> > zfield_;
  ROL::Ptr<std::vector<Real> >     zparam_;

  bool isAssembled_;

  void assemble(const ROL::Ptr<const Tpetra::MultiVector<> > &zf,
                const ROL::Ptr<const std::vector<Real> >     &zp) {
    if ( !isAssembled_ ) {
      ufield_ = assembler_->createStateVector();
      ufield_->putScalar(static_cast<Real>(0));
      if ( zf != ROL::nullPtr ) {
        zfield_ = assembler_->createControlVector();
        zfield_->putScalar(static_cast<Real>(0));
      }
      if ( zp != ROL::nullPtr ) {
       zparam_ = ROL::makePtr<std::vector<Real>>(zp->size(),static_cast<Real>(0));
      }
      shift_ = assembler_->assembleQoIValue(qoi_,ufield_,zfield_,zparam_);
      assembler_->assembleQoIGradient1(vecG1_,qoi_,ufield_,zfield_,zparam_);
      isAssembled_ = true;
    }
  }

public:
  IntegralAffineSimObjective(const ROL::Ptr<QoI<Real> > &qoi,
                             const ROL::Ptr<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler), isAssembled_(false) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective_SimOpt<Real>::setParameter(param);
    qoi_->setParameter(param);
    isAssembled_ = false;
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

    assemble(zf,zp);

    const size_t n = uf->getNumVectors();
    Teuchos::Array<Real> val(n,0);
    vecG1_->dot(*uf, val.view(0,n));
    
    return val[0] + shift_;
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u,
                  const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<Tpetra::MultiVector<> >       gf = getField(g);
    ROL::Ptr<const Tpetra::MultiVector<> > zf = getConstField(z);
    ROL::Ptr<const std::vector<Real> >     zp = getConstParameter(z);

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

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<> > xp;
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

  ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::StdVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
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

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::StdVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
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
};

#endif
