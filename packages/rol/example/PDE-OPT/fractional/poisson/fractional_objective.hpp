// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FRACTIONAL_OBJECTIVE_SIMOPT_H
#define ROL_FRACTIONAL_OBJECTIVE_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"

template <class Real>
class FractionalObjective : public ROL::Objective_SimOpt<Real> {
private:
  const ROL::Ptr<ROL::Objective_SimOpt<Real> > obj_;

public:
  FractionalObjective(const ROL::Ptr<ROL::Objective_SimOpt<Real> > &obj)
    : obj_(obj) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective_SimOpt<Real>::setParameter(param);
    obj_->setParameter(param);
  }

  void update( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    obj_->update(u,z,flag,iter);
  }

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    return obj_->value(*ur,z,tol);
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<Tpetra::MultiVector<> > gf = getField(g);
    ROL::Ptr<ROL::Vector<Real> > gr
      = ROL::makePtr<ROL::TpetraMultiVector<Real>(gf->subViewNonConst(cols>()));
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    g.zero();
    obj_->gradient_1(*gr,*ur,z,tol);
    //ROL::Ptr<Tpetra::MultiVector<Real> > grf = getField(*gr);
    //gf->getVectorNonConst(0)->scale(static_cast<Real>(1),*grf);
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    obj_->gradient_2(g,*ur,z,tol);
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<Tpetra::MultiVector<> > hvf = getField(hv);
    ROL::Ptr<ROL::Vector<Real> > hvr
      = ROL::makePtr<ROL::TpetraMultiVector<Real>(hvf->subViewNonConst(cols>()));
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<ROL::Vector<Real> > vr
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(vf)->subViewNonConst(cols()));
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    obj_->hessVec_11(*hvr,*vr,*ur,z,tol);
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<Tpetra::MultiVector<> > hvf = getField(hv);
    ROL::Ptr<ROL::Vector<Real> > hvr
      = ROL::makePtr<ROL::TpetraMultiVector<Real>(hvf->subViewNonConst(cols>()));
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    obj_->hessVec_12(*hvr,v,*ur,z,tol);
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<const Tpetra::MultiVector<> > vf = getConstField(v);
    ROL::Ptr<ROL::Vector<Real> > vr
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(vf)->subViewNonConst(cols()));
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    obj_->hessVec_21(hv,*vr,*ur,z,tol);
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
             const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    Teuchos::Array<size_t> cols(1,0);
    ROL::Ptr<const Tpetra::MultiVector<> > uf = getConstField(u);
    ROL::Ptr<ROL::Vector<Real> > ur
      = ROL::makePtr<ROL::TpetraMultiVector<Real>>(
          ROL::constPtrCast<Tpetra::MultiVector<> >(uf)->subViewNonConst(cols()));
    obj_->hessVec_22(hv,v,*ur,z,tol);
  }

private: // Vector accessor functions

  ROL::Ptr<const Tpetra::MultiVector<> > getConstField(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
  }

  ROL::Ptr<Tpetra::MultiVector<> > getField(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
  }

}; // class Objective_SimOpt

#endif
