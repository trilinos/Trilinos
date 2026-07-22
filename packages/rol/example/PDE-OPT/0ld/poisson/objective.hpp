// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  objective.hpp
    \brief Defines the SimOpt objective function for the 'poisson' example.
*/

#ifndef ROL_PDEOPT_POISSON_OBJECTIVE_H
#define ROL_PDEOPT_POISSON_OBJECTIVE_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "data.hpp"

template<class Real>
class Objective_PDEOPT_Poisson : public ROL::Objective_SimOpt<Real> {
private:

  ROL::Ptr<PoissonData<Real> > data_;
  Real alpha_;

public:

  Objective_PDEOPT_Poisson(const ROL::Ptr<PoissonData<Real> > &data,
                           const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
    data_ = data;
    alpha_ = parlist->sublist("Problem").get("Control penalty parameter", 1e-2);
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<> > up =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    Teuchos::Array<Real> dotvalU(1, 0);
    Teuchos::Array<Real> dotvalZ(1, 0);

    // Set difference vector diffp to up.
    ROL::Ptr<Tpetra::MultiVector<> > diffp =
      ROL::makePtr<Tpetra::MultiVector<>>(*up, Teuchos::Copy);
    // Temporary matvec vector.
    ROL::Ptr<Tpetra::MultiVector<> > matvecp =
      ROL::makePtr<Tpetra::MultiVector<>>(*up, Teuchos::Copy);

    // (u-ud)
    diffp->update(-1.0, *(data_->getVecUd()), 1.0);
    // M*(u-ud)
    data_->getMatM()->apply(*diffp, *matvecp);
    // (u-ud)'*M*(u-ud)
    diffp->dot(*matvecp, dotvalU);

    // R*z
    data_->getMatR()->apply(*zp, *matvecp);
    // z'*R*z
    zp->dot(*matvecp, dotvalZ);

    // 1/2 * (u-ud)'*M*(u-ud) + alpha/2 * z'*R*z
    return(0.5*dotvalU[0] + 0.5*alpha_*dotvalZ[0]);
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > gp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();

    // Set difference vector diffp to up.
    ROL::Ptr<Tpetra::MultiVector<> > diffp =
      ROL::makePtr<Tpetra::MultiVector<>>(*up, Teuchos::Copy);
    // (u-ud)
    diffp->update(-1.0, *(data_->getVecUd()), 1.0);
    // M*(u-ud)
    data_->getMatM()->apply(*diffp, *gp);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > gp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    // alpha * R*z
    data_->getMatR()->apply(*zp, *gp);
    gp->scale(alpha_);
  }

  void hessVec_11(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > hvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // M*v
    data_->getMatM()->apply(*vp, *hvp);
  }

  void hessVec_12(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > hvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();

    // zero
    hvp->scale(0);
  }

  void hessVec_21(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > hvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();

    // zero
    hvp->scale(0);
  }

  void hessVec_22(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > hvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // alpha * R*v
    data_->getMatR()->apply(*vp, *hvp);
    hvp->scale(alpha_);
  }

};

#endif
