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

#ifndef ROL_PDEOPT_STOCH_ADV_DIFF_OBJECTIVE_H
#define ROL_PDEOPT_STOCH_ADV_DIFF_OBJECTIVE_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "data.hpp"

template<class Real>
class Objective_PDEOPT_Poisson : public ROL::Objective_SimOpt<Real> {
private:

  const ROL::Ptr<PoissonData<Real> > data_;
  Real cost_control_;
  Real cost_state_;

public:

  Objective_PDEOPT_Poisson(const ROL::Ptr<PoissonData<Real> > &data,
                           const Teuchos::RCP<Teuchos::ParameterList> &parlist)
      : data_(data) {
    cost_control_ = parlist->sublist("Problem").get("Control Cost", 1e0);
    cost_state_   = parlist->sublist("Problem").get("State Cost", 1e5);
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<> > up =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    //ROL::Ptr<const Tpetra::MultiVector<> > zp =
    //  (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
    ROL::Ptr<const std::vector<Real> > zp =
      (dynamic_cast<const ROL::StdVector<Real>&>(z)).getVector();

    Teuchos::Array<Real> dotvalU(1, 0);
    //Teuchos::Array<Real> dotvalZ(1, 0);

    // Set difference vector diffp to up.
    ROL::Ptr<Tpetra::MultiVector<> > diffp =
      ROL::makePtr<Tpetra::MultiVector<>>(*up, Teuchos::Copy);
    // Temporary matvec vector.
    ROL::Ptr<Tpetra::MultiVector<> > matvecp =
      ROL::makePtr<Tpetra::MultiVector<>>(*up, Teuchos::Copy);

    // (u-ud)
    Real one(1);
    diffp->update(-one, *(data_->getVecUd()), one);
    // M*(u-ud)
    data_->getMatM()->apply(*diffp, *matvecp); // mass matrix product
    // (u-ud)'*M*(u-ud)
    diffp->dot(*matvecp, dotvalU);

    // R*z
//    data_->getMatR()->apply(*zp, *matvecp);
    // z'*R*z
//    zp->dot(*matvecp, dotvalZ);
    Real sumz(0);
    int nz = static_cast<int>(zp->size());
    for (int i = 0; i < nz; ++i) {
      sumz += (*zp)[i];
    }

    // 1/2 * (u-ud)'*M*(u-ud) + alpha/2 * z'*R*z
    Real half(0.5);
//    return half*(dotvalU[0] + alpha_*dotvalZ[0]);
    return half*cost_state_*dotvalU[0] + cost_control_*sumz;
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
    Real one(1);
    diffp->update(-one, *(data_->getVecUd()), one);
    // M*(u-ud)
    data_->getMatM()->apply(*diffp, *gp); // mass matrix product
    gp->scale(cost_state_);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > gp =
      (dynamic_cast<ROL::StdVector<Real>&>(g)).getVector();
    gp->assign(gp->size(),cost_control_);
//    ROL::Ptr<Tpetra::MultiVector<> > gp =
//      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
//    ROL::Ptr<const Tpetra::MultiVector<> > zp =
//      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
//
//    // alpha * R*z
//    data_->getMatR()->apply(*zp, *gp);
//    gp->scale(alpha_);
  }

  void hessVec_11(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > hvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // M*v
    data_->getMatM()->apply(*vp, *hvp); // mass matrix product
    hvp->scale(cost_state_);
  }

  void hessVec_12(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    hv.zero();
  }

  void hessVec_21(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    hv.zero();
  }

  void hessVec_22(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    hv.zero();
//    ROL::Ptr<Tpetra::MultiVector<> > hvp =
//      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
//    ROL::Ptr<const Tpetra::MultiVector<> > vp =
//      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
//
//    // alpha * R*v
//    data_->getMatR()->apply(*vp, *hvp);
//    hvp->scale(alpha_);
  }

};

#endif
