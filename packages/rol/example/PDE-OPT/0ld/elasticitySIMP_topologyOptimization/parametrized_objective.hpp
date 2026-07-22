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

#ifndef ROL_PDEOPT_ELASTICITYSIMP_OBJECTIVE_H
#define ROL_PDEOPT_ELASTICITYSIMP_OBJECTIVE_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_ScaledTpetraMultiVector.hpp"
#include "filter.hpp"
#include "data.hpp"

template<class Real>
class ParametrizedObjective_PDEOPT_ElasticitySIMP : public ROL::Objective_SimOpt<Real> {
private:

  const ROL::Ptr<ElasticitySIMPOperators<Real> > data_;
  const ROL::Ptr<DensityFilter<Real> > filter_;
  Real scale_;
  bool useFU_;
  bool print_;

public:

  ParametrizedObjective_PDEOPT_ElasticitySIMP(const ROL::Ptr<ElasticitySIMPOperators<Real> > &data,
                                              const ROL::Ptr<DensityFilter<Real> > &filter,
                                              const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                                              const Real scale)
    : data_(data), filter_(filter), scale_(scale) {
    useFU_ = parlist->sublist("ElasticitySIMP").get("Use Force Dot Displacement Objective",true);
    print_ = parlist->sublist("ElasticitySIMP").get("Print Density",false);
  }

  ParametrizedObjective_PDEOPT_ElasticitySIMP(const ROL::Ptr<ElasticitySIMPOperators<Real> > &data,
                                              const ROL::Ptr<DensityFilter<Real> > &filter,
                                              const Teuchos::RCP<Teuchos::ParameterList> &parlist)
    : data_(data), filter_(filter) {
    // Compute compliance scaling
    Teuchos::Array<Real> dotF(1, 0);
    (data_->getVecF())->dot(*(data_->getVecF()),dotF);
    Real minDensity = parlist->sublist("ElasticitySIMP").get<Real>("Minimum Density");
    scale_ = minDensity/dotF[0];
    useFU_ = parlist->sublist("ElasticitySIMP").get("Use Force Dot Displacement Objective",true);
    print_ = parlist->sublist("ElasticitySIMP").get("Print Density",false);
  }

  void update(const ROL::Vector<Real> &u,
              const ROL::Vector<Real> &z,
              bool flag = true,
              int iter = -1) {
    if ( !useFU_ ) {
      ROL::Ptr<const Tpetra::MultiVector<> > zp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
      ROL::Ptr<Tpetra::MultiVector<> > Fzp
        = ROL::makePtr<Tpetra::MultiVector<>>(zp->getMap(), 1);
      filter_->apply(Fzp, zp);
      data_->updateMaterialDensity(Fzp);
    }
    if ( (flag && iter >= 0) && print_ ) {
      ROL::Ptr<const Tpetra::MultiVector<> > zp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
      std::stringstream filename;
      filename << "density_" << iter << ".txt";
      data_->outputTpetraVector(zp,filename.str());
    }
  }

  Real value(const ROL::Vector<Real> &u,
             const ROL::Vector<Real> &z,
             Real &tol) {
    ROL::Ptr<const Tpetra::MultiVector<> > up
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    
   Teuchos::Array<Real> dotvalU(1, 0);
    if ( useFU_ ) {
      data_->updateF(ROL::Objective_SimOpt<Real>::getParameter());
      up->dot(*(data_->getVecF()), dotvalU.view(0,1));
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> > matvecp
        = ROL::makePtr<Tpetra::MultiVector<>>(*up, Teuchos::Copy);
      data_->ApplyJacobian1ToVec(matvecp, up);
      up->dot(*matvecp, dotvalU.view(0,1));
    }

    return scale_*dotvalU[0];
  }

  void gradient_1(ROL::Vector<Real> &g,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > gp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
    if ( useFU_ ) {
      data_->updateF(ROL::Objective_SimOpt<Real>::getParameter());
      Tpetra::deep_copy(*gp, *(data_->getVecF()));
    }
    else {
      ROL::Ptr<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      data_->ApplyJacobian1ToVec(gp, up);
      Real two(2);
      gp->scale(two);
    }
    gp->scale(scale_);
  }

  void gradient_2(ROL::Vector<Real> &g,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > gp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
    if ( useFU_ ) {
      g.zero();
    }
    else {
      ROL::Ptr<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::Ptr<Tpetra::MultiVector<> > tmp
        = ROL::makePtr<Tpetra::MultiVector<>>(gp->getMap(), 1);
      data_->ApplyAdjointJacobian2ToVec (tmp, up, up);
      filter_->apply(gp, tmp);
    }
  }

  void hessVec_11(ROL::Vector<Real> &hv,
            const ROL::Vector<Real> &v,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    if ( useFU_ ) {
      hv.zero();
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> > hvp
        = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
      ROL::Ptr<const Tpetra::MultiVector<> > vp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
      data_->ApplyJacobian1ToVec(hvp, vp);
      Real two(2);
      hvp->scale(two*scale_);
    }
  }

  void hessVec_12(ROL::Vector<Real> &hv,
            const ROL::Vector<Real> &v,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    if ( useFU_ ) {
      hv.zero();
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> > hvp
        = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
      ROL::Ptr<const Tpetra::MultiVector<> > vp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
      ROL::Ptr<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::Ptr<Tpetra::MultiVector<> > Fvp
        = ROL::makePtr<Tpetra::MultiVector<>>(vp->getMap(), 1);
      filter_->apply(Fvp, vp);
      data_->ApplyJacobian2ToVec (hvp, up, Fvp);
      Real two(2);
      hvp->scale(two*scale_);
    }
  }

  void hessVec_21(ROL::Vector<Real> &hv,
            const ROL::Vector<Real> &v,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    if ( useFU_ ) {
      hv.zero();
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<> > hvp
        = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
      ROL::Ptr<const Tpetra::MultiVector<> > vp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
      ROL::Ptr<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::Ptr<Tpetra::MultiVector<> > tmp
        = ROL::makePtr<Tpetra::MultiVector<>>(hvp->getMap(), 1);
      data_->ApplyAdjointJacobian2ToVec (tmp, up, vp);
      filter_->apply(hvp, tmp);
      Real two(2);
      hvp->scale(two*scale_);
    }
  }

  void hessVec_22(ROL::Vector<Real> &hv,
            const ROL::Vector<Real> &v,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > hvp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    if ( useFU_ ) {
      hv.zero();
    }
    else {
      ROL::Ptr<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::Ptr<Tpetra::MultiVector<> > Fvp
        = ROL::makePtr<Tpetra::MultiVector<>>(vp->getMap(), 1);
      filter_->apply(Fvp, vp);
      ROL::Ptr<Tpetra::MultiVector<> > tmp
        = ROL::makePtr<Tpetra::MultiVector<>>(hvp->getMap(), 1);
      data_->ApplyAdjointHessian22ToVec (tmp, up, Fvp, up);
      filter_->apply(hvp, tmp);
    }
  }

};

#endif
