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
class Objective_PDEOPT_ElasticitySIMP : public ROL::Objective_SimOpt<Real> {
private:

  const ROL::SharedPointer<ElasticitySIMPOperators<Real> > data_;
  const ROL::SharedPointer<DensityFilter<Real> > filter_;
  Real scale_;
  bool useFU_;
  bool print_;

public:

  Objective_PDEOPT_ElasticitySIMP(const ROL::SharedPointer<ElasticitySIMPOperators<Real> > &data,
                                  const ROL::SharedPointer<DensityFilter<Real> > &filter,
                                  const ROL::SharedPointer<Teuchos::ParameterList> &parlist)
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
      ROL::SharedPointer<const Tpetra::MultiVector<> > zp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
      ROL::SharedPointer<Tpetra::MultiVector<> > Fzp
        = ROL::makeShared<Tpetra::MultiVector<>>(zp->getMap(), 1);
      filter_->apply(Fzp, zp);
      data_->updateMaterialDensity(Fzp);
    }
    if ( (flag && iter >= 0) && print_ ) {
      ROL::SharedPointer<const Tpetra::MultiVector<> > zp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
      std::stringstream filename;
      filename << "density_" << iter << ".txt";
      data_->outputTpetraVector(zp,filename.str());
    }
  }

  Real value(const ROL::Vector<Real> &u,
             const ROL::Vector<Real> &z,
             Real &tol) {
    ROL::SharedPointer<const Tpetra::MultiVector<> > up
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::SharedPointer<const Tpetra::MultiVector<> > zp
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
    
   Teuchos::Array<Real> dotvalU(1, 0);
    if ( useFU_ ) {
      up->dot(*(data_->getVecF()), dotvalU.view(0,1));
    }
    else {
      ROL::SharedPointer<Tpetra::MultiVector<> > matvecp
        = ROL::makeShared<Tpetra::MultiVector<>>(*up, Teuchos::Copy);
      data_->ApplyJacobian1ToVec(matvecp, up);
      up->dot(*matvecp, dotvalU.view(0,1));
    }

    return scale_*dotvalU[0];
  }

  void gradient_1(ROL::Vector<Real> &g,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    ROL::SharedPointer<Tpetra::MultiVector<> > gp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
    if ( useFU_ ) {
      Tpetra::deep_copy(*gp, *(data_->getVecF()));
    }
    else {
      ROL::SharedPointer<const Tpetra::MultiVector<> > up
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
    ROL::SharedPointer<Tpetra::MultiVector<> > gp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(g)).getVector();
    ROL::SharedPointer<const Tpetra::MultiVector<> > zp
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
    if ( useFU_ ) {
      g.zero();
    }
    else {
      ROL::SharedPointer<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::SharedPointer<Tpetra::MultiVector<> > tmp
        = ROL::makeShared<Tpetra::MultiVector<>>(gp->getMap(), 1);
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
      ROL::SharedPointer<Tpetra::MultiVector<> > hvp
        = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
      ROL::SharedPointer<const Tpetra::MultiVector<> > vp
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
      ROL::SharedPointer<Tpetra::MultiVector<> > hvp
        = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
      ROL::SharedPointer<const Tpetra::MultiVector<> > vp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
      ROL::SharedPointer<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::SharedPointer<Tpetra::MultiVector<> > Fvp
        = ROL::makeShared<Tpetra::MultiVector<>>(vp->getMap(), 1);
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
      ROL::SharedPointer<Tpetra::MultiVector<> > hvp
        = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
      ROL::SharedPointer<const Tpetra::MultiVector<> > vp
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
      ROL::SharedPointer<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::SharedPointer<Tpetra::MultiVector<> > tmp
        = ROL::makeShared<Tpetra::MultiVector<>>(hvp->getMap(), 1);
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
    ROL::SharedPointer<Tpetra::MultiVector<> > hvp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(hv)).getVector();
    ROL::SharedPointer<const Tpetra::MultiVector<> > vp
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    if ( useFU_ ) {
      hv.zero();
    }
    else {
      ROL::SharedPointer<const Tpetra::MultiVector<> > up
        = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
      ROL::SharedPointer<Tpetra::MultiVector<> > Fvp
        = ROL::makeShared<Tpetra::MultiVector<>>(vp->getMap(), 1);
      filter_->apply(Fvp, vp);
      ROL::SharedPointer<Tpetra::MultiVector<> > tmp
        = ROL::makeShared<Tpetra::MultiVector<>>(hvp->getMap(), 1);
      data_->ApplyAdjointHessian22ToVec (tmp, up, Fvp, up);
      filter_->apply(hvp, tmp);
    }
  }

};

#endif
