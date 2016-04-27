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
#include "data.hpp"

template<class Real>
class Objective_PDEOPT_ElasticitySIMP : public ROL::Objective_SimOpt<Real> {
private:

  const Teuchos::RCP<ElasticitySIMPOperators<Real> > data_;
  Real alpha_;
  Real scale_;
  bool useFU_;

  Teuchos::RCP<const Tpetra::Vector<> > Zscale_;

public:

  Objective_PDEOPT_ElasticitySIMP(const Teuchos::RCP<ElasticitySIMPOperators<Real> > &data,
                                  const Teuchos::RCP<Teuchos::ParameterList> &parlist)
    : data_(data) {
    alpha_ = parlist->sublist("Problem").get<Real>("Penalty parameter");
    // Compute compliance scaling
    Teuchos::Array<Real> dotF(1, 0);
    (data_->getVecF())->dot(*(data_->getVecF()),dotF);
    Real minDensity = parlist->sublist("ElasticitySIMP").get<Real>("Minimum Density");
    scale_ = minDensity/dotF[0];
    useFU_ = parlist->sublist("ElasticitySIMP").get<bool>("Use Force Dot Displacement Objective");
    // Compute density scaling vector
    bool useZscale = parlist->sublist("Problem").get<bool>("Use Scaled Density Vectors");
    Teuchos::RCP<Tpetra::MultiVector<> > cellMeasures = data->getCellAreas();
    if ( !useZscale ) {
      Real one(1);
      cellMeasures->putScalar(one);
    }
    Zscale_ = cellMeasures->getVector(0);
  }

  void update(const ROL::Vector<Real> &u,
              const ROL::Vector<Real> &z,
              bool flag = true,
              int iter = -1) {
    Teuchos::RCP<const Tpetra::MultiVector<> > zp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    data_->updateMaterialDensity(zp);
  }

  Real value(const ROL::Vector<Real> &u,
             const ROL::Vector<Real> &z,
             Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > up
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
   Teuchos::Array<Real> dotvalU(1, 0);
    if ( useFU_ ) {
      up->dot(*(data_->getVecF()), dotvalU);
    }
    else {
      Teuchos::RCP<Tpetra::MultiVector<> > matvecp
        = Teuchos::rcp(new Tpetra::MultiVector<>(*up, Teuchos::Copy));
      data_->ApplyJacobian1ToVec(matvecp, up);
      up->dot(*matvecp, dotvalU);
    }

    //z z
    Real dotvalZ = z.dot(z);

    Real half(0.5);
    return scale_*dotvalU[0] + half*alpha_*dotvalZ;
  }

  void gradient_1(ROL::Vector<Real> &g,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > gp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
    if ( useFU_ ) {
      Tpetra::deep_copy(*gp, *(data_->getVecF()));
    }
    else {
      Teuchos::RCP<const Tpetra::MultiVector<> > up
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
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
    Teuchos::RCP<Tpetra::MultiVector<> > gp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    if ( useFU_ ) {
      g.zero();
    }
    else {
      Teuchos::RCP<const Tpetra::MultiVector<> > up
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      data_->ApplyAdjointJacobian2ToVec (gp, up, up);
    }

    Tpetra::MultiVector<> Wzp(*zp); Wzp.elementWiseMultiply(1,*Zscale_,*zp,0);
    gp->update(alpha_, Wzp, scale_);
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
      Teuchos::RCP<Tpetra::MultiVector<> > hvp
        = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
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
      Teuchos::RCP<Tpetra::MultiVector<> > hvp
        = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      data_->ApplyJacobian2ToVec (hvp, up, vp);
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
      Teuchos::RCP<Tpetra::MultiVector<> > hvp
        = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > vp
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
      Teuchos::RCP<const Tpetra::MultiVector<> > up
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      data_->ApplyAdjointJacobian2ToVec (hvp, up, vp);
      Real two(2);
      hvp->scale(two*scale_);
    }
  }

  void hessVec_22(ROL::Vector<Real> &hv,
            const ROL::Vector<Real> &v,
            const ROL::Vector<Real> &u,
            const ROL::Vector<Real> &z,
                  Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    if ( useFU_ ) {
      hv.zero();
    }
    else {
      Teuchos::RCP<const Tpetra::MultiVector<> > up
        = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
      data_->ApplyAdjointHessian22ToVec (hvp, up, vp, up);
    }

    Tpetra::MultiVector<> Wvp(*vp); Wvp.elementWiseMultiply(1,*Zscale_,*vp,0);
    hvp->update(alpha_, Wvp, scale_);
  }

};

#endif
