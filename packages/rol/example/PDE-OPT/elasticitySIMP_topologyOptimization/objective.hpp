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
#include "ROL_TpetraMultiVector.hpp"
#include "data.hpp"

template<class Real>
class Objective_PDEOPT_ElasticitySIMP : public ROL::Objective_SimOpt<Real> {
private:

  Teuchos::RCP<ElasticitySIMPOperators<Real> > data_;
  Real alpha_;

public:

  Objective_PDEOPT_ElasticitySIMP(const Teuchos::RCP<ElasticitySIMPOperators<Real> > &data,
                           const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
    data_ = data;
    alpha_ = parlist->sublist("Problem").get<Real>("Penalty parameter");
  }

  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > up = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    
    // K u
    Teuchos::RCP<Tpetra::MultiVector<> > matvecp = Teuchos::rcp(new Tpetra::MultiVector<>(*up, Teuchos::Copy));
    data_->ApplyJacobian1ToVec(matvecp, up);
    
    Teuchos::Array<Real> dotvalU(1, 0);
    Teuchos::Array<Real> dotvalZ(1, 0);
    //u K u
    up->dot(*matvecp, dotvalU);
    //z z
    zp->dot(*zp, dotvalZ);

    return(0.5*dotvalU[0] + 0.5*alpha_*dotvalZ[0]);
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > gp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    data_->ApplyJacobian1ToVec(gp, up);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > gp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    data_->ApplyAdjointJacobian2ToVec (gp, up, up);
    gp->update(alpha_, *zp, 0.5);
  }

  void hessVec_11(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    data_->ApplyJacobian1ToVec(hvp, vp);
  }

  void hessVec_12(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    data_->ApplyJacobian2ToVec (hvp, up, vp);
  }

  void hessVec_21(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    data_->ApplyAdjointJacobian2ToVec (hvp, up, vp);
  }

  void hessVec_22(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    
    data_->updateMaterialDensity (zp);
    data_->ApplyAdjointHessian22ToVec (hvp, up, vp, up);
    hvp->update(alpha_, *vp, 0.5);
  }

};

#endif
