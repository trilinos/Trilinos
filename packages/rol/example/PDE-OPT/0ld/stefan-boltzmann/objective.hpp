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
    \brief Defines the SimOpt objective function for the 'stefan-boltzmann' example.
*/

#ifndef ROL_PDEOPT_STEFANBOLTZMANN_OBJECTIVE_H
#define ROL_PDEOPT_STEFANBOLTZMANN_OBJECTIVE_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "data.hpp"

template<class Real>
class Objective_PDEOPT_StefanBoltzmann : public ROL::Objective_SimOpt<Real> {
private:

  Teuchos::RCP<StefanBoltzmannData<Real> > data_;
  Real alpha_;

public:

  Objective_PDEOPT_StefanBoltzmann(const Teuchos::RCP<StefanBoltzmannData<Real> > &data,
                                   const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
    data_ = data;
    alpha_ = parlist->sublist("Problem").get("Penalty parameter", 1e-2);
  }

  using ROL::Objective_SimOpt<Real>::value;
  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

    Teuchos::Array<Real> dotvalU(1, 0);
    Teuchos::Array<Real> dotvalZ(1, 0);

    // Set difference vector diffp to up.
    Teuchos::RCP<Tpetra::MultiVector<> > diffp =
      Teuchos::rcp(new Tpetra::MultiVector<>(*up, Teuchos::Copy));
    // Temporary matvec vector.
    Teuchos::RCP<Tpetra::MultiVector<> > matvecp =
      Teuchos::rcp(new Tpetra::MultiVector<>(*up, Teuchos::Copy));

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
    Teuchos::RCP<Tpetra::MultiVector<> > gp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();

    // Set difference vector diffp to up.
    Teuchos::RCP<Tpetra::MultiVector<> > diffp =
      Teuchos::rcp(new Tpetra::MultiVector<>(*up, Teuchos::Copy));
    // (u-ud)
    diffp->update(-1.0, *(data_->getVecUd()), 1.0);
    // M*(u-ud)
    data_->getMatM()->apply(*diffp, *gp);
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > gp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(g)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();

    // alpha * R*z
    data_->getMatR()->apply(*zp, *gp);
    gp->scale(alpha_);
  }

  void hessVec_11(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();

    // M*v
    data_->getMatM()->apply(*vp, *hvp);
  }

  void hessVec_12(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();

    // zero
    hvp->scale(0);
  }

  void hessVec_21(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();

    // zero
    hvp->scale(0);
  }

  void hessVec_22(ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
                  const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > hvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(hv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();

    // alpha * R*v
    data_->getMatR()->apply(*vp, *hvp);
    hvp->scale(alpha_);
  }

};

#endif
