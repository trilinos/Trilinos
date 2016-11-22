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


#ifndef ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_H
#define ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_H

#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "Amesos2.hpp"
#include "filter.hpp"
#include "data.hpp"

template<class Real>
class ParametrizedEqualityConstraint_PDEOPT_ElasticitySIMP : public ROL::EqualityConstraint_SimOpt<Real> {
private:

  const Teuchos::RCP<ElasticitySIMPOperators<Real> > data_;
  const Teuchos::RCP<DensityFilter<Real> > filter_;

public:

  ParametrizedEqualityConstraint_PDEOPT_ElasticitySIMP(const Teuchos::RCP<ElasticitySIMPOperators<Real> > &data,
                                                       const Teuchos::RCP<DensityFilter<Real> > &filter,
                                                       const Teuchos::RCP<Teuchos::ParameterList> &parlist)
    : data_(data), filter_(filter) {}

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  
  void value(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &u,
       const ROL::Vector<Real> &z,
             Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > cp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(c)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    
    data_->ApplyJacobian1ToVec(cp, up);
    Real one(1);
    data_->updateF(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    cp->update(-one, *(data_->getVecF()), one);
  }

  void solve(ROL::Vector<Real> &c,
             ROL::Vector<Real> &u,
       const ROL::Vector<Real> &z,
             Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > up
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(u)).getVector();
    // Solve PDE    
    data_->updateF(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    data_->ApplyInverseJacobian1ToVec(up, data_->getVecF(), false);
    // Compute residual
    ParametrizedEqualityConstraint_PDEOPT_ElasticitySIMP<Real>::value(c,u,z,tol);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z,
                       Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    
    data_->ApplyJacobian1ToVec(jvp, vp);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv,
                 const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u,
                 const ROL::Vector<Real> &z,
                       Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();

    Teuchos::RCP<Tpetra::MultiVector<> > Fvp
      = Teuchos::rcp(new Tpetra::MultiVector<>(vp->getMap(), 1));
    filter_->apply(Fvp, vp);
    data_->ApplyJacobian2ToVec (jvp, up, Fvp);
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ajvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    
    data_->ApplyJacobian1ToVec(ajvp, vp);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ajvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();

    Teuchos::RCP<Tpetra::MultiVector<> > tmp
      = Teuchos::rcp(new Tpetra::MultiVector<>(ajvp->getMap(), 1));
    data_->ApplyAdjointJacobian2ToVec (tmp, up, vp);
    filter_->apply(ajvp, tmp);
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    applyAdjointJacobian_2(ahwv, w, v, z, tol);
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    applyJacobian_2(ahwv, v, w, z, tol);
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
                        const ROL::Vector<Real> &w,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ahwvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ahwv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > wp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(w)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
    
    Teuchos::RCP<Tpetra::MultiVector<> > Fvp
      = Teuchos::rcp(new Tpetra::MultiVector<>(vp->getMap(), 1));
    filter_->apply(Fvp, vp);
    Teuchos::RCP<Tpetra::MultiVector<> > tmp
      = Teuchos::rcp(new Tpetra::MultiVector<>(ahwvp->getMap(), 1));
    data_->ApplyAdjointHessian22ToVec (tmp, up, Fvp, wp);
    filter_->apply(ahwvp, tmp);
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,
                        const ROL::Vector<Real> &z,
                              Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ijvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ijv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    
    data_->ApplyInverseJacobian1ToVec (ijvp, vp, false);
  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &u,
                               const ROL::Vector<Real> &z,
                                     Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > iajvp
      = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(iajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    
    data_->ApplyInverseJacobian1ToVec (iajvp, vp, true);
  }

  void update_2(const ROL::Vector<Real> &z, bool flag = true, int iter = -1) {
    Teuchos::RCP<const Tpetra::MultiVector<> > zp
      = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    Teuchos::RCP<Tpetra::MultiVector<> > Fzp
      = Teuchos::rcp(new Tpetra::MultiVector<>(zp->getMap(), 1));
    filter_->apply(Fzp, zp);
    data_->updateMaterialDensity(Fzp);
    data_->constructSolverWithNewMaterial();
  }

};

#endif
