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

/*! \file  constraint.hpp
    \brief Defines the SimOpt constraint for the 'poisson' example.
*/

#ifndef ROL_PDEOPT_STOCH_ADV_DIFF_CONSTRAINT_H
#define ROL_PDEOPT_STOCH_ADV_DIFF_CONSTRAINT_H

#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "Amesos2.hpp"
#include "data.hpp"

template<class Real>
class EqualityConstraint_PDEOPT_Poisson : public ROL::EqualityConstraint_SimOpt<Real> {
private:

  const Teuchos::RCP<PoissonData<Real> > data_;

public:

  EqualityConstraint_PDEOPT_Poisson(const Teuchos::RCP<PoissonData<Real> > &data,
                                    const Teuchos::RCP<Teuchos::ParameterList> &parlist)
    : data_(data) {}

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > cp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(c)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > up =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(u)).getVector();
//    Teuchos::RCP<const Tpetra::MultiVector<> > zp =
//      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(z)).getVector();


    Real one(1);

    data_->updateA(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    data_->updateF(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
 
    // A*u
    data_->getMatA()->apply(*up, *cp);

    // B*z + A*u
    data_->applyMatB(cp,*zp,true);
    //data_->getMatB()->apply(*zp, *cp, Teuchos::NO_TRANS, one, one);

    // A*u + B*z - f
    cp->update(-one, *(data_->getVecF()), one);
  }


  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();

    // A*v
    data_->updateA(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    data_->getMatA()->apply(*vp, *jvp);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
//    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
//      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const std::vector<Real> > vp =
      (Teuchos::dyn_cast<const ROL::StdVector<Real> >(v)).getVector();

    // B*v
    data_->applyMatB(jvp,*vp,false);
//    data_->getMatB()->apply(*vp, *jvp);
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ajvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();

    // A'*v
    bool transpose = true;
    data_->updateA(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    data_->getMatA(transpose)->apply(*vp, *ajvp);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
//    Teuchos::RCP<Tpetra::MultiVector<> > ajvp =
//      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<std::vector<Real> > ajvp =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(ajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();

    // B'*v
    data_->applyMatBtranspose(ajvp,*vp);
//    bool transpose = true;
//    data_->getMatB(transpose)->apply(*vp, *ajvp);
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > ijvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(ijv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();

    data_->updateA(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    data_->getSolver()->setX(ijvp);
    data_->getSolver()->setB(vp);
    data_->getSolver()->solve();

    /*    
    // Construct solver using Amesos2 factory.
    Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver;
    try{
      solver = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", data_->getMatA(), ijvp, vp);
    } catch (std::invalid_argument e) {
      std::cout << e.what() << std::endl;
    }
    solver->numericFactorization();
    solver->solve();
    */

  }


  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                                     const ROL::Vector<Real> &z, Real &tol) {
    Teuchos::RCP<Tpetra::MultiVector<> > iajvp =
      (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(iajv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp =
      (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    
    bool transpose = true;    
    data_->updateA(ROL::EqualityConstraint_SimOpt<Real>::getParameter());
    data_->getSolver(transpose)->setX(iajvp);
    data_->getSolver(transpose)->setB(vp);
    data_->getSolver(transpose)->solve();

    /*
    bool transpose = true;
    // Construct solver using Amesos2 factory.
    Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver;
    try{
      solver = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", data_->getMatA(transpose), iajvp, vp);
    } catch (std::invalid_argument e) {
      std::cout << e.what() << std::endl;
    }
    solver->numericFactorization();
    solver->solve();
    */

  }

};

#endif
