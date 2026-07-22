// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  constraint.hpp
    \brief Defines the SimOpt constraint for the 'poisson' example.
*/

#ifndef ROL_PDEOPT_POISSON_CONSTRAINT_H
#define ROL_PDEOPT_POISSON_CONSTRAINT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "Amesos2.hpp"
#include "data.hpp"

template<class Real>
class EqualityConstraint_PDEOPT_Poisson : public ROL::Constraint_SimOpt<Real> {
private:

  ROL::Ptr<PoissonData<Real> > data_;

public:

  EqualityConstraint_PDEOPT_Poisson(const ROL::Ptr<PoissonData<Real> > &data,
                                    const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
    data_ = data;
  }

  using ROL::Constraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > cp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(c)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    Real one(1);
 
    // A*u
    data_->getMatA()->apply(*up, *cp);

    // B*z + A*u
    data_->getMatB()->apply(*zp, *cp, Teuchos::NO_TRANS, one, one);

    // A*u + B*z - f
    cp->update(-one, *(data_->getVecF()), one);
  }


  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > jvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(jv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // A*v
    data_->getMatA()->apply(*vp, *jvp);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > jvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(jv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // B*v
    data_->getMatB()->apply(*vp, *jvp);
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > ajvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(ajv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // A'*v
    bool transpose = true;
    data_->getMatA(transpose)->apply(*vp, *ajvp);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<Tpetra::MultiVector<> > ajvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(ajv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    // B'*v
    bool transpose = true;
    data_->getMatB(transpose)->apply(*vp, *ajvp);
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
    ROL::Ptr<Tpetra::MultiVector<> > ijvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(ijv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();

    data_->getSolver()->setX(ijvp);
    data_->getSolver()->setB(vp);
    data_->getSolver()->solve();

    /*    
    // Construct solver using Amesos2 factory.
    ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver;
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
    ROL::Ptr<Tpetra::MultiVector<> > iajvp =
      (dynamic_cast<ROL::TpetraMultiVector<Real>&>(iajv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp =
      (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    
    bool transpose = true;    
    data_->getSolver(transpose)->setX(iajvp);
    data_->getSolver(transpose)->setB(vp);
    data_->getSolver(transpose)->solve();

    /*
    bool transpose = true;
    // Construct solver using Amesos2 factory.
    ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver;
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
