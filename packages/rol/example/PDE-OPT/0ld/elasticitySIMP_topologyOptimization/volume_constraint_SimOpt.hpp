// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_VOLUME_SimOpt_H
#define ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_VOLUME_SimOpt_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "Amesos2.hpp"

template<class Real>
class EqualityConstraint_PDEOPT_ElasticitySIMP_Volume_SimOpt : public ROL::Constraint_SimOpt<Real> {
private:

  Real volFrac_;

public:

  EqualityConstraint_PDEOPT_ElasticitySIMP_Volume_SimOpt(const ROL::Ptr<ElasticitySIMPOperators<Real> > &data,
                                                         const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
	
	  volFrac_     = parlist->sublist("ElasticityTopoOpt").get("volfrac", 0.5);
  }

  using ROL::Constraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<std::vector<Real> > cp = (dynamic_cast<ROL::StdVector<Real>&>(c)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();
    ROL::Ptr<Tpetra::MultiVector<> > unit = ROL::makePtr<Tpetra::MultiVector<>>(zp->getMap(), 1, true);
    unit->putScalar(1.0);
    Teuchos::Array<Real> sumZ(1, 0);
    zp->dot(*unit, sumZ);
    (*cp)[0] = sumZ[0] - 400*volFrac_;  
  }


  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) 
  {
	jv.zero();
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<std::vector<Real> > jvp = (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    ROL::Ptr<Tpetra::MultiVector<> > unit = ROL::makePtr<Tpetra::MultiVector<>>(vp->getMap(), 1, true);
    unit->putScalar(1.0);
    Teuchos::Array<Real> sumV(1, 0);
    vp->dot(*unit, sumV);
    (*jvp)[0] = sumV[0];  
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) 
  {
  	ajv.zero();
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<Tpetra::MultiVector<> > ajvp = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(ajv)).getVector();
    ROL::Ptr<const std::vector<Real> > vp = (dynamic_cast<const ROL::StdVector<Real>&>(v)).getVector();
    ajvp->putScalar(1.0);
    ajvp->scale((*vp)[0]);
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
	ahwv.zero();
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
	ahwv.zero();
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
	ahwv.zero();
  }


  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
	ahwv.zero();	
  }

};

#endif
