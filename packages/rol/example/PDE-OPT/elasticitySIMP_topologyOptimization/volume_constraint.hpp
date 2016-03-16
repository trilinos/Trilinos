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


#ifndef ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_VOLUMN_H
#define ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_VOLUMN_H

#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "Amesos2.hpp"

template<class Real>
class EqualityConstraint_PDEOPT_ElasticitySIMP_Volumn : public ROL::EqualityConstraint_SimOpt<Real> {
private:

  Real volFrac_;

public:

  EqualityConstraint_PDEOPT_ElasticitySIMP_Volumn(const Teuchos::RCP<ElasticitySIMPOperators<Real> > &data,
                                    const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
	
	volFrac_     = this->parlist_->sublist("ElasticityTopoOpt").get("volfrac", 0.5);
  }

  using ROL::EqualityConstraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
    Teuchos::RCP<Tpetra::MultiVector<> > cp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(c)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > zp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(z)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > unit = Teuchos::rcp(new Tpetra::MultiVector<>(zp->getMap(), 1, true));
    unit->PutScalar(1.0);
    Teuchos::Array<Real> sumZ(1, 0);
    zp->dot(*unit, sumZ);
    cp->PutScalar(sumZ[0]);  
  }


  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) 
  {
	jv.zero();
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) 
  {
    Teuchos::RCP<Tpetra::MultiVector<> > jvp = (Teuchos::dyn_cast<ROL::TpetraMultiVector<Real> >(jv)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > vp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(v)).getVector();
    Teuchos::RCP<const Tpetra::MultiVector<> > unit = Teuchos::rcp(new Tpetra::MultiVector<>(vp->getMap(), 1, true));
    unit->PutScalar(1.0);
    Teuchos::Array<Real> sumV(1, 0);
    vp->dot(*unit, sumV);
    jv->PutScalar(sumV[0]);  
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) 
  {
  	ajvp.zero();
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) 
  {
    Teuchos::RCP<const Tpetra::MultiVector<> > ajvp = (Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real> >(ajv)).getVector();
    Teuchos::RCP<std::vector<Real> > vp = (Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector();
    ajvp->PutScalar(1.0);
    ajvp->Scale((*vp)[0]);
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
