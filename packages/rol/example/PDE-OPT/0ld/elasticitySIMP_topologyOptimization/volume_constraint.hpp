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


#ifndef ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_VOLUME_H
#define ROL_PDEOPT_ELASTICITYSIMP_CONSTRAINT_VOLUME_H

#include "ROL_Constraint.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RiskVector.hpp"

template<class Real>
class EqualityConstraint_PDEOPT_ElasticitySIMP_Volume : public ROL::Constraint<Real> {
private:
  Real volFrac_;
  Real totalMeasure_;
  ROL::Ptr<Tpetra::MultiVector<> > cellMeasures_;

  ROL::Vector<Real> & castRiskVector(ROL::Vector<Real> &x) const {
    try {
      ROL::RiskVector<Real> & rx = dynamic_cast<ROL::RiskVector<Real>&>(x);
      return *(rx.getVector());
    }
    catch (std::exception &e) {
      return x;
    }
  }

  const ROL::Vector<Real> & castConstRiskVector(const ROL::Vector<Real> &x) const {
    try {
      const ROL::RiskVector<Real> & rx = dynamic_cast<const ROL::RiskVector<Real>&>(x);
      return *(rx.getVector());
    }
    catch (std::exception &e) {
      return x;
    }
  }

public:
  EqualityConstraint_PDEOPT_ElasticitySIMP_Volume(const ROL::Ptr<ElasticitySIMPOperators<Real> > &data,
                                                  const Teuchos::RCP<Teuchos::ParameterList> &parlist) {
    volFrac_ = parlist->sublist("ElasticityTopoOpt").get<Real>("Volume Fraction");
    cellMeasures_= data->getCellAreas();
    Teuchos::Array<Real> sumM(1, 0);
    cellMeasures_->norm1(sumM);
    totalMeasure_ = sumM[0];
//    std::cout << "Volume Constraints: volfrac_= " << volFrac_ 
//              << ", totalMeasure_= "              << totalMeasure_
//              << std::endl;
  }

  using ROL::Constraint<Real>::value;
  void value(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &z,
             Real &tol) {
    ROL::Ptr<std::vector<Real> > cp = (dynamic_cast<ROL::StdVector<Real>&>(c)).getVector();
    const ROL::Vector<Real> & zr = castConstRiskVector(z);
    ROL::Ptr<const Tpetra::MultiVector<> > zp
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(zr)).getVector();

//    ROL::Ptr<Tpetra::MultiVector<> > unit
//      = ROL::makePtr<Tpetra::MultiVector<>>(zp->getMap(), 1, true);
//    unit->putScalar(1.0);
//    Teuchos::Array<Real> sumZ(1, 0);
//    zp->dot(*unit, sumZ);
//    (*cp)[0] = sumZ[0] - static_cast<Real>(zp->getGlobalLength())*volFrac_;

    Teuchos::Array<Real> sumZ(1, 0);
    cellMeasures_->dot(*zp, sumZ);
    (*cp)[0] = sumZ[0] - totalMeasure_ * volFrac_;
//    std::cout << sumZ[0]       << ", "
//              << totalMeasure_ << ", "
//              << volFrac_      << std::endl;
  }

  void applyJacobian(ROL::Vector<Real> &jv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &z,
                     Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp = (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    const ROL::Vector<Real> & vr = castConstRiskVector(v);
    ROL::Ptr<const Tpetra::MultiVector<> > vp
      = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(vr)).getVector();

//    ROL::Ptr<Tpetra::MultiVector<> > unit
//      = ROL::makePtr<Tpetra::MultiVector<>>(vp->getMap(), 1, true);
//    unit->putScalar(1.0);
//    Teuchos::Array<Real> sumV(1, 0);
//    vp->dot(*unit, sumV);
//    (*jvp)[0] = sumV[0];

    Teuchos::Array<Real> sumV(1, 0);
    cellMeasures_->dot(*vp, sumV);
    (*jvp)[0] = sumV[0];
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                      const ROL::Vector<Real> &v,
                      const ROL::Vector<Real> &z,
                            Real &tol) {
    ROL::Vector<Real> & ajvr = castRiskVector(ajv);
    ROL::Ptr<Tpetra::MultiVector<> > ajvp
      = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(ajvr)).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = (dynamic_cast<const ROL::StdVector<Real>&>(v)).getVector();

//    ajvp->putScalar(1.0);
//    ajvp->scale((*vp)[0]);

    Tpetra::deep_copy(*ajvp, *cellMeasures_);	
    ajvp->scale((*vp)[0]);
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                     const ROL::Vector<Real> &w,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &z,
                           Real &tol) {
    ahwv.zero();
  }

};

#endif
