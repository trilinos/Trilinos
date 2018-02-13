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


#ifndef ROL_CoupledConstraint_SimOpt_H
#define ROL_CoupledConstraint_SimOpt_H

#include "ROL_Constraint_SimOpt.hpp"

#include "constraint.hpp"
#include "volume_constraint_SimOpt.hpp"

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"
#include "Amesos2.hpp"

template<class Real>
class EqualityConstraint_PDEOPT_ElasticitySIMP_Coupled : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<EqualityConstraint_PDEOPT_ElasticitySIMP<Real> > pde_;
  const ROL::Ptr<EqualityConstraint_PDEOPT_ElasticitySIMP_Volume_SimOpt<Real> > vol_;

public:
  EqualityConstraint_PDEOPT_ElasticitySIMP_Coupled(const ROL::Ptr<EqualityConstraint_PDEOPT_ElasticitySIMP<Real> > &pde,
                                                   const ROL::Ptr<EqualityConstraint_PDEOPT_ElasticitySIMP_Volume_SimOpt<Real> > &vol) : pde_(pde), vol_(vol) {}

  using ROL::Constraint_SimOpt<Real>::value;
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<std::vector<Real> > cp = (dynamic_cast<ROL::StdVector<Real>&>(c)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    pde_->value(*cp1,*up,*zp,tol);
    vol_->value(*cp2,*up,*zp,tol);
  }


  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<std::vector<Real> > jvp = (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    pde_->applyJacobian_1(*jvp1,*vp,*up,*zp,tol);
    vol_->applyJacobian_1(*jvp2,*vp,*up,*zp,tol);
  }


  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<std::vector<Real> > jvp = (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    pde_->applyJacobian_2(*jvp1,*vp,*up,*zp,tol);
    vol_->applyJacobian_2(*jvp2,*vp,*up,*zp,tol);
  }


  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<std::vector<Real> > ajvp = (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > vp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(v)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    pde_->applyAdjointJacobian_1(*ajvp1,*vp,*up,*zp,tol);
    vol_->applyAdjointJacobian_1(*ajvp2,*vp,*up,*zp,tol);
  }


  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) 
  {
    ROL::Ptr<Tpetra::MultiVector<> > ajvp = (dynamic_cast<ROL::TpetraMultiVector<Real>&>(ajv)).getVector();
    ROL::Ptr<const std::vector<Real> > vp = (dynamic_cast<const ROL::StdVector<Real>&>(v)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > up = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(u)).getVector();
    ROL::Ptr<const Tpetra::MultiVector<> > zp = (dynamic_cast<const ROL::TpetraMultiVector<Real>&>(z)).getVector();

    pde_->applyAdjointJacobian_2(*ajvp1,*vp,*up,*zp,tol);
    vol_->applyAdjointJacobian_2(*ajvp2,*vp,*up,*zp,tol);
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
