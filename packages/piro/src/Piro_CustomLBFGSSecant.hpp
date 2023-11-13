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

#ifndef PIRO_CUSTOM_LBFGS_SECANT_H
#define PIRO_CUSTOM_LBFGS_SECANT_H

/** \class ROL::CustomLBFGSSecant
    \brief Provides interface for and implements limited-memory secant operators.
*/

#include "ROL_Secant.hpp"

namespace Piro {

template<class Real>
class CustomLBFGSSecant : public ROL::Secant<Real> {
protected:

  Teuchos::RCP<const Thyra::LinearOpBase<Real> > hessianApprox_;
  Teuchos::RCP<const Thyra::LinearOpBase<Real> > invHessianApprox_;

public:

  virtual ~CustomLBFGSSecant() {}

  // Constructor
  CustomLBFGSSecant( const Teuchos::RCP<const Thyra::LinearOpBase<Real> >& hessianApprox,
                const Teuchos::RCP<const Thyra::LinearOpBase<Real> >& invHessianApprox, 
                int M = 10,
                Real scaling = Real(1)) : ROL::Secant<Real>(M, false, scaling),
                hessianApprox_(hessianApprox),
                invHessianApprox_(invHessianApprox) {}
  // Apply Secant Approximate Inverse Hessian
  virtual void applyH( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v ) const {
    const Real zero(0);

    auto tmp = v.clone();
    tmp->set(v); 
    std::vector<Real> alpha(this->state_->current+1,zero);
    for (int i = this->state_->current; i>=0; i--) {
      alpha[i]  = this->state_->iterDiff[i]->apply(*tmp);
      alpha[i] /= this->state_->product[i];
      tmp->axpy(-alpha[i],*this->state_->gradDiff[i]);
    }

    // Apply initial inverse Hessian approximation to v
    applyH0(Hv,*tmp);

    Real beta(0);
    for (int i = 0; i <= this->state_->current; i++) {
      //beta  = Hv.dot((state_->gradDiff[i])->dual());
      beta  = Hv.apply(*this->state_->gradDiff[i]);
      beta /= this->state_->product[i];
      Hv.axpy((alpha[i]-beta),*(this->state_->iterDiff[i]));
    }
  }

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0( ROL::Vector<Real> &Hv, const ROL::Vector<Real> &v ) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
          Teuchos::is_null(invHessianApprox_),
          std::logic_error,
          std::endl <<
          "Piro::CustomSecant::applyH0:  " <<
          "The approximate inverse Hessian is not defined." <<
          std::endl);

    ROL::ThyraVector<Real>  & thyra_Hv = dynamic_cast<ROL::ThyraVector<Real>&>(Hv);
    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    Thyra::apply(*invHessianApprox_, Thyra::NOTRANS, *thyra_v.getVector(), thyra_Hv.getVector().ptr(), static_cast<Real>(1)/this->Bscaling_, static_cast<Real>(0));
  }

  // Apply Secant Approximate Hessian
  virtual void applyB( ROL::Vector<Real> &Bv, const ROL::Vector<Real> &v ) const {
    const Real one(1);

    // Apply initial Hessian approximation to v
    applyB0(Bv,v);

    std::vector<ROL::Ptr<ROL::Vector<Real>>> a(this->state_->current+1);
    std::vector<ROL::Ptr<ROL::Vector<Real>>> b(this->state_->current+1);
    Real bv(0), av(0), bs(0), as(0);
    for (int i = 0; i <= this->state_->current; i++) {
      b[i] = Bv.clone();
      b[i]->set(*(this->state_->gradDiff[i]));
      b[i]->scale(one/sqrt(this->state_->product[i]));
      //bv = v.dot(b[i]->dual());
      bv = v.apply(*b[i]);
      Bv.axpy(bv,*b[i]);

      a[i] = Bv.clone();
      applyB0(*a[i],*(this->state_->iterDiff[i]));

      for (int j = 0; j < i; j++) {
        //bs = (state_->iterDiff[i])->dot(b[j]->dual());
        bs = (this->state_->iterDiff[i])->apply(*b[j]);
        a[i]->axpy(bs,*b[j]);
        //as = (state_->iterDiff[i])->dot(a[j]->dual());
        as = (this->state_->iterDiff[i])->apply(*a[j]);
        a[i]->axpy(-as,*a[j]);
      }
      //as = (state_->iterDiff[i])->dot(a[i]->dual());
      as = (this->state_->iterDiff[i])->apply(*a[i]);
      a[i]->scale(one/sqrt(as));
      //av = v.dot(a[i]->dual());
      av = v.apply(*a[i]);
      Bv.axpy(-av,*a[i]);
    }
  }

  // Apply Initial Secant Approximate Hessian 
  virtual void applyB0( ROL::Vector<Real> &Bv, const ROL::Vector<Real> &v ) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
        Teuchos::is_null(hessianApprox_ ),
        std::logic_error,
        std::endl <<
        "Piro::CustomSecant::apply0:  " <<
        "The approximate Hessian is not defined." <<
        std::endl);

    ROL::ThyraVector<Real>  & thyra_Bv = dynamic_cast<ROL::ThyraVector<Real>&>(Bv);
    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    Thyra::apply(*hessianApprox_, Thyra::NOTRANS, *thyra_v.getVector(), thyra_Bv.getVector().ptr(), this->Bscaling_, static_cast<Real>(0));
  }
};

}

#endif
