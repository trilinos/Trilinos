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

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_QOI_NONLINEAR_ELLIPTIC_HPP
#define PDEOPT_QOI_NONLINEAR_ELLIPTIC_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_nonlinear_elliptic.hpp"

template <class Real>
class QoI_StateTracking_Nonlinear_Elliptic : public QoI<Real> {
private:
  ROL::Ptr<FE<Real> > fe_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > target_;

protected:
  void computeTracking(void) {
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    target_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_)(i,j) = targetFunc(pt);
      }
    }
  }

public:
  QoI_StateTracking_Nonlinear_Elliptic(const ROL::Ptr<FE<Real> > &fe) : fe_(fe) {
    computeTracking();
  }

  virtual Real targetFunc(const std::vector<Real> & x) const {
    const int size = x.size();
    const Real zero(0), one(1);
    Real val(0);
    for (int i = 0; i < size; ++i) {
      val = std::max(val, std::min(x[i] - zero, one - x[i]));
    }
    return val;
    //const int size = x.size();
    //bool flag = true;
    //for (int i = 0; i < size; ++i) {
    //  flag *= (x[i] < static_cast<Real>(0.6))*(x[i] > static_cast<Real>(0.4));
    //}
    //return (flag ? static_cast<Real>(1) : static_cast<Real>(0));
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valU_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int f = fe_->N()->dimension(1);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real>
class QoI_ControlPenalty_Nonlinear_Elliptic : public QoI<Real> {
private:
  ROL::Ptr<FE<Real> > fe_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > ones_;

public:
  QoI_ControlPenalty_Nonlinear_Elliptic(const ROL::Ptr<FE<Real> > &fe) : fe_(fe) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Build field container of ones
    ones_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*ones_)(i,j) = static_cast<Real>(1);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Build local state tracking term
    ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    fe_->computeIntegral(val,ones_,valZ_eval);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *ones_,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_22 is zero.");
  }

}; // QoI_L2Penalty

#endif
