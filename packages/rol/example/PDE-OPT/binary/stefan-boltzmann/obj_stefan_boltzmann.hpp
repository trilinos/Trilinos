// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef BINARY_QOI_STEFAN_BOLTZMANN_HPP
#define BINARY_QOI_STEFAN_BOLTZMANN_HPP

#include "../../TOOLS/qoi.hpp"
#include "../../TOOLS/fe.hpp"

template <class Real>
class QoI_StateCost : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> target_;
  bool truncate_;

public:
  QoI_StateCost(const ROL::Ptr<FE<Real>> &fe,
                Teuchos::ParameterList   &parlist) : fe_(fe) {
    Real T = parlist.sublist("Problem").get("Desired Temperature",1000.0);
    // Nondimensionalize
    bool nondim = parlist.sublist("Problem").get("Nondimensionalize", true);
    if (nondim) {
      Real T0 = parlist.sublist("Problem").get("Reference Temperature",1000.0);
      T /= T0;
    }
    truncate_ = parlist.sublist("Problem").get("Penalize Negative Part",false);

    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    target_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    target_->initialize(T);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>>             & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval;
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute squared L2-norm of diff
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    if (truncate_) {
      const Real zero(0);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          Real val = (*valU_eval)(i,j);
          (*valU_eval)(i,j) = std::min(val,zero);
        }
      }
    }
    fe_->computeIntegral(val,valU_eval,valU_eval);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>>             & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval;
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute gradient of squared L2-norm of diff
    Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
    if (truncate_) {
      const Real zero(0);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          Real val = (*valU_eval)(i,j);
          (*valU_eval)(i,j) = std::min(val,zero);
        }
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valU_eval,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>>             & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval;
    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valV_eval, v_coeff);
    // Compute gradient of squared L2-norm of diff
    if (truncate_) {
      const Real zero(0), one(1);
      ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval;
      valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valU_eval, u_coeff);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          Real val = (*valU_eval)(i,j);
          (*valV_eval)(i,j) *= (val <= zero ? one : zero);
        }
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *fe_->NdetJ(),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

#endif
