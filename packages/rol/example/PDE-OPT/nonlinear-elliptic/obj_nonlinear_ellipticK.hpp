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

#ifndef PDEOPT_QOI_NONLINEAR_ELLIPTICK_HPP
#define PDEOPT_QOI_NONLINEAR_ELLIPTICK_HPP

#include "../TOOLS/qoiK.hpp"
#include "pde_nonlinear_ellipticK.hpp"

template <class Real, class DeviceType>
class QoI_StateTracking_Nonlinear_Elliptic : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;
  scalar_view target_;

protected:
  void computeTracking(void) {
    int c = fe_->cubPts().extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int d = fe_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    target_ = scalar_view("target_", c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          pt[k] = (fe_->cubPts())(i,j,k);
        target_(i,j) = targetFunc(pt);
      }
    }
  }

public:
  QoI_StateTracking_Nonlinear_Elliptic(const ROL::Ptr<fe_type> &fe) : fe_(fe) {
    computeTracking();
  }

  virtual Real targetFunc(const std::vector<Real> & x) const {
    const Real zero(0), one(1);
    Real val(0);
    for (auto xi : x)
      val = std::max(val, std::min(xi-zero, one-xi));
    return val;
    //const int size = x.size();
    //bool flag = true;
    //for (int i = 0; i < size; ++i) {
    //  flag *= (x[i] < static_cast<Real>(0.6))*(x[i] > static_cast<Real>(0.4));
    //}
    //return (flag ? static_cast<Real>(1) : static_cast<Real>(0));
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval,target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval,target_);
    // Compute gradient of squared L2-norm of diff
    fst::integrate(grad,valU_eval,fe_->NdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    int c = v_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    hess = scalar_view("hess", c, f);
    scalar_view valV_eval("valV_eval", c, p);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess,valV_eval,fe_->NdetJ(),false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Tracking_Nonlinear_Elliptic::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real, class DeviceType>
class QoI_ControlPenalty_Nonlinear_Elliptic : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;
  scalar_view ones_;

public:
  QoI_ControlPenalty_Nonlinear_Elliptic(const ROL::Ptr<fe_type> &fe) : fe_(fe) {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    // Build field container of ones
    ones_ = scalar_view("ones_",c,p);
    Kokkos::deep_copy(ones_,static_cast<Real>(1));
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Build local state tracking term
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    fe_->computeIntegral(val,ones_,valZ_eval);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Build local gradient of state tracking term
    fst::integrate(grad,ones_,fe_->NdetJ(),false);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_Nonlinear_Elliptic::HessVec_22 is zero.");
  }

}; // QoI_L2Penalty

#endif
