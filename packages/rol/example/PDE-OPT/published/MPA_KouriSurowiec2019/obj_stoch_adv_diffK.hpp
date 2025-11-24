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

#ifndef PDEOPT_QOI_STOCH_ADV_DIFFK_HPP
#define PDEOPT_QOI_STOCH_ADV_DIFFK_HPP

#include "../../TOOLS/qoiK.hpp"
#include "pde_stoch_adv_diffK.hpp"

template <class Real, class DeviceType>
class QoI_State_Cost_stoch_adv_diff : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;

public:
  QoI_State_Cost_stoch_adv_diff(const ROL::Ptr<fe_type> &fe) : fe_(fe) {
    const int d = fe_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->gradN().extent_int(0);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
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
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute gradient of squared L2-norm of diff
    fst::integrate(grad,valU_eval,fe_->NdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<scalar_view> & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::gradient_3 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c = fe_->gradN().extent_int(0);
    const int f = fe_->gradN().extent_int(1);
    const int p = fe_->gradN().extent_int(2);
    scalar_view valV_eval("valV_eval", c, p);
    hess = scalar_view("hess", c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess,valV_eval,fe_->NdetJ(),false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_12 is zero.");
  }

  void HessVec_13(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_13 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_22 is zero.");
  }

  void HessVec_23(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<scalar_view> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_33 is zero.");
  }

}; // QoI_State_Cost

template <class Real, class DeviceType>
class QoI_Control_Cost_stoch_adv_diff : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

public:
  QoI_Control_Cost_stoch_adv_diff(void) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    Real sum(0);
    for (const auto zi : *z_param) sum += zi;
    val = scalar_view();
    return sum;
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<scalar_view> & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int size = z_param->size();
    std::vector<Real> g(size,static_cast<Real>(1));
    return g;
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_12 is zero.");
  }

  void HessVec_13(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_13 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_22 is zero.");
  }

  void HessVec_23(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<scalar_view> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_33 is zero.");
  }

}; // QoI_Control_Cost

template <class Real>
class StdObjective_stoch_adv_diff : public ROL::StdObjective<Real> {
private:
  Real alpha1_, alpha2_;

public:
  StdObjective_stoch_adv_diff(ROL::ParameterList &parlist) {
    alpha1_ = parlist.sublist("Problem").get("State Cost",1.e5);
    alpha2_ = parlist.sublist("Problem").get("Control Cost",1.e0);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    return alpha1_*x[0] + alpha2_*x[1];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    g[0] = alpha1_;
    g[1] = alpha2_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;
  }

}; // OBJ_SCALAR

#endif
