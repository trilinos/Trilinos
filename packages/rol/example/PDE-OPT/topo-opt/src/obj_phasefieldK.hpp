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

#ifndef PDEOPT_QOI_PHASEFIELDK_HPP
#define PDEOPT_QOI_PHASEFIELDK_HPP

#include "../../TOOLS/qoiK.hpp"
#include "../../TOOLS/feK.hpp"

template <class Real, class DeviceType>
class QoI_ModicaMortolaEnergy_PhaseField : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  const Real scale_;

public:
  QoI_ModicaMortolaEnergy_PhaseField(const ROL::Ptr<fe_type> &fe,
                                     const Real scale = Real(1))
    : fe_(fe), scale_(scale) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1), cPhi(0.75);
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize output val
    val = scalar_view("val",c);
    // Evaluate on FE basis
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view gradZ_eval("gradZ_eval", c, p, d);
    scalar_view valPhi_eval("valPhi_eval", c, p);
    fe_->evaluateValue(valZ_eval,z_coeff);
    fe_->evaluateGradient(gradZ_eval,z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        Real zij = valZ_eval(i,j);
        valPhi_eval(i,j) = cPhi*(zij*zij - one)/scale_;
      }
    }
    // Integrate
    fe_->computeIntegral(val,gradZ_eval,gradZ_eval,false);
    fe_->computeIntegral(val,valPhi_eval,valPhi_eval,true);
    rst::scale(val, static_cast<Real>(0.5)*scale_);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_ModicaMortolaEnergy_PhaseField::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1), cPhi(9.0/8.0);
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize Gradients.
    grad = scalar_view("grad",c,f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view gradZ_eval("gradZ_eval", c, p, d);
    scalar_view valPhi_eval("valPhi_eval", c, p);
    fe_->evaluateValue(valZ_eval,z_coeff);
    fe_->evaluateGradient(gradZ_eval,z_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        Real zij = valZ_eval(i,j);
        valPhi_eval(i,j) = cPhi*zij*(zij*zij-one)/scale_;
      }
    }
    // Evaluate gradient of energy.
    fst::integrate(grad,gradZ_eval,fe_->gradNdetJ(),false);
    rst::scale(grad,scale_);
    fst::integrate(grad,valPhi_eval,fe_->NdetJ(),true);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_ModicaMortolaEnergy_PhaseField::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_ModicaMortolaEnergy_PhaseField::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_ModicaMortolaEnergy_PhaseField::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const Real one(1), three(3), cPhi(9.0/8.0);
    // Retrieve dimensions.
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize Hessians.
    hess = scalar_view("hess",c,f);
    // Evaluate/interpolate finite element fields on cells.
    scalar_view valZ_eval("valZ_eval", c, p);
    scalar_view valV_eval("valV_eval", c, p);
    scalar_view gradV_eval("gradV_eval", c, p, d);
    scalar_view valPhi_eval("valPhi_eval", c, p);
    fe_->evaluateValue(valZ_eval,z_coeff);
    fe_->evaluateValue(valV_eval,v_coeff);
    fe_->evaluateGradient(gradV_eval,v_coeff);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        Real zij = valZ_eval(i,j);
        valPhi_eval(i,j) = valV_eval(i,j)*cPhi*(three*zij*zij-one)/scale_;
      }
    }
    // Evaluate hessian-times-vector of energy.
    fst::integrate(hess,gradV_eval,fe_->gradNdetJ(),false);
    rst::scale(hess,scale_);
    fst::integrate(hess,valPhi_eval,fe_->NdetJ(),true);
  }

}; // QoI_ModicaMortolaEnergy_PhaseField

template <class Real, class DeviceType>
class QoI_Volume_PhaseField : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;

public:
  QoI_Volume_PhaseField(const ROL::Ptr<fe_type> &fe) : fe_(fe) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val",c);
    // Evaluate on FE basis
    scalar_view valZ_eval("valZ_eval",c,p);
    fe_->evaluateValue(valZ_eval,z_coeff);
    // Approximate characteristic function
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) += static_cast<Real>(1);
      }
    }
    // Compute volume
    fe_->computeIntegral(val,valZ_eval,valZ_eval);
    rst::scale(val,static_cast<Real>(0.25));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_PhaseField::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Evaluate on FE basis
    scalar_view valZ_eval("valZ_eval",c,p);
    fe_->evaluateValue(valZ_eval,z_coeff);
    // Approximate derivative of characteristic function
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        valZ_eval(i,j) += static_cast<Real>(1);
      }
    }
    // Compute gradient of volume
    fst::integrate(grad,valZ_eval,(fe_->NdetJ()),false);
    rst::scale(grad,static_cast<Real>(0.5));
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_PhaseField::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_PhaseField::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_PhaseField::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);
    // Initialize output grad
    hess = scalar_view("hess",c,f);
    // Evaluate on FE basis
    scalar_view valV_eval("valV_eval",c,p);
    fe_->evaluateValue(valV_eval,v_coeff);
    // Compute gradient of volume
    fst::integrate(hess,valV_eval,fe_->NdetJ(),false);
    rst::scale(hess,static_cast<Real>(0.5));
  }

}; // QoI_Volume_PhaseField

#endif
