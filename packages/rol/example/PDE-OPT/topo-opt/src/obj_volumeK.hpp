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

#ifndef PDEOPT_QOI_VOLUMEK_HPP
#define PDEOPT_QOI_VOLUMEK_HPP

#include "../../TOOLS/qoiK.hpp"
#include "../../TOOLS/feK.hpp"

template <class Real, class DeviceType>
class QoI_Volume_TopoOpt : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  scalar_view ones_, volFrac_;

public:
  QoI_Volume_TopoOpt(const ROL::Ptr<fe_type> &fe,
                     const Real v0 = Real(0))
  : fe_(fe) {
    // Get relevant dimensions
    int c = fe_->cubPts().extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    ones_ = scalar_view("ones_",c,p);
    volFrac_ = scalar_view("volFrac_",c,p);
    Kokkos::deep_copy(ones_,static_cast<Real>(1));
    Kokkos::deep_copy(volFrac_,v0);
    //for (int i = 0; i < c; ++i) {
    //  for (int j = 0; j < p; ++j) {
    //    ones_(i,j) = static_cast<Real>(1);
    //    volFrac_(i,j) = v0;
    //  }
    //}
  }

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
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    rst::subtract(valZ_eval,volFrac_);

    // Compute energy
    fe_->computeIntegral(val,valZ_eval,ones_);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::gradient_1 is zero.");
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

    // Compute gradient of energy
    fst::integrate(grad,ones_,fe_->NdetJ(),false);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_TopoOpt::HessVec_22 is zero.");
  }

}; // QoI_VolumeObj

#endif
