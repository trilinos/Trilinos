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

#ifndef PDEOPT_QOI_HELMHOLTZK_HPP
#define PDEOPT_QOI_HELMHOLTZK_HPP

#include "../../TOOLS/qoiK.hpp"
#include "pde_helmholtz_realK.hpp"

template <class Real, class DeviceType>
class QoI_Helmholtz_StateTracking : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;

  scalar_view weight_;
  std::vector<scalar_view> target_;

  Real RoiRadius_;
  Real waveNumber_;
  Real angle_;
  int example_;

  unsigned comp_;

protected:
  bool insideDomain(const std::vector<Real> &x) const {
    bool val = true;
    if (example_==1) {
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      Real xnorm(0);
      for (const auto xi : x) xnorm += xi*xi;
      xnorm = std::sqrt(xnorm);
      val = (xnorm <= RoiRadius_+eps);
    }
    return val;
  }

  void computeDomainWeight(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    weight_ = scalar_view("weight",c,p);
   
    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          x[k] = (fe_->cubPts())(i,j,k);
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j)
        weight_(i,j) = (inside ? one : zero);
    }
  }

  Real evaluateTarget(const std::vector<Real> &x, const int component) const {
    const Real arg = waveNumber_ * (std::cos(angle_)*x[0] + std::sin(angle_)*x[1]);
    return (component==0) ? std::cos(arg) : std::sin(arg);
  }

  void computeTarget(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    target_.clear(); target_.resize(2);
    target_[0] = scalar_view("target",c,p);
    target_[1] = scalar_view("target",c,p);
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          x[k] = (fe_->cubPts())(i,j,k);
        (target_[0])(i,j) = evaluateTarget(x,0);
        (target_[1])(i,j) = evaluateTarget(x,1);
      }
    }
  }

  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }
  
public:
  QoI_Helmholtz_StateTracking(const ROL::Ptr<fe_type> &fe,
                              ROL::ParameterList &parlist,
                              int comp = 0)
    : fe_(fe), RoiRadius_(2), comp_(comp) {
    waveNumber_ = parlist.sublist("Problem").get("Wave Number",10.0);
    angle_      = parlist.sublist("Problem").get("Target Angle",45.0);
    angle_      = DegreesToRadians(angle_);
    example_    = parlist.sublist("Problem").get("Example",1);
    computeDomainWeight();
    computeTarget();
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
    // Evaluate tracking term
    scalar_view diffU("diffU", c, p);
    scalar_view WdiffU("WdiffU", c, p);
    fe_->evaluateValue(diffU, u_coeff);
    rst::subtract(diffU,target_[comp_]);
    fst::scalarMultiplyDataData(WdiffU,weight_,diffU);
    fe_->computeIntegral(val,WdiffU,diffU,true);
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Evaluate tracking term
    scalar_view diffU("diffU", c, p);
    scalar_view WdiffU("WdiffU", c, p);
    fe_->evaluateValue(diffU, u_coeff);
    rst::subtract(diffU,target_[comp_]);
    fst::scalarMultiplyDataData(diffU,weight_,diffU);
    grad = scalar_view("grad", c, f);
    fst::integrate(grad,WdiffU,fe_->NdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Evaluate tracking term
    scalar_view valV("valV", c, p);
    scalar_view WvalV("WvalV", c, p);
    fe_->evaluateValue(valV, v_coeff);
    fst::scalarMultiplyDataData(WvalV,weight_,valV);
    hess = scalar_view("hess", c, f);
    fst::integrate(hess,WvalV,fe_->NdetJ(),false);
  }

  void Hessian_11(scalar_view & hess,
            const scalar_view u_coeff,
            const scalar_view z_coeff = scalar_view(),
            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Build force/control term
    scalar_view F("F", c, f, p);
    fst::scalarMultiplyDataField(F, weight_, fe_->N());
    hess = scalar_view("hess", c, f, f);
    fst::integrate(hess,F,fe_->NdetJ(),false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_22 is zero.");
  }

}; // QoI_Helmholtz


template <class Real, class DeviceType>
class QoI_Helmholtz_ControlPenalty : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  int example_;

  scalar_view weight_;

  bool insideDomain(const std::vector<Real> &x) const {
    bool val = true;
    if (example_==1) {
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      Real xnorm(0);
      for (const auto xi : x) xnorm += xi*xi;
      xnorm = std::sqrt(xnorm);
      val = (xnorm <= outerAnnulusRadius_+eps && xnorm >= innerAnnulusRadius_-eps);
    }
    return val;
  }

  void computeDomainWeight(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    weight_ = scalar_view("weigt",c,p);
   
    const Real one(1), zero(0);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          x[k] = (fe_->cubPts())(i,j,k);
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j)
        weight_(i,j) = (inside ? one : zero);
    }
  }

public:
  QoI_Helmholtz_ControlPenalty(const ROL::Ptr<fe_type> &fe,
                               ROL::ParameterList &parlist)
  : fe_(fe), innerAnnulusRadius_(2.5), outerAnnulusRadius_(2.6) {
    example_ = parlist.sublist("Problem").get("Example",1);
    computeDomainWeight();
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
    // Evaluate control penalty
    scalar_view valZ("valZ", c, p);
    scalar_view WvalZ("WvalZ", c, p);
    fe_->evaluateValue(valZ, z_coeff);
    fst::scalarMultiplyDataData(WvalZ,weight_,valZ);
    fe_->computeIntegral(val,WvalZ,valZ,true);
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Evaluate control penalty
    scalar_view valZ("valZ", c, p);
    scalar_view WvalZ("WvalZ", c, p);
    fe_->evaluateValue(valZ, z_coeff);
    fst::scalarMultiplyDataData(WvalZ,weight_,valZ);
    grad = scalar_view("grad", c, f);
    fst::integrate(grad,WvalZ,fe_->NdetJ(),false);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_21 is zero.");
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
    // Evaluate control penalty
    scalar_view valV("valV", c, p);
    scalar_view WvalV("WvalV", c, p);
    fe_->evaluateValue(valV, v_coeff);
    fst::scalarMultiplyDataData(WvalV,weight_,valV);
    hess = scalar_view("hess", c, f);
    fst::integrate(hess,WvalV,fe_->NdetJ(),false);
  }

  void Hessian_22(scalar_view & hess,
            const scalar_view u_coeff,
            const scalar_view z_coeff = scalar_view(),
            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Build force/control term
    scalar_view F("F", c, f, p);
    fst::scalarMultiplyDataField(F, weight_, fe_->N());
    hess = scalar_view("hess", c, f, f);
    fst::integrate(hess,F,fe_->NdetJ(),false);
  }

}; // QoI_Helmholtz_ControlPenalty

#endif
