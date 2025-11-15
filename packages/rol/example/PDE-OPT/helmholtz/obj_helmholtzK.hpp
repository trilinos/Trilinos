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

#include "../TOOLS/qoiK.hpp"
#include "pde_helmholtzK.hpp"

template <class Real, class DeviceType>
class QoI_Helmholtz_StateTracking : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  const ROL::Ptr<FieldHelper<Real, DeviceType>> fieldHelper_;

  scalar_view weight_;
  std::vector<scalar_view> target_;

  Real RoiRadius_;
  Real waveNumber_;
  Real angle_;

protected:
  void computeDomainWeight(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    weight_ = scalar_view("weight_", c, p);

    const Real zero(0), one(1);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = fe_->cubPts()(i, j, k);
        }
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        weight_(i, j) = (inside ? one : zero);
      }
    }
  }

  void computeTarget(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    target_.clear();
    target_.resize(2);
    target_[0] = scalar_view("target_[0]", c, p);
    target_[1] = scalar_view("target_[1]", c, p);

    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = fe_->cubPts()(i, j, k);
        }
        (target_[0])(i, j) = evaluateTarget(x, 0);
        (target_[1])(i, j) = evaluateTarget(x, 1);
      }
    }
  }

  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

public:
  QoI_Helmholtz_StateTracking(const ROL::Ptr<fe_type> &fe,
                              const ROL::Ptr<FieldHelper<Real, DeviceType>> &fieldHelper,
                              Teuchos::ParameterList &parlist)
    : fe_(fe), fieldHelper_(fieldHelper) {
    RoiRadius_  = parlist.sublist("Problem").get("ROI Radius", 2.0);
    waveNumber_ = parlist.sublist("Problem").get("Wave Number", 10.0);
    angle_      = parlist.sublist("Problem").get("Target Angle", 45.0);
    angle_      = DegreesToRadians(angle_);
    computeDomainWeight();
    computeTarget();
  }

  virtual bool insideDomain(const std::vector<Real> &x) const {
    Real xnorm(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    return (xnorm <= RoiRadius_);
  }

  virtual Real evaluateTarget(const std::vector<Real> &x, const int component) const {
    const Real arg = waveNumber_ * (std::cos(angle_)*x[0] + std::sin(angle_)*x[1]);
    return (component==0) ? std::cos(arg) : std::sin(arg);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate tracking term
    scalar_view diffU("diffU", c, p);
    scalar_view WdiffU("WdiffU", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(diffU, static_cast<Real>(0));
      Kokkos::deep_copy(WdiffU, static_cast<Real>(0));
      fe_->evaluateValue(diffU, U[i]);
      rst::subtract(diffU, target_[i]);
      fst::scalarMultiplyDataData(WdiffU, weight_, diffU);
      fe_->computeIntegral(val, WdiffU, diffU, true);
    }
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
    // Initialize output grad
    std::vector<scalar_view> G(2);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate tracking term
    scalar_view diffU("diffU", c, p);
    scalar_view WdiffU("WdiffU", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(diffU, static_cast<Real>(0));
      Kokkos::deep_copy(WdiffU, static_cast<Real>(0));
      fe_->evaluateValue(diffU, U[i]);
      rst::subtract(diffU, target_[i]);
      fst::scalarMultiplyDataData(WdiffU, weight_, diffU);
      G[i] = scalar_view("G", c, f);
      fst::integrate(G[i], WdiffU, fe_->NdetJ(), false);
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Initialize output hessvec
    std::vector<scalar_view> H(2);
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate tracking term
    scalar_view valV("valV", c, p);
    scalar_view WvalV("WvalV", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valV, static_cast<Real>(0));
      Kokkos::deep_copy(WvalV, static_cast<Real>(0));
      fe_->evaluateValue(valV, V[i]);
      fst::scalarMultiplyDataData(WvalV, weight_, valV);
      H[i] = scalar_view("H", c, f);
      fst::integrate(H[i], WvalV, fe_->NdetJ(), false);
    }
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_22 is zero.");
  }

}; // QoI_Helmholtz


template <class Real, class DeviceType>
class QoI_Helmholtz_ControlPenalty : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  const ROL::Ptr<FieldHelper<Real, DeviceType>> fieldHelper_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  Real RoiRadius_;
  Real alpha_;

  scalar_view weight_;

  void computeDomainWeight(void) {
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    int d = fe_->gradN().extent_int(3);

    weight_ = scalar_view("weight_", c, p);

    const Real zero(0);
    bool inside(false);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      inside = false;
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = fe_->cubPts()(i, j, k);
        }
        if ( insideDomain(x) ) {
          inside = true;
          break;
        }
      }
      for (int j = 0; j < p; ++j) {
        weight_(i,j) = (inside ? alpha_ : zero);
      }
    }
  }

public:
  QoI_Helmholtz_ControlPenalty(const ROL::Ptr<fe_type> &fe,
                               const ROL::Ptr<FieldHelper<Real, DeviceType>> &fieldHelper,
                               Teuchos::ParameterList &parlist)
  : fe_(fe), fieldHelper_(fieldHelper) {
    Real dist2annulus   = parlist.sublist("Problem").get("Distance to Control Annulus",0.5);
    Real annulusWidth   = parlist.sublist("Problem").get("Control Annulus Width",0.1);
    RoiRadius_          = parlist.sublist("Problem").get("ROI Radius",2.0);
    innerAnnulusRadius_ = RoiRadius_ + dist2annulus;
    outerAnnulusRadius_ = innerAnnulusRadius_ + annulusWidth;
    alpha_              = parlist.sublist("Problem").get("Control Penalty",1e-4);
    computeDomainWeight();
  }

  virtual bool insideDomain(const std::vector<Real> &x) const {
    Real xnorm(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    return (xnorm <= outerAnnulusRadius_ && xnorm >= innerAnnulusRadius_);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int p = fe_->gradN().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate control penalty
    scalar_view valZ("valZ", c, p);
    scalar_view WvalZ("WvalZ", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valZ, static_cast<Real>(0));
      Kokkos::deep_copy(WvalZ, static_cast<Real>(0));
      fe_->evaluateValue(valZ, Z[i]);
      fst::scalarMultiplyDataData(WvalZ, weight_, valZ);
      fe_->computeIntegral(val, WvalZ, valZ, true);
    }
    rst::scale(val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> G(2);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate control penalty
    scalar_view valZ("valZ", c, p);
    scalar_view WvalZ("WvalZ", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valZ, static_cast<Real>(0));
      Kokkos::deep_copy(WvalZ, static_cast<Real>(0));
      fe_->evaluateValue(valZ, Z[i]);
      fst::scalarMultiplyDataData(WvalZ, weight_, valZ);
      G[i] = scalar_view("G", c, f);
      fst::integrate(G[i], WvalZ, fe_->NdetJ(), false);
    }
    fieldHelper_->combineFieldCoeff(grad, G);
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
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = fe_->gradN().extent_int(0);
    int f = fe_->gradN().extent_int(1);
    int p = fe_->gradN().extent_int(2);
    // Initialize output hessvec
    std::vector<scalar_view> H(2);
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate control penalty
    scalar_view valV("valV", c, p);
    scalar_view WvalV("WvalV", c, p);
    for (int i=0; i<2; ++i) {
      Kokkos::deep_copy(valV, static_cast<Real>(0));
      Kokkos::deep_copy(WvalV, static_cast<Real>(0));
      fe_->evaluateValue(valV, V[i]);
      fst::scalarMultiplyDataData(WvalV, weight_, valV);
      H[i] = scalar_view("H", c, f);
      fst::integrate(H[i], WvalV, fe_->NdetJ(), false);
    }
    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_Helmholtz_ControlPenalty

#endif
