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

#ifndef PDEOPT_QOI_L2TRACKING_POISSON_BOLTZMANNK_HPP
#define PDEOPT_QOI_L2TRACKING_POISSON_BOLTZMANNK_HPP

#include "../TOOLS/qoiK.hpp"
#include "pde_poisson_boltzmannK.hpp"

template <class Real, class DeviceType>
class QoI_L2Tracking_Poisson_Boltzmann : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;

  scalar_view target_;

  Real targetFunc(const std::vector<Real> & x) const {
    int size = x.size();
    Real val(0);
    for (int i = 0; i < size; ++i) {
      val += x[i]*x[i];
    }
    return val;
  }

public:
  QoI_L2Tracking_Poisson_Boltzmann(const ROL::Ptr<fe_type> & fe) : fe_(fe) {
    int c = fe_->cubPts().extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int d = fe_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    target_ = scalar_view("target_", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (fe_->cubPts())(i, j, k);
        }
        (target_)(i,j) = targetFunc(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval,target_);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    rst::scale(val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    // Initialize output grad
    grad = scalar_view("gradient_1", c, f);
    // Evaluate state on FE basis
    scalar_view valU_eval("valU_eval", c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute difference between state and target
    rst::subtract(valU_eval, target_);
    // Compute gradient of squared L2-norm of diff
    fst::integrate(grad, valU_eval, fe_->NdetJ(), false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Poisson_Boltzmann::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    scalar_view valV_eval("valV_eval", c, p);
    hess = scalar_view("hess11", c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess, valV_eval, fe_->NdetJ(), false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Tracking_Poisson_Boltzmann::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real, class DeviceType>
class QoI_L2Penalty_Poisson_Boltzmann : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using fe_type = FE<Real, DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<fe_type> fe_;

public:
  QoI_L2Penalty_Poisson_Boltzmann(const ROL::Ptr<fe_type> & fe) : fe_(fe) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Build local state tracking term
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    fe_->computeIntegral(val, valZ_eval, valZ_eval);
    rst::scale(val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Poisson_Boltzmann::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = z_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Build local gradient of state tracking term
    scalar_view valZ_eval("valZ_eval", c, p);
    fe_->evaluateValue(valZ_eval, z_coeff);
    fst::integrate(grad, valZ_eval, fe_->NdetJ(), false);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Poisson_Boltzmann::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Poisson_Boltzmann::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_Poisson_Boltzmann::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    int c = v_coeff.extent_int(0);
    int p = fe_->cubPts().extent_int(1);
    int f = fe_->N().extent_int(1);
    scalar_view valV_eval("valV_eval", c, p);
    hess = scalar_view("hess22", c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    fst::integrate(hess, valV_eval, fe_->NdetJ(), false);
  }

}; // QoI_L2Penalty

template <class Real>
class StdObjective_Poisson_Boltzmann : public ROL::StdObjective<Real> {
private:
  Real alpha_;

public:
  StdObjective_Poisson_Boltzmann(Teuchos::ParameterList &parlist) {
    alpha_ = parlist.sublist("Poisson Boltzmann Objective Function").get("Control Penalty",1.e-4);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    return x[0] + alpha_*x[1];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    const Real one(1);
    g[0] = one;
    g[1] = alpha_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;
  }

}; // OBJ_SCALAR

#endif
