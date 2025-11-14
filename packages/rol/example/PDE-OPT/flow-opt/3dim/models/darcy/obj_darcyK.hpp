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

#ifndef OBJ_DARCYK_HPP
#define OBJ_DARCYK_HPP

#include "../../../../TOOLS/qoiK.hpp"
#include "pde_darcyK.hpp"
#include "permeabilityK.hpp"


template <class Real, class DeviceType>
class QoI_VelocityTracking_Darcy : public QoI<Real, DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type>            fePrs_, feCtrl_;
  const ROL::Ptr<Permeability<Real,DeviceType>> perm_;
  scalar_view                        target_, weight_;
  Real                               rad_, yvel_, frac_, twpow_;
  bool                               onlyAxial_;

  Real xTarget(const std::vector<Real> &x) const {
    const Real X = x[0], Y = x[1];
    //return xWeight(x) ? -X*Y/(rad_*rad_-Y*Y) : zero;
    //return xWeight(x) ? -X*Y/std::sqrt(rad_*rad_-Y*Y) : zero;
    //return polyWeight(x) * (-X*Y/std::sqrt(rad_*rad_-Y*Y));
    //return -X*Y/std::sqrt(rad_*rad_-Y*Y);
    return -X*Y/((rad_*rad_-Y*Y)*(rad_*rad_-Y*Y));
  }

  Real yTarget(const std::vector<Real> &x) const {
    const Real one(1), Y = x[1];
    //return yWeight(x) ? one : zero;
    //return yWeight(x) ? std::sqrt(rad_*rad_-Y*Y) : zero;
    //return polyWeight(x) * std::sqrt(rad_*rad_-Y*Y);
    //return std::sqrt(rad_*rad_-Y*Y);
    return one/(rad_*rad_-Y*Y);
  }

  Real xWeight(const std::vector<Real> &x) const {
    return yWeight(x);
  }

  Real yWeight(const std::vector<Real> &x) const {
    //const Real zero(0), one(1), Y = x[1];
    //return (std::abs(Y) <= frac_*rad_ ? one : zero);
    return polyWeight(x);
  }

  Real polyWeight(const std::vector<Real> &x) const {
    const Real zero(0), one(1), Y = x[1], p = twpow_;
    const Real yTOP = 9.976339196;
    const Real yBOT = -yTOP;
    Real val = 0, at = 0, bt = 0;
    at = one / std::pow(-yTOP,p);
    bt = one / std::pow(-yBOT,p);
    if (Y > zero) val = at*std::pow(Y-yTOP,p);
    else          val = bt*std::pow(Y-yBOT,p);
    //std::cout << Y << "  " << val << std::endl;
    return val;
  }

public:
  QoI_VelocityTracking_Darcy(ROL::ParameterList             &list,
                              const ROL::Ptr<fe_type>           &fePrs,
                              const ROL::Ptr<fe_type>           &feCtrl,
                              const ROL::Ptr<Permeability<Real,DeviceType>> &perm)
    : fePrs_(fePrs), feCtrl_(feCtrl), perm_(perm) {
    rad_         = list.sublist("Problem").get("Diffuser Radius",5.0);
    yvel_        = list.sublist("Problem").get("Target Axial Velocity",15.0);
    frac_        = list.sublist("Problem").get("Integration Domain Fraction",0.95);
    onlyAxial_   = list.sublist("Problem").get("Only Use Axial Velocity",false);
    twpow_       = list.sublist("Problem").get("Target Weighting Power",0.0);
    Real xWScal  = list.sublist("Problem").get("Radial Tracking Scale",1.0);
    Real yWScal  = list.sublist("Problem").get("Axial Tracking Scale",1.0);
    bool useNorm = list.sublist("Problem").get("Use Normalized Misfit",false);
    useNorm = onlyAxial_ ? false : useNorm;
    xWScal  = onlyAxial_ ? static_cast<Real>(0) : xWScal;
    const int c = fePrs_->gradN().extent_int(0);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    target_ = scalar_view("target",c,p,d);
    weight_ = scalar_view("weight",c,p,d);
    std::vector<Real> x(2);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        x[0] = (fePrs_->cubPts())(i,j,0);
        x[1] = (fePrs_->cubPts())(i,j,1);
        target_(i,j,0) = xTarget(x);
        target_(i,j,1) = yTarget(x);
        if (useNorm && yWeight(x)) {
          xWScal = static_cast<Real>(1)
                  /(std::pow(target_(i,j,0),2) + std::pow(target_(i,j,1),2));
          yWScal = xWScal;
        }
        weight_(i,j,0) = x[0] * xWScal * xWeight(x);
        weight_(i,j,1) = x[0] * yWScal * yWeight(x);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output val
    val = scalar_view("val",c);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view vel("vel",c,p,d);
    scalar_view wvel("wvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          vel(i,j,k)  = alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k);
          wvel(i,j,k) = weight_(i,j,k)*vel(i,j,k);
        }
      }
    }

    fePrs_->computeIntegral(val,vel,wvel);
    // Scale by one half
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    grad = scalar_view("grad", c, f);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view awvel("awvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          awvel(i,j,k) = weight_(i,j,k)*alpha(i,j)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k));
        }
      }
    }

    fst::integrate(grad,awvel,fePrs_->gradNdetJ(),false);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    grad = scalar_view("grad",c,fc);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    scalar_view deriv("deriv",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          deriv(i,j) += weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k))
                        *dalpha(i,j)*gradU(i,j,k);
        }
      }
    }

    fst::integrate(grad,deriv,feCtrl_->NdetJ(),false);
  }

  std::vector<Real> gradient_3(std::vector<scalar_view> & grad,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int i = 0; i < d; ++i)
        grad[i] = scalar_view("grad",c);
      // Compute cost integral
      scalar_view gradU("gradU",c,p,d);
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view alpha("alpha",c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k));
          }
        }
      }

      fePrs_->computeIntegral(grad[0],wvel,target_);

      return g_param;
    }
    else {
      throw Exception::Zero(">>> QoI_VelocityTracking_Darcy::gradient_3 is zero.");
    }
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    // Initialize output hess
    hess = scalar_view("hess",c,f);
    // Compute cost integral
    scalar_view gradV("gradV",c,p,d);
    scalar_view awvel("awvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    fePrs_->evaluateGradient(gradV, v_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          awvel(i,j,k) = alpha(i,j)*alpha(i,j)*weight_(i,j,k)*gradV(i,j,k);
        }
      }
    }

    fst::integrate(hess,awvel,fePrs_->gradNdetJ(),false);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fePrs_->gradN().extent_int(0);
    const int f = fePrs_->gradN().extent_int(1);
    const int p = fePrs_->gradN().extent_int(2);
    const int d = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    hess = scalar_view("hess",c,f);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view awvel("awvel",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view valV("valV",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    feCtrl_->evaluateValue(valV, v_coeff);
    perm_->compute( alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          awvel(i,j,k)  = static_cast<Real>(2)*alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k);
          awvel(i,j,k) *= weight_(i,j,k)*dalpha(i,j)*valV(i,j);
        }
      }
    }

    fst::integrate(hess,awvel,fePrs_->gradNdetJ(),false);
  }

  void HessVec_13(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int f = fePrs_->gradN().extent_int(1);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      hess = scalar_view("hess",c,f);
      // Compute cost integral
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view alpha("alpha",c,p);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*alpha(i,j)*target_(i,j,k)*(*v_param)[0];
          }
        }
      }

      fst::integrate(hess,wvel,fePrs_->gradNdetJ(),false);
    }
    else {
      throw Exception::NotImplemented(">>> HessVec_13 not implemented.");
    }
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    hess = scalar_view("hess",c,fc);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view gradV("gradV",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    scalar_view deriv("deriv",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    fePrs_->evaluateGradient(gradV, v_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          deriv(i,j) += weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k))
                        *dalpha(i,j)*gradV(i,j,k);
          deriv(i,j) += weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)*alpha(i,j)*gradV(i,j,k);
        }
      }
    }

    fst::integrate(hess,deriv,feCtrl_->NdetJ(),false);
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = fePrs_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = fePrs_->gradN().extent_int(2);
    const int d  = fePrs_->gradN().extent_int(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output hess
    hess = scalar_view("hess",c,fc);
    // Compute cost integral
    scalar_view gradU("gradU",c,p,d);
    scalar_view valZ("valZ",c,p);
    scalar_view valV("valV",c,p);
    scalar_view alpha("alpha",c,p);
    scalar_view dalpha("dalpha",c,p);
    scalar_view ddalpha("ddalpha",c,p);
    scalar_view deriv("deriv",c,p);
    fePrs_->evaluateGradient(gradU, u_coeff);
    feCtrl_->evaluateValue(valZ, z_coeff);
    feCtrl_->evaluateValue(valV, v_coeff);
    perm_->compute(  alpha, valZ, fePrs_->cubPts(), 0);
    perm_->compute( dalpha, valZ, fePrs_->cubPts(), 1);
    perm_->compute(ddalpha, valZ, fePrs_->cubPts(), 2);

    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          deriv(i,j) += weight_(i,j,k)*(alpha(i,j)*gradU(i,j,k)+yvel*target_(i,j,k))
                        *ddalpha(i,j)*valV(i,j)*gradU(i,j,k);
          deriv(i,j) += weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)
                        *dalpha(i,j)*valV(i,j)*gradU(i,j,k);
        }
      }
    }

    fst::integrate(hess,deriv,feCtrl_->NdetJ(),false);
  }

  void HessVec_23(scalar_view & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c  = fePrs_->gradN().extent_int(0);
      const int fc = feCtrl_->N().extent_int(1);
      const int p  = fePrs_->gradN().extent_int(2);
      const int d  = fePrs_->gradN().extent_int(3);
      // Initialize output val
      hess = scalar_view("hess",c,fc);
      // Compute cost integral
      scalar_view gradU("gradU",c,p,d);
      scalar_view wvel("wvel",c,p);
      scalar_view valZ("valZ",c,p);
      scalar_view dalpha("dalpha",c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j) += weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)*target_(i,j,k)*(*v_param)[0];
          }
        }
      }

      fst::integrate(hess,wvel,feCtrl_->NdetJ(),false);
    }
    else {
      throw Exception::Zero(">>> QoI_VelocityTracking_Darcy::HessVec_23 is zero.");
    }
  }

  std::vector<Real> HessVec_31(std::vector<scalar_view> & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess",c);
      // Compute cost integral
      scalar_view gradV("gradV",c,p,d);
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view alpha("alpha",c,p);
      fePrs_->evaluateGradient(gradV, v_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      perm_->compute(alpha, valZ, fePrs_->cubPts(), 0);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*alpha(i,j)*gradV(i,j,k);
          }
        }
      }

      fePrs_->computeIntegral(hess[0],wvel,target_);

      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_VelocityTracking_Darcy::HessVec_31 is zero.");
    }
  }

  std::vector<Real> HessVec_32(std::vector<scalar_view> & hess,
                               const scalar_view v_coeff,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess",c);
      // Compute cost integral
      scalar_view gradU("gradU",c,p,d);
      scalar_view wvel("wvel",c,p,d);
      scalar_view valZ("valZ",c,p);
      scalar_view valV("valV",c,p);
      scalar_view dalpha("dalpha",c,p);
      fePrs_->evaluateGradient(gradU, u_coeff);
      feCtrl_->evaluateValue(valZ, z_coeff);
      feCtrl_->evaluateValue(valV, v_coeff);
      perm_->compute(dalpha, valZ, fePrs_->cubPts(), 1);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wvel(i,j,k) = weight_(i,j,k)*dalpha(i,j)*gradU(i,j,k)*valV(i,j);
          }
        }
      }

      fePrs_->computeIntegral(hess[0],wvel,target_);

      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_VelocityTracking_Darcy::HessVec_32 is zero.");
    }
  }

  std::vector<Real> HessVec_33(std::vector<scalar_view> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const scalar_view u_coeff,
                               const scalar_view z_coeff = scalar_view(),
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN().extent_int(0);
      const int p = fePrs_->gradN().extent_int(2);
      const int d = fePrs_->gradN().extent_int(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = scalar_view("hess",c);
      // Compute cost integral
      scalar_view wtarget("wtarget",c,p,d);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            wtarget(i,j,k) = weight_(i,j,k)*target_(i,j,k);
          }
        }
      }
      fePrs_->computeIntegral(hess[0],wtarget,target_);
      rst::scale(hess[0],(*v_param)[0]);

      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_VelocityTracking_Darcy::HessVec_33 is zero.");
    }
  }

}; // QoI_VelocityTracking_Darcy

#endif
