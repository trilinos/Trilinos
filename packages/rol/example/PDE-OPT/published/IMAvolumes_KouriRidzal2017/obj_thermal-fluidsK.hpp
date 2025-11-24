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

#ifndef PDEOPT_QOI_THERMALFLUIDSK_HPP
#define PDEOPT_QOI_THERMALFLUIDSK_HPP

#include "../../TOOLS/qoiK.hpp"
#include "pde_thermal-fluidsK.hpp"

template <class Real, class DeviceType>
class QoI_Vorticity_ThermalFluids : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_, fePrs_, feThr_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  scalar_view weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Vorticity_ThermalFluids(const ROL::Ptr<fe_type> &feVel,
                              const ROL::Ptr<fe_type> &fePrs,
                              const ROL::Ptr<fe_type> &feThr,
                              const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper,
                              ROL::ParameterList &parlist)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
    int c = feVel_->cubPts().extent_int(0);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    weight_ = scalar_view("weight", c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          pt[k] = (feVel_->cubPts())(i,j,k);
        weight_(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<scalar_view> gradU_vec(d);
    for (int i = 0; i < d; ++i) {
      gradU_vec[i] = scalar_view("gradU_vec", c, p, d);
      feVel_->evaluateGradient(gradU_vec[i], U[i]);
    }
    // Compute weighted curl
    scalar_view curlU_eval;
    if (d==2)
      curlU_eval = scalar_view("curlU_eval", c, p);
    else if (d==3)
      curlU_eval = scalar_view"curlU_eval", (c, p, d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if (d==2)
          curlU_eval(i,j) = weight_(i,j) * ((gradU_vec[1])(i,j,0)-(gradU_vec[0])(i,j,1));
        else if (d==3) {
          for (int k = 0; k < d; ++k) {
            int i1 = (k+2)%d, i2 = (k+1)%d;
            curlU_eval(i,j,k) = weight_(i,j) * ((gradU_vec[i1])(i,j,i2)-(gradU_vec[i2])(i,j,i1));
          }
        }
      }
    }
    // Compute L2 norm squared
    feVel_->computeIntegral(val,curlU_eval,curlU_eval,false);
    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int fv = feVel_->N().extent_int(1);
    int fp = fePrs_->N().extent_int(1);
    int fh = feThr_->N().extent_int(1);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> G(d+2);
    for (int i = 0; i < d; ++i)
      G[i] = scalar_view("grad", c, fv);
    G[d] = scalar_view("grad", c, fp);
    G[d+1] = scalar_view("grad", c, fh);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<scalar_view> gradU_vec(d);
    for (int i = 0; i < d; ++i) {
      gradU_vec[i] = scalar_view("gradU_vec", c, p, d);
      feVel_->evaluateGradient(gradU_vec[i], U[i]);
    }
    // Compute weighted curl
    int size = (d==2) ? 1 : d;
    std::vector<scalar_view> curlU_vec(size);
    if (d==2) {
      curlU_vec[0] = scalar_view("curlU_vec", c, p);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j)
          (curlU_vec[0])(i,j) = weight_(i,j) * ((gradU_vec[1])(i,j,0)-(gradU_vec[0])(i,j,1));
      }
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        curlU_vec[i] = scalar_view("curlU_vec", c, p);
        for (int j = 0; j < c; ++j) {
          for (int k = 0; k < p; ++k) {
            int i1 = (i+2)%d, i2 = (i+1)%d;
            (curlU_vec[i])(j,k) = weight_(j,k) * ((gradU_vec[i1])(j,k,i2)-(gradU_vec[i2])(j,k,i1));
          }
        }
      }
    }
    // Build local gradient of state tracking term
    if (d==2) {
      fst::integrate(G[0],curlU_vec[0],feVel_->DNDdetJ(1),false);
      rst::scale(G[0],static_cast<Real>(-1));
      fst::integrate(G[1],curlU_vec[0],feVel_->DNDdetJ(0),false);
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        int i1 = (i+2)%d, i2 = (i+1)%d;
        fst::integrate(G[i],curlU_vec[i1],feVel_->DNDdetJ(i2),false);
        rst::scale(G[i],static_cast<Real>(-1));
        fst::integrate(G[i],curlU_vec[i2],feVel_->DNDdetJ(i1),false);
      }
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c  = z_coeff.extent_int(0);
    int p  = feVel_->cubPts().extent_int(1);
    int fv = feVel_->N().extent_int(1);
    int fp = fePrs_->N().extent_int(1);
    int fh = feThr_->N().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> H(d+2);
    for (int i = 0; i < d; ++i)
      H[i] = scalar_view("hess", c, fv);
    H[d] = scalar_view("hess", c, fp);
    H[d+1] = scalar_view("hess", c, fh);
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    std::vector<scalar_view> gradV_vec(d);
    for (int i = 0; i < d; ++i) {
      gradV_vec[i] = scalar_view("gradV_vec", c, p, d);
      feVel_->evaluateGradient(gradV_vec[i], V[i]);
    }
    // Compute weighted curl
    int size = (d==2) ? 1 : d;
    std::vector<scalar_view> curlV_vec(size);
    if (d==2) {
      curlV_vec[0] = scalar_view("curlV_vec", c, p);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j)
          (curlV_vec[0])(i,j) = weight_(i,j) * ((gradV_vec[1])(i,j,0)-(gradV_vec[0])(i,j,1));
      }
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        curlV_vec[i] = scalar_view("curlV_vec", c, p);
        for (int j = 0; j < c; ++j) {
          for (int k = 0; k < p; ++k) {
            int i1 = (i+2)%d, i2 = (i+1)%d;
            (curlV_vec[i])(j,k) = weight_(j,k) * ((gradV_vec[i1])(j,k,i2)-(gradV_vec[i2])(j,k,i1));
          }
        }
      }
    }
    // Build local gradient of state tracking term
    if (d==2) {
      fst::integrate(H[0],curlV_vec[0],feVel_->DNDdetJ(1),false);
      rst::scale(H[0],static_cast<Real>(-1));
      fst::integrate(H[1],curlV_vec[0],feVel_->DNDdetJ(0),false);
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        int i1 = (i+2)%d, i2 = (i+1)%d;
        fst::integrate(H[i],curlV_vec[i1],feVel_->DNDdetJ(i2),false);
        rst::scale(*H[i],static_cast<Real>(-1));
        fst::integrate(H[i],curlV_vec[i2],feVel_->DNDdetJ(i1),false);
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Vorticity_ThermalFluids

template <class Real>
class QoI_Circulation_ThermalFluids : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_, fePrs_, feThr_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  scalar_view weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Circulation_ThermalFluids(const ROL::Ptr<fe_type> &feVel,
                                const ROL::Ptr<fe_type> &fePrs,
                                const ROL::Ptr<fe_type> &feThr,
                                const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
    int c = feVel_->cubPts().extent_int(0);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    weight_ = scalar_view("weight",c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          pt[k] = (feVel_->cubPts())(i,j,k);
        weight_(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    scalar_view gradUX_eval("gradUX_eval", c, p, d);
    scalar_view gradUY_eval("gradUY_eval", c, p, d);
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    scalar_view curlU_eval("curlU_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j)
        curlU_eval(i,j) = gradUY_eval(i,j,0) - gradUX_eval(i,j,1);
    }
    // Compute circulation
    feVel_->computeIntegral(val,curlU_eval,weight_,false);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int fv = feVel_->N().extent_int(1);
    int fp = fePrs_->N().extent_int(1);
    int fh = feThr_->N().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> G(d+2);
    for (int i = 0; i < d; ++i)
      G[i] = scalar_view("grad", c, fv);
    G[d] = scalar_view("grad", c, fp);
    G[d+1] = scalar_view("grad", c, fh);
    // Build local gradient of state tracking term
    fst::integrate(G[0],weight_,feVel_->DNDdetJ(1),false);
    rst::scale(G[0],static_cast<Real>(-1));
    fst::integrate(G[1],weight_,feVel_->DNDdetJ(0),false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Circulation_ThermalFluids

template <class Real>
class QoI_Horizontal_ThermalFluids : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_, fePrs_, feThr_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  scalar_view weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Horizontal_ThermalFluids(const ROL::Ptr<fe_type> &feVel,
                               const ROL::Ptr<fe_type> &fePrs,
                               const ROL::Ptr<fe_type> &feThr,
                               const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
    int c = feVel_->cubPts().extent_int(0);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    std::vector<Real> pt(d);
    weight_ = scalar_view("weight",c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k)
          pt[k] = (feVel_->cubPts())(i,j,k);
        weight_(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int p = feVel_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    scalar_view valUX_eval("valUX_eval", c, p);
    scalar_view valUY_eval("valUY_eval", c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Compute negative part of x-velocity
    scalar_view minUX_eval("minUX_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j)
        minUX_eval(i,j) = std::min(static_cast<Real>(0),valUX_eval(i,j));
    }
    // Multiply by weight
    scalar_view weighted_minUX_eval("weighted_minUX_eval", c, p);
    fst::scalarMultiplyDataData(weighted_minUX_eval,weight_,minUX_eval);
    scalar_view weighted_valUY_eval("weighted_valUY_eval", c, p);
    fst::scalarMultiplyDataData(weighted_valUY_eval,weight_,valUY_eval);
    // Compute L2 norm squared
    feVel_->computeIntegral(val,minUX_eval,weighted_minUX_eval,false);
    feVel_->computeIntegral(val,valUY_eval,weighted_valUY_eval,true);

    rst::scale(val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c = u_coeff.extent_int(0);
    int fv = feVel_->N().extent_int(1);
    int fp = fePrs_->N().extent_int(1);
    int fh = feThr_->N().extent_int(1);
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> G(d+2);
    for (int i = 0; i < d; ++i)
      G[i] = scalar_view("grad", c, fv);
    G[d] = scalar_view("grad", c, fp);
    G[d+1] = scalar_view("grad", c, fh);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    scalar_view valUX_eval("valUX_eval", c, p);
    scalar_view valUY_eval("valUY_eval", c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Compute negative part of x-velocity
    scalar_view weighted_minUX_eval("weighted_minUX_eval", c, p);
    scalar_view weighted_valUY_eval("weighted_valUY_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        weighted_minUX_eval(i,j) = weight_(i,j) * std::min(static_cast<Real>(0),valUX_eval(i,j));
        weighted_valUY_eval(i,j) = weight_(i,j) * valUY_eval(i,j);
      }
    }
    // Build local gradient of state tracking term
    fst::integrate(G[0],weighted_minUX_eval,feVel_->NdetJ(),false);
    fst::integrate(G[1],weighted_valUY_eval,feVel_->NdetJ(),false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    int c  = z_coeff.extent_int(0);
    int p  = feVel_->cubPts().extent_int(1);
    int fv = feVel_->N().extent_int(1);
    int fp = fePrs_->N().extent_int(1);
    int fh = feThr_->N().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    std::vector<scalar_view> H(d+2);
    for (int i = 0; i < d; ++i)
      H[i] = scalar_view("hess", c, fv);
    H[d] = scalar_view("hess", c, fp);
    H[d+1] = scalar_view("hess", c, fh);
    // Get components of the control
    std::vector<scalar_view> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    scalar_view valUX_eval("valUX_eval", c, p);
    scalar_view valVX_eval("valVX_eval", c, p);
    scalar_view valVY_eval("valVY_eval", c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valVX_eval, V[0]);
    feVel_->evaluateValue(valVY_eval, V[1]);
    // Compute negative part of x-velocity
    scalar_view weighted_minVX_eval("weighted_minVX_eval", c, p);
    scalar_view weighted_valVY_eval("weighted_valVY_eval", c, p);
    const Real zero(0), one(1);
    Real scale(0), uij(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        uij = valUX_eval(i,j);
        scale = (uij < zero ? one : (uij > zero ? zero : one));
        weighted_minVX_eval(i,j) = scale * weight_(i,j) * valVX_eval(i,j);
        weighted_valVY_eval(i,j) = weight_(i,j) * valVY_eval(i,j);
      }
    }
    // Build local gradient of state tracking term
    fst::integrate(H[0],weighted_minVX_eval,feVel_->NdetJ(),false);
    fst::integrate(H[1],weighted_valVY_eval,feVel_->NdetJ(),false);

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Horizontal_ThermalFluids

template <class Real>
class QoI_State_ThermalFluids : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  ROL::Ptr<QoI<Real,DeviceType>> qoi_;

public:
  QoI_State_ThermalFluids(ROL::ParameterList &parlist,
                         const ROL::Ptr<fe_type> &feVel,
                         const ROL::Ptr<fe_type> &fePrs,
                         const ROL::Ptr<fe_type> &feThr,
                         const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper) {
    std::string stateObj = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj != "Vorticity" && stateObj != "Circulation" && stateObj != "Directional" )
      throw Exception::NotImplemented(">>> (QoI_State_ThermalFluids): Unknown objective type."); 
    if ( stateObj == "Vorticity" )
      qoi_ = ROL::makePtr<QoI_Vorticity_ThermalFluids<Real,DeviceType>>(feVel,fePrs,feThr,fieldHelper,parlist);
    else if ( stateObj == "Directional" )
      qoi_ = ROL::makePtr<QoI_Horizontal_ThermalFluids<Real,DeviceType>>(feVel,fePrs,feThr,fieldHelper);
    else
      qoi_ = ROL::makePtr<QoI_Circulation_ThermalFluids<Real,DeviceType>>(feVel,fePrs,feThr,fieldHelper);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    return qoi_->value(val, u_coeff, z_coeff, z_param);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    qoi_->gradient_1(grad, u_coeff, z_coeff, z_param);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    qoi_->gradient_2(grad, u_coeff, z_coeff, z_param);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    qoi_->HessVec_11(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    qoi_->HessVec_12(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    qoi_->HessVec_21(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    qoi_->HessVec_22(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

};

template <class Real>
class QoI_L2Penalty_ThermalFluids : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_, fePrs_, feThr_;
  const std::vector<std::vector<ROL::Ptr<fe_type>>> feThrBdry_;
  const std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  scalar_view getBoundaryCoeff(const scalar_view cell_coeff, int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feThr_->N().extent_int(1);
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j)
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
    }
    return bdry_coeff;
  }

public:
  QoI_L2Penalty_ThermalFluids(const ROL::Ptr<fe_type> &feVel,
                              const ROL::Ptr<fe_type> &fePrs,
                              const ROL::Ptr<fe_type> &feThr,
                              const std::vector<std::vector<ROL::Ptr<fe_type>>> &feThrBdry,
                              const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds,
                              const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), feThrBdry_(feThrBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c = feVel_->gradN().extent_int(0);
    const int d = feVel_->gradN().extent_int(3);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Initialize output val
    val = scalar_view("val", c);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ( i == 1 || i == 2 ) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts().extent_int(1);
              // Evaluate control on FE basis
              scalar_view z_coeff_bdry = getBoundaryCoeff(Z[d+1], i, j);
              scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
              // Integrate cell L2 norm squared
              scalar_view intVal("intVal", numCellsSide);
              feThrBdry_[i][j]->computeIntegral(intVal,valZ_eval,valZ_eval,false);
              // Add to integral value
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                val(cidx) += static_cast<Real>(0.5)*intVal(k);
              }
            }
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int fh = feThr_->gradN().extent_int(1);
    const int d  = feThr_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> G(d+2);
    for (int i = 0; i < d; ++i)
      G[i] = scalar_view("grad", c, fv);
    G[d] = scalar_view("grad", c, fp);
    G[d+1] = scalar_view("grad", c, fh);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ( i == 1 || i == 2 ) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts().extent_int(1);
              // Evaluate control on FE basis
              scalar_view z_coeff_bdry = getBoundaryCoeff(Z[d+1], i, j);
              scalar_view valZ_eval("valZ_eval", numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
              // Compute gradient of squared L2-norm of diff
              scalar_view intGrad("intGrad", numCellsSide, fh);
              fst::integrate(intGrad,valZ_eval,feThrBdry_[i][j]->NdetJ(),false);
              // Add to integral value
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l)
                  (G[d+1])(cidx,l) += intGrad(k,l);
              }
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int fh = feThr_->gradN().extent_int(1);
    const int d  = feThr_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> H(d+2);
    for (int i = 0; i < d; ++i)
      H[i] = scalar_view("hess", c, fv);
    H[d] = scalar_view("hess", c, fp);
    H[d+1] = scalar_view("hess", c, fh);
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ( i == 1 || i == 2 ) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts().extent_int(1);
              // Evaluate control on FE basis
              scalar_view v_coeff_bdry = getBoundaryCoeff(V[d+1], i, j);
              scalar_view valV_eval("valV_eval", numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valV_eval, v_coeff_bdry);
              // Compute gradient of squared L2-norm of diff
              scalar_view intHess("intHess", numCellsSide, fh);
              fst::integrate(intHess,valV_eval,feThrBdry_[i][j]->NdetJ(),false);
              // Add to integral value
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l)
                  (H[d+1])(cidx,l) += intHess(k,l);
              }
            }
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_L2Penalty_ThermalFluids

template <class Real>
class StdObjective_ThermalFluids : public ROL::StdObjective<Real> {
private:
  Real alpha_;
  std::string stateObj_;

public:
  StdObjective_ThermalFluids(ROL::ParameterList &parlist) {
    alpha_    = parlist.sublist("Problem").get("Control penalty parameter",1.e-4);
    stateObj_ = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj_ != "Vorticity" && stateObj_ != "Circulation" && stateObj_ != "Directional") {
      throw Exception::NotImplemented(">>> (StdObjective_ThermalFluids): Unknown objective type."); 
    }
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    Real val = alpha_*x[1];
    if ( stateObj_ == "Vorticity" || stateObj_ == "Directional" ) {
      val += x[0];
    }
    else {
      val += static_cast<Real>(0.5)*x[0]*x[0];
    }
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    const Real one(1);
    if ( stateObj_ == "Vorticity" || stateObj_ == "Directional" ) {
      g[0] = one;
    }
    else {
      g[0] = x[0];
    }
    g[1] = alpha_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    if ( stateObj_ == "Vorticity" || stateObj_ == "Directional" ) {
      hv[0] = zero;
    }
    else {
      hv[0] = v[0];
    }
    hv[1] = zero;
  }

}; // OBJ_SCALAR

#endif
