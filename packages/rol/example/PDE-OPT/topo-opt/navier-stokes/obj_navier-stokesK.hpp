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

#ifndef PDEOPT_QOI_NAVIERSTOKESK_HPP
#define PDEOPT_QOI_NAVIERSTOKESK_HPP

#include "../../TOOLS/qoiK.hpp"
#include "pde_navier-stokesK.hpp"
#include "impermiabilityK.hpp"

template <class Real, class DeviceType>
class QoI_Power_NavierStokes : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_, fePrs_, feCtrl_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_, fieldInfoCtrl_;

  ROL::Ptr<Impermiability<Real,DeviceType>> imp_;
  Real viscosity_;

public:
  QoI_Power_NavierStokes(ROL::ParameterList                          &list,
                         const ROL::Ptr<fe_type>                     &feVel,
                         const ROL::Ptr<fe_type>                     &fePrs,
                         const ROL::Ptr<fe_type>                     &feCtrl,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfoCtrl)
    : feVel_(feVel), fePrs_(fePrs), feCtrl_(feCtrl), fieldInfo_(fieldInfo), fieldInfoCtrl_(fieldInfoCtrl) {
    viscosity_ = list.sublist("Problem").get("Viscosity",5e-3);
    imp_ = ROL::makePtr<Impermiability<Real,DeviceType>>(list);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = feVel_->gradN().extent_int(0);
    const int p = feVel_->gradN().extent_int(2);
    const int d = feVel_->gradN().extent_int(3);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate on FE basis
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Uval("Uval", c, p);
    scalar_view aUval("aUval", c, p);
    scalar_view Ugrad("Ugrad", c, p, d);
    scalar_view nUgrad("nUgrad", c, p, d);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 0); 
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      feVel_->evaluateGradient(Ugrad, U[i]);
      // Scale
      fst::scalarMultiplyDataData(aUval, alpha, Uval);
      rst::scale(nUgrad, Ugrad, viscosity_);
      // Integrate
      feVel_->computeIntegral(val, nUgrad, Ugrad, true);
      feVel_->computeIntegral(val, aUval, Uval, true);
    }
    rst::scale(val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> G;
    G.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      G[i] = scalar_view("grad1", c, fv);
    }
    G[d] = scalar_view("grad1", c, fp);
    // Get components of the control
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Uval("Uval", c, p);
    scalar_view aUval("aUval", c, p);
    scalar_view Ugrad("Ugrad", c, p, d);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 0); 
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      feVel_->evaluateGradient(Ugrad, U[i]);
      // Scale
      fst::scalarMultiplyDataData(aUval, alpha, Uval);
      rst::scale(Ugrad, viscosity_);
      // Integrate
      fst::integrate(G[i],Ugrad,feVel_->gradNdetJ(),false);
      fst::integrate(G[i],aUval,feVel_->NdetJ(),true);
    }
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = feVel_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> G;
    G.resize(fieldInfoCtrl_->numFields);
    G[0] = scalar_view("grad2", c, fc);
    // Get components of the control
    std::vector<scalar_view> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Uval("Uval", c, p);
    scalar_view aUval("aUval", c, p);
    scalar_view aU2val("aU2val", c, p);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 1); 
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      fst::scalarMultiplyDataData(aUval, alpha, Uval);
      fst::scalarMultiplyDataData(aU2val, aUval, Uval);
      // Integrate
      fst::integrate(G[0],aU2val,feCtrl_->NdetJ(),true);
    }
    rst::scale(G[0], static_cast<Real>(0.5));
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfoCtrl_);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      H[i] = scalar_view("hess11", c, fv);
    }
    H[d] = scalar_view("hess11", c, fp);
    // Get components of the control
    std::vector<scalar_view> V, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Vval("Vval", c, p);
    scalar_view aVval("aVval", c, p);
    scalar_view Vgrad("Vgrad", c, p, d);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 0); 
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Vval, V[i]);
      feVel_->evaluateGradient(Vgrad, V[i]);
      // Scale
      fst::scalarMultiplyDataData(aVval, alpha, Vval);
      rst::scale(Vgrad, viscosity_);
      // Integrate
      fst::integrate(H[i],Vgrad,feVel_->gradNdetJ(),false);
      fst::integrate(H[i],aVval,feVel_->NdetJ(),true);
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      H[i] = scalar_view("hess12", c, fv);
    }
    H[d] = scalar_view("hess12", c, fp);
    // Get components of the control
    std::vector<scalar_view> V, U, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfoCtrl_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Uval("Uval", c, p);
    scalar_view aUval("aUval", c, p);
    scalar_view Vval("Vval", c, p);
    scalar_view aVval("aVval", c, p);
    feCtrl_->evaluateValue(Vval, V[0]);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 1); 
    fst::scalarMultiplyDataData(aVval, alpha, Vval);
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      fst::scalarMultiplyDataData(aUval, aVval, Uval);
      // Integrate
      fst::integrate(H[i],aUval,feVel_->NdetJ(),false);
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = feVel_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> H;
    H.resize(fieldInfoCtrl_->numFields);
    H[0] = scalar_view("hess21", c, fc);
    // Get components of the control
    std::vector<scalar_view> V, U, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Uval("Uval", c, p);
    scalar_view aUval("aUval", c, p);
    scalar_view aUVval("aUVval", c, p);
    scalar_view Vval("Vval", c, p);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 1); 
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Vval, V[i]);
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      fst::scalarMultiplyDataData(aUval, alpha, Uval);
      fst::scalarMultiplyDataData(aUVval, aUval, Vval);
      // Integrate
      fst::integrate(H[0],aUVval,feCtrl_->NdetJ(), true);
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_);
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c  = feVel_->gradN().extent_int(0);
    const int fc = feCtrl_->gradN().extent_int(1);
    const int p  = feVel_->gradN().extent_int(2);
    const int d  = feVel_->gradN().extent_int(3);
    // Initialize output grad
    std::vector<scalar_view> H;
    H.resize(fieldInfoCtrl_->numFields);
    H[0] = scalar_view("hess22", c, fc);
    // Get components of the control
    std::vector<scalar_view> V, U, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfoCtrl_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    scalar_view Zval("Zval", c, p);
    scalar_view alpha("alpha", c, p);
    scalar_view Uval("Uval", c, p);
    scalar_view aUval("aUval", c, p);
    scalar_view aU2val("aU2val", c, p);
    scalar_view Vval("Vval", c, p);
    scalar_view aVval("aVval", c, p);
    feCtrl_->evaluateValue(Vval, V[0]);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 2); 
    fst::scalarMultiplyDataData(aVval, alpha, Vval);
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      fst::scalarMultiplyDataData(aUval, aVval, Uval);
      fst::scalarMultiplyDataData(aU2val, aUval, Uval);
      // Integrate
      fst::integrate(H[0],aU2val,feCtrl_->NdetJ(),true);
    }
    rst::scale(H[0], static_cast<Real>(0.5));
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_);
  }

}; // QoI_Power_NavierStokes

template <class Real, class DeviceType>
class QoI_Volume_NavierStokes : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> fe_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;

  scalar_view weight_;
  Real volFraction_;

public:
  QoI_Volume_NavierStokes(ROL::ParameterList &list,
                          const ROL::Ptr<fe_type> &fe,
                          const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : fe_(fe), fieldInfo_(fieldInfo) {
    const Real one(1);
    volFraction_   = list.sublist("Problem").get("Volume Fraction",0.5);
    const int c = fe_->cubPts().extent_int(0);
    const int p = fe_->cubPts().extent_int(1);
    weight_ = scalar_view("weight_", c, p);
    Kokkos::deep_copy(weight_,one);
  }

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->cubPts().extent_int(0);
    const int p = fe_->cubPts().extent_int(1);
    // Initialize output val
    val = scalar_view("val", c);
    // Get components of the control
    std::vector<scalar_view> Z;
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfo_);
    // Evaluate on FE basis
    scalar_view Z0("Z0", c, p);
    fe_->evaluateValue(Z0, Z[0]);
    // Integrate the density minus volume fraction
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        Z0(i,j) -= volFraction_;
      }
    }
    fe_->computeIntegral(val,weight_,Z0,true);
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    // Get relevant dimensions
    const int c = fe_->cubPts().extent_int(0);
    const int f = fe_->N().extent_int(1);
    // Initialize output grad
    std::vector<scalar_view> G;
    G.resize(fieldInfo_->numFields);
    G[0] = scalar_view("grad2", c, f);
    // Integrate density
    fst::integrate(G[0],weight_,fe_->NdetJ(),false);
    FieldUtils::combineFieldCoeff(grad, G, fieldInfo_);
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Volume_NavierStokes

#endif
