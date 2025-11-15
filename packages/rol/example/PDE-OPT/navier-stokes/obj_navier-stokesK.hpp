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

#include "../TOOLS/qoiK.hpp"
#include "pde_navier-stokesK.hpp"

template <class Real, class DeviceType>
class QoI_Vorticity_NavierStokes : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_, fePrs_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  scalar_view weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return (((x[1] <= 0.5+eps_)&&(x[0]>= 1.0-eps_)&&(x[0] <= 4.0+eps_)) ?
                static_cast<Real>(1) : static_cast<Real>(0));
  }

public:
  QoI_Vorticity_NavierStokes(const ROL::Ptr<fe_type> &feVel,
                             const ROL::Ptr<fe_type> &fePrs,
                             const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
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
        curlU_eval(i,j)   = gradUY_eval(i,j,0) - gradUX_eval(i,j,1);
    }
    // Multiply by weight
    scalar_view weighted_curlU_eval("curlU_eval", c, p);
    fst::scalarMultiplyDataData(weighted_curlU_eval,weight_,curlU_eval);
    // Compute L2 norm squared
    feVel_->computeIntegral(val,curlU_eval,weighted_curlU_eval,false);
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
    int p = feVel_->cubPts().extent_int(1);
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    scalar_view velUX_grad("velUX_grad", c, fv);
    scalar_view velUY_grad("velUY_grad", c, fv);
    scalar_view presU_grad("presU_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velUX_grad;
    G[1] = velUY_grad;
    G[2] = presU_grad;
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
        curlU_eval(i,j) = weight_(i,j)*(gradUY_eval(i,j,0) - gradUX_eval(i,j,1));
    }
    // Build local gradient of state tracking term
    fst::integrate(velUX_grad,curlU_eval,feVel_->DNDdetJ(1),false);
    rst::scale(velUX_grad,static_cast<Real>(-1));
    fst::integrate(velUY_grad,curlU_eval,feVel_->DNDdetJ(0),false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::gradient_2 is zero.");
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
    int d = feVel_->cubPts().extent_int(2);
    // Initialize output grad
    scalar_view velVX_grad("velVX_grad", c, fv);
    scalar_view velVY_grad("velVY_grad", c, fv);
    scalar_view presV_grad("presV_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velVX_grad;
    G[1] = velVY_grad;
    G[2] = presV_grad;
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    scalar_view gradVX_eval("gradVX_eval", c, p, d);
    scalar_view gradVY_eval("gradVY_eval", c, p, d);
    feVel_->evaluateGradient(gradVX_eval, V[0]);
    feVel_->evaluateGradient(gradVY_eval, V[1]);
    // Compute curl
    scalar_view curlV_eval("curlV_eval", c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        curlV_eval(i,j) = weight_(i,j)*(gradVY_eval(i,j,0) - gradVX_eval(i,j,1));
      }
    }
    // Build local gradient of state tracking term
    fst::integrate(velVX_grad,curlV_eval,feVel_->DNDdetJ(1),false);
    rst::scale(velVX_grad,static_cast<Real>(-1));
    fst::integrate(velVY_grad,curlV_eval,feVel_->DNDdetJ(0),false);

    fieldHelper_->combineFieldCoeff(hess, G);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Vorticity_NavierStokes

template <class Real, class DeviceType>
class QoI_Circulation_NavierStokes : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_;
  const ROL::Ptr<fe_type> fePrs_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  scalar_view weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return (((x[1] <= 1.0+eps_)&&(x[0]>= 1.0-eps_)&&(x[0] <= 4.0+eps_)) ?
                static_cast<Real>(1) : static_cast<Real>(0));
  }

public:
  QoI_Circulation_NavierStokes(const ROL::Ptr<fe_type> &feVel,
                               const ROL::Ptr<fe_type> &fePrs,
                               const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
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
    val = scalar_view("val",c);
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
    // Initialize output grad
    scalar_view velUX_grad("velUX_grad", c, fv);
    scalar_view velUY_grad("velUY_grad", c, fv);
    scalar_view presU_grad("presU_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velUX_grad;
    G[1] = velUY_grad;
    G[2] = presU_grad;
    // Build local gradient of state tracking term
    fst::integrate(velUX_grad,weight_,feVel_->DNDdetJ(1),false);
    rst::scale(velUX_grad,static_cast<Real>(-1));
    fst::integrate(velUY_grad,weight_,feVel_->DNDdetJ(0),false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Circulation_NavierStokes

template <class Real, class DeviceType>
class QoI_Horizontal_NavierStokes : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_;
  const ROL::Ptr<fe_type> fePrs_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;
  scalar_view weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return (((x[1] <= 0.5+eps_)&&(x[0]>= 1.0-eps_)&&(x[0] <= 4.0+eps_)) ?
                static_cast<Real>(1) : static_cast<Real>(0));
  }

public:
  QoI_Horizontal_NavierStokes(const ROL::Ptr<fe_type> &feVel,
                              const ROL::Ptr<fe_type> &fePrs,
                              const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
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
    scalar_view weighted_minUX_eval("minUX_eval", c, p);
    fst::scalarMultiplyDataData(weighted_minUX_eval,weight_,minUX_eval);
    scalar_view weighted_valUY_eval("valUY_eval", c, p);
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
    int p = feVel_->cubPts().extent_int(1);
    // Initialize output grad
    scalar_view velUX_grad("velUX_grad", c, fv);
    scalar_view velUY_grad("velUY_grad", c, fv);
    scalar_view presU_grad("presU_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velUX_grad;
    G[1] = velUY_grad;
    G[2] = presU_grad;
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
        weighted_minUX_eval(i,j)
          = weight_(i,j) * std::min(static_cast<Real>(0),valUX_eval(i,j));
        weighted_valUY_eval(i,j)
          = weight_(i,j) * valUY_eval(i,j);
      }
    }
    // Build local gradient of state tracking term
    fst::integrate(velUX_grad,weighted_minUX_eval,feVel_->NdetJ(),false);
    fst::integrate(velUY_grad,weighted_valUY_eval,feVel_->NdetJ(),false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::gradient_2 is zero.");
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
    // Initialize output grad
    scalar_view velVX_grad("velVX_grad", c, fv);
    scalar_view velVY_grad("velVY_grad", c, fv);
    scalar_view presV_grad("presV_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velVX_grad;
    G[1] = velVY_grad;
    G[2] = presV_grad;
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
    Real scale(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if ( valUX_eval(i,j) < static_cast<Real>(0) ) {
          scale = static_cast<Real>(1);
        }
        else if ( valUX_eval(i,j) > static_cast<Real>(0) ) {
          scale = static_cast<Real>(0);
        }
        else {
          //scale = static_cast<Real>(0);
          //scale = static_cast<Real>(0.5);
          scale = static_cast<Real>(1);
        }
        weighted_minVX_eval(i,j) = scale * weight_(i,j) * valVX_eval(i,j);
        weighted_valVY_eval(i,j) = weight_(i,j) * valVY_eval(i,j);
      }
    }
    // Build local gradient of state tracking term
    fst::integrate(velVX_grad,weighted_minVX_eval,feVel_->NdetJ(),false);
    fst::integrate(velVY_grad,weighted_valVY_eval,feVel_->NdetJ(),false);

    fieldHelper_->combineFieldCoeff(hess, G);
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Horizontal_NavierStokes

template <class Real, class DeviceType>
class QoI_State_NavierStokes : public QoI<Real,DeviceType> {
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
  QoI_State_NavierStokes(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<fe_type> &feVel,
                         const ROL::Ptr<fe_type> &fePrs,
                         const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper) {
    std::string stateObj = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj != "Vorticity" && stateObj != "Circulation" && stateObj != "Directional" ) {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Unknown objective type."); 
    }
    if ( stateObj == "Vorticity" ) {
      qoi_ = ROL::makePtr<QoI_Vorticity_NavierStokes<Real,DeviceType>>(feVel,fePrs,fieldHelper);
    }
    else if ( stateObj == "Directional" ) {
      qoi_ = ROL::makePtr<QoI_Horizontal_NavierStokes<Real,DeviceType>>(feVel,fePrs,fieldHelper);
    }
    else {
      qoi_ = ROL::makePtr<QoI_Circulation_NavierStokes<Real,DeviceType>>(feVel,fePrs,fieldHelper);
    }
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

template <class Real, class DeviceType>
class QoI_L2Penalty_NavierStokes : public QoI<Real,DeviceType> {
public:
  using scalar_view = Kokkos::DynRankView<Real, DeviceType>;
  using basis_ptr   = Intrepid2::BasisPtr<DeviceType, Real, Real>;
  using fe_type     = FE<Real,DeviceType>;
  using fst         = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct          = Intrepid2::CellTools<DeviceType>;
  using rst         = Intrepid2::RealSpaceTools<DeviceType>;

private:
  const ROL::Ptr<fe_type> feVel_;
  const ROL::Ptr<fe_type> fePrs_;
  const std::vector<ROL::Ptr<fe_type>> feVelBdry_;
  const std::vector<std::vector<int>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real,DeviceType>> fieldHelper_;

  scalar_view getBoundaryCoeff(
      const scalar_view & cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feVel_->N().extent_int(1);
    
    scalar_view bdry_coeff("bdry_coeff", numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        bdry_coeff(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_L2Penalty_NavierStokes(const ROL::Ptr<fe_type> &feVel,
                             const ROL::Ptr<fe_type> &fePrs,
                             const std::vector<ROL::Ptr<fe_type>> &feVelBdry,
                             const std::vector<std::vector<int>> &bdryCellLocIds,
                             const ROL::Ptr<FieldHelper<Real,DeviceType>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feVelBdry_(feVelBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {}

  Real value(scalar_view & val,
             const scalar_view u_coeff,
             const scalar_view z_coeff = scalar_view(),
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c = feVel_->gradN().extent_int(0);
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Initialize output val
    val = scalar_view("val", c);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts().extent_int(1);
        // Evaluate x-component of control on FE basis
        scalar_view zx_coeff_bdry = getBoundaryCoeff(Z[0], l);
        scalar_view valZX_eval("valZX_eval", numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZX_eval, zx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        scalar_view zy_coeff_bdry = getBoundaryCoeff(Z[1], l);
        scalar_view valZY_eval("valZY_eval", numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZY_eval, zy_coeff_bdry);
        // Integrate cell L2 norm squared
        scalar_view intVal("intVal", numCellsSide);
        feVelBdry_[l]->computeIntegral(intVal,valZX_eval,valZX_eval,false);
        feVelBdry_[l]->computeIntegral(intVal,valZY_eval,valZY_eval,true);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          val(cidx) += static_cast<Real>(0.5)*intVal(i);
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::gradient_1 is zero.");
  }

  void gradient_2(scalar_view & grad,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    // Initialize output grad
    scalar_view velZX_grad("velZX_grad", c, fv);
    scalar_view velZY_grad("velZY_grad", c, fv);
    scalar_view presZ_grad("presZ_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velZX_grad;
    G[1] = velZY_grad;
    G[2] = presZ_grad;
    // Get components of the control
    std::vector<scalar_view> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts().extent_int(1);
        // Evaluate x-component of control on FE basis
        scalar_view zx_coeff_bdry = getBoundaryCoeff(Z[0], l);
        scalar_view valZX_eval("valZX_eval", numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZX_eval, zx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        scalar_view zy_coeff_bdry = getBoundaryCoeff(Z[1], l);
        scalar_view valZY_eval("valZY_eval", numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZY_eval, zy_coeff_bdry);
        // Compute gradient of squared L2-norm of diff
        scalar_view intGradX("intGradX", numCellsSide, fv);
        fst::integrate(intGradX,valZX_eval,feVelBdry_[l]->NdetJ(),false);
        scalar_view intGradY("intGradY", numCellsSide, fv);
        fst::integrate(intGradY,valZY_eval,feVelBdry_[l]->NdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fv; ++j) {
            (G[0])(cidx,j) += intGradX(i,j);
            (G[1])(cidx,j) += intGradY(i,j);
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
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(scalar_view & hess,
                  const scalar_view v_coeff,
                  const scalar_view u_coeff,
                  const scalar_view z_coeff = scalar_view(),
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) override {
    const int c  = feVel_->gradN().extent_int(0);
    const int fv = feVel_->gradN().extent_int(1);
    const int fp = fePrs_->gradN().extent_int(1);
    // Initialize output grad
    scalar_view velVX_grad("velVX_grad", c, fv);
    scalar_view velVY_grad("velVY_grad", c, fv);
    scalar_view presV_grad("presV_grad", c, fp);
    std::vector<scalar_view> G;
    G.resize(fieldHelper_->numFields());
    G[0] = velVX_grad;
    G[1] = velVY_grad;
    G[2] = presV_grad;
    // Get components of the control
    std::vector<scalar_view> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts().extent_int(1);
        // Evaluate x-component of control on FE basis
        scalar_view vx_coeff_bdry = getBoundaryCoeff(V[0], l);
        scalar_view valVX_eval("valVX_eval", numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valVX_eval, vx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        scalar_view vy_coeff_bdry = getBoundaryCoeff(V[1], l);
        scalar_view valVY_eval("valVY_eval", numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valVY_eval, vy_coeff_bdry);
        // Compute gradient of squared L2-norm of diff
        scalar_view intHessX("intHessX", numCellsSide, fv);
        fst::integrate(intHessX,valVX_eval,feVelBdry_[l]->NdetJ(),false);
        scalar_view intHessY("intHessY", numCellsSide, fv);
        fst::integrate(intHessY,valVY_eval,feVelBdry_[l]->NdetJ(),false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fv; ++j) {
            (G[0])(cidx,j) += intHessX(i,j);
            (G[1])(cidx,j) += intHessY(i,j);
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, G);
  }

}; // QoI_L2Penalty_NavierStokes

template <class Real>
class StdObjective_NavierStokes : public ROL::StdObjective<Real> {
private:
  Real alpha_;
  std::string stateObj_;

public:
  StdObjective_NavierStokes(Teuchos::ParameterList &parlist) {
    alpha_    = parlist.sublist("Problem").get("Control penalty parameter",1.e-4);
    stateObj_ = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj_ != "Vorticity" && stateObj_ != "Circulation" && stateObj_ != "Directional") {
      throw Exception::NotImplemented(">>> (StdObjective_NavierStokes): Unknown objective type."); 
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
