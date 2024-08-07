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

#ifndef PDEOPT_QOI_NAVIERSTOKES_HPP
#define PDEOPT_QOI_NAVIERSTOKES_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_navier-stokes.hpp"
#include "impermiability.hpp"

template <class Real>
class QoI_Power_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>>                    feVel_;
  const ROL::Ptr<FE<Real>>                    fePrs_;
  const ROL::Ptr<FE<Real>>                    feCtrl_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfoCtrl_;

  ROL::Ptr<Impermiability<Real>> imp_;
  Real viscosity_;

public:
  QoI_Power_NavierStokes(Teuchos::ParameterList                      &list,
                         const ROL::Ptr<FE<Real>>                    &feVel,
                         const ROL::Ptr<FE<Real>>                    &fePrs,
                         const ROL::Ptr<FE<Real>>                    &feCtrl,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfoCtrl)
    : feVel_(feVel), fePrs_(fePrs), feCtrl_(feCtrl), fieldInfo_(fieldInfo), fieldInfoCtrl_(fieldInfoCtrl) {
    viscosity_ = list.sublist("Problem").get("Viscosity",5e-3);
    imp_ = ROL::makePtr<Impermiability<Real>>(list);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = feVel_->gradN()->dimension(0);
    const int p = feVel_->gradN()->dimension(2);
    const int d = feVel_->gradN()->dimension(3);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> alpha, Zval, Uval, aUval, Ugrad, nUgrad;
    Zval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Ugrad  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    nUgrad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 0); 
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      Ugrad->initialize();
      aUval->initialize();
      nUgrad->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      feVel_->evaluateGradient(Ugrad, U[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUval, *alpha, *Uval);
      Intrepid::RealSpaceTools<Real>::scale(*nUgrad, *Ugrad, viscosity_);
      // Integrate
      feVel_->computeIntegral(val, nUgrad, Ugrad, true);
      feVel_->computeIntegral(val, aUval, Uval, true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val, static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> Zval, alpha, Uval, aUval, Ugrad;
    Zval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Ugrad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 0); 
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      Ugrad->initialize();
      aUval->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      feVel_->evaluateGradient(Ugrad, U[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUval, *alpha, *Uval);
      Intrepid::RealSpaceTools<Real>::scale(*Ugrad, viscosity_);
      // Integrate
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *Ugrad,
                                                    *(feVel_->gradNdetJ()),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *aUval,
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfoCtrl_->numFields);
    G[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, Z;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> Zval, alpha, Uval, aUval, aU2val;
    Zval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aU2val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 1); 
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      aUval->initialize();
      aU2val->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUval, *alpha, *Uval);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aU2val, *aUval, *Uval);
      // Integrate
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                    *aU2val,
                                                    *(feCtrl_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*G[0], static_cast<Real>(0.5));
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfoCtrl_);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> Zval, alpha, Vval, Vgrad, aVval;
    Zval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Vval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aVval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Vgrad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 0); 
    for (int i = 0; i < d; ++i) {
      Vval->initialize();
      Vgrad->initialize();
      aVval->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Vval, V[i]);
      feVel_->evaluateGradient(Vgrad, V[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aVval, *alpha, *Vval);
      Intrepid::RealSpaceTools<Real>::scale(*Vgrad, viscosity_);
      // Integrate
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *Vgrad,
                                                    *(feVel_->gradNdetJ()),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *aVval,
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V, U, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfoCtrl_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> Zval, alpha, Uval, aUval, Vval, aVval;
    Zval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Vval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aVval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCtrl_->evaluateValue(Vval, V[0]);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 1); 
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aVval, *alpha, *Vval);
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      aUval->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUval, *aVval, *Uval);
      // Integrate
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *aUval,
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfoCtrl_->numFields);
    H[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V, U, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> Zval, alpha, Uval, aUval, aUVval, Vval;
    Zval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUVval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Vval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 1); 
    for (int i = 0; i < d; ++i) {
      Vval->initialize();
      Uval->initialize();
      aUval->initialize();
      aUVval->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Vval, V[i]);
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUval, *alpha, *Uval);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUVval, *aUval, *Vval);
      // Integrate
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                    *aUVval,
                                                    *(feCtrl_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_);
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fc = feCtrl_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfoCtrl_->numFields);
    H[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fc);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V, U, Z;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfoCtrl_);
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfoCtrl_);
    // Create storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> Zval, alpha, Uval, aUval, aU2val, Vval, aVval;
    Zval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    alpha  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aUval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aU2val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Vval   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    aVval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feCtrl_->evaluateValue(Vval, V[0]);
    feCtrl_->evaluateValue(Zval, Z[0]);
    imp_->compute(alpha, Zval, feVel_->cubPts(), 2); 
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aVval, *alpha, *Vval);
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      aUval->initialize();
      aU2val->initialize();
      // Evaluate on FE basis
      feVel_->evaluateValue(Uval, U[i]);
      // Scale
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aUval, *aVval, *Uval);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*aU2val, *aUval, *Uval);
      // Integrate
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                    *aU2val,
                                                    *(feCtrl_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*H[0], static_cast<Real>(0.5));
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfoCtrl_);
  }

}; // QoI_Power_NavierStokes

template <class Real>
class QoI_Volume_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;
  Real volFraction_;

public:
  QoI_Volume_NavierStokes(ROL::ParameterList &list,
                          const ROL::Ptr<FE<Real>> &fe,
                          const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : fe_(fe), fieldInfo_(fieldInfo) {
    const Real one(1);
    volFraction_   = list.sublist("Problem").get("Volume Fraction",0.5);
    const int c = fe_->cubPts()->dimension(0);
    const int p = fe_->cubPts()->dimension(1);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    weight_->initialize(one);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->cubPts()->dimension(0);
    const int p = fe_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    FieldUtils::splitFieldCoeff<Real>(Z, z_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> Z0;
    Z0  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(Z0, Z[0]);
    // Integrate the density minus volume fraction
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*Z0)(i,j) -= volFraction_;
      }
    }
    fe_->computeIntegral(val,weight_,Z0,true);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->cubPts()->dimension(0);
    const int f = fe_->N()->dimension(1);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    G[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Integrate density
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *weight_,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    FieldUtils::combineFieldCoeff(grad, G, fieldInfo_);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Volume_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Volume_NavierStokes

#endif
