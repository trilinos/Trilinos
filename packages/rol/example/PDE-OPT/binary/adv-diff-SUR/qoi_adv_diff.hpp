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

#ifndef PDEOPT_QOI_ADV_DIFF_SUR_HPP
#define PDEOPT_QOI_ADV_DIFF_SUR_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_adv_diff.hpp"

template <class Real>
class QoI_State_Cost_adv_diff : public QoI<Real> {
private:
  ROL::Ptr<FE<Real>> fe_;

public:
  QoI_State_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe) : fe_(fe) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,valU_eval,valU_eval);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    // Evaluate state on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *valU_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::gradient_3 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *valV_eval,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_12 is zero.");
  }

  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_22 is zero.");
  }

  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_33 is zero.");
  }

}; // QoI_State_Cost

template <class Real>
class QoI_Control_Cost_adv_diff : public QoI<Real> {
private:
  ROL::Ptr<FE<Real>> fe_;
  Real vol_;

public:
  QoI_Control_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                            ROL::ParameterList & parlist) : fe_(fe) {
    int order = parlist.sublist("Problem").get("Hilbert Curve Order", 2);
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    int n = std::pow(2,order);
    vol_ = (XU-XL)/static_cast<Real>(n) * (YU-YL)/static_cast<Real>(n);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    if (z_coeff != ROL::nullPtr) {
      // Evaluate state on FE basis
      ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, ones;
      valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      ones      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valZ_eval, z_coeff);
      ones->initialize(static_cast<Real>(1));
      // Compute squared L1-norm of diff
      fe_->computeIntegral(val,valZ_eval,ones);
    }
    else {
      const int size = z_param->size();
      Real sum(0);
      for (int i = 0; i < size; ++i) {
        sum += (*z_param)[i];
      }
      val = ROL::nullPtr;
      return vol_*sum;
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int f = fe_->gradN()->dimension(1);
      const int p = fe_->gradN()->dimension(2);
      // Initialize output grad
      grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Evaluate state on FE basis
      ROL::Ptr<Intrepid::FieldContainer<Real>> ones;
      ones = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      ones->initialize(static_cast<Real>(1));
      // Compute gradient of L1-norm
      Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                    *ones,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
    }
    else {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_2 is zero.");
    }
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_3 is zero.");
    }
    else {
      const int size = z_param->size();
      std::vector<Real> g(size,vol_);
      return g;
    }
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_12 is zero.");
  }

  void HessVec_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_22 is zero.");
  }

  void HessVec_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_33 is zero.");
  }

}; // QoI_Control_Cost

template <class Real>
class QoI_Control_Cost_L2_adv_diff : public QoI<Real> {
private:
  ROL::Ptr<FE<Real>> fe_;
  Real vol_;

public:
  QoI_Control_Cost_L2_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                               ROL::ParameterList & parlist) : fe_(fe) {
    int order = parlist.sublist("Problem").get("Hilbert Curve Order", 2);
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    int n = std::pow(2,order);
    vol_ = (XU-XL)/static_cast<Real>(n) * (YU-YL)/static_cast<Real>(n);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int p = fe_->gradN()->dimension(2);
      // Initialize output val
      val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Evaluate state on FE basis
      ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval;
      valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valZ_eval, z_coeff);
      // Compute squared L2-norm of diff
      fe_->computeIntegral(val,valZ_eval,valZ_eval);
      Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    }
    else {
      const int size = z_param->size();
      Real sum(0);
      for (int i = 0; i < size; ++i) {
        sum += std::pow((*z_param)[i],static_cast<Real>(2));
      }
      val = ROL::nullPtr;
      return static_cast<Real>(0.5)*vol_*sum;
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int f = fe_->gradN()->dimension(1);
      const int p = fe_->gradN()->dimension(2);
      // Initialize output grad
      grad = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      // Evaluate state on FE basis
      ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval;
      valZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      fe_->evaluateValue(valZ_eval, z_coeff);
      // Compute gradient of L2-norm
      Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                    *valZ_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
    }
    else {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_2 is zero.");
    }
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_3 is zero.");
    }
    else {
      const int size = z_param->size();
      std::vector<Real> g(size);
      for (int i = 0; i < size; ++i) {
        g[i] = vol_ * (*z_param)[i];
      }
      return g;
    }
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_12 is zero.");
  }

  void HessVec_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      const int c = fe_->gradN()->dimension(0);
      const int f = fe_->gradN()->dimension(1);
      const int p = fe_->gradN()->dimension(2);
      ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval =
        ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      hess = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      fe_->evaluateValue(valV_eval, v_coeff);
      Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                    *valV_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
    }
    else {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_22 is zero.");
    }
  }

  void HessVec_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_coeff != ROL::nullPtr) {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::HessVec_33 is zero.");
    }
    else {
      const int size = z_param->size();
      std::vector<Real> h(size);
      for (int i = 0; i < size; ++i) {
        h[i] = vol_ * (*v_param)[i];
      }
      return h;
    }
  }

}; // QoI_Control_Cost

#endif
