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

#ifndef PDEOPT_QOI_DYNAMIC_NAVIERSTOKES_HPP
#define PDEOPT_QOI_DYNAMIC_NAVIERSTOKES_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_navier-stokes.hpp"

template <class Real>
class QoI_Vorticity_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Vorticity_NavierStokes(const ROL::Ptr<FE<Real>> &feVel,
                             const ROL::Ptr<FE<Real>> &fePrs,
                             const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    ROL::Ptr<Intrepid::FieldContainer<Real>> curlU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlU_eval)(i,j)   = (*gradUY_eval)(i,j,0) - (*gradUX_eval)(i,j,1);
      }
    }
    // Multiply by weight
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_curlU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_curlU_eval,
                                                               *weight_,
                                                               *curlU_eval);
    // Compute L2 norm squared
    feVel_->computeIntegral(val,curlU_eval,weighted_curlU_eval,false);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = feVel_->cubPts()->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velUX_grad);
    G[1] = ROL::makePtrFromRef(velUY_grad);
    G[2] = ROL::makePtrFromRef(presU_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    ROL::Ptr<Intrepid::FieldContainer<Real>> curlU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlU_eval)(i,j) = (*weight_)(i,j)*((*gradUY_eval)(i,j,0) - (*gradUX_eval)(i,j,1));
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velUX_grad,
                                                  *curlU_eval,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(velUX_grad,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(velUY_grad,
                                                  *curlU_eval,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = feVel_->cubPts()->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int d  = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradVX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradVY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feVel_->evaluateGradient(gradVX_eval, V[0]);
    feVel_->evaluateGradient(gradVY_eval, V[1]);
    // Compute curl
    ROL::Ptr<Intrepid::FieldContainer<Real>> curlV_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlV_eval)(i,j) = (*weight_)(i,j)*((*gradVY_eval)(i,j,0) - (*gradVX_eval)(i,j,1));
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velVX_grad,
                                                  *curlV_eval,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(velVX_grad,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(velVY_grad,
                                                  *curlV_eval,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(hess, G);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Vorticity_NavierStokes

template <class Real>
class QoI_Circulation_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Circulation_NavierStokes(const ROL::Ptr<FE<Real>> &feVel,
                               const ROL::Ptr<FE<Real>> &fePrs,
                               const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> gradUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    ROL::Ptr<Intrepid::FieldContainer<Real>> curlU_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*curlU_eval)(i,j)   = (*gradUY_eval)(i,j,0) - (*gradUX_eval)(i,j,1);
      }
    }
    // Compute circulation
    feVel_->computeIntegral(val,curlU_eval,weight_,false);
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = feVel_->cubPts()->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velUX_grad);
    G[1] = ROL::makePtrFromRef(velUY_grad);
    G[2] = ROL::makePtrFromRef(presU_grad);
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velUX_grad,
                                                  *weight_,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(velUX_grad,static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(velUY_grad,
                                                  *weight_,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Circulation_NavierStokes

template <class Real>
class QoI_Horizontal_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Horizontal_NavierStokes(const ROL::Ptr<FE<Real>> &feVel,
                              const ROL::Ptr<FE<Real>> &fePrs,
                              const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Compute negative part of x-velocity
    ROL::Ptr<Intrepid::FieldContainer<Real>> minUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*minUX_eval)(i,j) = std::min(static_cast<Real>(0),(*valUX_eval)(i,j));
      }
    }
    // Multiply by weight
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_minUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_minUX_eval,
                                                               *weight_,
                                                               *minUX_eval);
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_valUY_eval,
                                                               *weight_,
                                                               *valUY_eval);
    // Compute L2 norm squared
    feVel_->computeIntegral(val,minUX_eval,weighted_minUX_eval,false);
    feVel_->computeIntegral(val,valUY_eval,weighted_valUY_eval,true);

    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = feVel_->cubPts()->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int p  = feVel_->cubPts()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velUX_grad);
    G[1] = ROL::makePtrFromRef(velUY_grad);
    G[2] = ROL::makePtrFromRef(presU_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Compute negative part of x-velocity
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_minUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*weighted_minUX_eval)(i,j)
          = (*weight_)(i,j) * std::min(static_cast<Real>(0),(*valUX_eval)(i,j));
        (*weighted_valUY_eval)(i,j)
          = (*weight_)(i,j) * (*valUY_eval)(i,j);
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velUX_grad,
                                                  *weighted_minUX_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velUY_grad,
                                                  *weighted_valUY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = feVel_->cubPts()->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> valUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valVX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valVY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valVX_eval, V[0]);
    feVel_->evaluateValue(valVY_eval, V[1]);
    // Compute negative part of x-velocity
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_minVX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real>> weighted_valVY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Real scale(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if ( (*valUX_eval)(i,j) < static_cast<Real>(0) ) {
          scale = static_cast<Real>(1);
        }
        else if ( (*valUX_eval)(i,j) > static_cast<Real>(0) ) {
          scale = static_cast<Real>(0);
        }
        else {
          //scale = static_cast<Real>(0);
          //scale = static_cast<Real>(0.5);
          scale = static_cast<Real>(1);
        }
        (*weighted_minVX_eval)(i,j)
          = scale * (*weight_)(i,j) * (*valVX_eval)(i,j);
        (*weighted_valVY_eval)(i,j)
          = (*weight_)(i,j) * (*valVY_eval)(i,j);
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velVX_grad,
                                                  *weighted_minVX_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velVY_grad,
                                                  *weighted_valVY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(hess, G);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Horizontal_NavierStokes

template <class Real>
class QoI_Tracking_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> target_;
  bool useParabolicInflow_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    const Real one(1);
    return one;
  }

  Real target(const std::vector<Real> &x, const int dir) const {
    const Real one(1), two(2);
    Real val(0);
    if (dir == 0) {
     val = (useParabolicInflow_ ? (two + x[1])*(two - x[1]) : one);
    }
    return val;
  }

public:
  QoI_Tracking_NavierStokes(Teuchos::ParameterList &parlist,
                            const ROL::Ptr<FE<Real>> &feVel,
                            const ROL::Ptr<FE<Real>> &fePrs,
                            const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    useParabolicInflow_ = parlist.sublist("Problem").get("Use Parabolic Inflow", true);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    target_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i, j, k);
        }
        (*weight_)(i, j)  = weightFunc(pt);
        for (int k = 0; k < d; ++k) {
          (*target_)(i, j, k) = target(pt, k);
        }
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = feVel_->gradN()->dimension(0);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Compute tracking term and integrate
    ROL::Ptr<Intrepid::FieldContainer<Real>> U_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      // Evaluate u value on FE basis
      U_eval->initialize();
      feVel_->evaluateValue(U_eval, U[i]);
      // Compute distance to target
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*U_eval)(j, k) -= (*target_)(j, k, i);
          (*U_eval)(j, k) *= std::sqrt((*weight_)(j, k));
        }
      }
      // Compute L2 norm squared
      feVel_->computeIntegral(val,U_eval,U_eval,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Compute weighted u value and integrate
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+1, ROL::nullPtr);
    ROL::Ptr<Intrepid::FieldContainer<Real>> U_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      // Evaluate on FE basis
      U_eval->initialize();
      feVel_->evaluateValue(U_eval, U[i]);
      // Compute tracking term
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*U_eval)(j, k) -= (*target_)(j, k, i);
          (*U_eval)(j, k) *= (*weight_)(j, k);
        }
      }
      // Build local gradient of state tracking term
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *U_eval,
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute weighted v value and integrate
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+1, ROL::nullPtr);
    ROL::Ptr<Intrepid::FieldContainer<Real>> V_eval
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      // Evaluate v value on FE basis
      V_eval->initialize();
      feVel_->evaluateValue(V_eval, V[i]);
      // Weight v value
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*V_eval)(j, k) *= (*weight_)(j, k);
        }
      }
      // Integrate v value
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *V_eval,
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Tracking_NavierStokes

template <class Real>
class QoI_Dissipation_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  Real nu_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Dissipation_NavierStokes(Teuchos::ParameterList &parlist,
                               const ROL::Ptr<FE<Real>> &feVel,
                               const ROL::Ptr<FE<Real>> &fePrs,
                               const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), fieldHelper_(fieldHelper) {
    Real Re = parlist.sublist("Problem").get("Reynolds Number", 40.0);
    nu_ = static_cast<Real>(1.0)/Re;
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const Real half(0.5);
    // Get relevant dimensions
    const int c = feVel_->gradN()->dimension(0);
    const int p = feVel_->gradN()->dimension(2);
    const int d = feVel_->gradN()->dimension(3);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradU(d);
    for (int i = 0; i < d; ++i) {
      gradU[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradU[i], U[i]);
    }
    // Compute energy dissipation rate
    ROL::Ptr<Intrepid::FieldContainer<Real>> sigmaU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < c; ++k) {
          for (int l = 0; l < p; ++l) {
            (*sigmaU)(k,l) = (*weight_)(k,l)*((*gradU[i])(k,l,j)+(*gradU[j])(k,l,i));
          }
        }
        feVel_->computeIntegral(val,sigmaU,sigmaU,true);
      }
    }
    Intrepid::RealSpaceTools<Real>::scale(*val, half*nu_);
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+1);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradU(d);
    for (int i = 0; i < d; ++i) {
      gradU[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradU[i], U[i]);
    }
    // Compute energy dissipation gradient
    ROL::Ptr<Intrepid::FieldContainer<Real>> sigmaU
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < c; ++k) {
          for (int l = 0; l < p; ++l) {
            (*sigmaU)(k,l) = std::pow((*weight_)(k,l),2)*((*gradU[i])(k,l,j)+(*gradU[j])(k,l,i));
          }
        }
        Intrepid::FunctionSpaceTools::integrate<Real>(*G[j],
                                                      *sigmaU,
                                                      *(feVel_->DNDdetJ(i)),
                                                      Intrepid::COMP_CPP,
                                                      true);
        Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                      *sigmaU,
                                                      *(feVel_->DNDdetJ(j)),
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
    }
    fieldHelper_->combineFieldCoeff(grad, G);
    Intrepid::RealSpaceTools<Real>::scale(*grad, nu_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Dissipation_NavierStokes::gradient_2 is zero.");
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+1);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradV(d);
    for (int i = 0; i < d; ++i) {
      gradV[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradV[i], V[i]);
    }
    // Compute energy dissipation gradient
    ROL::Ptr<Intrepid::FieldContainer<Real>> sigmaV
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        for (int k = 0; k < c; ++k) {
          for (int l = 0; l < p; ++l) {
            (*sigmaV)(k,l) = std::pow((*weight_)(k,l),2)*((*gradV[i])(k,l,j)+(*gradV[j])(k,l,i));
          }
        }
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[j],
                                                      *sigmaV,
                                                      *(feVel_->DNDdetJ(i)),
                                                      Intrepid::COMP_CPP,
                                                      true);
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                      *sigmaV,
                                                      *(feVel_->DNDdetJ(j)),
                                                      Intrepid::COMP_CPP,
                                                      true);
      }
    }
    fieldHelper_->combineFieldCoeff(hess, H);
    Intrepid::RealSpaceTools<Real>::scale(*hess, nu_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Dissipation_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Dissipation_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Dissipation_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Dissipation_NavierStokes

template <class Real>
class QoI_State_NavierStokes : public QoI<Real> {
private:
  ROL::Ptr<QoI<Real>> qoi_;

  void initialize(const std::string &stateObj,
                  Teuchos::ParameterList &parlist,
                  const ROL::Ptr<FE<Real>> &feVel,
                  const ROL::Ptr<FE<Real>> &fePrs,
                  const ROL::Ptr<FieldHelper<Real>> &fieldHelper) {
    if ( stateObj == "Circulation" ) {
      qoi_ = ROL::makePtr<QoI_Circulation_NavierStokes<Real>>(feVel,fePrs,fieldHelper);
    }
    else if ( stateObj == "Vorticity" ) {
      qoi_ = ROL::makePtr<QoI_Vorticity_NavierStokes<Real>>(feVel,fePrs,fieldHelper);
    }
    else if ( stateObj == "Directional" ) {
      qoi_ = ROL::makePtr<QoI_Horizontal_NavierStokes<Real>>(feVel,fePrs,fieldHelper);
    }
    else if ( stateObj == "Tracking" ) {
      qoi_ = ROL::makePtr<QoI_Tracking_NavierStokes<Real>>(parlist,feVel,fePrs,fieldHelper);
    }
    else if ( stateObj == "Dissipation" ) {
      qoi_ = ROL::makePtr<QoI_Dissipation_NavierStokes<Real>>(parlist,feVel,fePrs,fieldHelper);
    }
    else {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Unknown objective type."); 
    }
  }

public:
  QoI_State_NavierStokes(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const ROL::Ptr<FieldHelper<Real>> &fieldHelper) {
    std::string stateObj = parlist.sublist("Problem").get("Objective type","Vorticity");
    initialize(stateObj,parlist,feVel,fePrs,fieldHelper);
  }

  QoI_State_NavierStokes(const std::string &stateObj,
                         Teuchos::ParameterList &parlist,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const ROL::Ptr<FieldHelper<Real>> &fieldHelper) {
    initialize(stateObj,parlist,feVel,fePrs,fieldHelper);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    return qoi_->value(val, u_coeff, z_coeff, z_param);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    qoi_->gradient_1(grad, u_coeff, z_coeff, z_param);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    qoi_->gradient_2(grad, u_coeff, z_coeff, z_param);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    qoi_->HessVec_11(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    qoi_->HessVec_12(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    qoi_->HessVec_21(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    qoi_->HessVec_22(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

};

template <class Real>
class QoI_DownStreamPower_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const std::vector<ROL::Ptr<FE<Real>>> feVelBdry_;
  const std::vector<std::vector<int>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  const std::vector<Real> target_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> &cell_coeff,
      const int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feVel_->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real >> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_DownStreamPower_NavierStokes(const ROL::Ptr<FE<Real>>              &feVel,
                                   const ROL::Ptr<FE<Real>>              &fePrs,
                                   const std::vector<ROL::Ptr<FE<Real>>> &feVelBdry,
                                   const std::vector<std::vector<int>>   &bdryCellLocIds,
                                   const ROL::Ptr<FieldHelper<Real>>     &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feVelBdry_(feVelBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper),
      target_({1, 0, 0}) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>>             &val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    const Real half(0.5);
    const int c = feVel_->gradN()->dimension(0);
    const int d = feVel_->gradN()->dimension(3);
    const int numLocSides = bdryCellLocIds_.size();
    // Initialize storage
    int numCellsSide(0), numCubPerSide(0), cidx(0);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, u(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> tmp, u_diff, intVal;
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Compute cost integral
    for (int l = 0; l < numLocSides; ++l) {
      numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        u_diff = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intVal = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
        // Evaluate objective function value
        for (int i = 0; i < d; ++i) {
          u[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          tmp  = getBoundaryCoeff(*U[i], l);
          feVelBdry_[l]->evaluateValue(u[i], tmp);
          for (int j = 0; j < numCellsSide; ++j) {
            for (int k = 0; k < numCubPerSide; ++k) {
              (*u_diff)(j,k) += std::pow((*u[i])(j,k)-target_[i],2);
            }
          }
        }
        // Integrate objective function
        feVelBdry_[l]->computeIntegral(intVal,u_diff,u[0],false);
        // Add to volume integral value
        for (int i = 0; i < numCellsSide; ++i) {
          cidx = bdryCellLocIds_[l][i];
          (*val)(cidx) += half*(*intVal)(i);
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>>             &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    const Real half(0.5);
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
    const int numLocSides = bdryCellLocIds_.size();
    // Initialize output grad
    int numCellsSide(0), numCubPerSide(0), cidx(0);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+1), u(d), U;
    ROL::Ptr<Intrepid::FieldContainer<Real>> tmp, du, intGrad;
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Compute cost integral
    for (int l = 0; l < numLocSides; ++l) {
      numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate objective function gradient
        du      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intGrad = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        for (int i = 0; i < d; ++i) {
          tmp  = getBoundaryCoeff(*U[i], l);
          u[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feVelBdry_[l]->evaluateValue(u[i], tmp);
        }
        for (int i = 0; i < d; ++i) {
          for (int j = 0; j < numCellsSide; ++j) {
            for (int k = 0; k < numCubPerSide; ++k) {
              (*du)(j,k) = ((*u[i])(j,k)-target_[i])*(*u[0])(j,k);
              if (i==0) {
                for (int m = 0; m < d; ++m) {
                  (*du)(j,k) += half*std::pow((*u[m])(j,k)-target_[m],2);
                }
              }
            }
          }
          // Integrate gradient
          Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                        *du,
                                                        *(feVelBdry_[l]->NdetJ()),
                                                        Intrepid::COMP_CPP,
                                                        false);
          for (int j = 0; j < numCellsSide; ++j) {
            cidx = bdryCellLocIds_[l][j];
            for (int k = 0; k < fv; ++k) {
              (*G[i])(cidx,k) += (*intGrad)(j,k);
            }
          }
        }
      }
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>>             &grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_DownStreamPower_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
    const int numLocSides = bdryCellLocIds_.size();
    // Initialize output hess
    int numCellsSide(0), numCubPerSide(0), cidx(0);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+1), U, V, u(d), v(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> tmp, du, intHess;
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    for (int l = 0; l < numLocSides; ++l) {
      numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate state and direction on FE basis
        du      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intHess = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        for (int i = 0; i < d; ++i) {
          tmp  = getBoundaryCoeff(*U[i], l);
          u[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feVelBdry_[l]->evaluateValue(u[i], tmp);
          tmp  = getBoundaryCoeff(*V[i], l);
          v[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          feVelBdry_[l]->evaluateValue(v[i], tmp);
        }
        // Compute hessian times a vector
        for (int i = 0; i < d; ++i) {
          for (int j = 0; j < numCellsSide; ++j) {
            for (int k = 0; k < numCubPerSide; ++k) {
              if (i==0) {
                (*du)(j,k) = (((*u[0])(j,k)-target_[i])+(*u[0])(j,k))*(*v[0])(j,k);
                for (int m = 0; m < d; ++m) {
                  (*du)(j,k) += ((*u[m])(j,k)-target_[m])*(*v[m])(j,k);
                }
              }
              else {
                (*du)(j,k) = (*u[i])(j,k)*(*v[0])(j,k) + (*u[0])(j,k)*(*v[i])(j,k);
              }
            }
          }
          // Integrate hessian
          Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                        *du,
                                                        *(feVelBdry_[l]->NdetJ()),
                                                        Intrepid::COMP_CPP,
                                                        false);
          for (int j = 0; j < numCellsSide; ++j) {
            cidx = bdryCellLocIds_[l][j];
            for (int k = 0; k < fv; ++k) {
              (*H[i])(cidx,k) += (*intHess)(j,k);
            }
          }
        }
      }
    }
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_DownStreamPower_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_DownStreamPower_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_DownStreamPower_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_DownStreamPower_NavierStokes

template <class Real>
class QoI_L2Penalty_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const std::vector<ROL::Ptr<FE<Real>>> feVelBdry_;
  const std::vector<std::vector<int>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feVel_->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real >> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_L2Penalty_NavierStokes(const ROL::Ptr<FE<Real>>              &feVel,
                             const ROL::Ptr<FE<Real>>              &fePrs,
                             const std::vector<ROL::Ptr<FE<Real>>> &feVelBdry,
                             const std::vector<std::vector<int>>   &bdryCellLocIds,
                             const ROL::Ptr<FieldHelper<Real>>     &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feVelBdry_(feVelBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = feVel_->gradN()->dimension(0);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate x-component of control on FE basis
        ROL::Ptr<Intrepid::FieldContainer<Real>> zx_coeff_bdry
          = getBoundaryCoeff(*Z[0], l);
        ROL::Ptr<Intrepid::FieldContainer<Real>> valZX_eval
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZX_eval, zx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        ROL::Ptr<Intrepid::FieldContainer<Real>> zy_coeff_bdry
          = getBoundaryCoeff(*Z[1], l);
        ROL::Ptr<Intrepid::FieldContainer<Real>> valZY_eval
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZY_eval, zy_coeff_bdry);
        // Integrate cell L2 norm squared
        ROL::Ptr<Intrepid::FieldContainer<Real>> intVal
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
        feVelBdry_[l]->computeIntegral(intVal,valZX_eval,valZX_eval,false);
        feVelBdry_[l]->computeIntegral(intVal,valZY_eval,valZY_eval,true);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          (*val)(cidx) += static_cast<Real>(0.5)*(*intVal)(i);
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velZX_grad(c, fv);
    Intrepid::FieldContainer<Real> velZY_grad(c, fv);
    Intrepid::FieldContainer<Real> presZ_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velZX_grad);
    G[1] = ROL::makePtrFromRef(velZY_grad);
    G[2] = ROL::makePtrFromRef(presZ_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate x-component of control on FE basis
        ROL::Ptr<Intrepid::FieldContainer<Real>> zx_coeff_bdry
          = getBoundaryCoeff(*Z[0], l);
        ROL::Ptr<Intrepid::FieldContainer<Real>> valZX_eval
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZX_eval, zx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        ROL::Ptr<Intrepid::FieldContainer<Real>> zy_coeff_bdry
          = getBoundaryCoeff(*Z[1], l);
        ROL::Ptr<Intrepid::FieldContainer<Real>> valZY_eval
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valZY_eval, zy_coeff_bdry);
        // Compute gradient of squared L2-norm of diff
        ROL::Ptr<Intrepid::FieldContainer<Real>> intGradX
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        Intrepid::FunctionSpaceTools::integrate<Real>(*intGradX,
                                                      *valZX_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        ROL::Ptr<Intrepid::FieldContainer<Real>> intGradY
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        Intrepid::FunctionSpaceTools::integrate<Real>(*intGradY,
                                                      *valZY_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fv; ++j) {
            (*G[0])(cidx,j) += (*intGradX)(i,j);
            (*G[1])(cidx,j) += (*intGradY)(i,j);
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldHelper_->numFields());
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        // Evaluate x-component of control on FE basis
        ROL::Ptr<Intrepid::FieldContainer<Real>> vx_coeff_bdry
          = getBoundaryCoeff(*V[0], l);
        ROL::Ptr<Intrepid::FieldContainer<Real>> valVX_eval
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valVX_eval, vx_coeff_bdry);
        // Evaluate y-component of control on FE basis
        ROL::Ptr<Intrepid::FieldContainer<Real>> vy_coeff_bdry
          = getBoundaryCoeff(*V[1], l);
        ROL::Ptr<Intrepid::FieldContainer<Real>> valVY_eval
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        feVelBdry_[l]->evaluateValue(valVY_eval, vy_coeff_bdry);
        // Compute gradient of squared L2-norm of diff
        ROL::Ptr<Intrepid::FieldContainer<Real>> intHessX
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHessX,
                                                      *valVX_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        ROL::Ptr<Intrepid::FieldContainer<Real>> intHessY
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        Intrepid::FunctionSpaceTools::integrate<Real>(*intHessY,
                                                      *valVY_eval,
                                                      *(feVelBdry_[l]->NdetJ()),
                                                      Intrepid::COMP_CPP, false);
        // Add to integral value
        for (int i = 0; i < numCellsSide; ++i) {
          int cidx = bdryCellLocIds_[l][i];
          for (int j = 0; j < fv; ++j) {
            (*G[0])(cidx,j) += (*intHessX)(i,j);
            (*G[1])(cidx,j) += (*intHessY)(i,j);
          }
        }
      }
    }

    fieldHelper_->combineFieldCoeff(hess, G);
  }

}; // QoI_L2Penalty_NavierStokes

template <class Real>
class QoI_RotationControl_NavierStokes : public QoI<Real> {
public:
  QoI_RotationControl_NavierStokes(void) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int size = z_param->size();
    Real half(0.5), sum(0);
    for (int i = 0; i < size; ++i) {
      sum += std::pow((*z_param)[i],2);
    }
    val = ROL::nullPtr;
    return half * sum;
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    std::vector<Real> g; g.assign(z_param->begin(),z_param->end());
    return g;
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_22 is zero.");
  }

  void HessVec_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_NavierStokes::HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    std::vector<Real> h; h.assign(v_param->begin(),v_param->end());
    return h;
  }

}; // QoI_RotationControl_NavierStoke

#endif
