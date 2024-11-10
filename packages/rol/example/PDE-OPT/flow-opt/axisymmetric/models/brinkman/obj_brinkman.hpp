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

#ifndef OBJ_BRINKMAN_HPP
#define OBJ_BRINKMAN_HPP

#include "../../../../TOOLS/qoi.hpp"
#include "pde_brinkman.hpp"
#include "impermeability.hpp"

template <class Real>
class QoI_Vorticity_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Vorticity_NavierStokes(const ROL::Ptr<FE<Real>> &feVel,
                             const ROL::Ptr<FE<Real>> &fePrs,
                             const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : feVel_(feVel), fePrs_(fePrs), fieldInfo_(fieldInfo), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
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
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
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
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    G[0] = ROL::makePtrFromRef(velUX_grad);
    G[1] = ROL::makePtrFromRef(velUY_grad);
    G[2] = ROL::makePtrFromRef(presU_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
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

    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
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
    int c  = z_coeff->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
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

    FieldUtils::combineFieldCoeff<Real>(hess, G, fieldInfo_, fieldInfo_);
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
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Circulation_NavierStokes(const ROL::Ptr<FE<Real>> &feVel,
                               const ROL::Ptr<FE<Real>> &fePrs,
                               const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : feVel_(feVel), fePrs_(fePrs), fieldInfo_(fieldInfo), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
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
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
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
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velUX_grad(c, fv);
    Intrepid::FieldContainer<Real> velUY_grad(c, fv);
    Intrepid::FieldContainer<Real> presU_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
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

    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
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
class QoI_Directional_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real eps_;
  std::vector<Real> dir_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Directional_NavierStokes(Teuchos::ParameterList &list,
                               const ROL::Ptr<FE<Real>> &feVel,
                               const ROL::Ptr<FE<Real>> &fePrs,
                               const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : feVel_(feVel), fePrs_(fePrs), fieldInfo_(fieldInfo), eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    const int c = feVel_->gradN()->dimension(0);
    const int p = feVel_->gradN()->dimension(2);
    const int d = feVel_->gradN()->dimension(3);
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
    Teuchos::Array<Real> dir = Teuchos::getArrayFromStringParameter<Real>(list, "Direction");
    dir_ = dir.toVector();
    Real norm(0);
    for (int i = 0; i < d; ++i) {
      norm += std::pow(dir_[i],static_cast<Real>(2));
    }
    norm = std::sqrt(norm);
    for (int i = 0; i < d; ++i) {
      dir_[i] /= norm;
    }
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> Uval, dot, Upar, Uort;
    Uval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    dot   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Upar  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    Uort  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      feVel_->evaluateValue(Uval, U[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*dot)(j,k)   += dir_[i] * (*Uval)(j,k);
          (*Upar)(j,k,i) = std::sqrt((*weight_)(j,k)) * dir_[i];
          (*Uort)(j,k,i) = std::sqrt((*weight_)(j,k)) * (*Uval)(j,k);
        }
      }
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*Upar)(j,k,i) *= (*dot)(j,k);
          (*Uort)(j,k,i) -= (*Upar)(j,k,i);
          (*Upar)(j,k,i)  = std::min(static_cast<Real>(0),(*Upar)(j,k,i));
        }
      }
    }
    // Compute L2 norm squared
    feVel_->computeIntegral(val,Upar,Upar,false);
    feVel_->computeIntegral(val,Uort,Uort,true);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
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
    // Get components of the state
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Upar(d), Uort(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> Uval, dot;
    Uval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    dot  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < d; ++i) {
      Uval->initialize();
      Upar[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Uort[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      feVel_->evaluateValue(Uval, U[i]);
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*dot)(j,k)    += dir_[i] * (*Uval)(j,k);
          (*Upar[i])(j,k) = std::sqrt((*weight_)(j,k)) * dir_[i];
          (*Uort[i])(j,k) = std::sqrt((*weight_)(j,k)) * (*Uval)(j,k);
        }
      }
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < c; ++j) {
        for (int k = 0; k < p; ++k) {
          (*Upar[i])(j,k) *= (*dot)(j,k);
          (*Uort[i])(j,k) -= (*Upar[i])(j,k);
          (*Upar[i])(j,k)  = std::min(static_cast<Real>(0),(*Upar[i])(j,k));
        }
      }
    }
    // Build local gradient of state tracking term
    for (int i = 0; i < d; ++i) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *Upar[i],
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *Uort[i],
                                                    *(feVel_->NdetJ()),
                                                    Intrepid::COMP_CPP, true);
    }
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Directional_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = z_coeff->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    // Initialize output grad
    Intrepid::FieldContainer<Real> velVX_grad(c, fv);
    Intrepid::FieldContainer<Real> velVY_grad(c, fv);
    Intrepid::FieldContainer<Real> presV_grad(c, fp);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    G[0] = ROL::makePtrFromRef(velVX_grad);
    G[1] = ROL::makePtrFromRef(velVY_grad);
    G[2] = ROL::makePtrFromRef(presV_grad);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U, V;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> UX, UY, VX, VY, WVX, WVY;
    UX  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    UY  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    VX  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    VY  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WVX = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WVY = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(UX, U[0]);
    feVel_->evaluateValue(UY, U[1]);
    feVel_->evaluateValue(VX, V[0]);
    feVel_->evaluateValue(VY, V[1]);
    // Compute negative part of x-velocity
    const Real zero(0), one(1);
    Real scale(0);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        scale = ((*UY)(i,j) >= zero) ? one : zero;
        (*WVX)(i,j) =         (*weight_)(i,j) * (*VX)(i,j);
        (*WVY)(i,j) = scale * (*weight_)(i,j) * (*VY)(i,j);
      }
    }
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(velVX_grad,
                                                  *WVX,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(velVY_grad,
                                                  *WVY,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    FieldUtils::combineFieldCoeff<Real>(hess, G, fieldInfo_, fieldInfo_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Directional_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Directional_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Directional_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Directional_NavierStokes

template <class Real>
class QoI_Pressure_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  const Real width_;
  const Real outletHeight_;
  const Real eps_;

  Real weightFunc(const std::vector<Real> & x) const {
    //return (x[1]<=outletHeight_+eps_ && std::abs(x[0]-width_)<=width_-0.6) ? static_cast<Real>(-1) : static_cast<Real>(0);
    return static_cast<Real>(1);
  }

public:
  QoI_Pressure_NavierStokes(const ROL::Ptr<FE<Real>> &feVel,
                            const ROL::Ptr<FE<Real>> &fePrs,
                            const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : feVel_(feVel), fePrs_(fePrs), fieldInfo_(fieldInfo),
      width_(12.7), outletHeight_(0.1*1.6936),
      eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) {
    const int c = feVel_->cubPts()->dimension(0);
    const int p = feVel_->cubPts()->dimension(1);
    const int d = feVel_->cubPts()->dimension(2);
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
    const int c = fePrs_->gradN()->dimension(0);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> P, WP;
    P  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WP = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fePrs_->evaluateValue(P, U[d]);
    // Scale the pressure
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WP,*weight_,*P);
    feVel_->computeIntegral(val,P,WP,false);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> P, WP;
    P  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WP = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fePrs_->evaluateValue(P, U[d]);
    // Scale pressure
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WP,*weight_,*P);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[d],
                                                  *WP,
                                                  *(fePrs_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Pressure_NavierStokes::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    FieldUtils::splitFieldCoeff(V, v_coeff, fieldInfo_);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> P, WP;
    P  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    WP = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fePrs_->evaluateValue(P, V[d]);
    // Scale pressure
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*WP,*weight_,*P);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[d],
                                                  *WP,
                                                  *(fePrs_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_, fieldInfo_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Pressure_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Pressure_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Pressure_NavierStokes::HessVec_22 is zero.");
  }

}; // QoI_Pressure_NavierStokes

template <class Real>
class QoI_Velocity_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>>                    feVel_;
  const ROL::Ptr<FE<Real>>                    fePrs_;
  const std::vector<ROL::Ptr<FE<Real>>>       feVelBdry_;
  const std::vector<std::vector<int>>         bdryCellLocIds_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  std::vector<Real>                           target_;
  bool                                        onlyAxial_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int locSideId) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[locSideId];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feVel_->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real>> bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_Velocity_NavierStokes(Teuchos::ParameterList   &list,
                            const ROL::Ptr<FE<Real>> &feVel,
                            const ROL::Ptr<FE<Real>> &fePrs,
                            const std::vector<ROL::Ptr<FE<Real>>> &feVelBdry,
                            const std::vector<std::vector<int>>   &bdryCellLocIds,
                            const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : feVel_(feVel), fePrs_(fePrs),
      feVelBdry_(feVelBdry), bdryCellLocIds_(bdryCellLocIds),
      fieldInfo_(fieldInfo) {
    target_.clear(); target_.resize(2);
    target_[0] = list.sublist("Problem").get("Target Radial Velocity",0.0);
    target_[1] = list.sublist("Problem").get("Target Axial Velocity",-15.0);
    onlyAxial_ = list.sublist("Problem").get("Only Use Axial Velocity",false);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = feVel_->gradN()->dimension(0);
    const int d = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, valU_eval, intVal;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intVal = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
        for (int k = 0; k < d; ++k) {
          if ((k==0 && !onlyAxial_) || k==1) {
            u_coeff_bdry = getBoundaryCoeff(*U[k], l);
            valU_eval->initialize();
            feVelBdry_[l]->evaluateValue(valU_eval, u_coeff_bdry);
            for (int i = 0; i < numCellsSide; ++i) {
              for (int j = 0; j < numCubPerSide; ++j) {
                (*valU_eval)(i,j) -= target[k];
              }
            }
            feVelBdry_[l]->computeIntegral(intVal,valU_eval,valU_eval,true);
          }
        }
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
    // Get relevant dimensions
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i)
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the state
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Compute cost integral
    std::vector<Real> target = (z_param == ROL::nullPtr) ? target_ : *z_param;
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, valU_eval, intGrad;
        valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intGrad = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        for (int k = 0; k < d; ++k) {
          if ((k==0 && !onlyAxial_) || k==1) {
            u_coeff_bdry = getBoundaryCoeff(*U[k], l);
            valU_eval->initialize();
            feVelBdry_[l]->evaluateValue(valU_eval, u_coeff_bdry);
            for (int i = 0; i < numCellsSide; ++i) {
              for (int j = 0; j < numCubPerSide; ++j) {
                (*valU_eval)(i,j) -= target[k];
              }
            }
            // Compute gradient of squared L2-norm of diff
            Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                          *valU_eval,
                                                          *(feVelBdry_[l]->NdetJ()),
                                                          Intrepid::COMP_CPP, false);
            // Add to integral value
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              for (int j = 0; j < fv; ++j) {
                (*G[k])(cidx,j) += (*intGrad)(i,j);
              }
            }
          }
        }
      }
    }
    FieldUtils::combineFieldCoeff(grad, G, fieldInfo_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Velocity_NavierStokes::gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = feVel_->gradN()->dimension(0);
      const int d = feVel_->gradN()->dimension(3);
      // Get components of the control
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
      FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int k = 0; k < d; ++k)
        grad[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> intVal(d);
          ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff_bdry, valU_eval, weight;
          valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              u_coeff_bdry = getBoundaryCoeff(*U[k], l);
              valU_eval->initialize();
              feVelBdry_[l]->evaluateValue(valU_eval, u_coeff_bdry);
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j) = static_cast<Real>(1);
                  (*valU_eval)(i,j) *= static_cast<Real>(-1);
                  (*valU_eval)(i,j) += (*z_param)[k];
                }
              }
              feVelBdry_[l]->computeIntegral(intVal[k],weight,valU_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (*grad[k])(cidx) += (*intVal[k])(i);
            }
          }
        }
      }
      return g_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_NavierStokes::gradient_3 is zero.");
    }
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
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i)
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the direction
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    // Compute cost integral
    const int numLocSides = bdryCellLocIds_.size();
    for (int l = 0; l < numLocSides; ++l) {
      const int numCellsSide  = bdryCellLocIds_[l].size();
      if ( numCellsSide ) {
        const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
        ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff_bdry, valV_eval, intHess;
        valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
        intHess = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
        for (int k = 0; k < d; ++k) {
          if ((k==0 && !onlyAxial_) || k==1) {
            v_coeff_bdry = getBoundaryCoeff(*V[k], l);
            valV_eval->initialize();
            feVelBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
            // Compute hessian of squared L2-norm of diff
            Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                          *valV_eval,
                                                          *(feVelBdry_[l]->NdetJ()),
                                                          Intrepid::COMP_CPP, false);
            // Add to integral value
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              for (int j = 0; j < fv; ++j) {
                (*H[k])(cidx,j) += (*intHess)(i,j);
              }
            }
          }
        }
      }
    }
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Velocity_NavierStokes::HessVec_12 is zero.");
  }

  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c  = feVel_->gradN()->dimension(0);
      const int fv = feVel_->gradN()->dimension(1);
      const int fp = fePrs_->gradN()->dimension(1);
      const int d  = feVel_->gradN()->dimension(3);
      // Initialize output grad
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
      H.resize(fieldInfo_->numFields);
      for (int i = 0; i < d; ++i)
        H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
      H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
          ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, intHess;
          valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          intHess = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fv);
          for (int k = 0; k < d; ++k) {
            if ((k==0 && !onlyAxial_) || k==1) {
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*valU_eval)(i,j) = -(*v_param)[k];
                }
              }
              // Compute gradient of squared L2-norm of diff
              Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                            *valU_eval,
                                                            *(feVelBdry_[l]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              // Add to integral value
              for (int i = 0; i < numCellsSide; ++i) {
                int cidx = bdryCellLocIds_[l][i];
                for (int j = 0; j < fv; ++j) {
                  (*H[k])(cidx,j) += (*intHess)(i,j);
                }
              }
            }
          }
        }
      }
      FieldUtils::combineFieldCoeff(hess, H, fieldInfo_);
    }
    else {
      throw Exception::NotImplemented(">>> HessVec_13 not implemented.");
    }
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Velocity_NavierStokes::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Velocity_NavierStokes::HessVec_22 is zero.");
  }

  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Velocity_NavierStokes::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = feVel_->gradN()->dimension(0);
      const int d = feVel_->gradN()->dimension(3);
      // Get components of the state
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
      FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> intVal(d);
          ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff_bdry, valV_eval, weight;
          valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              v_coeff_bdry = getBoundaryCoeff(*V[k], l);
              valV_eval->initialize();
              feVelBdry_[l]->evaluateValue(valV_eval, v_coeff_bdry);
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j) = static_cast<Real>(1);
                  (*valV_eval)(i,j) *= static_cast<Real>(-1);
                }
              }
              feVelBdry_[l]->computeIntegral(intVal[k],weight,valV_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (*hess[k])(cidx) += (*intVal[k])(i);
            }
          }
        }
      }
      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_NavierStokes::HessVec_31 is zero.");
    }
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = feVel_->gradN()->dimension(0);
      const int d = feVel_->gradN()->dimension(3);
      // Get components of the control
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
      FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      const int numLocSides = bdryCellLocIds_.size();
      for (int l = 0; l < numLocSides; ++l) {
        const int numCellsSide  = bdryCellLocIds_[l].size();
        if ( numCellsSide ) {
          const int numCubPerSide = feVelBdry_[l]->cubPts()->dimension(1);
          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> intVal(d);
          ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, weight;
          valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          weight    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
          for (int k = 0; k < d; ++k) {
            intVal[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
            if ((k==0 && !onlyAxial_) || k==1) {
              valU_eval->initialize();
              for (int i = 0; i < numCellsSide; ++i) {
                for (int j = 0; j < numCubPerSide; ++j) {
                  (*weight)(i,j) = static_cast<Real>(1);
                  (*valU_eval)(i,j) = (*v_param)[k];
                }
              }
              feVelBdry_[l]->computeIntegral(intVal[k],weight,valU_eval,true);
            }
          }
          // Add to integral value
          for (int k = 0; k < d; ++k) {
            for (int i = 0; i < numCellsSide; ++i) {
              int cidx = bdryCellLocIds_[l][i];
              (*hess[k])(cidx) += (*intVal[k])(i);
            }
          }
        }
      }
      return h_param;
    }
    else {
      throw Exception::Zero(">>> QoI_Velocity_NavierStokes::HessVec_33 is zero.");
    }
  }

}; // QoI_Velocity_NavierStokes


template <class Real>
class QoI_VelocityTracking_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>>                    feVel_, fePrs_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  ROL::Ptr<Intrepid::FieldContainer<Real>>    target_, weight_;
  Real                                        rad_, yvel_, frac_, twpow_;
  bool                                        onlyAxial_;

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
    if (Y > zero) {
      val = at*std::pow(Y-yTOP,p);
    } else {
      val = bt*std::pow(Y-yBOT,p);
    }
    //std::cout << Y << "  " << val << std::endl;
    return val;
  }

public:
  QoI_VelocityTracking_NavierStokes(Teuchos::ParameterList                      &list,
                                    const ROL::Ptr<FE<Real>>                    &feVel,
                                    const ROL::Ptr<FE<Real>>                    &fePrs,
                                    const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo)
    : feVel_(feVel), fePrs_(fePrs), fieldInfo_(fieldInfo) {
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
    const int c = fePrs_->gradN()->dimension(0);
    const int p = fePrs_->gradN()->dimension(2);
    target_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,2);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,2);
    std::vector<Real> x(2);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        x[0] = (*feVel_->cubPts())(i,j,0);
        x[1] = (*feVel_->cubPts())(i,j,1);
        (*target_)(i,j,0) = xTarget(x);
        (*target_)(i,j,1) = yTarget(x);
        if (useNorm && yWeight(x)) {
          xWScal = static_cast<Real>(1)
                  /(std::pow((*target_)(i,j,0),2) + std::pow((*target_)(i,j,1),2));
          yWScal = xWScal;
        }
        (*weight_)(i,j,0) = x[0] * xWScal * xWeight(x);
        (*weight_)(i,j,1) = x[0] * yWScal * yWeight(x);
      }
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fePrs_->gradN()->dimension(0);
    const int p = fePrs_->gradN()->dimension(2);
    const int d = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the velocity
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> vel, wvel, velx, vely;
    vel  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    wvel = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    velx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    vely = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    feVel_->evaluateValue(velx, U[0]);
    feVel_->evaluateValue(vely, U[1]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*vel)(i,j,0)  = (*velx)(i,j)-yvel*(*target_)(i,j,0);
        (*vel)(i,j,1)  = (*vely)(i,j)-yvel*(*target_)(i,j,1);
        (*wvel)(i,j,0) = (*weight_)(i,j,0)*(*vel)(i,j,0);
        (*wvel)(i,j,1) = (*weight_)(i,j,1)*(*vel)(i,j,1);
      }
    }

    feVel_->computeIntegral(val,vel,wvel);
    // Scale by one half
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G;
    G.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i)
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the velocity
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> velx, vely, wvelx, wvely;
    velx  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    vely  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    wvelx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    wvely = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    feVel_->evaluateValue(velx, U[0]);
    feVel_->evaluateValue(vely, U[1]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*wvelx)(i,j) = (*weight_)(i,j,0)*((*velx)(i,j)-yvel*(*target_)(i,j,0));
        (*wvely)(i,j) = (*weight_)(i,j,1)*((*vely)(i,j)-yvel*(*target_)(i,j,1));
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *wvelx,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *wvely,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    FieldUtils::combineFieldCoeff<Real>(grad, G, fieldInfo_);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> gradient_2 is zero.");
  }

  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      Real yvel = (z_param == ROL::nullPtr) ? yvel_ : (*z_param)[0];
      // Initialize output val
      std::vector<Real> g_param(d,static_cast<Real>(0));
      grad.clear(); grad.resize(d);
      for (int i = 0; i < d; ++i)
        grad[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Get components of the velocity
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
      FieldUtils::splitFieldCoeff<Real>(U, u_coeff, fieldInfo_);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> vel, wvel, velx, vely;
      vel  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      velx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      vely = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      feVel_->evaluateValue(velx, U[0]);
      feVel_->evaluateValue(vely, U[1]);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          (*vel)(i,j,0)  = -(*weight_)(i,j,0)*((*velx)(i,j)-yvel*(*target_)(i,j,0));
          (*vel)(i,j,1)  = -(*weight_)(i,j,1)*((*vely)(i,j)-yvel*(*target_)(i,j,1));
        }
      }

      feVel_->computeIntegral(grad[0],vel,target_);
      
      return g_param;
    }
    else {
      throw Exception::Zero(">>> gradient_3 is zero.");
    }
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c  = fePrs_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int p  = fePrs_->gradN()->dimension(2);
    const int d  = fePrs_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
    H.resize(fieldInfo_->numFields);
    for (int i = 0; i < d; ++i)
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    // Get components of the velocity
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
    // Compute cost integral
    ROL::Ptr<Intrepid::FieldContainer<Real>> velx, vely, wvelx, wvely;
    velx  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    vely  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    wvelx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    wvely = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
    feVel_->evaluateValue(velx, V[0]);
    feVel_->evaluateValue(vely, V[1]);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*wvelx)(i,j)  = (*weight_)(i,j,0)*(*velx)(i,j);
        (*wvely)(i,j)  = (*weight_)(i,j,1)*(*vely)(i,j);
      }
    }

    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *wvelx,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[1],
                                                  *wvely,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_12 is zero.");
  }

  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c  = fePrs_->gradN()->dimension(0);
      const int fv = feVel_->gradN()->dimension(1);
      const int fp = fePrs_->gradN()->dimension(1);
      const int p  = fePrs_->gradN()->dimension(2);
      const int d  = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H;
      H.resize(fieldInfo_->numFields);
      for (int i = 0; i < d; ++i)
        H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
      H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> wvelx, wvely;
      wvelx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      wvely = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);

      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          (*wvelx)(i,j) = -(*weight_)(i,j,0)*(*target_)(i,j,0)*(*v_param)[0];
          (*wvely)(i,j) = -(*weight_)(i,j,1)*(*target_)(i,j,1)*(*v_param)[0];
        }
      }

      Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                    *wvelx,
                                                    *feVel_->NdetJ(),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[1],
                                                    *wvely,
                                                    *feVel_->NdetJ(),
                                                    Intrepid::COMP_CPP, false);
      FieldUtils::combineFieldCoeff<Real>(hess, H, fieldInfo_);
    }
    else {
      throw Exception::Zero(">>> HessVec_13 is zero.");
    }
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_22 is zero.");
  }

  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const std::vector<Real>> & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Get components of the velocity
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
      FieldUtils::splitFieldCoeff<Real>(V, v_coeff, fieldInfo_);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> velx, vely, wvel;
      velx = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      vely = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p);
      wvel = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      feVel_->evaluateValue(velx, V[0]);
      feVel_->evaluateValue(vely, V[1]);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          (*wvel)(i,j,0) = -(*weight_)(i,j,0)*(*velx)(i,j);
          (*wvel)(i,j,1) = -(*weight_)(i,j,1)*(*vely)(i,j);
        }
      }

      feVel_->computeIntegral(hess[0],wvel,target_);

      return h_param;
    }
    else {
      throw Exception::Zero(">>> HessVec_31 is zero.");
    }
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> HessVec_32 is zero.");
  }

  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
                               const ROL::Ptr<const std::vector<Real>> & v_param,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                               const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                               const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    if (z_param != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fePrs_->gradN()->dimension(0);
      const int p = fePrs_->gradN()->dimension(2);
      const int d = fePrs_->gradN()->dimension(3);
      // Initialize output val
      std::vector<Real> h_param(d,static_cast<Real>(0));
      hess.clear(); hess.resize(d);
      for (int k = 0; k < d; ++k)
        hess[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
      // Compute cost integral
      ROL::Ptr<Intrepid::FieldContainer<Real>> wtarget;
      wtarget = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          for (int k = 0; k < d; ++k) {
            (*wtarget)(i,j,k) = (*weight_)(i,j,k)*(*target_)(i,j,k);
          }
        }
      }
      feVel_->computeIntegral(hess[0],wtarget,target_);
      Intrepid::RealSpaceTools<Real>::scale(*hess[0],(*v_param)[0]);
      
      return h_param;
    }
    else {
      throw Exception::Zero(">>> HessVec_33 is zero.");
    }
  }

}; // QoI_Velocity_Darcy2


template <class Real>
class QoI_Power_NavierStokes : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>>                    feVel_;
  const ROL::Ptr<FE<Real>>                    fePrs_;
  const ROL::Ptr<FE<Real>>                    feCtrl_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfo_;
  const ROL::Ptr<const FieldUtils::FieldInfo> fieldInfoCtrl_;

  ROL::Ptr<Impermeability<Real>> imp_;
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
    imp_ = ROL::makePtr<Impermeability<Real>>(list);
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
class QoI_State_NavierStokes : public QoI<Real> {
private:
  ROL::Ptr<QoI<Real>> qoi_;

public:
  QoI_State_NavierStokes(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo) {
    std::string stateObj = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj != "Vorticity"
      && stateObj != "Circulation"
      && stateObj != "Directional"
      && stateObj != "Pressure"
      && stateObj != "Velocity"
      && stateObj != "VelocityTracking"
      && stateObj != "Power" ) {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Unknown objective type."); 
    }
    if ( stateObj == "Vorticity" ) {
      qoi_ = ROL::makePtr<QoI_Vorticity_NavierStokes<Real>>(feVel,fePrs,fieldInfo);
    }
    else if ( stateObj == "Directional" ) {
      qoi_ = ROL::makePtr<QoI_Directional_NavierStokes<Real>>(parlist,feVel,fePrs,fieldInfo);
    }
    else if ( stateObj == "Pressure" ) {
      qoi_ = ROL::makePtr<QoI_Pressure_NavierStokes<Real>>(feVel,fePrs,fieldInfo);
    }
    else if ( stateObj == "Velocity" ) {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Incorrect Constructor.");
    }
    else if ( stateObj == "VelocityTracking" ) {
      qoi_ = ROL::makePtr<QoI_VelocityTracking_NavierStokes<Real>>(parlist,feVel,fePrs,fieldInfo);
    }
    else if ( stateObj == "Power" ) {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Incorrect Constructor.");
    }
    else {
      qoi_ = ROL::makePtr<QoI_Circulation_NavierStokes<Real>>(feVel,fePrs,fieldInfo);
    }
  }
  QoI_State_NavierStokes(Teuchos::ParameterList   &list,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const std::vector<ROL::Ptr<FE<Real>>> &feVelBdry,
                            const std::vector<std::vector<int>>   &bdryCellLocIds,
                            const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo) {
    std::string stateObj = list.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj != "Velocity" ) {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Incorrect Constructor.");
    }
    qoi_ = ROL::makePtr<QoI_Velocity_NavierStokes<Real>>(list,feVel,fePrs,feVelBdry,bdryCellLocIds,fieldInfo);
  }
  QoI_State_NavierStokes(Teuchos::ParameterList &list,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const ROL::Ptr<FE<Real>> &feCtrl,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfo,
                         const ROL::Ptr<const FieldUtils::FieldInfo> &fieldInfoCtrl) {
    std::string stateObj = list.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj != "Power" ) {
      throw Exception::NotImplemented(">>> (QoI_State_NavierStokes): Incorrect Constructor.");
    }
    qoi_ = ROL::makePtr<QoI_Power_NavierStokes<Real>>(feVel,fePrs,feCtrl,fieldInfo,fieldInfoCtrl);
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
    const Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real zero(0), one(1);
    volFraction_   = list.sublist("Problem").get("Volume Fraction",0.5);
    Real outHeight = list.sublist("Problem").get("Outlet Height",1.5*1.6936);
    const int c = fe_->cubPts()->dimension(0);
    const int p = fe_->cubPts()->dimension(1);
    const int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = (pt[1] <= outHeight+tol ? one : zero);
      }
    }
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

template <class Real>
class StdObjective_NavierStokes : public ROL::StdObjective<Real> {
private:
  std::string stateObj_;

public:
  StdObjective_NavierStokes(Teuchos::ParameterList &parlist) {
    stateObj_ = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj_ != "Vorticity"
      && stateObj_ != "Circulation"
      && stateObj_ != "Directional"
      && stateObj_ != "Pressure"
      && stateObj_ != "Velocity"
      && stateObj_ != "VelocityTracking"
      && stateObj_ != "Power" ) {
      throw Exception::NotImplemented(">>> (StdObjective_NavierStokes): Unknown objective type."); 
    }
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    Real val(0);
    if ( stateObj_ == "Vorticity"
      || stateObj_ == "Directional"
      || stateObj_ == "Pressure"
      || stateObj_ == "Velocity"
      || stateObj_ == "VelocityTracking"
      || stateObj_ == "Power" ) {
      val = x[0];
    }
    else {
      val = static_cast<Real>(0.5)*x[0]*x[0];
    }
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    const Real one(1);
    if ( stateObj_ == "Vorticity"
      || stateObj_ == "Directional"
      || stateObj_ == "Pressure"
      || stateObj_ == "Velocity"
      || stateObj_ == "VelocityTracking"
      || stateObj_ == "Power" ) {
      g[0] = one;
    }
    else {
      g[0] = x[0];
    }
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    if ( stateObj_ == "Vorticity"
      || stateObj_ == "Directional"
      || stateObj_ == "Pressure"
      || stateObj_ == "Velocity"
      || stateObj_ == "VelocityTracking"
      || stateObj_ == "Power" ) {
      hv[0] = zero;
    }
    else {
      hv[0] = v[0];
    }
  }

}; // OBJ_SCALAR

#endif
