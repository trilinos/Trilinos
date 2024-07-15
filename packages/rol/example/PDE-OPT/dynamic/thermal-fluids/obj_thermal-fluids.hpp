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

#ifndef DYNAMIC_PDEOPT_QOI_THERMALFLUIDS_HPP
#define DYNAMIC_PDEOPT_QOI_THERMALFLUIDS_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_thermal-fluids.hpp"

template <class Real>
class QoI_Vorticity_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Vorticity_ThermalFluids(const ROL::Ptr<FE<Real>> &feVel,
                              const ROL::Ptr<FE<Real>> &fePrs,
                              const ROL::Ptr<FE<Real>> &feThr,
                              const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
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
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradU_vec(d);
    for (int i = 0; i < d; ++i) {
      gradU_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradU_vec[i], U[i]);
    }
    // Compute weighted curl
    ROL::Ptr<Intrepid::FieldContainer<Real>> curlU_eval;
    if (d==2) {
      curlU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    }
    else if (d==3) {
      curlU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    }
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        if (d==2) {
          (*curlU_eval)(i,j) = (*weight_)(i,j)
                               * ((*gradU_vec[1])(i,j,0) - (*gradU_vec[0])(i,j,1));
        }
        else if (d==3) {
          for (int k = 0; k < d; ++k) {
            int i1 = (k+2)%d, i2 = (k+1)%d;
            (*curlU_eval)(i,j,k) = (*weight_)(i,j)
                                   * ((*gradU_vec[i1])(i,j,i2) - (*gradU_vec[i2])(i,j,i1));
          }
        }
      }
    }
    // Compute L2 norm squared
    feVel_->computeIntegral(val,curlU_eval,curlU_eval,false);
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradU_vec(d);
    for (int i = 0; i < d; ++i) {
      gradU_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradU_vec[i], U[i]);
    }
    // Compute weighted curl
    int size = (d==2) ? 1 : d;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> curlU_vec(size);
    if (d==2) {
      curlU_vec[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          (*curlU_vec[0])(i,j) = (*weight_)(i,j)
                               * ((*gradU_vec[1])(i,j,0) - (*gradU_vec[0])(i,j,1));
        }
      }
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        curlU_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
        for (int j = 0; j < c; ++j) {
          for (int k = 0; k < p; ++k) {
            int i1 = (i+2)%d, i2 = (i+1)%d;
            (*curlU_vec[i])(j,k) = (*weight_)(j,k)
                                   * ((*gradU_vec[i1])(j,k,i2) - (*gradU_vec[i2])(j,k,i1));
          }
        }
      }
    }
    // Build local gradient of state tracking term
    if (d==2) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                    *curlU_vec[0],
                                                    *(feVel_->DNDdetJ(1)),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*G[0],static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                    *curlU_vec[0],
                                                    *(feVel_->DNDdetJ(0)),
                                                    Intrepid::COMP_CPP, false);
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        int i1 = (i+2)%d, i2 = (i+1)%d;
        Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                      *curlU_vec[i1],
                                                      *(feVel_->DNDdetJ(i2)),
                                                      Intrepid::COMP_CPP, false);
        Intrepid::RealSpaceTools<Real>::scale(*G[i],static_cast<Real>(-1));
        Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                      *curlU_vec[i2],
                                                      *(feVel_->DNDdetJ(i1)),
                                                      Intrepid::COMP_CPP, true);
      }
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::gradient_2 is zero.");
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> gradV_vec(d);
    for (int i = 0; i < d; ++i) {
      gradV_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradV_vec[i], V[i]);
    }
    // Compute weighted curl
    int size = (d==2) ? 1 : d;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> curlV_vec(size);
    if (d==2) {
      curlV_vec[0] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          (*curlV_vec[0])(i,j) = (*weight_)(i,j)
                               * ((*gradV_vec[1])(i,j,0) - (*gradV_vec[0])(i,j,1));
        }
      }
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        curlV_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
        for (int j = 0; j < c; ++j) {
          for (int k = 0; k < p; ++k) {
            int i1 = (i+2)%d, i2 = (i+1)%d;
            (*curlV_vec[i])(j,k) = (*weight_)(j,k)
                                   * ((*gradV_vec[i1])(j,k,i2) - (*gradV_vec[i2])(j,k,i1));
          }
        }
      }
    }
    // Build local gradient of state tracking term
    if (d==2) {
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                    *curlV_vec[0],
                                                    *(feVel_->DNDdetJ(1)),
                                                    Intrepid::COMP_CPP, false);
      Intrepid::RealSpaceTools<Real>::scale(*H[0],static_cast<Real>(-1));
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[1],
                                                    *curlV_vec[0],
                                                    *(feVel_->DNDdetJ(0)),
                                                    Intrepid::COMP_CPP, false);
    }
    else if (d==3) {
      for (int i = 0; i < d; ++i) {
        int i1 = (i+2)%d, i2 = (i+1)%d;
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                      *curlV_vec[i1],
                                                      *(feVel_->DNDdetJ(i2)),
                                                      Intrepid::COMP_CPP, false);
        Intrepid::RealSpaceTools<Real>::scale(*H[i],static_cast<Real>(-1));
        Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                      *curlV_vec[i2],
                                                      *(feVel_->DNDdetJ(i1)),
                                                      Intrepid::COMP_CPP, true);
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Vorticity_ThermalFluids

template <class Real>
class QoI_Circulation_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Circulation_ThermalFluids(const ROL::Ptr<FE<Real>> &feVel,
                                const ROL::Ptr<FE<Real>> &fePrs,
                                const ROL::Ptr<FE<Real>> &feThr,
                                const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
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
    int c  = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int fh = feThr_->gradN()->dimension(1);
    int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Build local gradient of state tracking term
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *weight_,
                                                  *(feVel_->DNDdetJ(1)),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::RealSpaceTools<Real>::scale(*G[0],static_cast<Real>(-1));
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *weight_,
                                                  *(feVel_->DNDdetJ(0)),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Circulation_ThermalFluids

template <class Real>
class QoI_Horizontal_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Horizontal_ThermalFluids(const ROL::Ptr<FE<Real>> &feVel,
                               const ROL::Ptr<FE<Real>> &fePrs,
                               const ROL::Ptr<FE<Real>> &feThr,
                               const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
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
    int c  = feVel_->gradN()->dimension(0);
    int fv = feVel_->gradN()->dimension(1);
    int fp = fePrs_->gradN()->dimension(1);
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
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
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[0],
                                                  *weighted_minUX_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*G[1],
                                                  *weighted_valUY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::gradient_2 is zero.");
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
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
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[0],
                                                  *weighted_minVX_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
    Intrepid::FunctionSpaceTools::integrate<Real>(*H[1],
                                                  *weighted_valVY_eval,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Horizontal_ThermalFluids

template <class Real>
class QoI_Tracking_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> target_;

  Real weightFunc(const std::vector<Real> & x) const {
    const Real one(1);
    return one;
  }

  Real target(const std::vector<Real> &x, const int dir) const {
    const Real zero(0), one(1);
    return (dir == 0 ? one : zero);
  }

public:
  QoI_Tracking_ThermalFluids(const ROL::Ptr<FE<Real>> &feVel,
                             const ROL::Ptr<FE<Real>> &fePrs,
                             const ROL::Ptr<FE<Real>> &feThr,
                             const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Compute weighted u value and integrate
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2, ROL::nullPtr);
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
    G[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::gradient_2 is zero.");
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute weighted v value and integrate
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2, ROL::nullPtr);
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
    H[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Tracking_ThermalFluids

template <class Real>
class QoI_Dissipation_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;

  Real nu_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Dissipation_ThermalFluids(Teuchos::ParameterList &parlist,
                               const ROL::Ptr<FE<Real>> &feVel,
                               const ROL::Ptr<FE<Real>> &fePrs,
                               const ROL::Ptr<FE<Real>> &feThr,
                               const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
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
    const int fh = feThr_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
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
    throw Exception::Zero(">>> QoI_Dissipation_ThermalFluids::gradient_2 is zero.");
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
    const int fh = feThr_->gradN()->dimension(1);
    const int p  = feVel_->gradN()->dimension(2);
    const int d  = feVel_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
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
    throw Exception::Zero(">>> QoI_Dissipation_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Dissipation_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Dissipation_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Dissipation_ThermalFluids

template <class Real>
class QoI_Bouyancy_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real>> weight_;
  Real alpha_;

  Real weightFunc(const std::vector<Real> & x) const {
    const Real one(1);
    return one;
  }

public:
  QoI_Bouyancy_ThermalFluids(Teuchos::ParameterList &parlist,
                             const ROL::Ptr<FE<Real>> &feVel,
                             const ROL::Ptr<FE<Real>> &fePrs,
                             const ROL::Ptr<FE<Real>> &feThr,
                             const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), fieldHelper_(fieldHelper) {
    Real Re = parlist.sublist("Problem").get("Reynolds Number", 200.0);
    Real Gr = parlist.sublist("Problem").get("Grashof Number",  40000.0);
    alpha_ = static_cast<Real>(-1)*Gr/(Re*Re);
    int c = feVel_->cubPts()->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    weight_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    std::vector<Real> pt(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*feVel_->cubPts())(i, j, k);
        }
        (*weight_)(i, j)  = alpha_*weightFunc(pt);
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
    ROL::Ptr<Intrepid::FieldContainer<Real>> Uv, Ut, Uvt;
    Uv  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Ut  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Uvt = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(Uv, U[d-1]);
    feThr_->evaluateValue(Ut, U[d+1]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*Uvt, *Uv, *Ut);
    feVel_->computeIntegral(val,weight_,Uvt,false);
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Compute weighted u value and integrate
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2, ROL::nullPtr);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);

    ROL::Ptr<Intrepid::FieldContainer<Real>> Uv, Ut, wUv, wUt;
    Uv  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Ut  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    wUv = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    wUt = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(Uv, U[d-1]);
    feThr_->evaluateValue(Ut, U[d+1]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wUv, *weight_, *Uv);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wUt, *weight_, *Ut);

    Intrepid::FunctionSpaceTools::integrate<Real>(*G[d-1],
                                                  *wUt,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);

    Intrepid::FunctionSpaceTools::integrate<Real>(*G[d+1],
                                                  *wUv,
                                                  *(feThr_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Bouyancy_ThermalFluids::gradient_2 is zero.");
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
    int fh = feThr_->gradN()->dimension(1);
    int p  = feVel_->gradN()->dimension(2);
    int d  = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute weighted v value and integrate
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2, ROL::nullPtr);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);

    ROL::Ptr<Intrepid::FieldContainer<Real>> Vv, Vt, wVv, wVt;
    Vv  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Vt  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    wVv = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    wVt = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(Vv, V[d-1]);
    feThr_->evaluateValue(Vt, V[d+1]);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wVv, *weight_, *Vv);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wVt, *weight_, *Vt);

    Intrepid::FunctionSpaceTools::integrate<Real>(*H[d-1],
                                                  *wVt,
                                                  *(feVel_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);

    Intrepid::FunctionSpaceTools::integrate<Real>(*H[d+1],
                                                  *wVv,
                                                  *(feThr_->NdetJ()),
                                                  Intrepid::COMP_CPP,
                                                  false);

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Tracking_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Tracking_ThermalFluids

template <class Real>
class QoI_State_ThermalFluids : public QoI<Real> {
private:
  ROL::Ptr<QoI<Real>> qoi_;

  void initialize(const std::string &stateObj,
                  Teuchos::ParameterList &parlist,
                  const ROL::Ptr<FE<Real>> &feVel,
                  const ROL::Ptr<FE<Real>> &fePrs,
                  const ROL::Ptr<FE<Real>> &feThr,
                  const ROL::Ptr<FieldHelper<Real>> &fieldHelper) {
    if ( stateObj == "Circulation" ) {
      qoi_ = ROL::makePtr<QoI_Circulation_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper);
    }
    else if ( stateObj == "Vorticity" ) {
      qoi_ = ROL::makePtr<QoI_Vorticity_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper);
    }
    else if ( stateObj == "Directional" ) {
      qoi_ = ROL::makePtr<QoI_Horizontal_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper);
    }
    else if ( stateObj == "Tracking" ) {
      qoi_ = ROL::makePtr<QoI_Tracking_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper);
    }
    else if ( stateObj == "Dissipation" ) {
      qoi_ = ROL::makePtr<QoI_Dissipation_ThermalFluids<Real>>(parlist,feVel,fePrs,feThr,fieldHelper);
    }
    else if ( stateObj == "Bouyancy" ) {
      qoi_ = ROL::makePtr<QoI_Bouyancy_ThermalFluids<Real>>(parlist,feVel,fePrs,feThr,fieldHelper);
    }
    else {
      throw Exception::NotImplemented(">>> (QoI_State_ThermalFluids): Unknown objective type."); 
    }
  }

public:
  QoI_State_ThermalFluids(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const ROL::Ptr<FE<Real>> &feThr,
                         const ROL::Ptr<FieldHelper<Real>> &fieldHelper) {
    std::string stateObj = parlist.sublist("Problem").get("Objective type","Vorticity");
    initialize(stateObj,parlist,feVel,fePrs,feThr,fieldHelper);
  }

  QoI_State_ThermalFluids(const std::string &stateObj,
                         Teuchos::ParameterList &parlist,
                         const ROL::Ptr<FE<Real>> &feVel,
                         const ROL::Ptr<FE<Real>> &fePrs,
                         const ROL::Ptr<FE<Real>> &feThr,
                         const ROL::Ptr<FieldHelper<Real>> &fieldHelper) {
    initialize(stateObj,parlist,feVel,fePrs,feThr,fieldHelper);
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
class QoI_DownStreamPower_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
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
  QoI_DownStreamPower_ThermalFluids(const ROL::Ptr<FE<Real>>              &feVel,
                                    const ROL::Ptr<FE<Real>>              &fePrs,
                                    const ROL::Ptr<FE<Real>>              &feThr,
                                    const std::vector<ROL::Ptr<FE<Real>>> &feVelBdry,
                                    const std::vector<std::vector<int>>   &bdryCellLocIds,
                                    const ROL::Ptr<FieldHelper<Real>>     &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), feVelBdry_(feVelBdry),
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
    const int fh = feThr_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
    const int numLocSides = bdryCellLocIds_.size();
    // Initialize output grad
    int numCellsSide(0), numCubPerSide(0), cidx(0);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2), u(d), U;
    ROL::Ptr<Intrepid::FieldContainer<Real>> tmp, du, intGrad;
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
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
    throw Exception::Zero(">>> QoI_DownStreamPower_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int d  = feVel_->gradN()->dimension(3);
    const int numLocSides = bdryCellLocIds_.size();
    // Initialize output hess
    int numCellsSide(0), numCubPerSide(0), cidx(0);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2), U, V, u(d), v(d);
    ROL::Ptr<Intrepid::FieldContainer<Real>> tmp, du, intHess;
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d]   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
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
    throw Exception::Zero(">>> QoI_DownStreamPower_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_DownStreamPower_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>>             &hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> &z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>>              &z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_DownStreamPower_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_DownStreamPower_ThermalFluids

template <class Real>
class QoI_L2Penalty_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> feVel_;
  const ROL::Ptr<FE<Real>> fePrs_;
  const ROL::Ptr<FE<Real>> feThr_;
  const std::vector<std::vector<ROL::Ptr<FE<Real>>>> feThrBdry_;
  const std::vector<std::vector<std::vector<int>>> bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real>> fieldHelper_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feThr_->N()->dimension(1);
    
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
  QoI_L2Penalty_ThermalFluids(const ROL::Ptr<FE<Real>> &feVel,
                              const ROL::Ptr<FE<Real>> &fePrs,
                              const ROL::Ptr<FE<Real>> &feThr,
                              const std::vector<std::vector<ROL::Ptr<FE<Real>>>> &feThrBdry,
                              const std::vector<std::vector<std::vector<int>>> &bdryCellLocIds,
                              const ROL::Ptr<FieldHelper<Real>> &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), feThrBdry_(feThrBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = feVel_->gradN()->dimension(0);
    const int d = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==4) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts()->dimension(1);
              // Evaluate control on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real>> z_coeff_bdry
                = getBoundaryCoeff(*Z[d+1], i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
              // Integrate cell L2 norm squared
              ROL::Ptr<Intrepid::FieldContainer<Real>> intVal
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide);
              feThrBdry_[i][j]->computeIntegral(intVal,valZ_eval,valZ_eval,false);
              // Add to integral value
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                (*val)(cidx) += static_cast<Real>(0.5)*(*intVal)(k);
              }
            }
          }
        }
      }
    }
    return static_cast<Real>(0);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int d  = feThr_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==4) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts()->dimension(1);
              // Evaluate control on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real>> z_coeff_bdry
                = getBoundaryCoeff(*Z[d+1], i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
              // Compute gradient of squared L2-norm of diff
              ROL::Ptr<Intrepid::FieldContainer<Real>> intGrad
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fh);
              Intrepid::FunctionSpaceTools::integrate<Real>(*intGrad,
                                                            *valZ_eval,
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              // Add to integral value
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  (*G[d+1])(cidx,l) += (*intGrad)(k,l);
                }
              }
            }
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
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int d  = feThr_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if (i==4) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts()->dimension(1);
              // Evaluate control on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff_bdry
                = getBoundaryCoeff(*V[d+1], i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valV_eval, v_coeff_bdry);
              // Compute gradient of squared L2-norm of diff
              ROL::Ptr<Intrepid::FieldContainer<Real>> intHess
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, fh);
              Intrepid::FunctionSpaceTools::integrate<Real>(*intHess,
                                                            *valV_eval,
                                                            *(feThrBdry_[i][j]->NdetJ()),
                                                            Intrepid::COMP_CPP, false);
              // Add to integral value
              for (int k = 0; k < numCellsSide; ++k) {
                int cidx = bdryCellLocIds_[i][j][k];
                for (int l = 0; l < fh; ++l) {
                  (*H[d+1])(cidx,l) += (*intHess)(k,l);
                }
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
class QoI_RotationControl_ThermalFluids : public QoI<Real> {
public:
  QoI_RotationControl_ThermalFluids(void) {}

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
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::gradient_2 is zero.");
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
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_13(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_13 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_22 is zero.");
  }

  void HessVec_23(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const std::vector<Real> > & v_param,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_23 is zero.");
  }

  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_31 is zero.");
  }

  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_RotationControl_ThermalFluids::HessVec_32 is zero.");
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
