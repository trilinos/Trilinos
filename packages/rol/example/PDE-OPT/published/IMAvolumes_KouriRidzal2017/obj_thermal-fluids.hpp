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

#ifndef PDEOPT_QOI_THERMALFLUIDS_HPP
#define PDEOPT_QOI_THERMALFLUIDS_HPP

#include "../../TOOLS/qoi.hpp"
#include "pde_thermal-fluids.hpp"

template <class Real>
class QoI_Vorticity_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > feVel_;
  const ROL::Ptr<FE<Real> > fePrs_;
  const ROL::Ptr<FE<Real> > feThr_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Vorticity_ThermalFluids(const ROL::Ptr<FE<Real> > &feVel,
                              const ROL::Ptr<FE<Real> > &fePrs,
                              const ROL::Ptr<FE<Real> > &feThr,
                              const ROL::Ptr<FieldHelper<Real> > &fieldHelper,
                              Teuchos::ParameterList &parlist)
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

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_vec(d);
    for (int i = 0; i < d; ++i) {
      gradU_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradU_vec[i], U[i]);
    }
    // Compute weighted curl
    ROL::Ptr<Intrepid::FieldContainer<Real> > curlU_eval;
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

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int fh = feThr_->N()->dimension(1);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradU_vec(d);
    for (int i = 0; i < d; ++i) {
      gradU_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradU_vec[i], U[i]);
    }
    // Compute weighted curl
    int size = (d==2) ? 1 : d;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > curlU_vec(size);
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
                                                      Intrepid::COMP_CPP, false);
      }
    }
    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = z_coeff->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int fh = feThr_->N()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > gradV_vec(d);
    for (int i = 0; i < d; ++i) {
      gradV_vec[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      feVel_->evaluateGradient(gradV_vec[i], V[i]);
    }
    // Compute weighted curl
    int size = (d==2) ? 1 : d;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > curlV_vec(size);
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
                                                      Intrepid::COMP_CPP, false);
      }
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Vorticity_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Vorticity_ThermalFluids

template <class Real>
class QoI_Circulation_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > feVel_;
  const ROL::Ptr<FE<Real> > fePrs_;
  const ROL::Ptr<FE<Real> > feThr_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Circulation_ThermalFluids(const ROL::Ptr<FE<Real> > &feVel,
                                const ROL::Ptr<FE<Real> > &fePrs,
                                const ROL::Ptr<FE<Real> > &feThr,
                                const ROL::Ptr<FieldHelper<Real> > &fieldHelper)
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

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > gradUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
    feVel_->evaluateGradient(gradUX_eval, U[0]);
    feVel_->evaluateGradient(gradUY_eval, U[1]);
    // Compute curl
    ROL::Ptr<Intrepid::FieldContainer<Real> > curlU_eval =
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

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int fh = feThr_->N()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d+2);
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

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Circulation_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Circulation_ThermalFluids

template <class Real>
class QoI_Horizontal_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > feVel_;
  const ROL::Ptr<FE<Real> > fePrs_;
  const ROL::Ptr<FE<Real> > feThr_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weight_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

public:
  QoI_Horizontal_ThermalFluids(const ROL::Ptr<FE<Real> > &feVel,
                               const ROL::Ptr<FE<Real> > &fePrs,
                               const ROL::Ptr<FE<Real> > &feThr,
                               const ROL::Ptr<FieldHelper<Real> > &fieldHelper)
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

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int p = feVel_->cubPts()->dimension(1);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Compute negative part of x-velocity
    ROL::Ptr<Intrepid::FieldContainer<Real> > minUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*minUX_eval)(i,j) = std::min(static_cast<Real>(0),(*valUX_eval)(i,j));
      }
    }
    // Multiply by weight
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_minUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*weighted_minUX_eval,
                                                               *weight_,
                                                               *minUX_eval);
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_valUY_eval =
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

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c = u_coeff->dimension(0);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int fh = feThr_->N()->dimension(1);
    int p = feVel_->cubPts()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valUY_eval, U[1]);
    // Compute negative part of x-velocity
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_minUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_valUY_eval =
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

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::gradient_2 is zero.");
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    int c  = z_coeff->dimension(0);
    int p  = feVel_->cubPts()->dimension(1);
    int fv = feVel_->N()->dimension(1);
    int fp = fePrs_->N()->dimension(1);
    int fh = feThr_->N()->dimension(1);
    int d = feVel_->cubPts()->dimension(2);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    ROL::Ptr<Intrepid::FieldContainer<Real> > valUX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valVY_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    feVel_->evaluateValue(valUX_eval, U[0]);
    feVel_->evaluateValue(valVX_eval, V[0]);
    feVel_->evaluateValue(valVY_eval, V[1]);
    // Compute negative part of x-velocity
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_minVX_eval =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_valVY_eval =
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

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_Horizontal_ThermalFluids::HessVec_22 is zero.");
  }

}; // QoI_Horizontal_ThermalFluids

template <class Real>
class QoI_State_ThermalFluids : public QoI<Real> {
private:
  ROL::Ptr<QoI<Real> > qoi_;

public:
  QoI_State_ThermalFluids(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<FE<Real> > &feVel,
                         const ROL::Ptr<FE<Real> > &fePrs,
                         const ROL::Ptr<FE<Real> > &feThr,
                         const ROL::Ptr<FieldHelper<Real> > &fieldHelper) {
    std::string stateObj = parlist.sublist("Problem").get("Objective type","Vorticity");
    if ( stateObj != "Vorticity" && stateObj != "Circulation" && stateObj != "Directional" ) {
      throw Exception::NotImplemented(">>> (QoI_State_ThermalFluids): Unknown objective type."); 
    }
    if ( stateObj == "Vorticity" ) {
      qoi_ = ROL::makePtr<QoI_Vorticity_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper,parlist);
    }
    else if ( stateObj == "Directional" ) {
      qoi_ = ROL::makePtr<QoI_Horizontal_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper);
    }
    else {
      qoi_ = ROL::makePtr<QoI_Circulation_ThermalFluids<Real>>(feVel,fePrs,feThr,fieldHelper);
    }
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    return qoi_->value(val, u_coeff, z_coeff, z_param);
  }

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    qoi_->gradient_1(grad, u_coeff, z_coeff, z_param);
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    qoi_->gradient_2(grad, u_coeff, z_coeff, z_param);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    qoi_->HessVec_11(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    qoi_->HessVec_12(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    qoi_->HessVec_21(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    qoi_->HessVec_22(hess, v_coeff, u_coeff, z_coeff, z_param);
  }

};

template <class Real>
class QoI_L2Penalty_ThermalFluids : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real> > feVel_;
  const ROL::Ptr<FE<Real> > fePrs_;
  const ROL::Ptr<FE<Real> > feThr_;
  const std::vector<std::vector<ROL::Ptr<FE<Real> > > > feThrBdry_;
  const std::vector<std::vector<std::vector<int> > > bdryCellLocIds_;
  const ROL::Ptr<FieldHelper<Real> > fieldHelper_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > getBoundaryCoeff(
      const Intrepid::FieldContainer<Real> & cell_coeff,
      int sideSet, int cell) const {
    std::vector<int> bdryCellLocId = bdryCellLocIds_[sideSet][cell];
    const int numCellsSide = bdryCellLocId.size();
    const int f = feThr_->N()->dimension(1);
    
    ROL::Ptr<Intrepid::FieldContainer<Real > > bdry_coeff = 
      ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, f);
    for (int i = 0; i < numCellsSide; ++i) {
      for (int j = 0; j < f; ++j) {
        (*bdry_coeff)(i, j) = cell_coeff(bdryCellLocId[i], j);
      }
    }
    return bdry_coeff;
  }

public:
  QoI_L2Penalty_ThermalFluids(const ROL::Ptr<FE<Real> > &feVel,
                              const ROL::Ptr<FE<Real> > &fePrs,
                              const ROL::Ptr<FE<Real> > &feThr,
                              const std::vector<std::vector<ROL::Ptr<FE<Real> > > > &feThrBdry,
                              const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds,
                              const ROL::Ptr<FieldHelper<Real> > &fieldHelper)
    : feVel_(feVel), fePrs_(fePrs), feThr_(feThr), feThrBdry_(feThrBdry),
      bdryCellLocIds_(bdryCellLocIds), fieldHelper_(fieldHelper) {}

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real> > & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c = feVel_->gradN()->dimension(0);
    const int d = feVel_->gradN()->dimension(3);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    // Compute cost integral
    const int numSideSets = bdryCellLocIds_.size();
    if (numSideSets > 0) {
      for (int i = 0; i < numSideSets; ++i) {
        if ( i == 1 || i == 2 ) {
          const int numLocalSideIds = bdryCellLocIds_[i].size();
          for (int j = 0; j < numLocalSideIds; ++j) {
            const int numCellsSide = bdryCellLocIds_[i][j].size();
            if ( numCellsSide ) {
              const int numCubPerSide = feThrBdry_[i][j]->cubPts()->dimension(1);
              // Evaluate control on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
                = getBoundaryCoeff(*Z[d+1], i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
              // Integrate cell L2 norm squared
              ROL::Ptr<Intrepid::FieldContainer<Real> > intVal
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

  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::gradient_1 is zero.");
  }

  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real> > & grad,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int d  = feThr_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > G(d+2);
    for (int i = 0; i < d; ++i) {
      G[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    G[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    G[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > Z;
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
              const int numCubPerSide = feThrBdry_[i][j]->cubPts()->dimension(1);
              // Evaluate control on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff_bdry
                = getBoundaryCoeff(*Z[d+1], i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real> > valZ_eval
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valZ_eval, z_coeff_bdry);
              // Compute gradient of squared L2-norm of diff
              ROL::Ptr<Intrepid::FieldContainer<Real> > intGrad
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

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_11 is zero.");
  }

  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_12 is zero.");
  }

  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    throw Exception::Zero(">>> QoI_L2Penalty_ThermalFluids::HessVec_21 is zero.");
  }

  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real> > & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real> > & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    const int c  = feVel_->gradN()->dimension(0);
    const int fv = feVel_->gradN()->dimension(1);
    const int fp = fePrs_->gradN()->dimension(1);
    const int fh = feThr_->gradN()->dimension(1);
    const int d  = feThr_->gradN()->dimension(3);
    // Initialize output grad
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > H(d+2);
    for (int i = 0; i < d; ++i) {
      H[i] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fv);
    }
    H[d] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fp);
    H[d+1] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, fh);
    // Get components of the control
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > V;
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
              const int numCubPerSide = feThrBdry_[i][j]->cubPts()->dimension(1);
              // Evaluate control on FE basis
              ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff_bdry
                = getBoundaryCoeff(*V[d+1], i, j);
              ROL::Ptr<Intrepid::FieldContainer<Real> > valV_eval
                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCellsSide, numCubPerSide);
              feThrBdry_[i][j]->evaluateValue(valV_eval, v_coeff_bdry);
              // Compute gradient of squared L2-norm of diff
              ROL::Ptr<Intrepid::FieldContainer<Real> > intHess
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
class StdObjective_ThermalFluids : public ROL::StdObjective<Real> {
private:
  Real alpha_;
  std::string stateObj_;

public:
  StdObjective_ThermalFluids(Teuchos::ParameterList &parlist) {
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
