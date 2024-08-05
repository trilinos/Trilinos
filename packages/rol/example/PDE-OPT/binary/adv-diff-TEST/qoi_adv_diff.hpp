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
  int nsens_;
  Real sig_;
  std::vector<std::vector<Real>> pts_;
  std::vector<Real> vals_;

  ROL::Ptr<Intrepid::FieldContainer<Real>> phi_, dphi_, dphi1_;

  Real evaluateSensor(const std::vector<Real> &pt, const int isens) const {
    const int d = pt.size();
    const Real half(0.5), two(2), pi(M_PI), dr(d);
    Real dist(0), s2 = sig_*sig_, nc = std::pow(two*pi*s2,half*dr);
    for (int i = 0; i < d; ++i) {
      dist += std::pow(pt[i]-pts_[isens][i],two);
    }
    return std::exp(-half*dist/s2)/nc;
  }

  void computeSensor(ROL::Ptr<Intrepid::FieldContainer<Real>> &phi,
                     ROL::Ptr<Intrepid::FieldContainer<Real>> &dphi,
                     ROL::Ptr<Intrepid::FieldContainer<Real>> &dphi1) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    Real phil(0);
    std::vector<Real> pt(d);
    phi->initialize(static_cast<Real>(0));
    dphi->initialize(static_cast<Real>(0));
    dphi1->initialize(static_cast<Real>(0));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        // Compute diffusivity kappa
        for (int l = 0; l < nsens_; ++l) {
          phil = evaluateSensor(pt,l);
          (*phi)(i,j)   += phil;
          (*dphi)(i,j)  += vals_[l]*phil;
          (*dphi1)(i,j) += vals_[l]*std::sqrt(phil);
        }
      }
    }
  }

public:
  QoI_State_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                          ROL::ParameterList &list) : fe_(fe) {
    std::string filename = list.sublist("Problem").get("Sensor File Name", "sensor.txt");
    nsens_ = list.sublist("Problem").get("Number of Sensors",200);
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    const int d = fe_->gradN()->dimension(3);
    sig_   = list.sublist("Problem").get("Sensor Width",1e-2);
    std::fstream file;
    file.open(filename.c_str(),std::ios::in);
    if (file.is_open()) {
      std::vector<Real> pt(d);
      pts_.clear(); pts_.resize(nsens_,pt);
      vals_.clear(); vals_.resize(nsens_);
      for (int i = 0; i < nsens_; ++i) {
        for (int j = 0; j < d; ++j) {
          file >> pts_[i][j];
        }
        file >> vals_[i];
      }
    }
    file.close();
    phi_   = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    dphi_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    dphi1_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    computeSensor(phi_,dphi_,dphi1_);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, val0, Usqr;
    val        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    val0       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    valU_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    Usqr       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // Evaluate state on FE basis
    fe_->evaluateValue(valU_eval, u_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*Usqr,*valU_eval,*valU_eval);
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,phi_,Usqr,false);
    fe_->computeIntegral(val,dphi1_,dphi1_,true);
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    fe_->computeIntegral(val0,dphi_,valU_eval,false);
    Intrepid::RealSpaceTools<Real>::subtract(*val,*val0);
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
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, dU;
    grad      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    dU        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // Evaluate state on FE basis
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*dU, *phi_, *valU_eval);
    Intrepid::RealSpaceTools<Real>::subtract(*dU,*dphi_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *dU,
                                                  *fe_->NdetJ(),
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

  void Hessian_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> H;
    hess      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, f);
    H         = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*H, *phi_, *fe_->N());
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *H,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval, dV;
    hess      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    dV        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valV_eval, v_coeff);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*dV, *phi_, *valV_eval);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *dV,
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

//template <class Real>
//class QoI_State_Cost_adv_diff : public QoI<Real> {
//private:
//  ROL::Ptr<FE<Real>> fe_;
//  int nsens_;
//  Real sig_;
//  std::vector<std::vector<Real>> pts_;
//  std::vector<Real> vals_;
//
//  Real evaluateSensor(const Real val, const std::vector<Real> &pt, const int isens, const int deriv = 0) const {
//    const Real two(2), pi(M_PI);
//    const int d = pt.size();
//    Real dist(0), s2 = sig_*sig_, nc = std::sqrt(std::pow(two*pi*s2,static_cast<Real>(d)));
//    for (int i = 0; i < d; ++i) {
//      dist += std::pow(pt[i]-pts_[isens][i],two);
//    }
//    Real gd = std::exp(-dist/(two*s2))/nc;
//    if (deriv==0) {
//      return gd * std::pow(val-vals_[isens],two);
//    }
//    else if (deriv==1) {
//      return gd * (val-vals_[isens]);
//    }
//    else if (deriv==2) {
//      return gd * val;
//    }
//    else {
//      throw ROL::Exception::NotImplemented(">>> deriv must be 0, 1, or 2!");
//    }
//  }
//
//  void computeSensor(ROL::Ptr<Intrepid::FieldContainer<Real>> &data,
//                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u,
//                     const int deriv = 0) const {
//    // GET DIMENSIONS
//    int c = fe_->gradN()->dimension(0);
//    int p = fe_->gradN()->dimension(2);
//    int d = fe_->gradN()->dimension(3);
//    std::vector<Real> pt(d);
//    data->initialize(static_cast<Real>(0));
//    for (int i = 0; i < c; ++i) {
//      for (int j = 0; j < p; ++j) {
//        for ( int k = 0; k < d; ++k) {
//          pt[k] = (*fe_->cubPts())(i,j,k);
//        }
//        // Compute diffusivity kappa
//        for (int l = 0; l < nsens_; ++l) {
//          (*data)(i,j) += evaluateSensor((*u)(i,j),pt,l,deriv);
//        }
//      }
//    }
//  }
//
//public:
//  QoI_State_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
//                          ROL::ParameterList &list) : fe_(fe) {
//    std::string filename = list.sublist("Problem").get("Sensor File Name", "sensor.txt");
//    nsens_ = list.sublist("Problem").get("Number of Sensors",200);
//    const int d = fe_->gradN()->dimension(3);
//    sig_   = list.sublist("Problem").get("Sensor Width",1e-2);
//    std::fstream file;
//    file.open(filename.c_str(),std::ios::in);
//    if (file.is_open()) {
//      std::vector<Real> pt(d);
//      pts_.clear(); pts_.resize(nsens_,pt);
//      vals_.clear(); vals_.resize(nsens_);
//      for (int i = 0; i < nsens_; ++i) {
//        for (int j = 0; j < d; ++j) {
//          file >> pts_[i][j];
//        }
//        file >> vals_[i];
//      }
//    }
//    file.close();
//  }
//
//  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
//             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    // Get relevant dimensions
//    const int c = fe_->gradN()->dimension(0);
//    const int p = fe_->gradN()->dimension(2);
//    // Initialize storage
//    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, data, ones;
//    val       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
//    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    data      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    ones      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    // Evaluate state on FE basis
//    fe_->evaluateValue(valU_eval, u_coeff);
//    // Compute sensor data
//    computeSensor(data,valU_eval,0);
//    // Compute squared L2-norm of diff
//    ones->initialize(static_cast<Real>(1));
//    fe_->computeIntegral(val,data,ones);
//    // Scale by one half
//    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
//    return static_cast<Real>(0);
//  }
//
//  void gradient_1(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    // Get relevant dimensions
//    const int c = fe_->gradN()->dimension(0);
//    const int f = fe_->gradN()->dimension(1);
//    const int p = fe_->gradN()->dimension(2);
//    // Initialize output grad
//    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, data;
//    grad      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
//    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    data      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    // Evaluate state on FE basis
//    fe_->evaluateValue(valU_eval, u_coeff);
//    // Compute sensor data
//    computeSensor(data,valU_eval,1);
//    // Compute gradient of squared L2-norm of diff
//    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
//                                                  *data,
//                                                  *fe_->NdetJ(),
//                                                  Intrepid::COMP_CPP, false);
//  }
//
//  void gradient_2(ROL::Ptr<Intrepid::FieldContainer<Real>> & grad,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::gradient_2 is zero.");
//  }
//
//  std::vector<Real> gradient_3(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & grad,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::gradient_3 is zero.");
//  }
//
//  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    const int c = fe_->gradN()->dimension(0);
//    const int f = fe_->gradN()->dimension(1);
//    const int p = fe_->gradN()->dimension(2);
//    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval, data;
//    hess      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
//    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    data      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
//    fe_->evaluateValue(valV_eval, v_coeff);
//    // Compute sensor data
//    computeSensor(data,valV_eval,2);
//    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
//                                                  *data,
//                                                  *(fe_->NdetJ()),
//                                                  Intrepid::COMP_CPP, false);
//  }
//
//  void HessVec_12(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_12 is zero.");
//  }
//
//  void HessVec_13(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
//                  const ROL::Ptr<const std::vector<Real>> & v_param,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_13 is zero.");
//  }
//
//  void HessVec_21(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_21 is zero.");
//  }
//
//  void HessVec_22(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_22 is zero.");
//  }
//
//  void HessVec_23(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
//                  const ROL::Ptr<const std::vector<Real>> & v_param,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_23 is zero.");
//  }
//
//  std::vector<Real> HessVec_31(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_31 is zero.");
//  }
//
//  std::vector<Real> HessVec_32(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_32 is zero.");
//  }
//
//  std::vector<Real> HessVec_33(std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> & hess,
//                  const ROL::Ptr<const std::vector<Real>> & v_param,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
//                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
//                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
//    throw Exception::Zero(">>> QoI_State_Cost_stoch_adv_diff::HessVec_33 is zero.");
//  }
//
//}; // QoI_State_Cost

template <class Real>
class QoI_Control_Cost_adv_diff : public QoI<Real> {
private:
  ROL::Ptr<FE<Real>> fe_;
  Real vol_;

public:
  QoI_Control_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                            ROL::ParameterList & parlist) : fe_(fe) {
    int nx  = parlist.sublist("Problem").get("Number of X-Cells", 4);
    int ny  = parlist.sublist("Problem").get("Number of Y-Cells", 2);
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    vol_ = (XU-XL)/static_cast<Real>(nx) * (YU-YL)/static_cast<Real>(ny);
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
    int nx  = parlist.sublist("Problem").get("Number of X-Cells", 4);
    int ny  = parlist.sublist("Problem").get("Number of Y-Cells", 2);
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    vol_ = (XU-XL)/static_cast<Real>(nx) * (YU-YL)/static_cast<Real>(ny);
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

template <class Real>
class QoI_TVControl_Cost_adv_diff : public QoI<Real> {
private:
  ROL::Ptr<FE<Real>> fe_;
  Real volx_,voly_, vol_, eps_;
  int nx_, ny_;

public:
  QoI_TVControl_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                              ROL::ParameterList & parlist) : fe_(fe) {
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    nx_ = parlist.sublist("Problem").get("Number of X-Cells", 4);
    ny_ = parlist.sublist("Problem").get("Number of Y-Cells", 2);
    volx_ = (XU-XL)/static_cast<Real>(nx_);
    voly_ = (YU-YL)/static_cast<Real>(ny_);
    vol_  = volx_*voly_;
    eps_  = parlist.sublist("Problem").get("TV Smoothing Parameter",1e-3);
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    if (z_coeff != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int p = fe_->gradN()->dimension(2);
      const int d = fe_->gradN()->dimension(3);
      // Initialize storage
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradZ_eval, dotZ, ones;
      gradZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      dotZ        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      ones        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      ones->initialize(static_cast<Real>(1));
      // Evaluate state on FE basis
      fe_->evaluateGradient(gradZ_eval, z_coeff);
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*dotZ,*gradZ_eval,*gradZ_eval);
      Real vij(0);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          vij = (*dotZ)(i,j);
          (*dotZ)(i,j) = std::sqrt(vij + eps_);
        }
      }
      // Compute integral
      fe_->computeIntegral(val,dotZ,ones,false);
      return static_cast<Real>(0);
    }
    else {
      const Real half(0.5), two(2);
      Real sum(0);
      std::vector<Real> tmpx(nx_), tmpy(ny_);
      std::vector<std::vector<Real>> Dx(nx_+1,tmpy), Dy(ny_+1,tmpx);
      for (int i = 0; i < nx_+1; ++i) {
        for (int j = 0; j < ny_; ++j) {
          if (i==0) {
            Dx[i][j] = std::pow(-(*z_param)[i+j*nx_]/volx_,two);
          }
          else if (i > 0 && i < nx_) {
            Dx[i][j] = std::pow(((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])/volx_,two);
          }
          else {
            Dx[i][j] = std::pow((*z_param)[(i-1)+j*nx_]/volx_,two);
          }
        }
      }
      for (int i = 0; i < ny_+1; ++i) {
        for (int j = 0; j < nx_; ++j) {
          if (i==0) {
            Dy[i][j] = std::pow(-(*z_param)[j+i*nx_]/voly_,two);
          }
          else if (i > 0 && i < ny_) {
            Dy[i][j] = std::pow(((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])/voly_,two);
          }
          else {
            Dy[i][j] = std::pow((*z_param)[j+(i-1)*nx_]/voly_,two);
          }
        }
      }
      for (int i = 0; i < nx_+1; ++i) {
        for (int j = 0; j < ny_; ++j) {
          if (i==0) {
            sum += std::sqrt(Dx[i][j]+eps_);
          }
          else if (i > 0 && i < nx_) {
            sum += std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
          }
          else {
            sum += std::sqrt(half*Dx[i][j]+eps_);
          }
        }
      }
      for (int i = 0; i < ny_+1; ++i) {
        for (int j = 0; j < nx_; ++j) {
          if (i==0) {
            sum += std::sqrt(Dy[i][j]+eps_);
          }
          else if (i > 0 && i < ny_) {
            sum += std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
          }
          else {
            sum += std::sqrt(half*Dy[i][j]+eps_);
          }
        }
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
      const int d = fe_->gradN()->dimension(3);
      // Initialize storage
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradZ_eval, dotZ, dtv;
      grad       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      gradZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      dotZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      dtv        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      // Evaluate state on FE basis
      fe_->evaluateGradient(gradZ_eval, z_coeff);
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*dotZ,*gradZ_eval,*gradZ_eval);
      Real vij(0);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          vij = (*dotZ)(i,j);
          (*dotZ)(i,j) = static_cast<Real>(1)/std::sqrt(vij + eps_);
        }
      }
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*dtv,*dotZ,*gradZ_eval);
      // Compute integral
      Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                    *dtv,
                                                    *fe_->gradNdetJ(),
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
      std::vector<Real> g(size,static_cast<Real>(0));
      const Real half(0.5), two(2), volx2 = volx_*volx_, voly2 = voly_*voly_;
      std::vector<Real> tmpx(nx_), tmpy(ny_);
      std::vector<std::vector<Real>> Dx(nx_+1,tmpy), Dy(ny_+1,tmpx);
      for (int i = 0; i < nx_+1; ++i) {
        for (int j = 0; j < ny_; ++j) {
          if (i==0) {
            Dx[i][j] = std::pow(-(*z_param)[i+j*nx_]/volx_,two);
          }
          else if (i > 0 && i < nx_) {
            Dx[i][j] = std::pow(((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])/volx_,two);
          }
          else {
            Dx[i][j] = std::pow((*z_param)[(i-1)+j*nx_]/volx_,two);
          }
        }
      }
      for (int i = 0; i < ny_+1; ++i) {
        for (int j = 0; j < nx_; ++j) {
          if (i==0) {
            Dy[i][j] = std::pow(-(*z_param)[j+i*nx_]/voly_,two);
          }
          else if (i > 0 && i < ny_) {
            Dy[i][j] = std::pow(((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])/voly_,two);
          }
          else {
            Dy[i][j] = std::pow((*z_param)[j+(i-1)*nx_]/voly_,two);
          }
        }
      }
      for (int i = 0; i < nx_+1; ++i) {
        for (int j = 0; j < ny_; ++j) {
          if (i==0) {
            Real cx = (vol_/volx2)/std::sqrt(Dx[i][j]+eps_);
            g[i+j*nx_] += cx*(*z_param)[i+j*nx_];
          }
          else if (i > 0 && i < nx_-1) {
            Real cx = half*vol_/(volx2*std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_));
            g[(i-1)+j*nx_] -= cx*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_]);
            g[i+j*nx_]     += cx*(two*(*z_param)[i+j*nx_]-(*z_param)[(i+1)+j*nx_]-(*z_param)[(i-1)+j*nx_]);
            g[(i+1)+j*nx_] += cx*((*z_param)[(i+1)+j*nx_]-(*z_param)[i+j*nx_]);
          }
          else if (i == nx_-1) {
            Real cx = half*vol_/(volx2*std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_));
            g[(i-1)+j*nx_] -= cx*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_]);
            g[i+j*nx_]     += cx*(two*(*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_]);
          }
          else if (i==nx_) {
            Real cx = vol_*half/(volx2*std::sqrt(half*Dx[i][j]+eps_));
            g[(i-1)+j*nx_] += cx*(*z_param)[(i-1)+j*nx_];
          }
        }
      }
      for (int i = 0; i < ny_+1; ++i) {
        for (int j = 0; j < nx_; ++j) {
          if (i==0) {
            Real cy = (vol_/voly2)/std::sqrt(Dy[i][j]+eps_);
            g[j+i*nx_] += cy*(*z_param)[j+i*nx_];
          }
          else if (i > 0 && i < ny_-1) {
            Real cy = half*vol_/(voly2*std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_));
            g[j+(i-1)*nx_] -= cy*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_]);
            g[j+i*nx_]     += cy*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i+1)*nx_]-(*z_param)[j+(i-1)*nx_]);
            g[j+(i+1)*nx_] += cy*((*z_param)[j+(i+1)*nx_]-(*z_param)[j+i*nx_]);
          }
          else if (i == ny_-1) {
            Real cy = half*vol_/(voly2*std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_));
            g[j+(i-1)*nx_] -= cy*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_]);
            g[j+i*nx_]     += cy*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_]);
          }
          else if (i==ny_) {
            Real cy = vol_*half/(voly2*std::sqrt(half*Dy[i][j]+eps_));
            g[j+(i-1)*nx_] += cy*(*z_param)[j+(i-1)*nx_];
          }
        }
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
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int f = fe_->gradN()->dimension(1);
      const int p = fe_->gradN()->dimension(2);
      const int d = fe_->gradN()->dimension(3);
      // Initialize storage
      ROL::Ptr<Intrepid::FieldContainer<Real>> gradZ_eval, gradV_eval, dotZ, dotZ3, dotZV, dtv, dtv2;
      hess       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      gradZ_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      gradV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      dotZ       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      dotZ3      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      dotZV      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      dtv        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      dtv2       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p, d);
      // Evaluate state on FE basis
      fe_->evaluateGradient(gradZ_eval, z_coeff);
      fe_->evaluateGradient(gradV_eval, v_coeff);
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*dotZ,*gradZ_eval,*gradZ_eval);
      Intrepid::FunctionSpaceTools::dotMultiplyDataData<Real>(*dotZV,*gradZ_eval,*gradV_eval);
      Real vij(0);
      for (int i = 0; i < c; ++i) {
        for (int j = 0; j < p; ++j) {
          vij = (*dotZ)(i,j);
          (*dotZ)(i,j)  = static_cast<Real>(1)/std::sqrt(vij + eps_);
          (*dotZ3)(i,j) = -(*dotZV)(i,j)/std::pow(vij + eps_,static_cast<Real>(1.5));
        }
      }
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*dtv,*dotZ,*gradV_eval);
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*dtv2,*dotZ3,*gradZ_eval);
      Intrepid::RealSpaceTools<Real>::add(*dtv,*dtv2);
      // Compute integral
      Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                    *dtv,
                                                    *fe_->gradNdetJ(),
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
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::hessVec_33 is zero.");
    }
    else {
      const int size = z_param->size();
      std::vector<Real> h(size,static_cast<Real>(0));
      const Real half(0.5), two(2), three(3), volx2 = volx_*volx_, voly2 = voly_*voly_;
      std::vector<Real> tmpx(nx_), tmpy(ny_);
      std::vector<std::vector<Real>> Dx(nx_+1,tmpy), Dy(ny_+1,tmpx);
      for (int i = 0; i < nx_+1; ++i) {
        for (int j = 0; j < ny_; ++j) {
          if (i==0) {
            Dx[i][j] = std::pow(-(*z_param)[i+j*nx_]/volx_,two);
          }
          else if (i > 0 && i < nx_) {
            Dx[i][j] = std::pow(((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])/volx_,two);
          }
          else {
            Dx[i][j] = std::pow((*z_param)[(i-1)+j*nx_]/volx_,two);
          }
        }
      }
      for (int i = 0; i < ny_+1; ++i) {
        for (int j = 0; j < nx_; ++j) {
          if (i==0) {
            Dy[i][j] = std::pow(-(*z_param)[j+i*nx_]/voly_,two);
          }
          else if (i > 0 && i < ny_) {
            Dy[i][j] = std::pow(((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])/voly_,two);
          }
          else {
            Dy[i][j] = std::pow((*z_param)[j+(i-1)*nx_]/voly_,two);
          }
        }
      }
      for (int i = 0; i < nx_+1; ++i) {
        for (int j = 0; j < ny_; ++j) {
          if (i==0) {
            Real cx1 = (vol_/volx2)/std::sqrt(Dx[i][j]+eps_);
            Real cx2 = (vol_/volx2)*(two/volx2)*(-half/std::pow(std::sqrt(Dx[i][j]+eps_),three));
            h[i+j*nx_] += (cx1+cx2*std::pow((*z_param)[i+j*nx_],two))*(*v_param)[i+j*nx_];
          }
          else if (i > 0 && i < nx_-1) {
            Real cx1 = (half*vol_/volx2)/std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
            Real cx2 = (half*vol_/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_),three));
            h[(i-1)+j*nx_] += cx2*std::pow((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_],two)*(*v_param)[(i-1)+j*nx_]
                              -cx2*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*(two*(*z_param)[i+j*nx_]-(*z_param)[(i+1)+j*nx_]-(*z_param)[(i-1)+j*nx_])*(*v_param)[i+j*nx_]
                              -cx2*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*((*z_param)[(i+1)+j*nx_]-(*z_param)[i+j*nx_])*(*v_param)[(i+1)+j*nx_]
                              -cx1*((*v_param)[i+j*nx_]-(*v_param)[(i-1)+j*nx_]);
            h[i+j*nx_]     += cx2*std::pow(two*(*z_param)[i+j*nx_]-(*z_param)[(i+1)+j*nx_]-(*z_param)[(i-1)+j*nx_],two)*(*v_param)[i+j*nx_]
                              -cx2*(two*(*z_param)[i+j*nx_]-(*z_param)[(i+1)+j*nx_]-(*z_param)[(i-1)+j*nx_])*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*(*v_param)[(i-1)+j*nx_]
                              +cx2*(two*(*z_param)[i+j*nx_]-(*z_param)[(i+1)+j*nx_]-(*z_param)[(i-1)+j*nx_])*((*z_param)[(i+1)+j*nx_]-(*z_param)[i+j*nx_])*(*v_param)[(i+1)+j*nx_]
                              +cx1*(two*(*v_param)[i+j*nx_]-(*v_param)[(i+1)+j*nx_]-(*v_param)[(i-1)+j*nx_]);
            h[(i+1)+j*nx_] += cx2*std::pow((*z_param)[(i+1)+j*nx_]-(*z_param)[i+j*nx_],two)*(*v_param)[(i+1)+j*nx_]
                              -cx2*((*z_param)[(i+1)+j*nx_]-(*z_param)[i+j*nx_])*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*(*v_param)[(i-1)+j*nx_]
                              +cx2*((*z_param)[(i+1)+j*nx_]-(*z_param)[i+j*nx_])*(two*(*z_param)[i+j*nx_]-(*z_param)[(i+1)+j*nx_]-(*z_param)[(i-1)+j*nx_])*(*v_param)[i+j*nx_]
                              +cx1*((*v_param)[(i+1)+j*nx_]-(*v_param)[i+j*nx_]);
          }
          else if (i == nx_-1) {
            Real cx1 = (half*vol_/volx2)/std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
            Real cx2 = (half*vol_/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_),three));
            h[(i-1)+j*nx_] += cx2*std::pow((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_],two)*(*v_param)[(i-1)+j*nx_]
                              -cx2*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*(two*(*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*(*v_param)[i+j*nx_]
                              -cx1*((*v_param)[i+j*nx_]-(*v_param)[(i-1)+j*nx_]);
            h[i+j*nx_]     += cx2*std::pow(two*(*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_],two)*(*v_param)[i+j*nx_]
                              -cx2*(two*(*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*((*z_param)[i+j*nx_]-(*z_param)[(i-1)+j*nx_])*(*v_param)[(i-1)+j*nx_]
                              +cx1*(two*(*v_param)[i+j*nx_]-(*v_param)[(i-1)+j*nx_]);
          }
          else if (i==nx_) {
            Real cx1 = (vol_*half/volx2)/std::sqrt(half*Dx[i][j]+eps_);
            Real cx2 = (vol_*half/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*Dx[i][j]+eps_),three));
            h[(i-1)+j*nx_] += (cx1+cx2*std::pow((*z_param)[(i-1)+j*nx_],two))*(*v_param)[(i-1)+j*nx_];
          }
        }
      }
      for (int i = 0; i < ny_+1; ++i) {
        for (int j = 0; j < nx_; ++j) {
          if (i==0) {
            Real cy1 = (vol_/voly2)/std::sqrt(Dy[i][j]+eps_);
            Real cy2 = (vol_/voly2)*(two/voly2)*(-half/std::pow(std::sqrt(Dy[i][j]+eps_),three));
            h[j+i*nx_] += (cy1+cy2*std::pow((*z_param)[j+i*nx_],two))*(*v_param)[j+i*nx_];
          }
          else if (i > 0 && i < ny_-1) {
            Real cy1 = (half*vol_/voly2)/std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
            Real cy2 = (half*vol_/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_),three));
            h[j+(i-1)*nx_] += cy2*std::pow((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_],two)*(*v_param)[j+(i-1)*nx_]
                              -cy2*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i+1)*nx_]-(*z_param)[j+(i-1)*nx_])*(*v_param)[j+i*nx_]
                              -cy2*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*((*z_param)[j+(i+1)*nx_]-(*z_param)[j+i*nx_])*(*v_param)[j+(i+1)*nx_]
                              -cy1*((*v_param)[j+i*nx_]-(*v_param)[j+(i-1)*nx_]);
            h[j+i*nx_]     += cy2*std::pow(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i+1)*nx_]-(*z_param)[j+(i-1)*nx_],two)*(*v_param)[j+i*nx_]
                              -cy2*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i+1)*nx_]-(*z_param)[j+(i-1)*nx_])*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*(*v_param)[j+(i-1)*nx_]
                              +cy2*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i+1)*nx_]-(*z_param)[j+(i-1)*nx_])*((*z_param)[j+(i+1)*nx_]-(*z_param)[j+i*nx_])*(*v_param)[j+(i+1)*nx_]
                              +cy1*(two*(*v_param)[j+i*nx_]-(*v_param)[j+(i+1)*nx_]-(*v_param)[j+(i-1)*nx_]);
            h[j+(i+1)*nx_] += cy2*std::pow((*z_param)[j+(i+1)*nx_]-(*z_param)[j+i*nx_],two)*(*v_param)[j+(i+1)*nx_]
                              -cy2*((*z_param)[j+(i+1)*nx_]-(*z_param)[j+i*nx_])*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*(*v_param)[j+(i-1)*nx_]
                              +cy2*((*z_param)[j+(i+1)*nx_]-(*z_param)[j+i*nx_])*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i+1)*nx_]-(*z_param)[j+(i-1)*nx_])*(*v_param)[j+i*nx_]
                              +cy1*((*v_param)[j+(i+1)*nx_]-(*v_param)[j+i*nx_]);
          }
          else if (i == ny_-1) {
            Real cy1 = (half*vol_/voly2)/std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
            Real cy2 = (half*vol_/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_),three));
            h[j+(i-1)*nx_] += cy2*std::pow((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_],two)*(*v_param)[j+(i-1)*nx_]
                              -cy2*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*(*v_param)[j+i*nx_]
                              -cy1*((*v_param)[j+i*nx_]-(*v_param)[j+(i-1)*nx_]);
            h[j+i*nx_]     += cy2*std::pow(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_],two)*(*v_param)[j+i*nx_]
                              -cy2*(two*(*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*((*z_param)[j+i*nx_]-(*z_param)[j+(i-1)*nx_])*(*v_param)[j+(i-1)*nx_]
                              +cy1*(two*(*v_param)[j+i*nx_]-(*v_param)[j+(i-1)*nx_]);
          }
          else if (i==ny_) {
            Real cy1 = (vol_*half/voly2)/std::sqrt(half*Dy[i][j]+eps_);
            Real cy2 = (vol_*half/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*Dy[i][j]+eps_),three));
            h[j+(i-1)*nx_] += (cy1+cy2*std::pow((*z_param)[j+(i-1)*nx_],two))*(*v_param)[j+(i-1)*nx_];
          }
        }
      }
      return h;
    }
  }

}; // QoI_Control_Cost

template <class Real>
class QoI_IntegralityControl_Cost_adv_diff : public QoI<Real> {
private:
  const ROL::Ptr<FE<Real>> fe_;
  Real volx_,voly_, vol_;
  int nx_, ny_;

public:
  QoI_IntegralityControl_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                                       ROL::ParameterList & parlist) : fe_(fe) {
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    nx_ = parlist.sublist("Problem").get("Number of X-Cells", 4);
    ny_ = parlist.sublist("Problem").get("Number of Y-Cells", 2);
    volx_ = (XU-XL)/static_cast<Real>(nx_);
    voly_ = (YU-YL)/static_cast<Real>(ny_);
    vol_  = volx_*voly_;
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    // Initialize output val
    val = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    if (z_coeff != ROL::nullPtr) {
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int p = fe_->gradN()->dimension(2);
      // Initialize storage
      ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, Zdiff;
      valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Zdiff      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Zdiff->initialize(static_cast<Real>(1));
      // Evaluate state on FE basis
      fe_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::subtract(*Zdiff,*valZ_eval);
      // Compute integral
      fe_->computeIntegral(val,valZ_eval,Zdiff,false);
      return static_cast<Real>(0);
    }
    else {
      const Real one(1);
      Real sum(0);
      for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
          sum += (*z_param)[i+j*nx_]*(one-(*z_param)[i+j*nx_]);
        }
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
      // Initialize storage
      ROL::Ptr<Intrepid::FieldContainer<Real>> valZ_eval, Zdiff;
      grad       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      valZ_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Zdiff      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      Zdiff->initialize(static_cast<Real>(1));
      // Evaluate state on FE basis
      fe_->evaluateValue(valZ_eval, z_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valZ_eval,static_cast<Real>(2));
      Intrepid::RealSpaceTools<Real>::subtract(*Zdiff,*valZ_eval);
      // Compute integral
      Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                    *Zdiff,
                                                    *fe_->NdetJ(),
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
      std::vector<Real> g(size,static_cast<Real>(0));
      const Real one(1), two(2);
      for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
          g[i+j*nx_] = vol_*(one - two*(*z_param)[i+j*nx_]);
        }
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
      // Get relevant dimensions
      const int c = fe_->gradN()->dimension(0);
      const int f = fe_->gradN()->dimension(1);
      const int p = fe_->gradN()->dimension(2);
      // Initialize storage
      ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval;
      hess       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
      valV_eval  = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
      // Evaluate state on FE basis
      fe_->evaluateValue(valV_eval, v_coeff);
      Intrepid::RealSpaceTools<Real>::scale(*valV_eval,static_cast<Real>(-2));
      // Compute integral
      Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                    *valV_eval,
                                                    *fe_->NdetJ(),
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
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::hessVec_33 is zero.");
    }
    else {
      const int size = z_param->size();
      std::vector<Real> h(size,static_cast<Real>(0));
      const Real two(2);
      for (int i = 0; i < nx_; ++i) {
        for (int j = 0; j < ny_; ++j) {
          h[i+j*nx_] = -vol_*two*(*v_param)[i+j*nx_];
        }
      }
      return h;
    }
  }

}; // QoI_Control_Cost
#endif
