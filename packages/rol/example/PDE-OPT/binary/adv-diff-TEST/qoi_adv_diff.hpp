// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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

  Real evaluateSensor(const Real val, const std::vector<Real> &pt, const int isens, const int deriv = 0) const {
    const Real two(2), pi(M_PI);
    Real dist(0), s2 = sig_*sig_, nc = std::sqrt(two*pi)*sig_;
    int d = pt.size();
    for (int i = 0; i < d; ++i) {
      dist += std::pow(pt[i]-pts_[isens][i],two);
    }
    Real gd = std::exp(-dist/(two*s2))/nc;
    if (deriv==0) {
      return gd * std::pow(val-vals_[isens],two);
    }
    else if (deriv==1) {
      return gd * (val-vals_[isens]);
    }
    else if (deriv==2) {
      return gd * val;
    }
    else {
      throw ROL::Exception::NotImplemented(">>> deriv must be 0, 1, or 2!");
    }
  }

  void computeSensor(ROL::Ptr<Intrepid::FieldContainer<Real>> &data,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real>> &u,
                     const int deriv = 0) const {
    // GET DIMENSIONS
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
    std::vector<Real> pt(d);
    data->initialize(static_cast<Real>(0));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for ( int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        // Compute diffusivity kappa
        for (int l = 0; l < nsens_; ++l) {
          (*data)(i,j) += evaluateSensor((*u)(i,j),pt,l,deriv);
        }
      }
    }
  }

public:
  QoI_State_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                          ROL::ParameterList &list) : fe_(fe) {
    std::string filename = list.sublist("Problem").get("Sensor File Name", "sensor.txt");
    nsens_ = list.sublist("Problem").get("Number of Sensors",200);
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
  }

  Real value(ROL::Ptr<Intrepid::FieldContainer<Real>> & val,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
             const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
             const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize storage
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, data, ones;
    val       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    data      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    ones      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // Evaluate state on FE basis
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute sensor data
    computeSensor(data,valU_eval,0);
    // Compute squared L2-norm of diff
    ones->initialize(static_cast<Real>(1));
    fe_->computeIntegral(val,data,ones);
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
    ROL::Ptr<Intrepid::FieldContainer<Real>> valU_eval, data;
    grad      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    valU_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    data      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    // Evaluate state on FE basis
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute sensor data
    computeSensor(data,valU_eval,1);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  *data,
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

  void HessVec_11(ROL::Ptr<Intrepid::FieldContainer<Real>> & hess,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & v_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & u_coeff,
                  const ROL::Ptr<const Intrepid::FieldContainer<Real>> & z_coeff = ROL::nullPtr,
                  const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr) {
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    ROL::Ptr<Intrepid::FieldContainer<Real>> valV_eval, data;
    hess      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, f);
    valV_eval = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    data      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe_->evaluateValue(valV_eval, v_coeff);
    // Compute sensor data
    computeSensor(data,valV_eval,2);
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  *data,
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

template <class Real>
class QoI_TVControl_Cost_adv_diff : public QoI<Real> {
private:
  ROL::Ptr<FE<Real>> fe_;
  Real volx_,voly_, vol_, eps_;
  int n_;

public:
  QoI_TVControl_Cost_adv_diff(const ROL::Ptr<FE<Real>> &fe,
                              ROL::ParameterList & parlist) : fe_(fe) {
    int order = parlist.sublist("Problem").get("Hilbert Curve Order", 2);
    Real XL = parlist.sublist("Geometry").get("X0", 0.0);
    Real YL = parlist.sublist("Geometry").get("Y0", 0.0);
    Real XU = XL + parlist.sublist("Geometry").get("Width",  1.0);
    Real YU = YL + parlist.sublist("Geometry").get("Height", 1.0);
    n_    = std::pow(2,order);
    volx_ = (XU-XL)/static_cast<Real>(n_);
    voly_ = (YU-YL)/static_cast<Real>(n_);
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
    }
    else {
      const Real half(0.5), two(2);
      Real sum(0);
      std::vector<Real> tmp(n_);
      std::vector<std::vector<Real>> Dx(n_+1,tmp), Dy(n_+1,tmp);
      for (int i = 0; i < n_+1; ++i) {
        for (int j = 0; j < n_; ++j) {
          if (i==0) {
            Dx[i][j] = std::pow(-(*z_param)[i+j*n_]/volx_,two);
            Dy[i][j] = std::pow(-(*z_param)[j+i*n_]/voly_,two);
          }
          else if (i > 0 && i < n_) {
            Dx[i][j] = std::pow(((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])/volx_,two);
            Dy[i][j] = std::pow(((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])/voly_,two);
          }
          else {
            Dx[i][j] = std::pow((*z_param)[(i-1)+j*n_]/volx_,two);
            Dy[i][j] = std::pow((*z_param)[j+(i-1)*n_]/voly_,two);
          }
        }
      }
      for (int i = 0; i < n_+1; ++i) {
        for (int j = 0; j < n_; ++j) {
          if (i==0) {
            sum += std::sqrt(Dx[i][j]+eps_);
            sum += std::sqrt(Dy[i][j]+eps_);
          }
          else if (i > 0 && i < n_) {
            sum += std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
            sum += std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
          }
          else {
            sum += std::sqrt(half*Dx[i][j]+eps_);
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
    throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::gradient_2 is zero.");
    if (z_coeff != ROL::nullPtr) {
    }
    else {
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
      std::vector<Real> tmp(n_);
      std::vector<std::vector<Real>> Dx(n_+1,tmp), Dy(n_+1,tmp);
      for (int i = 0; i < n_+1; ++i) {
        for (int j = 0; j < n_; ++j) {
          if (i==0) {
            Dx[i][j] = std::pow(-(*z_param)[i+j*n_]/volx_,two);
            Dy[i][j] = std::pow(-(*z_param)[j+i*n_]/voly_,two);
          }
          else if (i > 0 && i < n_) {
            Dx[i][j] = std::pow(((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])/volx_,two);
            Dy[i][j] = std::pow(((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])/voly_,two);
          }
          else {
            Dx[i][j] = std::pow((*z_param)[(i-1)+j*n_]/volx_,two);
            Dy[i][j] = std::pow((*z_param)[j+(i-1)*n_]/voly_,two);
          }
        }
      }
      for (int i = 0; i < n_+1; ++i) {
        for (int j = 0; j < n_; ++j) {
          if (i==0) {
            Real cx = (vol_/volx2)/std::sqrt(Dx[i][j]+eps_);
            g[i+j*n_] += cx*(*z_param)[i+j*n_];

            Real cy = (vol_/voly2)/std::sqrt(Dy[i][j]+eps_);
            g[j+i*n_] += cy*(*z_param)[j+i*n_];
          }
          else if (i > 0 && i < n_-1) {
            Real cx = half*vol_/(volx2*std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_));
            g[(i-1)+j*n_] -= cx*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_]);
            g[i+j*n_]     += cx*(two*(*z_param)[i+j*n_]-(*z_param)[(i+1)+j*n_]-(*z_param)[(i-1)+j*n_]);
            g[(i+1)+j*n_] += cx*((*z_param)[(i+1)+j*n_]-(*z_param)[i+j*n_]);

            Real cy = half*vol_/(voly2*std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_));
            g[j+(i-1)*n_] -= cy*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_]);
            g[j+i*n_]     += cy*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i+1)*n_]-(*z_param)[j+(i-1)*n_]);
            g[j+(i+1)*n_] += cy*((*z_param)[j+(i+1)*n_]-(*z_param)[j+i*n_]);
          }
          else if (i == n_-1) {
            Real cx = half*vol_/(volx2*std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_));
            g[(i-1)+j*n_] -= cx*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_]);
            g[i+j*n_]     += cx*(two*(*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_]);

            Real cy = half*vol_/(voly2*std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_));
            g[j+(i-1)*n_] -= cy*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_]);
            g[j+i*n_]     += cy*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_]);
          }
          else if (i==n_) {
            Real cx = vol_*half/(volx2*std::sqrt(half*Dx[i][j]+eps_));
            g[(i-1)+j*n_] += cx*(*z_param)[(i-1)+j*n_];

            Real cy = vol_*half/(voly2*std::sqrt(half*Dy[i][j]+eps_));
            g[j+(i-1)*n_] += cy*(*z_param)[j+(i-1)*n_];
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
    if (z_coeff != ROL::nullPtr) {
      throw Exception::Zero(">>> QoI_Control_Cost_stoch_adv_diff::hessVec_33 is zero.");
    }
    else {
      const int size = z_param->size();
      std::vector<Real> h(size,static_cast<Real>(0));
      const Real half(0.5), two(2), three(3), volx2 = volx_*volx_, voly2 = voly_*voly_;
      std::vector<Real> tmp(n_);
      std::vector<std::vector<Real>> Dx(n_+1,tmp), Dy(n_+1,tmp);
      for (int i = 0; i < n_+1; ++i) {
        for (int j = 0; j < n_; ++j) {
          if (i==0) {
            Dx[i][j] = std::pow(-(*z_param)[i+j*n_]/volx_,two);
            Dy[i][j] = std::pow(-(*z_param)[j+i*n_]/voly_,two);
          }
          else if (i > 0 && i < n_) {
            Dx[i][j] = std::pow(((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])/volx_,two);
            Dy[i][j] = std::pow(((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])/voly_,two);
          }
          else {
            Dx[i][j] = std::pow((*z_param)[(i-1)+j*n_]/volx_,two);
            Dy[i][j] = std::pow((*z_param)[j+(i-1)*n_]/voly_,two);
          }
        }
      }
      for (int i = 0; i < n_+1; ++i) {
        for (int j = 0; j < n_; ++j) {
          if (i==0) {
            Real cx1 = (vol_/volx2)/std::sqrt(Dx[i][j]+eps_);
            Real cx2 = (vol_/volx2)*(two/volx2)*(-half/std::pow(std::sqrt(Dx[i][j]+eps_),three));
            h[i+j*n_] += (cx1+cx2*std::pow((*z_param)[i+j*n_],two))*(*v_param)[i+j*n_];

            Real cy1 = (vol_/voly2)/std::sqrt(Dy[i][j]+eps_);
            Real cy2 = (vol_/voly2)*(two/voly2)*(-half/std::pow(std::sqrt(Dy[i][j]+eps_),three));
            h[j+i*n_] += (cy1+cy2*std::pow((*z_param)[j+i*n_],two))*(*v_param)[j+i*n_];
          }
          else if (i > 0 && i < n_-1) {
            Real cx1 = (half*vol_/volx2)/std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
            Real cx2 = (half*vol_/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_),three));
            h[(i-1)+j*n_] += cx2*std::pow((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_],two)*(*v_param)[(i-1)+j*n_]
                              -cx2*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*(two*(*z_param)[i+j*n_]-(*z_param)[(i+1)+j*n_]-(*z_param)[(i-1)+j*n_])*(*v_param)[i+j*n_]
                              -cx2*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*((*z_param)[(i+1)+j*n_]-(*z_param)[i+j*n_])*(*v_param)[(i+1)+j*n_]
                              -cx1*((*v_param)[i+j*n_]-(*v_param)[(i-1)+j*n_]);
            h[i+j*n_]     += cx2*std::pow(two*(*z_param)[i+j*n_]-(*z_param)[(i+1)+j*n_]-(*z_param)[(i-1)+j*n_],two)*(*v_param)[i+j*n_]
                              -cx2*(two*(*z_param)[i+j*n_]-(*z_param)[(i+1)+j*n_]-(*z_param)[(i-1)+j*n_])*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*(*v_param)[(i-1)+j*n_]
                              +cx2*(two*(*z_param)[i+j*n_]-(*z_param)[(i+1)+j*n_]-(*z_param)[(i-1)+j*n_])*((*z_param)[(i+1)+j*n_]-(*z_param)[i+j*n_])*(*v_param)[(i+1)+j*n_]
                              +cx1*(two*(*v_param)[i+j*n_]-(*v_param)[(i+1)+j*n_]-(*v_param)[(i-1)+j*n_]);
            h[(i+1)+j*n_] += cx2*std::pow((*z_param)[(i+1)+j*n_]-(*z_param)[i+j*n_],two)*(*v_param)[(i+1)+j*n_]
                              -cx2*((*z_param)[(i+1)+j*n_]-(*z_param)[i+j*n_])*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*(*v_param)[(i-1)+j*n_]
                              +cx2*((*z_param)[(i+1)+j*n_]-(*z_param)[i+j*n_])*(two*(*z_param)[i+j*n_]-(*z_param)[(i+1)+j*n_]-(*z_param)[(i-1)+j*n_])*(*v_param)[i+j*n_]
                              +cx1*((*v_param)[(i+1)+j*n_]-(*v_param)[i+j*n_]);

            Real cy1 = (half*vol_/voly2)/std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
            Real cy2 = (half*vol_/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_),three));
            h[j+(i-1)*n_] += cy2*std::pow((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_],two)*(*v_param)[j+(i-1)*n_]
                              -cy2*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i+1)*n_]-(*z_param)[j+(i-1)*n_])*(*v_param)[j+i*n_]
                              -cy2*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*((*z_param)[j+(i+1)*n_]-(*z_param)[j+i*n_])*(*v_param)[j+(i+1)*n_]
                              -cy1*((*v_param)[j+i*n_]-(*v_param)[j+(i-1)*n_]);
            h[j+i*n_]     += cy2*std::pow(two*(*z_param)[j+i*n_]-(*z_param)[j+(i+1)*n_]-(*z_param)[j+(i-1)*n_],two)*(*v_param)[j+i*n_]
                              -cy2*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i+1)*n_]-(*z_param)[j+(i-1)*n_])*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*(*v_param)[j+(i-1)*n_]
                              +cy2*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i+1)*n_]-(*z_param)[j+(i-1)*n_])*((*z_param)[j+(i+1)*n_]-(*z_param)[j+i*n_])*(*v_param)[j+(i+1)*n_]
                              +cy1*(two*(*v_param)[j+i*n_]-(*v_param)[j+(i+1)*n_]-(*v_param)[j+(i-1)*n_]);
            h[j+(i+1)*n_] += cy2*std::pow((*z_param)[j+(i+1)*n_]-(*z_param)[j+i*n_],two)*(*v_param)[j+(i+1)*n_]
                              -cy2*((*z_param)[j+(i+1)*n_]-(*z_param)[j+i*n_])*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*(*v_param)[j+(i-1)*n_]
                              +cy2*((*z_param)[j+(i+1)*n_]-(*z_param)[j+i*n_])*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i+1)*n_]-(*z_param)[j+(i-1)*n_])*(*v_param)[j+i*n_]
                              +cy1*((*v_param)[j+(i+1)*n_]-(*v_param)[j+i*n_]);
          }
          else if (i == n_-1) {
            Real cx1 = (half*vol_/volx2)/std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_);
            Real cx2 = (half*vol_/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*(Dx[i][j]+Dx[i+1][j])+eps_),three));
            h[(i-1)+j*n_] += cx2*std::pow((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_],two)*(*v_param)[(i-1)+j*n_]
                              -cx2*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*(two*(*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*(*v_param)[i+j*n_]
                              -cx1*((*v_param)[i+j*n_]-(*v_param)[(i-1)+j*n_]);
            h[i+j*n_]     += cx2*std::pow(two*(*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_],two)*(*v_param)[i+j*n_]
                              -cx2*(two*(*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*((*z_param)[i+j*n_]-(*z_param)[(i-1)+j*n_])*(*v_param)[(i-1)+j*n_]
                              +cx1*(two*(*v_param)[i+j*n_]-(*v_param)[(i-1)+j*n_]);

            Real cy1 = (half*vol_/voly2)/std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_);
            Real cy2 = (half*vol_/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*(Dy[i][j]+Dy[i+1][j])+eps_),three));
            h[j+(i-1)*n_] += cy2*std::pow((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_],two)*(*v_param)[j+(i-1)*n_]
                              -cy2*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*(*v_param)[j+i*n_]
                              -cy1*((*v_param)[j+i*n_]-(*v_param)[j+(i-1)*n_]);
            h[j+i*n_]     += cy2*std::pow(two*(*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_],two)*(*v_param)[j+i*n_]
                              -cy2*(two*(*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*((*z_param)[j+i*n_]-(*z_param)[j+(i-1)*n_])*(*v_param)[j+(i-1)*n_]
                              +cy1*(two*(*v_param)[j+i*n_]-(*v_param)[j+(i-1)*n_]);
          }
          else if (i==n_) {
            Real cx1 = (vol_*half/volx2)/std::sqrt(half*Dx[i][j]+eps_);
            Real cx2 = (vol_*half/volx2)*(two/volx2)*(-half*half/std::pow(std::sqrt(half*Dx[i][j]+eps_),three));
            h[(i-1)+j*n_] += (cx1+cx2*std::pow((*z_param)[(i-1)+j*n_],two))*(*v_param)[(i-1)+j*n_];

            Real cy1 = (vol_*half/voly2)/std::sqrt(half*Dy[i][j]+eps_);
            Real cy2 = (vol_*half/voly2)*(two/voly2)*(-half*half/std::pow(std::sqrt(half*Dy[i][j]+eps_),three));
            h[j+(i-1)*n_] += (cy1+cy2*std::pow((*z_param)[j+(i-1)*n_],two))*(*v_param)[j+(i-1)*n_];
          }
        }
      }
      return h;
    }
  }

}; // QoI_Control_Cost
#endif
