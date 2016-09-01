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

#ifndef PDEOPT_QOI_STOCHASTIC_STEFAN_BOLTZMANN_HPP
#define PDEOPT_QOI_STOCHASTIC_STEFAN_BOLTZMANN_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_stoch_stefan_boltzmann.hpp"

template <class Real>
class QoI_StateCost : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > weight_;
  Real xmid_;
  Real T_;

  Real weightFunc(const std::vector<Real> & x) const {
    return ((x[1] < xmid_) ? static_cast<Real>(0) : static_cast<Real>(1));
  }

  Real cost(const Real u, const int deriv = 0) const {
    if ( u > T_ ) {
      if ( deriv == 0 ) {
        return static_cast<Real>(0.5)*(u-T_)*(u-T_);
      }
      if ( deriv == 1 ) {
        return (u-T_);
      }
      if ( deriv == 2 ) {
        return static_cast<Real>(1);
      }
    }
    return static_cast<Real>(0);
  }

public:
  QoI_StateCost(const Teuchos::RCP<FE<Real> > &fe,
                Teuchos::ParameterList &parlist) : fe_(fe) {
    xmid_ = parlist.sublist("Geometry").get<Real>("Step height");
    T_    = parlist.sublist("Problem").get("Desired engine temperature",373.0);
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute cost
    Teuchos::RCP<Intrepid::FieldContainer<Real> > costU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costU)(i,j) = cost((*valU_eval)(i,j),0);
      }
    }
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,costU,weight_);
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    grad = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valU_eval, u_coeff);
    // Compute cost
    Teuchos::RCP<Intrepid::FieldContainer<Real> > costU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costU)(i,j) = cost((*valU_eval)(i,j),1);
      }
    }
    // Multiply cost and weight
    Intrepid::FieldContainer<Real> WcostU(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostU,
                                                               *weight_,
                                                               *costU);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  WcostU,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_StateCost::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valU_eval, u_coeff);
    // Evaluate direction on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valV_eval, v_coeff);
    // Compute cost
    Teuchos::RCP<Intrepid::FieldContainer<Real> > costU
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costU)(i,j) = cost((*valU_eval)(i,j),2);
      }
    }
    // Multiply cost and weight
    Intrepid::FieldContainer<Real> WcostU(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostU,
                                                               *weight_,
                                                               *costU);
    // Multiply weighted cost and direction
    Intrepid::FieldContainer<Real> WcostUV(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostUV,
                                                               WcostU,
                                                               *valV_eval);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  WcostUV,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_StateCost::HessVec_22 is zero.");
  }

}; // QoI_L2Tracking

template <class Real>
class QoI_ControlCost : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > weight_;
  Real T_;

  Real weightFunc(const std::vector<Real> & x) const {
    return static_cast<Real>(1);
  }

  Real cost(const Real z, const int deriv = 0) const {
    if ( deriv == 0 ) {
      return static_cast<Real>(0.5)*(z-T_)*(z-T_);
    }
    if ( deriv == 1 ) {
      return (z-T_);
    }
    if ( deriv == 2 ) {
      return static_cast<Real>(1);
    }
    return static_cast<Real>(0);
  }

public:
  QoI_ControlCost(const Teuchos::RCP<FE<Real> > &fe,
                     Teuchos::ParameterList &parlist) : fe_(fe) {
    T_ = parlist.sublist("Problem").get("Desired control temperature",293.0);
    int c = fe_->cubPts()->dimension(0);
    int p = fe_->cubPts()->dimension(1);
    int d = fe_->cubPts()->dimension(2);
    std::vector<Real> pt(d);
    weight_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          pt[k] = (*fe_->cubPts())(i,j,k);
        }
        (*weight_)(i,j) = weightFunc(pt);
      }
    }
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Evaluate control on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Compute cost
    Teuchos::RCP<Intrepid::FieldContainer<Real> > costZ
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costZ)(i,j) = cost((*valZ_eval)(i,j),0);
      }
    }
    // Compute squared L2-norm of diff
    fe_->computeIntegral(val,costZ,weight_);
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    grad = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Compute cost
    Teuchos::RCP<Intrepid::FieldContainer<Real> > costZ
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costZ)(i,j) = cost((*valZ_eval)(i,j),1);
      }
    }
    // Multiply cost and weight
    Intrepid::FieldContainer<Real> WcostZ(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostZ,
                                                               *weight_,
                                                               *costZ);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*grad,
                                                  WcostZ,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    const int c = fe_->gradN()->dimension(0);
    const int f = fe_->gradN()->dimension(1);
    const int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    hess = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
    // Evaluate state on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valZ_eval, z_coeff);
    // Evaluate direction on FE basis
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    fe_->evaluateValue(valV_eval, v_coeff);
    // Compute cost
    Teuchos::RCP<Intrepid::FieldContainer<Real> > costZ
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        (*costZ)(i,j) = cost((*valZ_eval)(i,j),2);
      }
    }
    // Multiply cost and weight
    Intrepid::FieldContainer<Real> WcostZ(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostZ,
                                                               *weight_,
                                                               *costZ);
    // Multiply weighted cost and direction
    Intrepid::FieldContainer<Real> WcostZV(c, p);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(WcostZV,
                                                               WcostZ,
                                                               *valV_eval);
    // Compute gradient of squared L2-norm of diff
    Intrepid::FunctionSpaceTools::integrate<Real>(*hess,
                                                  WcostZV,
                                                  *(fe_->NdetJ()),
                                                  Intrepid::COMP_CPP, false);
  }

}; // QoI_L2Penalty

template <class Real>
class StochasticStefanBoltzmannStdObjective : public ROL::ParametrizedStdObjective<Real> {
private:
  Real alpha0_;
  Real alpha1_;

public:
  StochasticStefanBoltzmannStdObjective(Teuchos::ParameterList &parlist) {
    alpha0_ = parlist.sublist("Problem").get("State Cost",1.0);
    alpha1_ = parlist.sublist("Problem").get("Control Cost",1.0);
  }

  Real value(const std::vector<Real> &x, Real &tol) {
    return alpha0_*x[0] + alpha1_*x[1];
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &x, Real &tol) {
    g[0] = alpha0_;
    g[1] = alpha1_;
  }

  void hessVec(std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol) {
    const Real zero(0);
    hv[0] = zero;
    hv[1] = zero;
  }

}; // OBJ_SCALAR

#endif
