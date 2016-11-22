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

#ifndef PDEOPT_QOI_HELMHOLTZ_HPP
#define PDEOPT_QOI_HELMHOLTZ_HPP

#include "../TOOLS/qoi.hpp"
#include "pde_helmholtz.hpp"

template <class Real>
class QoI_Helmholtz_StateTracking : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > target_;

  Real RoiRadius_;
  Real waveNumber_;
  Real angle_;

protected:
  void computeTarget(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);

    target_.clear(); target_.resize(2);
    target_[0] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
    target_[1] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p));
   
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*target_[0])(i,j) = evaluateRealTarget(x);
        (*target_[1])(i,j) = evaluateImagTarget(x);
      }
    } 
  }

  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }
  
public:
  QoI_Helmholtz_StateTracking(const Teuchos::RCP<FE<Real> > &fe,
                                   const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
                                   Teuchos::ParameterList &parlist)
    : fe_(fe), fieldHelper_(fieldHelper) {
    RoiRadius_  = parlist.sublist("Problem").get("ROI Radius",2.0);
    waveNumber_ = parlist.sublist("Problem").get("Wave Number",10.0);
    angle_      = parlist.sublist("Problem").get("Target Angle",45.0);
    angle_      = DegreesToRadians(angle_);
    computeTarget();
  }

  virtual Real evaluateRealTarget(const std::vector<Real> &x) const {
    const Real arg = waveNumber_ * (std::cos(angle_)*x[0] + std::sin(angle_)*x[1]);
    Real xnorm(0), val(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    val   = (xnorm <= RoiRadius_) ? std::cos(arg) : static_cast<Real>(0);
    return val;
  }

  virtual Real evaluateImagTarget(const std::vector<Real> &x) const {
    const Real arg = waveNumber_ * (std::cos(angle_)*x[0] + std::sin(angle_)*x[1]);
    Real xnorm(0), val(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    val   = (xnorm <= RoiRadius_) ? std::sin(arg) : static_cast<Real>(0);
    return val;
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_[i]);
      fe_->computeIntegral(val,valU_eval,valU_eval,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G(2);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper_->splitFieldCoeff(U, u_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valU_eval, U[i]);
      Intrepid::RealSpaceTools<Real>::subtract(*valU_eval,*target_[i]);
      G[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *valU_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz_StateTracking::gradient_2 is zero.");
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output hessvec
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > H(2);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valV_eval, V[i]);
      H[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *valV_eval,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz::HessVec_22 is zero.");
  }

}; // QoI_Helmholtz


template <class Real>
class QoI_Helmholtz_ControlPenalty : public QoI<Real> {
private:
  const Teuchos::RCP<FE<Real> > fe_;
  const Teuchos::RCP<FieldHelper<Real> > fieldHelper_;

  Real innerAnnulusRadius_;
  Real outerAnnulusRadius_;
  Real RoiRadius_;
  Real alpha_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > ctrlWeight_;

  void computeControlWeight(void) {
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    int d = fe_->gradN()->dimension(3);
   
    ctrlWeight_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));

    const Real one(1);
    std::vector<Real> x(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < d; ++k) {
          x[k] = (*fe_->cubPts())(i,j,k);
        }
        (*ctrlWeight_)(i,j) = (alpha_-one) * evaluateControlWeight(x) + one;
      }
    }
  }

public:
  QoI_Helmholtz_ControlPenalty(const Teuchos::RCP<FE<Real> > &fe,
                               const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
                               Teuchos::ParameterList &parlist)
  : fe_(fe), fieldHelper_(fieldHelper) {
    Real dist2annulus   = parlist.sublist("Problem").get("Distance to Control Annulus",0.5);
    Real annulusWidth   = parlist.sublist("Problem").get("Control Annulus Width",0.1);
    RoiRadius_          = parlist.sublist("Problem").get("ROI Radius",2.0);
    innerAnnulusRadius_ = RoiRadius_ + dist2annulus;
    outerAnnulusRadius_ = innerAnnulusRadius_ + annulusWidth;
    alpha_              = parlist.sublist("Problem").get("Control Penalty",1e-4);
    computeControlWeight();
  }

  virtual Real evaluateControlWeight(const std::vector<Real> &x) const {
    const Real one(1), zero(0);
    Real xnorm(0), val(0);
    const int d = x.size();
    for (int i = 0; i < d; ++i) {
      xnorm += x[i]*x[i];
    }
    xnorm = std::sqrt(xnorm);
    val   = (xnorm <= outerAnnulusRadius_ && xnorm >= innerAnnulusRadius_) ? one : zero;
    return val;
  }

  Real value(Teuchos::RCP<Intrepid::FieldContainer<Real> > & val,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
             const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
             const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int p = fe_->gradN()->dimension(2);
    // Initialize output val
    val = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c));
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valZ_eval, Z[i]);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > wZ
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wZ,*ctrlWeight_,*valZ_eval);
      fe_->computeIntegral(val,wZ,valZ_eval,true);
    }
    Intrepid::RealSpaceTools<Real>::scale(*val,static_cast<Real>(0.5));
    return static_cast<Real>(0);
  }

  void gradient_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::gradient_1 is zero.");
  }

  void gradient_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & grad,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output grad
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > G(2);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > Z;
    fieldHelper_->splitFieldCoeff(Z, z_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valZ_eval, Z[i]);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > wZ
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wZ,*ctrlWeight_,*valZ_eval);
      G[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      Intrepid::FunctionSpaceTools::integrate<Real>(*G[i],
                                                    *wZ,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    fieldHelper_->combineFieldCoeff(grad, G);
  }

  void HessVec_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_11 is zero.");
  }

  void HessVec_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_12 is zero.");
  }

  void HessVec_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    throw Exception::Zero(">>> QoI_Helmholtz_ControlCost::HessVec_21 is zero.");
  }

  void HessVec_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & v_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff = Teuchos::null,
                  const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Get relevant dimensions
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int p = fe_->gradN()->dimension(2);
    // Initialize output hessvec
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > H(2);
    // Get components of the control
    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > V;
    fieldHelper_->splitFieldCoeff(V, v_coeff);
    // Evaluate on FE basis
    for (int i=0; i<2; ++i) {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > valV_eval
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      fe_->evaluateValue(valV_eval, V[i]);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > wV
        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, p));
      Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(*wV,*ctrlWeight_,*valV_eval);
      H[i] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c, f));
      Intrepid::FunctionSpaceTools::integrate<Real>(*H[i],
                                                    *wV,
                                                    *(fe_->NdetJ()),
                                                    Intrepid::COMP_CPP,
                                                    false);
    }

    fieldHelper_->combineFieldCoeff(hess, H);
  }

}; // QoI_Helmholtz_ControlPenalty

#endif
